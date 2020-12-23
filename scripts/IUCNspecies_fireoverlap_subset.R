## IUCN priority species subset


## Set working environment ####
rm(list = ls())
gc()
# system("ps")
# system("pkill -f R")
x <- c("data.table", "stringr", 
       "sp", "raster", "rgdal", "gdalUtils", "rgeos", 
       "alphahull", "ConR", "rnaturalearthdata", 
       "doMC", "foreach",
       "future", "future.apply", "parallel")
lapply(x, require, character.only = TRUE)
rm(x)

## File paths and folders
bugs_data = "~/gsdms_r_vol/tempdata/research-cifs/uom_data/nesp_bugs_data"
iucn_data = "~/gsdms_r_vol/tempdata/research-cifs/uom_data/iucn_bugs_data"

nesp_output_dir = file.path(bugs_data, "outputs")
iucn_output_dir = file.path(iucn_data, "outputs")


## I. Find data for species ####
## >> IUCN species list (see updates in IUCN_specieslist.xlsx in Dropbox) ####
iucn_species <- fread(file.path(iucn_data, "IUCN_specieslist.csv"))
message(cat("Total number of species in updated IUCN list: "),
        nrow(iucn_species))

## >> Add ScientificName column to IUCN list ####
iucn_species$ScientificName <- gsub("  ", " ", tolower(paste(iucn_species$Genus, 
                                                           iucn_species$SubGenus,
                                                           iucn_species$Species)))

## >> Find and remove incomplete and improper names ####
source("/tempdata/workdir/nesp_bugs/scripts/remove_improper_names.R")
species_record <- remove_improper_names(as.character(iucn_species$ScientificName),
                                        allow.higher.taxa = FALSE,
                                        allow.subspecies = TRUE)
message(cat("# Improper species names found in ALA data: "),
        length(species_record$improper_species))
message(cat("# Incomplete species names found in ALA data: "),
        length(species_record$incomplete_species))
message(cat("Duplicates in cleaned ALA species list: "),
        length(species_record$updated_list[duplicated(species_record$updated_list)]))

## >> Subset list if improper or incomplete names found ####
if(!length(species_record$updated_list) == length(iucn_species$ScientificName)){
  iucn_species <- iucn_species[ScientificName == species_record$updated_list]
}
iucn_species$ScientificName <- gsub(" ", "_", iucn_species$ScientificName)

## >> Load cleaned and masked ALA data for IUCN species ####
spmasked_dir <- file.path(nesp_output_dir, "ala_data" ,"spdata_masked")

datfiles <- list.files(spmasked_dir, pattern = "_masked.rds$",
                        full.names = TRUE, all.files = TRUE)
x <- as.character(sort(iucn_species$ScientificName))
datfiles <- datfiles[grep(paste(x,collapse="|"), datfiles)]

y <- basename(tools::file_path_sans_ext(datfiles))
y <- gsub("_masked", "", y)
message(cat("Number of IUCN species found in cleaned/masked ALA data: "),
        length(x[x %in% y]))

message(cat("Number of IUCN species not found in cleaned/masked ALA data: "),
        length(x[!x %in% y]))
message(cat("IUCN species not found in cleaned/masked ALA data: "))
x[!x %in% y]


## II. Species polygons ####
## Using ConR package - minimum convex polygon or alpha hulls
## https://cran.r-project.org/web/packages/ConR/index.html
## Notes: 
##  WGS84 required for data
##  data format: latitude, longitude (in decimal degrees), and taxon name
##  outputs for EOO and AOO same as from {red} when method.range = "convex.hull"
## IUCN.eval() requires working directory to be set

polygons_dir = file.path(iucn_output_dir,"polygons")
if(!dir.exists(polygons_dir)){dir.create(polygons_dir)}
working_dir <- polygons_dir

source("~/gsdms_r_vol/tempdata/workdir/nesp_bugs/scripts/conr_iucn_eval.R")

spfiles <- datfiles
basemap_file <- file.path(nesp_output_dir, "masks", "auslands_1poly_wgs84.shp")

## Package: future - for catching errors
plan(multiprocess, workers = future::availableCores()-2)
options(future.globals.maxSize = +Inf) ## CAUTION
errorlog <- paste0(iucn_output_dir, "/errorlog_polygons_", gsub("-", "", Sys.Date()), ".txt")
# if(file.exists(errorlog)){unlink(errorlog)}
writeLines(c(""), errorlog)

system.time(
  suppressWarnings(
    future.apply::future_lapply(
      spfiles,
      function(x){
        tmp <- tryCatch(expr = conr_iucn_eval(species_filename = x,
                                              basemap_path = basemap_file,
                                              working_dir = working_dir,
                                              iucn_outpath = polygons_dir),
                        error = function(e) {
                          cat(
                            paste(as.character(x), "\n"),
                            file = errorlog,
                            append = TRUE)
                        }
        )
      }, future.seed = TRUE)))

## >> Error checking ####
message(cat("Species showing errors: "),
        trimws(readLines(errorlog)))

## >> Check files ####
csvfiles <- list.files(polygons_dir, pattern = ".csv$",
                       full.names = TRUE, all.files = TRUE)
message(cat("Number of .csv output files created from IUCN.eval(): "),
        length(csvfiles))
rdsfiles <- list.files(polygons_dir, pattern = ".rds$", 
                       full.names = TRUE, all.files = TRUE)
message(cat("Number of .rds output files created from IUCN.eval(): "),
        length(rdsfiles))
IUCNshpfiles <- list.files(file.path(polygons_dir, "shapesIUCN"), 
                           pattern = ".shp$", 
                           full.names = TRUE, all.files = TRUE)
message(cat("Number of .shp output files created from IUCN.eval(): "),
        length(IUCNshpfiles))
pngfiles <- list.files(polygons_dir, pattern = "png$", recursive = TRUE, 
                       full.names = TRUE, all.files = TRUE)
message(cat("Number of .png map files created from IUCN.eval(): "),
        length(pngfiles))

## >> Save outputs ####
## Table with EOO and AOO area estimates
csvfiles <- list.files(polygons_dir, pattern = ".csv$",
                       full.names = TRUE, all.files = TRUE)
out <- do.call("rbind", lapply(csvfiles , read.csv))
dim(out)
setorder(out, EOO, AOO)
out <- as.data.table(out)
write.csv(out, file = file.path(iucn_output_dir, "polygon_areas.csv"), 
          row.names = FALSE)

message(cat("#species with EOOs: "),
        nrow(out[!is.na(EOO)]))
message(cat("#species without EOOs: "),
        nrow(out[is.na(EOO)]))
message(cat("max #records for species without EOOs: "),
        max(out[is.na(EOO)]$Nbe_unique_occ.))

## List with polygons
rdsfiles <- list.files(polygons_dir, pattern = ".rds$", 
                       full.names = TRUE, all.files = TRUE)
polynames <- basename(tools::file_path_sans_ext(rdsfiles))
polylist <- lapply(rdsfiles, readRDS)
names(polylist) <- polynames
length(polylist)
polylist <- lapply(polylist, "[[", 1)
polylist <- lapply(polylist, "[[", 2)
length(polylist)
saveRDS(polylist, file = file.path(iucn_output_dir, "polygons.rds"))

## >> Move png files to maps folder
maps_dir <- file.path(iucn_output_dir, "EOO_AOO_maps")
if(!dir.exists(maps_dir)){dir.create(maps_dir)}

pngfiles <- list.files(polygons_dir, pattern = "png$", recursive = TRUE, 
                       full.names = TRUE, all.files = TRUE)
x <- unlist(str_extract_all(pngfiles, "polygons/\\w+"))
x <- gsub("polygons/", "", x)
x <- gsub("results_map$", ".png", x)
duplicated(x)
x <- file.path(maps_dir, x)
file.copy(pngfiles, x)
dir.create(file.path(maps_dir, "shapefiles"))
system(paste0("cp -R ", file.path(polygons_dir, "shapesIUCN/"), " ", maps_dir, "/shapefiles"))

## III. Fire overlap analysis - EOOs ####
shapefile_dir = file.path(iucn_output_dir, "species_shapefiles")
if(!dir.exists(shapefile_dir)){dir.create(shapefile_dir)}

overlap_dir = file.path(iucn_output_dir, "overlaps")
if(!dir.exists(overlap_dir)){dir.create(overlap_dir)}

source("/tempdata/workdir/nesp_bugs/scripts/overlap.R")

## >> Load species data ####
# polylist <- readRDS(file.path(iucn_output_dir, "polygons.rds"))

## Find species with EOOs in list
id <- which(!sapply(polylist, length) == 0)
message(cat("# Number of species with EOOs: "),
        length(id))

## Load species polygons
species_maps <- polylist[id]
polygon_list <- names(species_maps)

## >> Load in fire severity raster (re-classed) and get unique classes ####
fire_severity <- raster(file.path(nesp_output_dir, "fire", "severity3_eqar250.tif"))
fire_vals <- fire_severity[]
fire_classes <- sort(unique(na.omit(fire_vals)))

## >> Run overlap analysis in parallel: doMC ####
registerDoMC(future::availableCores()-2)
system.time(log <- foreach(polys = polygon_list, 
                           .combine = rbind,
                           .errorhandling = "pass",
                           .packages = c('sp', 'raster',
                                         'rgdal', 'data.table')) %dopar%{
                             
                             overlap(species_name = polys,
                                     species_poly = species_maps[[polys]],
                                     shapefile_dir = shapefile_dir,
                                     fire_vals = fire_vals,
                                     fire_classes = fire_classes, 
                                     outdir = overlap_dir)
                           })


## >> Error checking ####
log

csvfiles <- list.files(overlap_dir, pattern = ".csv$",
                       full.names = TRUE, all.files = TRUE)
message(cat("Number of input species: "),
        length(polygon_list))
message(cat("Number of output files: "),
        length(csvfiles))
out <- do.call("rbind", lapply(csvfiles , fread))
dim(out)
setorder(out, Species)
out <- as.data.table(out)

## Find missing species from outputs
error1_list <- error_list <- polygon_list[!polygon_list %in% out$Species]

## Reruns 
## Repeat this till most of the errors are fixed
## Errors seem to be an artefact of the system rather than problem with data/code
rm(log)
if(length(error_list) <= future::availableCores()-2) {
  registerDoMC(length(error_list))
} else{
  registerDoMC(future::availableCores()-2)
}
system.time(log <- foreach(polys = error_list, 
                           .combine = rbind,
                           .errorhandling = "pass",
                           .packages = c('sp', 'raster',
                                         'rgdal', 'data.table')) %dopar%{
                                           
                                           overlap(species_name = polys,
                                                   species_poly = species_maps[[polys]],
                                                   shapefile_dir = shapefile_dir,
                                                   fire_vals = fire_vals,
                                                   fire_classes = fire_classes, 
                                                   outdir = overlap_dir)
                                         })

csvfiles <- list.files(overlap_dir, pattern = ".csv$",
                       full.names = TRUE, all.files = TRUE)
message(cat("Number of input species: "),
        length(polygon_list))
message(cat("Number of output files: "),
        length(csvfiles))
temp <- basename(tools::file_path_sans_ext(csvfiles))
error_list[!error_list %in% temp]

## >> Save overlap areas output table ####
csvfiles <- list.files(overlap_dir, pattern = ".csv$",
                       full.names = TRUE, all.files = TRUE)
out <- do.call("rbind", lapply(csvfiles , fread))
setorder(out, Species)
out <- as.data.table(out)
write.csv(out, file = file.path(iucn_output_dir, "EOO_fireoverlap.csv"), row.names = FALSE)
file.remove(file.path(overlap_dir, dir(path = overlap_dir)))


## III. Fire overlap analysis - point data ####

## >> Load in fire severity raster (re-classed) and get unique classes ####
fire_severity <- raster(file.path(nesp_output_dir, "fire", "severity3_eqar250.tif"))
fire_vals <- fire_severity[]
fire_classes <- sort(unique(na.omit(fire_vals)))

## >> Load species data ####
## Species names without EOOs  
polylist <- readRDS(file.path(iucn_output_dir, "polygons.rds"))
id <- which(sapply(polylist, length) == 0)
message(cat("# Number of species without EOOs: "),
        length(id))

## Find & load cleand/masked species data files
points_list <- names(polylist[id]) ## species names
x <- basename(tools::file_path_sans_ext(datfiles))
x <- gsub("_masked", "", x)
y <- datfiles[x %in% points_list] ## subset datfiles for all IUCN species

points_dat <- do.call("rbind", lapply(y , readRDS))
message(cat("Number of species in data: "),
        length(unique(points_dat$scientificName)))

## >> Find and remove duplicates from data ####
## by name-lat-long (same as ConR)
points_dat <- points_dat[!duplicated(points_dat[ , c("scientificName", "longitude", "latitude")]), ]
length(unique(points_dat$scientificName))

## Check if text modified names are the same as scientifcNames
length(unique(points_dat$scientificName)) == length(unique(points_dat$spfile))
cbind(c(unique(points_dat$scientificName)), unique(points_dat$spfile))

## >> Extract fire severity information for species data points ###
## Subset species data
## Note: Here, include fields required to identify senstivie status, uncertainty, other...
points_datsub <- as.data.frame(points_dat[, c("scientificName", "longitude", "latitude",
                                            "phylum", "class", "order", "family", 
                                            "genus" , "eventDate", 
                                            "coordinateUncertaintyInMetres", "assertions", 
                                            "dataGeneralizationsOriginal", "sensitive", 
                                            "spfile", "id")])

## Project data points in equal area projection
wgs_crs  <-  "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
eqarea_crs <- "+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=134 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
spdf <- SpatialPointsDataFrame(coords = points_datsub[c("longitude", "latitude")], data =points_datsub, proj4string = CRS(wgs_crs))
spdf <- spTransform(spdf, CRSobj = CRS(eqarea_crs))

## Extract fire severity values for data points
# spdf$FireClass <- extract(fire_severity, spdf)
points_datsub$FireClass <- extract(fire_severity, spdf)
  # temp <- points_datsub[, .N, by = .(spfile , FireClass)]
  # temp <- na.omit(temp)

## >> Create output table ####
points_datsub <- as.data.table(points_datsub)
spnames <- unique(points_datsub$spfile)

df <- data.frame(matrix(nrow = length(spnames), 
                        ncol = length(fire_classes) + 3))
df[ , 1] <- spnames
colnames(df) <- c("Species", paste0("Fire_Class_", fire_classes), 
                  "Total_Overlap", "Occurrence_Points")
df[ , -c(1, ncol(df)-1, ncol(df))] <- t(sapply(spnames, FUN = function(y) sapply(fire_classes, FUN = function(x) points_datsub[spfile == y & FireClass == x, length(FireClass)])))
df[ , ncol(df)-1] <- rowSums(df[, -c(1, ncol(df)-1, ncol(df))])
df[, ncol(df)] <- points_datsub[, .N, by = spfile]$N

## >> Error checking ####
all(df$Total_Overlap <= df$Occurrence_Points)

## Save output
write.csv(df, file = file.path(iucn_output_dir, "Points_fireoverlap.csv"), row.names = FALSE)



## EXTRAS ------------------------------ ####
## To check spdf
# quickPlot::clearPlot()
# quickPlot::Plot(fire_severity,
#                 title = "",
#                 axes = FALSE,
#                 legend = FALSE,
#                 addTo = "fire_severity",
#                 new = TRUE)
# quickPlot::Plot(spdf, pch = 19,
#                 col = "hotpink",
#                 title = "",
#                 addTo = "fire_severity")

  
