## IUCN priority species subset


## Set working environment ####
rm(list = ls())
gc()
# system("ps")
# system("pkill -f R")
x <- c("data.table", "stringr", "sp", "raster", "rgdal", 
       "gdalUtils", "rgeos", "doMC", "foreach")
lapply(x, require, character.only = TRUE)
rm(x)

## File paths and folders
bugs_data = "~/gsdms_r_vol/tempdata/research-cifs/uom_data/nesp_bugs_data"
output_dir = file.path(bugs_data, "outputs")
polygons_dir = file.path(output_dir,"polygons")
iucn_dir = file.path(output_dir,"IUCN_maps")
if(!dir.exists(iucn_dir)){dir.create(iucn_dir)}

shapefile_dir = file.path(output_dir, "species_shapefiles")
if(!dir.exists(shapefile_dir)){
  dir.create(shapefile_dir)
}

overlap_dir = file.path(output_dir, "overlaps")
if(!dir.exists(overlap_dir)){dir.create(overlap_dir)}

source("/tempdata/workdir/nesp_bugs/scripts/overlap.R")

## Load data ####
## ALA speices polygon information
## IUCN.eval output table for EOO and AOO estimates, inlcuding NA
## Gives list of ALA species (with at least 1 record, masked by land)
## See 'Species polygons (range): ALA data' in data processing document
ala_polys <- fread(file.path(output_dir, "ala_polygons_areas.csv"))


## Processing IUCN priority species list (provided by JM) ####
# iucn_species <- fread(file.path(bugs_data, "IUCN_species/IUCN_specieslist.csv"))
# 
# ## Add columns to IUCN list
# iucn_species$ScientificName <- paste(iucn_species$Genus, iucn_species$Species)
# iucn_species$simple_name <- stringr::word(iucn_species$ScientificName, 1,2)
# 
# ## Name matching I - only done once
# ## Match ScientificName against ala_polys
# dim(iucn_species[iucn_species$ScientificName  %in% ala_polys$X])
# dim(iucn_species[!iucn_species$ScientificName  %in% ala_polys$X])
# 
# ## Match simple_name against ala_polys 
# dim(iucn_species[iucn_species$simple_name %in% ala_polys$X])
# dim(iucn_species[!(iucn_species$simple_name %in% ala_polys$X)])
# 
# ## Write nesp_name column to IUCN list 
# y <- iucn_species[iucn_species$simple_name %in% ala_polys$X]
# x <- iucn_species[!(iucn_species$simple_name %in% ala_polys$X)]
# dat <- rbind(y, x)
# dat$nesp_data <- c(rep(1, nrow(y)), rep(0, nrow(x)))
# dat$nesp_names <- rep(character(), nrow(dat))
# dat[nesp_data == 1]$nesp_names = dat[nesp_data == 1]$simple_name
# 
# ## Write updated IUCN list as csv
#   # write.csv(dat, file = file.path(bugs_data, "IUCN_specieslist_updated.csv"), 
#   #           row.names = FALSE)
# 
# 
# ## Name matching II - only done once
# ## Manually check for names not found in ala_polys and update 'nesp_names' in IUCN_specieslist_updated.xlsx
# ## Save file as IUCN_specieslist_updated.csv
# x <- iucn_species[!(iucn_species$simple_name %in% ala_polys$X)]
# for (i in 121:124){
#   print(i)
#   
#   message(cat("Looking for Genus: "),
#           x[i,]$ScientificName)
#   print(x[i,])
#   
#   message("Macthes found in ALA data: ")
#   ## Subset by genus matches
#   z <- ala_polys[grep(stringr::word(x[i,]$ScientificName, 1,1), ala_polys$X)]$X
#   ## Find species matches in above
#   print(z[grep(stringr::word(x[i,]$ScientificName, 2,2), z)])
#   
#   print("------------------------------------------------------------")
#   print("------------------------------------------------------------")
# }
# 
# 
# z <- ala_polys[grep("Lasioglossum", ala_polys$X)]$X
# z[grep("grumiculum", z)]
# 
# ## Reload updated IUCN list
# iucn_species <- fread(file.path(bugs_data, "IUCN_species/IUCN_specieslist_updated.csv")) 
# 
# 
# ## Find and fix duplicates in IUCN_specieslist_updated.xlsx ####
# message(cat("Duplicates found in updated list: "),
#         sum(duplicated(iucn_speciessub)))
# iucn_speciessub[duplicated(iucn_speciessub)]
# # iucn_speciessub[which(duplicated(iucn_speciessub) | duplicated(iucn_speciessub[length(iucn_speciessub):1]) [length(iucn_speciessub):1])]
# 
# ala_polys[grep("Oreixenica orichora", ala_polys$X)]


## Load updated list ####
iucn_species <- fread(file.path(bugs_data, "IUCN_species/IUCN_specieslist_updated.csv"))
message(cat("Total number of species in updated IUCN list: "),
        nrow(iucn_species))


## Subset list of species found in ala_polys #### 
## NOTE: Rerun IUCN.eval() analysis for IUCN species to have within one IUCN project
iucn_speciessub <- iucn_species[iucn_species$nesp_names %in% ala_polys$X]$nesp_names

message(cat("Number of IUCN species found in our analysis (polygons list) - subset list: "),
        length(iucn_speciessub))
message(cat("Duplicates found in subset list : "),
        sum(duplicated(iucn_speciessub)))

message(cat("Number of IUCN species NOT found in subset list: "),
        sum(!(iucn_species$nesp_names %in% ala_polys$X)))


## Find & copy png files for spcies in subset list ####
pngfiles <- read.table(file.path(output_dir, "png_filenames.csv"), 
                       as.is = TRUE,
                       header = FALSE)
pngfiles <- pngfiles[-1,]
length(pngfiles)
# x1 <- grep("[A-Z]", str_split_fixed(pngfiles, " ", 3)[,2], value = TRUE)
# x1 <- unlist(str_extract_all(x1, "[A-Z]\\w+"))
# x2 <- str_split_fixed(pngfiles, " ", 3)[,3]
# pngnames <- paste0(x1, " ", x2)
# pngnames <- gsub(".png$", "", pngnames)

iucn_files <- list()
for (i in 1:length(iucn_speciessub)){
  print(paste0("Looking for...", iucn_speciessub[i]))
  temp <- list(pngfiles[grep(iucn_speciessub[i], pngfiles)])
  if(length(temp[[1]]) == 0){
    iucn_files[[i]] <- NA }
  else {
    iucn_files[[i]] <- temp
  }
}
message(cat("Number of species in subset IUCN list with png files: "),
        sum(!is.na(iucn_files)))
message(cat("Species with files: "))
iucn_speciessub[!is.na(iucn_files)]

## Make list of files to extract/copy
x <- unlist(iucn_files)
x <- x[!is.na(x)]
duplicated(x)
x <- unique(x)
x2 <- unlist(str_extract_all(x, "polygons/\\w+"))
x2 <- gsub("polygons/", "", x2)
x2 <- gsub("results_map$", ".png", x2)
duplicated(x2)
x2 <- file.path(iucn_dir, x2)
file.copy(x, x2)



## Overlap analysis ####
## Load in fire severity raster (re-classed) and get unique classes
fire_severity <- raster(file.path(output_dir, "fire", "severity3_eqar250.tif"))
fire_vals <- fire_severity[]
fire_classes <- sort(unique(na.omit(fire_vals)))

## Load in species rds
species_maps <- readRDS(file.path(output_dir, "ala_EOO.rds"))

## Subset list species - text modification
## Polygons information for subset list ####
iucn_polys <- ala_polys[ala_polys$X %in% iucn_speciessub]
message(cat("Number of species in subset list without EOOs: "),
        nrow(iucn_polys[is.na(EOO)]))
message(cat("Number of species in subset list with EOOs: "),
        nrow(iucn_polys[!is.na(EOO)]))

## Species with EOOs in subset list
iucn_polys_eoo <- iucn_polys[!is.na(EOO)]
polygon_list <- iucn_polys_eoo$X
polygon_list <- stringr::str_replace_all(polygon_list, " ", "00xx00")
message(cat("Check for duplicates: "),
        sum(duplicated(polygon_list)))
polygon_list <- str_replace_all(polygon_list, "[^[:alnum:]]", "")

# polygon_list2 <- gsub("\\s*\\([^\\)]+\\)", "", polygon_list)
#  ## remove everythig between round brakets
message(cat("Check for duplicates: "),
        sum(duplicated(polygon_list)))
polygon_list <- tolower(gsub("00xx00", "_", polygon_list))
message(cat("Check for duplicates: "),
        sum(duplicated(polygon_list)))
polygon_list

## Check polygons exist
for (i in 1:length(polygon_list)){
  print(i)
  print(species_maps[[polygon_list[i]]])
  print("---------------------------------")
}


## Run overlap analysis in parallel: doMC ####
registerDoMC(future::availableCores()-2)
system.time(foreach(polys = polygon_list, 
                    .errorhandling = "remove",
                    .packages = c('sp', 'raster', 'rgdal', 'data.table')) %dopar%{
                      
                      ## write spatial data to disk if coming from RDS file
                      writeOGR(species_maps[[polys]], dsn = shapefile_dir, 
                               layer = polys, driver = "ESRI Shapefile", 
                               overwrite_layer = TRUE)
                      
                      ## reproject species boundaries to match fire severity raster
                      system(paste0("gdal_rasterize -at -burn 1 -ot Byte -tr .0025 .0025 -l ",
                                    polys, " ",
                                    file.path(shapefile_dir, polys), ".shp ",
                                    file.path(shapefile_dir, polys), ".tif"))
                      
                      system(paste0("gdalwarp -overwrite -ot Byte -te -2214250 -4876750 2187750 -1110750 -tr 250 250 -s_srs 'EPSG:4326' -t_srs '+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=134 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs' ",
                                    file.path(shapefile_dir, polys), ".tif ",
                                    file.path(shapefile_dir, polys), "_p.tif"))
                      
                      ## Create table of areas within each fire class
                      species_map <- raster(paste0(file.path(shapefile_dir, 
                                                             polys), "_p.tif"))
                      
                      dt <- data.table("species_map" = species_map[],
                                       "fire_severity" = fire_vals)
                      
                      df <- data.frame(matrix(ncol = length(fire_classes) + 3))
                      df[ , 1] <- polys
                      df[ , -c(1, ncol(df)-1, ncol(df))] <- 
                        sapply(fire_classes, FUN = function(x) dt[species_map == 1 & fire_severity == x, length(fire_severity) * 250 * 250 / 1000000])
                      df[ , ncol(df)-1] <- rowSums(df[, -c(1,ncol(df)-1, ncol(df))])
                      df[ , ncol(df)] <- dt[species_map == 1, length(species_map)* 250 * 250 / 1000000]
                      colnames(df) <- c("Species", 
                                        paste0("Fire_Class_", fire_classes), 
                                        "Total_Overlap", "Species_Polygon")
                      
                      ## Remove files
                      file.remove(file.path(shapefile_dir, dir(path = shapefile_dir)))
                      
                      ## Save output file
                      write.csv(df, file = paste0(file.path(overlap_dir, polys), ".csv"), row.names = FALSE)
                      
                    })


## Check output files ####
csvfiles <- list.files(overlap_dir, pattern = ".csv$",
                       full.names = TRUE, all.files = TRUE)
message(cat("Number of input species: "),
        length(polygon_list))
message(cat("Number of output files: "),
        length(csvfiles))


## Combine outputs & save as csv ####
out <- do.call("rbind", lapply(csvfiles , fread))
dim(out)
setorder(out, Species)
out <- as.data.table(out)
write.csv(out, file = file.path(output_dir, "IUCNsp_EOO_fireoverlap.csv"), 
          row.names = FALSE)
file.remove(file.path(overlap_dir, dir(path = overlap_dir)))

## Errors ####
## Find missing species from outputs
error1_list <- error_list <- polygon_list[!polygon_list %in% out$Species]

## Reruns - I
if(length(error_list) <= future::availableCores()-2) {
  registerDoMC(length(error_list))
} else{
  registerDoMC(future::availableCores()-2)
}

system.time(foreach(polys = error_list, 
                    .errorhandling = "remove",
                    .packages = c('sp', 'raster', 'rgdal', 'data.table')) %dopar%{
                      
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
        length(error_list))
message(cat("Number of output files: "),
        length(csvfiles))

temp <- basename(tools::file_path_sans_ext(csvfiles))
error_list <- error_list[!error_list %in% temp]


csvfiles <- list.files(overlap_dir, pattern = ".csv$",
                       full.names = TRUE, all.files = TRUE)
out <- do.call("rbind", lapply(csvfiles , fread))
setorder(out, Species)
out <- as.data.table(out)
write.csv(areas, file = file.path(output_dir, "IUCNsp_EOO_fireoverlap_error1.csv"), row.names = FALSE)
file.remove(file.path(overlap_dir, dir(path = overlap_dir)))



## Combine outputs & save as csv ####

## Species without EOOs in subset list
iucn_points <- iucn_polys[is.na(EOO)]
points_list <- iucn_points$X
points_list <- stringr::str_replace_all(points_list, " ", "00xx00")
message(cat("Check for duplicates: "),
        sum(duplicated(points_list)))
points_list <- str_replace_all(points_list, "[^[:alnum:]]", "")
message(cat("Check for duplicates: "),
        sum(duplicated(polygon_list)))
points_list <- tolower(gsub("00xx00", "_", points_list))
message(cat("Check for duplicates: "),
        sum(duplicated(points_list)))
points_list

spmasked_dir <- file.path(output_dir, "ala_data" ,"spdata_masked")
dat_files <- paste0(file.path(output_dir, "ala_data" ,"spdata_masked"), points_list, ".rds")
file.exists(dat_files)




# x1 <- grep("[A-Z]", str_split_fixed(pngfiles, " ", 3)[,2], value = TRUE)
# x1 <- unlist(str_extract_all(x1, "[A-Z]\\w+"))
# x2 <- str_split_fixed(pngfiles, " ", 3)[,3]
# pngnames <- paste0(x1, " ", x2)
# pngnames <- gsub(".png$", "", pngnames)


## Find files not found
names(iucn_files) <- iucn_speciessub
message(cat("Files not found for... "),
        length(names(iucn_files)[is.na(iucn_files)]))
names(iucn_files)[is.na(iucn_files)]

pngfiles[grep("Brachyhesma", pngfiles)]



## Text modification
spfilename <- stringr::str_replace_all(names(iucn_files), " ", "00xx00")
message(cat("Duplicates found in spfilename: "),
        sum(duplicated(spfilename)))

spfilename <- str_replace_all(spfilename, "[^[:alnum:]]", "")
message(cat("Duplicates found in spfilename: "),
        sum(duplicated(spfilename)))

spfilename <- tolower(gsub("00xx00", "_", spfilename))
message(cat("Duplicates found in spfilename: "),
        sum(duplicated(spfilename)))
spfilename <- unique(spfilename)
length(unique(spfilename))

message(cat("Files found for... "),
        length(names(iucn_files)[!is.na(iucn_files)]))

## Find species with more than one files allocated
temp <- lapply(iucn_files, "[[", 1)
temp2 <- sapply(temp, length)
which(temp2 > 1)

