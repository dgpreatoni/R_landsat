################################################################################
# plots band statistics from selected rasters
# 
# version 1.0
# created prea 20170127
# revised prea 
#
################################################################################

rm(list=ls())

library(rgdal)
library(raster)
library(rgeos)
library(ggplot2)

# GIS data root for this project
dataRoot <- '/lan/gis/Satellite_Belviso'

# Area of Interest polygon shapefile, used to clip out LANDSAT imagery. Must be in BASE directory
AOIShape <- 'area_studio_alpi'

# Polygon layer used to create a "water mask" (in BASE directory)
maskShape <- 'LG_CTR'

# arenas shapefile. Must be in BASE directory
arenaShape <- 'arene_VBVS3'
# field name in arenas shapefile containing arena unique identifier
arenaID <- 'UID'
# field name in arenas shapefile containing first presence year for an arena
arenaFromYear <- 'DA_ANNO'
# field name in arenas shapefile containing most recent presence year for an arena, NULL if still active at present
arenaToYear <- 'A_ANNO'

# buffer widths to use
widths <- c(50, 100, 250, 300, 600, 750, 1000, 1500, 2000) # meters, radius

#### band-by-sensor lookup table
whichBands <- data.frame(sensor=c('LT4', 'LT5', 'LE7', 'LC8'), 
                          Aband=c('B7',  'B7',  'B7',  'B7'), 
                          Bband=c('B5',  'B5',  'B5',  'B6'), stringsAsFactors=FALSE)

whichNDVIBands <- data.frame(sensor=c('LT4', 'LT5', 'LE7', 'LC8'), 
                         Aband=c('B3',  'B3',  'B3',  'B4'), 
                         Bband=c('B4',  'B4',  'B4',  'B5'), stringsAsFactors=FALSE)


#### prepare both LANDSAT data and arenas data
baseDir <- paste(dataRoot, 'BASE', sep='/')
archiveDir <- paste(dataRoot, 'ARCHIVIO', sep='/')
landsatDir <- paste(archiveDir, 'LANDSAT', sep='/')

#### get AOI cookie cutter
AOIClipper <- readOGR(dsn=baseDir, layer=AOIShape)

#### get a list of all available LANDSAT scenes
imageryDirs <- list.dirs(path=landsatDir, full.names=FALSE)[-1]

#### process LANDSAT imagery, one by one
oldwd <- getwd()
indexList <- list()
for(img in imageryDirs) {
  cat('Processing scene', img)
  # get the sensor type from image directory name
  sensor <- substr(img,1,3) 
  # get the year
  year <- substr(img, 10,13)
  # get swath data
  rowpath <- substr(img,4,9)
  cat(' generating', paste(year, rowpath, sep='-'), '\n')
  # calculate an NDVI-based mask
  NDVIbands <- as.character(whichNDVIBands[whichNDVIBands$sensor==sensor,-1])
  NDVIrasterFiles <- paste(paste(img, NDVIbands, sep='_'), 'TIF', sep='.')
  setwd(paste(landsatDir, img, sep='/'))
  NDVIrasters <- stack(NDVIrasterFiles)
  # check for projection and reproject the cookie cutter, if need be
  if(projection(AOIClipper)!=projection(NDVIrasters)) {
    cutter <- spTransform(AOIClipper, CRS(projection(NDVIrasters)))
  } else {
    cutter <- AOIClipper
  }
  NDVIrasters <- crop(NDVIrasters, cutter)
  # set correct datatypes and NA
  if(sensor=="LC8") {
    typ <- "INT2U" # OLI are Uint32 typed
    dataType(NDVIrasters) <- typ
  }
  NAvalue(NDVIrasters) <- 0
  NDVI <- (NDVIrasters[[2]] - NDVIrasters[[1]]) / (NDVIrasters[[2]] + NDVIrasters[[1]])
  NDVImask <- NDVI > 0
  ## calcluate our index
  # which bands do I need?
  bands <- as.character(whichBands[whichBands$sensor==sensor,-1])
  rasterFiles <- paste(paste(img, bands, sep='_'), 'TIF', sep='.')
  rasters <- stack(rasterFiles)
  # check for projection and reproject the cookie cutter, if need be
  if(projection(AOIClipper)!=projection(rasters)) {
    cutter <- spTransform(AOIClipper, CRS(projection(rasters)))
  } else {
    cutter <- AOIClipper
  }
  rasters <- crop(rasters, cutter)
  # set correct datatypes and NA
  if(sensor=="LC8") {
    typ <- "INT2U" # OLI are Uint32 typed
    dataType(rasters) <- typ
  }
  NAvalue(rasters) <- 0
  # apply mask
  rasters <- rasters * NDVImask
  NAvalue(rasters) <- 0
  # basic index
  # indexList[[year]] <- rasters[[1]] / rasters[[2]]
  # alternate form, range constrained between -1 and 1
  indexList[[paste(year, rowpath, sep='-')]] <- (rasters[[2]] - rasters[[1]]) / (rasters[[2]] + rasters[[1]])
  setwd(oldwd)
}
# tidy up
rm(NDVI)
rm(NDVImask)
rm(NDVIrasters)
rm(rasters)
gc()

## mosaic, if any
indexMosaic <- list()
indexNames <- names(indexList)
for(n in 1:(length(indexNames)-1)) {
  year1 <- substr(indexNames[n],1,4)
  cat("using year", year1, '\n')
  for (m in (n+1):length(indexNames)) {
    year2 <- substr(indexNames[m],1,4)
    cat("\t test:", year2)
    if(year1==year2) {
      cat(' match found, merging\n')
      indexMosaic[[year1]] <- mosaic(indexList[[n]], indexList[[m]], fun=max)
    } else {
      cat(' no matches\n')
    }
  }
}
## store index rasters
setwd(baseDir)
for(n in names(indexMosaic)) {
  writeRaster(indexMosaic[[n]], paste("index_", n, ".img", sep=''), format='HFA')
}
## end of file