################################################################################
# mosaicking of indices from LANDSAT imagery
# 
# version 1.0
# created prea 20170317 happy birthday!
# revised 
#
# to be run on eira
#
# Assumes that index files are stored in localwksp.
# identifies each scene by row/path and timestamp.
# scenes having the same row/path and belonging to the seme time span (by years) are averaged.
# averaged rasters are then mosaicked, yielding a single raster.

# Note that if re-projecting and mosaicking is really a large part of your project, you may want to consider using the gdalwarp command line utility (gdalwarp) directly. The gdalUtils R package provides utilities to run GDAL commands from R, including gdalwarp, for reprojection, resampling and mosaicking.
#
################################################################################

rm(list=ls())

source('/home/gis/R_scripts/raster_utilities.R')
#library(raster)  # already invoked by raster_utilities.R
#library(rgeos)  # already invoked by raster_utilities.R
#library(RStoolbox) # satellite stuff...
library(rgdal)
#library(readxl)
#library(ggplot2)
library(raster)
library(car) # recode!


# how many CPUs to use when going parallel
numCores <- 7

# Storage directory for raw uncompressed LANDSAT files
landsatDir <- '/home/ARCHIVIO/ARCHIVIO/cartografia/Europa/LANDSAT_Alpi'

# GIS data root for this project
dataRoot <- '/home/gis/Satellite_Forcello'

# local directory to be used as workspace (WARNING: a temporary directory is also creatred here!)
localWksp <- '/srv/data/w_Satellite_Forcello'

# Area of Interest polygon shapefile, used to clip out LANDSAT imagery. Must be in BASE directory
AOIShape <- 'AOI_Alps_Convention_buffer10k_extent'

# Polygon layer used to create a mask (in BASE directory)
maskShape <- 'AOI_Alps_Convention_buffer10k'

## Indices lookup table, an Excel fie with Index	Sensor	BandA	BandB (in data directory)
#indexLookupFile <- 'Landsat indices.xls'
#indexLookupSheet <- 'Indices'

#### prepare standard paths
baseDir <- paste(dataRoot, 'BASE', sep='/')
archiveDir <- paste(dataRoot, 'ARCHIVIO', sep='/')
dataDir <- paste(dataRoot, 'data', sep='/')

#### get AOI cookie cutter
AOIClipper <- readOGR(dsn=baseDir, layer=AOIShape)

#### get mask vector
maskPoly <- readOGR(dsn=baseDir, layer=maskShape)
  

#### carefully set work directory _and_ raster temporary directory: this obviously EATS disk space!
oldwd <- getwd()
oldtmpdir <- rasterOptions()$tmpdir
tmpDir <- tempfile(pattern='raster', tmpdir=localWksp)
setwd(localWksp)
rasterOptions(tmpdir=tmpDir)


#### get a list of all available LANDSAT index raster
imageryFiles <- list.files(path=localWksp, pattern="\\.img$", full.names=FALSE)

# parse index rasters and set up a processing schedule dataframe
schedule <- data.frame()
metaData <- lapply(imageryFiles, parseLandsatName)
names(metaData) <- imageryFiles
for(i in names(metaData)) {
  idx <- strsplit(i, "[\\._]")[[1]][[2]]
  year <- metaData[[i]]['year']
  month <- metaData[[i]]['month']
  pathrow <- paste(metaData[[i]]['path'], metaData[[i]]['row'], sep="")
  schedule <- rbind(schedule, data.frame(file=i, index=idx, year=year, month=month, pathrow=pathrow))
}

# create indexes to aggregate rasters
schedule$epoch <- recode(schedule$year, "1984:1987='past'; 2014:2016='present'")

# aggregate schedule
schedule.split <- split(schedule, list(index=schedule$index, epoch=schedule$epoch, pathrow=schedule$pathrow))

# set up parallel processing
beginCluster(numCores)
for(s in names(schedule.split)) {
  cat("aggregating", s, '\n')
  fname <- paste(gsub("\\.", "_", s), 'grd', sep='.')
  if(length(list.files(pattern=fname))==1) { # that scene has been already processed
    cat("\tscene has alreday been processed, skipping.\n")  
  } else { # process it...
    sceneList <- schedule.split[[s]]
    if(nrow(sceneList)>0) { # we have something to chew for that pathrow*epoch
      files <- as.character(sceneList$file)
      # protect against different extents, resolutions and projections
      # here we use minimum common extent to avoid boundery conditions
      extents <- lapply(files, function(x) extent(raster(x)))
      xmin <- max(unlist(lapply(extents, function(x) x@xmin)))
      xmax <- min(unlist(lapply(extents, function(x) x@xmax)))
      ymin <- max(unlist(lapply(extents, function(x) x@ymin)))
      ymax <- min(unlist(lapply(extents, function(x) x@ymax)))
      bestextent <- extent(xmin, xmax, ymin, ymax)
      bestrows <- min(unlist(lapply(files, function(x) nrow(raster(x)))))
      bestcols <- min(unlist(lapply(files, function(x) ncol(raster(x)))))
      bestres <- min(unlist(lapply(files, function(x) res(raster(x)))))
      bestCRS <- unlist(lapply(files, function(x) projection(raster(x))))[1]
      #baseRaster <- raster(bestextent, nrows=bestrows, ncols=bestcols, crs=CRS(bestCRS), resolution=bestres)
      indexStack <- stack()
      for(r in files) {
        rr <- raster(r)
        # this takes ages
        # rr <- extend(rr, baseRaster)
        # faster?
        rr <- crop(rr, bestextent)
        indexStack <- addLayer(indexStack, rr)
      }
      avg <- mean(indexStack)
      writeRaster(avg, filename=fname, overwrite=TRUE)
    } else {
      cat("\tno scenes for that row/path and epoch combination, skipping.\n")  
    }
  }
}

# stop computation: exit from parallel
cat("closing parallel mode...\n")
endCluster()
gc()

# ensure projection is the same
#rawTiles <- list.files(path=localWksp, pattern="\\.grd$", full.names=FALSE)
# time-averaged tiles must follow the pattern below 
rawTiles <- list.files(path=localWksp, pattern="ND[TV]I_(past|present)_[0-9]{6}\\.grd$", full.names=FALSE)
theCRS <- CRS("+init=epsg:3035")
beginCluster(numCores)
for(tile in rawTiles) {
  tiffName <- paste(substr(tile, 1, nchar(tile)-4), '3035', sep='_')
  if(length(list.files(pattern=paste(tiffName, '\\.tif', sep='')))>0){ # already projected
    cat('\tscene', tile, 'has already been projected. Skipping.\n')
  } else {
    t <- raster(tile)
    cat("projecting", tile, '\n')
    if(!compareCRS(t, theCRS)) {
      t <- projectRaster(t, crs=theCRS, res=res(t))
      #writeRaster(t, filename=paste(substr(tile, 1, nchar(tile)-4), '3035', sep='_'))
      writeRaster(t, filename=tiffName, format='GTiff', overwrite=TRUE) # better than native raster, can be read by other software
    }
  }
}
endCluster()
cat("\nDone.\n")  
gc()

# do the actual mosaicking
# actually we have to do four mosaicks...
# note that to use raster::mosaic all rasters _must_ have same extent and origin, so we use an insanely high value for tolerance poarameter
# see http://gis.stackexchange.com/questions/226351/combine-multiple-partially-overlapping-rasters-into-a-single-raster-in-r
beginCluster(numCores)
for(i in c('NDVI', 'NDTI')) {
  for(t in c('past', 'present')) {
    cat("mosaicking", i, t, '\n')
    # Not Run:
    # i <- 'NDVI'; t <- 'past'
    mosaicTilesNames <- list.files(path=localWksp, pattern=paste(i, "_", t, "_.*_3035\\.tif$", sep=''), full.names=FALSE)
    mosaicTiles <- sapply(mosaicTilesNames, raster) # NOTE: sapply instead of lapply!
    names(mosaicTiles) <- NULL # clear names to please do.call
    # add to list named elements with actual parameters to mosaic function
    mosaicTiles$fun <- mean
    mosaicTiles$na.rm <- TRUE
    mosaicTiles$tolerance <- 10
    #run do call to implement mosaic over the list of raster objects.
    mos <- do.call(raster::mosaic, mosaicTiles)
    writeRaster(mos, filename=paste(paste(i, t, '3035', sep='_'), 'tif', sep='.'), format='GTiff', overwrite=TRUE)
  }
}
endCluster()
cat("\nDone.\n")  
gc()

# mosaics can have flaws tue to nodata or simply tiling.
# try out now some smoothing/filtering


stop()







# for this exercise GTiffs generated on eira have been copied locally on apodemus
mosaicDir <- '/data/w_Satellite_Forcello/landsat mosaic/NDVI'
# list tiles
#mosaicTilesNames <- list.files(path=mosaicDir, pattern=paste(i, "_", t, "_.*_3035\\.grd$", sep=''), full.names=FALSE)
mosaicTilesNames <- list.files(path=mosaicDir, pattern=paste('NDVI', "_", 'past', "_.*_3035\\.tif$", sep=''), full.names=FALSE)
# load tiles
setwd(mosaicDir)
mosaicTiles <- sapply(mosaicTilesNames, raster) # NOTE: sapply instead of lapply!
# note that for mosaic to behave, we need that:
# - we have NA pixels properly encoded (looks like, no fully values to code for NA, a NA is a NA)
# - ensure resolution is the same (again, all OK, they're all lansdat, so 30 30 for all)
# - ensure projection is the same (that was the tricky part, all GTiffs have been projected from whatever into 3035)

tiles <- mosaicTiles
















## Not Run:
## brute force format conversion from native raster to GTiff
#idxTiles <- list.files(path=localWksp, pattern="ND[TV]I.*_3035\\.grd$", full.names=FALSE)
#convertTiff <- function(aRaster) {
#  r <- raster(aRaster)
#  n <- substr(aRaster,1,nchar(aRaster)-4)
#  writeRaster(r, filename=paste(aRaster, 'tif', sep='.'), format='GTiff') # better than native raster, can be read by other software)}
#}
#lapply(idxTiles, function(x) convertTiff(x)) 







#@TODO mosaicking in R is actually a pain. Try usong gdalwarp directly. Problem is that raster .gri/.grd format is readable from gdal v 2.1 onwards.








setwd(oldwd)
rasterOptions(tmpdir=oldtmpdir)
unlink(tmpDir, recursive=TRUE, force=TRUE)
## end of file


