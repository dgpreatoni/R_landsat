################################################################################
# calculates indices from LANDSAT imagery
# 
# version 4.0
# created prea 20161128
# revised prea 20170314 happy pi day, modified to use RStoolbox functions 
#                       and parallel computing (snow,MPI via raster::beginCluster())
#
# to be run on eira
#
# Assumes that raw LANDSAT .tar.gz are stored in LandsatDir.
# Processes seoarately each scene.
# Processing involves extracting from a scene file (.tar.gz) the bands needed, calculating indices and cropping with a polygon shapefile.
# Processed indices are saved as GeoTIFF in OutputDir
#
# Base workflow:
# (to be completed)
################################################################################

rm(list=ls())

source('/home/gis/R_scripts/raster_utilities.R')
#library(raster)  # already invoked by raster_utilities.R
#library(rgeos)  # already invoked by raster_utilities.R
library(RStoolbox) # satellite stuff...
library(rgdal)
library(readxl)
library(ggplot2)

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

# Indices lookup table, an Excel fie with Index	Sensor	BandA	BandB (in data directory)
indexLookupFile <- 'Landsat indices.xls'
indexLookupSheet <- 'Indices'

#### prepare standard paths
baseDir <- paste(dataRoot, 'BASE', sep='/')
archiveDir <- paste(dataRoot, 'ARCHIVIO', sep='/')
dataDir <- paste(dataRoot, 'data', sep='/')

#### get AOI cookie cutter
AOIClipper <- readOGR(dsn=baseDir, layer=AOIShape)

#### get mask vector
maskPoly <- readOGR(dsn=baseDir, layer=maskShape)
  
#### get a list of all available LANDSAT scenes
imageryFiles <- list.files(path=landsatDir, pattern="\\.tar\\.gz$", full.names=FALSE)

#### read indices lookup table
indexLookup <- read_excel(path=paste(dataDir, indexLookupFile, sep='/'), sheet=indexLookupSheet)

#### carefully set work directory _and_ raster temporary directory: this obviously EATS disk space!
oldwd <- getwd()
oldtmpdir <- rasterOptions()$tmpdir
tmpDir <- tempfile(pattern='raster', tmpdir=localWksp)
setwd(localWksp)
rasterOptions(tmpdir=tmpDir)

## process scenes
calcIndex <- function(x, A, B) (x[[A]]-x[[B]])/(x[[A]]+x[[B]])
# prepare for computation: go parallel
cat("going parallel...\n")
beginCluster(numCores)
for(sceneFile in imageryFiles) {
  cat('Processing scene', sceneFile, '\n')
  # scene must be uncompressed to work with it: create a scratch directory and untar there
  sceneBasename <- gsub("\\.tar\\.gz$", "", sceneFile)
  if(length(list.files(pattern=sceneBasename))==4) { # that scene has been already processed
    cat("\tscene has alreday been processed, skipping.\n")  
  } else { # process it...
    cat("\textracting bands...\n")
    tmpRoot <- tempfile(pattern=sceneBasename, tmpdir=localWksp)
    untar(tarfile=paste(landsatDir, sceneFile, sep='/'), exdir=tmpRoot)
    ## read in metadata
    metaData <- readMeta(paste(tmpRoot, '/', sceneBasename, '_MTL.txt', sep=''))
    # calculate which bands do we need
    indicesTable <- indexLookup[indexLookup$SensorName==metaData$SENSOR,]
    bands <- unique(c(indicesTable$BandA, indicesTable$BandB))
      # read in bands
    scene <- stackMeta(metaData)
    cat("\tradiometric correction (SR) for bands", sort(bands), "...\n")
    # radiometric correction (SR surface radiance) for the needed bands only  
    scene.sr <- radCor(scene, metaData=metaData, method="rad", bandSet=paste(bands, 'dn', sep='_'))
    # calculate indices
    for(i in indicesTable$Index) {
      cat("\tcalculating", i, "\n")
      A <- paste(as.character(indicesTable[indicesTable$Index==i,'BandA']), 'tra', sep='_')
      B <- paste(as.character(indicesTable[indicesTable$Index==i,'BandB']), 'tra', sep='_')
      idx <- clusterR(scene.sr, calcIndex, args=list(A=A, B=B))
      names(idx) <- i
      writeRaster(idx, filename=paste(localWksp, paste(sceneBasename, '_', i, '.img', sep=''), sep='/'), format="HFA", overwrite=TRUE)
    }
    # clean up
    rm(idx, scene.sr, scene, metaData)
    cat("\tcleaning up...\n")
    unlink(tmpRoot, recursive=TRUE, force=TRUE)
  }
}
# stop computation: exit from parallel
cat("closing parallel mode...\n")
endCluster()
gc()
cat("\nDone.\n")  
setwd(oldwd)
rasterOptions(tmpdir=oldtmpdir)
unlink(tmpDir, recursive=TRUE, force=TRUE)
## end of file

