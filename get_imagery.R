###############################################################################
# UAGRA R Scripts                                                get_imagery.R
###############################################################################
# Landsat squirrel project, get Landsat imagery
###############################################################################
#
# Version 1.0
#
# Description: get Landsat imagery
#
# Usage: <example syntax>
#
# Requires: <list any other package that is require()d by this script>
#
# Returns: <if the script returns some kind of object, explain it here>
#
# References: <if the procedure is based on any kind of literature, place here
#             the relevant bibliographic references>
#
###############################################################################
# created prea 20191202
# updated <name of the responsible of the last update> <timestamp>
# revision history:
#   <author name> <timestamp> - <concise description of the changes made, one per line>
#
###############################################################################
# 
#    Copyright (C) 2019 D.G. Preatoni
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    A copy of the GNU General Public License is available at 
#    http://www.gnu.org/licenses/gpl.txt.
#
###############################################################################

library(getSpatialData)
library(rgdal)
library(car)
library(lubridate)
library(parallel)
library(doParallel) # also loads foreach
library(foreach)

username <- 'prea'

# product names can be has with getLandsat_names()
products <- c("LANDSAT_MSS_C1", "LANDSAT_TM_C1", "LANDSAT_ETM_C1", "LANDSAT_8_C1")

maxCloud <- 10 # percent, land cloud coverage

#### study areas characteristics
# - US: BASE/AOI-MtGraham.kml, path 35 row 37 
#       BASE/AOI-WhiteMountains.kml, path 35 row 36
# - IT: BASE/AOI-Lombardia.kml, path 193 row 28
pathrow <- list(MTG=c(35,37), WHM=c(35,36), LOM=c(193,28))
#### assemble years & seasons time frame
timerange <- list(MTG=c("1992-12-01", "2019-12-01"), WHM=c("1992-12-01", "2019-12-01"), LOM=c("1999-12-01", "2019-12-02"))
# define seasons as months (use numbers, not month names)
# winter: Dec-Mar; spring: Apr-Jun; fall: Sep-Nov
seasonrange <- list(winter=c(1,2,3,12), spring=c(4,5,6), fall=c(9,10,11))
# make a time frame (shortcut, use either MTG or WHM areas since is the widest time span)
timeframe <- expand.grid(year=seq(min(lubridate::year(timerange$MTG)), max(lubridate::year(timerange$MTG))), month=unlist(seasonrange))
timeframe$season <- car::recode(timeframe$month, "c(1,2,3,12)='winter'; c(4,5,6)='spring'; c(9,10,11)='fall'")
# recode december as 'winter of followig year'
timeframe$seasonyear <- timeframe$year
timeframe[timeframe$month==12,]$seasonyear <- timeframe[timeframe$month==12,]$year+1


#### read in AOI shapefiles
AOI.MTG <- readOGR(dsn='BASE', layer='AOI-MtGraham_4326')
AOI.WHM <- readOGR(dsn='BASE', layer='AOI-WhiteMountains_4326')
AOI.LOM <- readOGR(dsn='BASE', layer='AOI-Lombardia_4326')


#### start up things
login_USGS(username=username)
set_archive('ARCHIVIO/Landsat/')


#### make queries: note that getLandsat_query calls take ages to complete, thus results will be stored, to save time
if(file.exists('data/MTG.query.rds')==FALSE){
  set_aoi(AOI.MTG)
  MTG.query <- getLandsat_query(time_range=timerange$MTG, name="all", maxCloudLand=maxCloud)
  saveRDS(MTG.query, file='data/MTG.query.rds')
} else {
  MTG.query <- readRDS('data/MTG.query.rds')
}

if(file.exists('data/WHM.query.rds')==FALSE){
  set_aoi(AOI.WHM)
  WHM.query <- getLandsat_query(time_range=timerange$WHM, name="all", maxCloudLand=maxCloud)
  saveRDS(WHM.query, file='data/WHM.query.rds')
} else {
  WHM.query <- readRDS('data/WHM.query.rds')
}

if(file.exists('data/LOM.query.rds')==FALSE){
  set_aoi(AOI.LOM)
  LOM.query <- getLandsat_query(time_range=timerange$LOM, name="all", maxCloudLand=maxCloud)
  saveRDS(LOM.query, file='data/LOM.query.rds')
} else {
  LOM.query <- readRDS('data/WHM.query.rds')
}


# filter out by path/row and date (land cloud cover already done at query level)
MTG.query <- MTG.query[MTG.query$WRSPath==pathrow$MTG[1] & MTG.query$WRSRow==pathrow$MTG[2] & lubridate::month(MTG.query$acquisitionDate) %in% unlist(seasonrange), ]
WHM.query <- WHM.query[WHM.query$WRSPath==pathrow$WHM[1] & WHM.query$WRSRow==pathrow$WHM[2] & lubridate::month(WHM.query$acquisitionDate) %in% unlist(seasonrange), ]
LOM.query <- LOM.query[LOM.query$WRSPath==pathrow$LOM[1] & LOM.query$WRSRow==pathrow$LOM[2] & lubridate::month(LOM.query$acquisitionDate) %in% unlist(seasonrange), ]

# tag with proper season and AOI ID
tagSeason <- function(dateString, seasonsTable, collapsed=FALSE) {
  yr <- lubridate::year(dateString)
  mo <- lubridate::month(dateString)
  sdf <- data.frame(year=yr, month=mo)
  sdf <- merge(sdf, seasonsTable, by=c('year', 'month'), all.x=TRUE)
  if(collapsed==FALSE) {
    return(sdf)
  } else {
    return(paste(sdf$seasonyear, sdf$season))
  }
} 
MTG.query$season <- tagSeason(MTG.query$acquisitionDate, timeframe, collapse=TRUE)
WHM.query$season <- tagSeason(WHM.query$acquisitionDate, timeframe, collapse=TRUE)
LOM.query$season <- tagSeason(LOM.query$acquisitionDate, timeframe, collapse=TRUE)
MTG.query$AOI <- "MTG"
WHM.query$AOI <- "WHM"
LOM.query$AOI <- "LOM"

# merge all
Landsat.query <- rbind(MTG.query, WHM.query, LOM.query)

#### ready to download (parrrrrallllllel!)
downloadDir <- './ARCHIVIO/incoming/'
nCores <- detectCores() - 1


files <- foreach(i=1:nrow(Landsat.query[]), .combine=c, .packages='getSpatialData') %dopar% {
                 getLandsat_data(Landsat.query[i,], dir_out=downloadDir)     
         }
stopCluster(cl)
