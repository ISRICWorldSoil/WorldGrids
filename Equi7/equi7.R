## Processing Equi7 system for the purpose of global soil mapping
## Equi7 (see Fig. 6 in http://www.sciencedirect.com/science/article/pii/S0098300414001629)

library(rgdal)
library(sp)
library(maptools)
library(plotKML)
library(GSIF)
library(snowfall)
library(rgeos)
setwd("E:\\EQUI7")
if(.Platform$OS.type == "windows"){
  gdal.dir <- shortPathName("C:/Program files/GDAL")
  gdal_translate <- paste0(gdal.dir, "/gdal_translate.exe")
  gdalwarp <- paste0(gdal.dir, "/gdalwarp.exe") 
} else {
  gdal_translate = "gdal_translate"
  gdalwarp = "gdalwarp"
}

lst <- list.files(pattern="*_PROJ_ZONE.shp$", full.names = TRUE, recursive = TRUE)
equi7 <- list(NULL)
for(i in 1:length(lst)){
  equi7[[i]] <- readOGR(lst[i], layer=strsplit(basename(lst[i]),"\\.")[[1]][1])
}

for(i in 1:length(lst)){
  names(equi7)[i] <- strsplit(strsplit(basename(lst[i]), "\\.")[[1]][1], "_")[[1]][3]
}
save(equi7, file="equi7.rda", compress="xz") 

lst.ll <- list.files(pattern="*_GEOG_ZONE.shp$", full.names = TRUE, recursive = TRUE)
equi7.ll <- list(NULL)
for(i in 1:length(lst.ll)){
  equi7.ll[[i]] <- readOGR(lst.ll[i], layer=strsplit(basename(lst.ll[i]),"\\.")[[1]][1])
}

library(plotKML)
kml_open("equi7.kml")
for(i in 1:length(lst)){
  kml_layer(equi7.ll[[i]], subfolder.name=strsplit(basename(lst.ll[i]),"\\.")[[1]][1], colour=ZONE, colour_scale=rep("#FFFF00", 2), alpha=.4)
}
kml_close("equi7.kml")

## TILING SYSTEMS:
c.lst <- c("AF", "AN", "AS", "EU", "NA", "OC", "SA")
t.lst <- list.files(path="EQUI7_V13_GRIDS", pattern="*_PROJ_TILE_T3.shp$", full.names = TRUE, recursive = TRUE)
equi7t3 <- list(NULL)
for(i in 1:length(t.lst)){
  equi7t3[[i]] <- readOGR(t.lst[i], layer=strsplit(basename(t.lst[i]),"\\.")[[1]][1])
  ## subset to tiles with land:
  equi7t3[[i]] <- equi7t3[[i]][equi7t3[[i]]$COVERSLAND==1,]
}
## Few tiles needed manual fixing i.e. they show no-land but should be land
save(equi7t3, file="equi7t3.rda")

world <- raster("ne_10m_land.tif")
world <- as(as(world, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
t1.lst <- list.files(path="EQUI7_V13_GRIDS", pattern="*_PROJ_TILE_T1.shp$", full.names = TRUE, recursive = TRUE)
list_equi7 <- function(i){
  equi7t1 <- readOGR(t1.lst[i], layer=strsplit(basename(t1.lst[i]),"\\.")[[1]][1])
  ## subset to tiles with land:
  geog <- readOGR(gsub("PROJ", "GEOG", t1.lst[i]), layer=strsplit(basename(gsub("PROJ", "GEOG", t1.lst[i])),"\\.")[[1]][1])
  world@proj4string = geog@proj4string
  sel <- over(y=world, x=geog)
  equi7t1 <- equi7t1[which(!is.na(sel[,1])),]
  return(equi7t1)
}
#plot(equi7t1)

## TAKES 30 mins:
sfInit(parallel=TRUE, cpus=7)
sfLibrary(rgdal)
sfLibrary(sp)
sfLibrary(raster)
sfLibrary(rgeos)
sfExport("world", "list_equi7", "t1.lst")
equi7t1 <- sfLapply(1:length(t1.lst), list_equi7)
sfStop()
names(equi7t1) <- c.lst
save(equi7t1, file="equi7t1.rda")

## land polys:
land.lst <- list.files(path="EQUI7_V13_GRIDS", pattern="*_PROJ_LAND.shp$", full.names = TRUE, recursive = TRUE)
equi7land <- list(NULL)
for(i in 1:length(land.lst)){
  equi7land[[i]] <- readOGR(land.lst[i], layer=strsplit(basename(land.lst[i]),"\\.")[[1]][1])
}
plot(equi7land[[7]])
lines(as(equi7t3[[7]], "SpatialLines"))
plot(equi7land[[4]])
lines(as(equi7t3[[4]], "SpatialLines"))
names(equi7t3) <- c.lst <- c("AF", "AN", "AS", "EU", "NA", "OC", "SA")
names(equi7land) <- c.lst <- c("AF", "AN", "AS", "EU", "NA", "OC", "SA")
save(equi7land, file="equi7land.rda")

## write to Shapefiles:
for(i in 1:length(equi7t3)){
  writeOGR(equi7t3[[i]], paste0(names(equi7t3)[i], "_t3_tiles.shp"), paste0(names(equi7t3)[i], "_t3_tiles"), "ESRI Shapefile")
}

for(i in 1:length(equi7t1)){
  writeOGR(equi7t1[[i]], paste0(names(equi7t1)[i], "_t1_tiles.shp"), paste0(names(equi7t1)[i], "_t1_tiles"), "ESRI Shapefile")
}

## 8 representative areas for testing:
s8 <- list(Zone=c("AF", "OC", "AS", "AS", "EU", "NA", "SA", "SA"), TILE=c("072_048", "087_063", "072_087", "048_003", "051_012", "060_036", "072_066", "090_048"))
s8.equi7t3 <- list(NULL)
for(j in 1:length(s8$Zone)){
  s8.equi7t3[[j]] <- equi7t3[[s8$Zone[j]]][equi7t3[[s8$Zone[j]]]$TILE==s8$TILE[j],]
}
names(s8.equi7t3) <- paste0(s8$Zone, "_", s8$TILE)
save(s8.equi7t3, file="s8.equi7t3.rda")

s8.equi7t3.ll <- lapply(s8.equi7t3, spTransform, CRS("+proj=longlat +datum=WGS84"))
s8.equi7t3.ll <- do.call(rbind, s8.equi7t3.ll)
plotKML(s8.equi7t3.ll, filename="s8.equi7t3.kml", folder.name="S8")

## plot in R:
data(landmask)
gridded(landmask) <- ~x+y
proj4string(landmask) <- "+proj=longlat +datum=WGS84"
library(maps)
country.m = map('world', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(country.m$names, ":"), function(x) x[1])
library(maptools)
country <- as(map2SpatialPolygons(country.m, IDs=IDs), "SpatialLines")
mask <- landmask[landmask$mask>10,"mask"]
plot(raster(mask), col="grey")
lines(as(s8.equi7t3.ll, "SpatialLines"))
spplot(s8.equi7t3.ll["COVERSLAND"], col.regions="red", sp.layout=list("sp.lines", country), scales=list(draw=TRUE), xlim=c(-180,180), ylim=c(-90,90), colorkey=FALSE)
