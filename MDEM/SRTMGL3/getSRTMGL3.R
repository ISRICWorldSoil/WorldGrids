## Download and mosaic SRTMGL3
## Requires: wget, GDAL and EQUI7 installation

library(XML)
library(RCurl)
library(gdalUtils)
library(rgdal)
library(utils)
library(snowfall)
library(raster)
library(RSAGA)
gdal.dir <- shortPathName("C:/Program files/GDAL")
gdal_setInstallation(search_path=gdal.dir, rescan=TRUE)

## Download all files in the folder:
system("wget -r --no-parent http://e4ftl01.cr.usgs.gov/SRTM/SRTMGL3.003/2000.02.11/")
## Downloaded: 57,129 files, 16G in 19h 14m 19s (236 KB/s)

## list of files:
sel <- list.files(path="./e4ftl01.cr.usgs.gov/SRTM/SRTMGL3.003/2000.02.11", pattern=glob2rx("*.hgt.zip$"), full.names=TRUE, recursive=TRUE)
## 14,280 tiles
for(i in 1:length(sel)){
 system(paste("7za x", sel[i]))
}

## create a MEGA!! mosaic:
tmp.lst <- list.files(pattern=glob2rx("*.hgt$"))
unlink("my_liste.txt")
cat(tmp.lst, sep="\n", file="my_liste.txt")
gdalbuildvrt(input_file_list="my_liste.txt", output.vrt="SRTMGL3.vrt")

## Mosaics (per continent):
load("../../Equi7/equi7t3.rda")
for(j in 1:length(equi7t3)){
  dstfile <- paste0("../SRTMGL3_", names(equi7t3)[j], "_250m.sdat")
  if(!file.exists(dstfile)){
    te <- extent(raster(paste0("../DEM_", names(equi7t3)[j], "_250m.sdat")))
    gdalwarp("SRTMGL3.vrt", dstfile=dstfile, tr=c(250,250), r="average", t_srs=proj4string(equi7t3[[j]]), te=as.vector(te)[c(1,3,2,4)], dstnodata=-32767, ot="Int16", of="SAGA")
  }
}