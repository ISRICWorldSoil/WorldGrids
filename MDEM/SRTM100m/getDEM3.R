## Download and prepare DEM3 from viewfinderpanoramas.org
## http://www.viewfinderpanoramas.org/Coverage%20map%20viewfinderpanoramas_org3.htm

library(rgdal)
library(raster)
library(sp)
library(snowfall)
library(gdalUtils)
library(R.utils)
library(RSAGA)
gdal.dir = shortPathName("C:\\Program Files\\GDAL")
gdal_setInstallation(search_path=gdal.dir, rescan=TRUE)

lst <- read.table("list_of_tiles.txt")
zip.lst <- sapply(lst$V7, function(x){strsplit(strsplit(paste(x), "href=\"", fixed=TRUE)[[1]][2], "\"")[[1]][1]})
## 1140 tiles!
for(j in 1:length(zip.lst)){
  gz.file <- paste(td, strsplit(zip.lst[j], "/")[[1]][5], sep="\\")
  if(!file.exists(gz.file)){
    download.file(zip.lst[j], gz.file, method="wget", quiet=TRUE, extra="--no-check-certificate")
  }
}

sel <- list.files(pattern=glob2rx("*.zip$"), full.names=TRUE)
## 1127 tiles
## unzip
for(i in 1:length(sel)){
 system(paste("7za x", sel[i]))
}

## create a mosaic:
#tmp.lst <- unlist(sapply(sel, function(x){list.files(path=strsplit(x, ".zip")[[1]][1], pattern=glob2rx("*.hgt"), full.names=TRUE)}))
tmp.lst <- list.files(pattern=glob2rx("*.hgt$"), full.names=TRUE, recursive=TRUE)
## 26,721 tiles!!!
GDALinfo(tmp.lst[100])
GDALinfo(tmp.lst[1000])
GDALinfo(tmp.lst[5000])
GDALinfo(tmp.lst[15000])
unlink("my_liste.txt")
cat(tmp.lst, sep="\n", file="my_liste.txt")
gdalbuildvrt(input_file_list="my_liste.txt", output.vrt="globe.vrt")
## takes 20 mins!!
#gdalinfo("globe.vrt")

## GLOBAL MOSAIC HUGE!!!
#gdalwarp("globe.vrt", dstfile="DEMSRP5a_250m.tif", t_srs="+proj=longlat +datum=WGS84", tr=c(1/480,1/480), r="near", srcnodata=-32768, dstnodata=-32768, te=c(-180,-90,180,90), ot="Int16", overwrite=TRUE)
unlink("DEMSRP3a.tif")
gdalwarp("globe.vrt", dstfile="DEMSRP3a.tif", tr=c(1/120,1/120), r="average", srcnodata=-32768, dstnodata=-32768, te=c(-180,-90,180,90), ot="Int16", overwrite=TRUE)
## TAKES CA 20 mins!

## per continent:
load("G:\\Equi7\\equi7t3.rda")
tile.tif <- function(t, t_srs = proj4string(t)){
  for(j in 1:nrow(t)){
    nfile <- paste0("tiled/DEM_", strsplit(paste(t@data[j,"SHORTNAME"]), " ")[[1]][2], "_", t@data[j,"TILE"], ".tif")
    te <- as.vector(bbox(t[j,]))
    gdalwarp("globe.vrt", dstfile=nfile, t_srs=t_srs, tr=c(250,250), r="average", srcnodata=-32768, dstnodata=-32768, te=te, ot="Int16")
  }
}

sfInit(parallel=TRUE, cpus=6)
sfLibrary(gdalUtils)
sfLibrary(sp)
sfExport("equi7t3", "tile.tif")
x <- sfLapply(equi7t3, tile.tif)
sfStop() ## end of parallel processing
## 2,594 TILES!

## Mosaics (per continent):
for(j in 1:length(equi7t3)){
  if(!file.exists(paste0("DEM_", names(equi7t3)[j], "_250m.sdat"))){
    t.lst <- list.files(path="tiled", pattern=glob2rx(paste0("DEM_", names(equi7t3)[j], "_*_*.tif$")), full.names=TRUE, recursive=TRUE)
    if(!file.exists(paste0(names(equi7t3)[j], ".vrt"))){
      unlink("my_liste.txt")
      cat(t.lst, sep="\n", file="my_liste.txt")
      gdalbuildvrt(input_file_list="my_liste.txt", output.vrt=paste0(names(equi7t3)[j], ".vrt"))
    }
    gdalwarp(paste0(names(equi7t3)[j], ".vrt"), dstfile=paste0("DEM_", names(equi7t3)[j], "_250m.sdat"), r="near", srcnodata=-32768, dstnodata=-32767, ot="Int16", of="SAGA")
  }
}

## single country only:
sel2 <- c(paste0("M", 33:35, ".zip"), paste0("L", 33:35, ".zip"))
lu.csy <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
## create a mosaic:
tmp2.lst <- unlist(sapply(sel2, function(x){list.files(path=strsplit(x, ".zip")[[1]][1], pattern=glob2rx("*.hgt"), full.names=TRUE)}))
unlink("my_liste.txt")
cat(tmp2.lst, sep="\n", file="my_liste.txt")
gdalbuildvrt(input_file_list="my_liste.txt", output.vrt="hu.vrt")
gdalwarp("hu.vrt", "HU_DEMSRE6a_250m.tif", t_srs=lu.csy, tr=c(250,250), r="bilinear", ot="Int16", overwrite=TRUE, te=c(4787000, 2549000, 5280000, 2896000))
