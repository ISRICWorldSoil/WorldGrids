## Overview of Precipitation data sets: https://climatedataguide.ucar.edu/climate-data/precipitation-data-sets-overview-comparison-table
## Merged WorldClim estimated precipitation and GPCP Version 2.2 Combined Precipitation Dataset (http://www.esrl.noaa.gov/psd/data/gridded/data.gpcp.html)
## compare also with: http://trmm.gsfc.nasa.gov/trmm_rain/Events/trmm_climatology_3B43.html
## Bookhagen B. prepared this one: http://www.geog.ucsb.edu/~bodo/TRMM/

library(rgdal)
library(RSAGA)
gdal.dir <- shortPathName("C:/OSGeo4W64/bin")
gdal_translate <- paste0(gdal.dir, "/gdal_translate.exe")
gdalwarp <- paste0(gdal.dir, "/gdalwarp.exe")
gdaladdo <- paste0(gdal.dir, "/gdaladdo.exe")
gdalinfo <- paste0(gdal.dir, "/gdalinfo.exe")
t_srs <- "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs"
te <- c(-20037508, -8283275, 20037508, 18428920)

download.file("ftp://ftp.cdc.noaa.gov/Datasets/gpcc/full_v6/precip.mon.1981-2010.ltm.v6.nc", "precip.mon.1981-2010.ltm.v6.nc")
system(paste(gdalinfo, "precip.mon.1981-2010.ltm.v6.nc"))

lst.gpcp <- paste0("prec_", 1:12, "_10km.sgrd")
lst.wc <- paste0("E:/WORLDGRIDS/worldclim/", "prec_", 1:12, ".sgrd")
mon <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
cod <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
lst.trmm <- paste0("trmm2b31_", tolower(mon), "_mm_per_month.sgrd")
lst.GSMap <- list.files(path="E:\\WORLDGRIDS\\GSMaP", pattern=glob2rx("PRE200*.sgrd"), full.names =TRUE)

## Produce monthly estimates from GSMap 10 km images:
#library(raster)
#grids <- readGDAL("trmm2b31_dec_mm_per_month.sdat")

for(i in 1:length(cod)){
  if(!file.exists(paste0("PREm_", i, ".sgrd"))){
    sel <- lst.GSMap[grep(paste0(cod[i], ".sgrd"), x=lst.GSMap)]
    outm <- paste0("PRE_GSMaP_", cod[i])
    rsaga.geoprocessor(lib="statistics_grid", module=4, param=list(GRIDS=paste(sel, collapse=";", sep=""), MEAN="tmp.sgrd"), check.module.exists=FALSE)
    rsaga.geoprocessor(lib="pj_proj4", module=4, param=list(SOURCE="tmp.sgrd", CRS_PROJ4="+proj=longlat +datum=WGS84", TARGET_DEFINITION=0, TARGET_USER_XMIN=-180, TARGET_USER_XMAX=180, TARGET_USER_YMIN=-90, TARGET_USER_YMAX=90, TARGET_USER_SIZE=0.1, TARGET_GRID="tmp2.sgrd", INTERPOLATION=0), check.module.exists=FALSE)
    ## Scale to mm / month
    rsaga.geoprocessor(lib="grid_calculus", module=1, param=list(GRIDS="tmp2.sgrd", RESULT="tmp3.sgrd", USE_NODATA=0, FORMULA="g1*24*30", TYPE=3), check.module.exists=FALSE)
    rsaga.geoprocessor(lib="statistics_grid", module=4, param=list(GRIDS=paste0(lst.trmm[i], ";", "tmp3.sgrd;", lst.gpcp[i]), MEAN=paste0("PREm_", i, ".sgrd")), check.module.exists=FALSE)
  }
}

## Merge two maps - WorldClim and satellite based PREC:
for(j in 1:length(lst.gpcp)){
  #outtif <- paste0("PREC_", mon[j], "_1km_merc.tif")
  outtif <- paste0("P", cod[j], "MRG3a.tif")
  if(!file.exists(outtif)){
    if(!file.exists(paste0("PREm_", j, "_1km_sum.sgrd"))){
      ## resample to 1km using bicubic spline:
      rsaga.geoprocessor(lib="grid_tools", module=0, param=list(INPUT=paste0("PREm_", j, ".sgrd"), SCALE_DOWN_METHOD=3, TARGET_DEFINITION=0, TARGET_USER_XMIN=-179.9958333333, TARGET_USER_XMAX=179.9958318934, TARGET_USER_YMIN=-59.9958333333, TARGET_USER_YMAX=89.9958327334, TARGET_USER_SIZE=0.0083333333, TARGET_OUT_GRID=paste0("PREm_", j, "_1km.sgrd")), check.module.exists = FALSE)
      ## calculate average:
      rsaga.geoprocessor(lib="grid_calculus", module=1, param=list(GRIDS=paste0(lst.wc[j], ";", "PREm_", j, "_1km.sgrd"), RESULT=paste0("PREm_", j, "_1km_sum.sgrd"), USE_NODATA=0, FORMULA="(g1+g2)/2", TYPE=3), check.module.exists=FALSE)
    }
    ## resample to Mercator proj
    #system(paste0(gdalwarp, ' ', paste0("PREm_", j, "_1km_sum.sdat"), ' ', outtif, ' -r \"bilinear\" -te ', paste(te, collapse=" "), ' -s_srs \"+proj=longlat +datum=WGS84\" -t_srs \"', t_srs, '\" -tr 1000 1000'))
    system(paste0(gdalwarp, ' ', paste0("PREm_", j, "_1km_sum.sdat"), ' ', outtif, ' -r \"near\" -te -180 -90 180 90 -s_srs \"+proj=longlat +datum=WGS84\" -t_srs \"+proj=longlat +datum=WGS84\" -tr 0.008333333 0.008333333'))
    system(paste0(gdaladdo, ' ', outtif, ' 2 4 8 16 32 64'))
    #unlink(gsub("_25km.sgrd", "_1km_sum.****", lst.gpcp[j]))
    #unlink(gsub("_25km.sgrd", "_1km.****", lst.gpcp[j]))
  }
}



