## Derive mean Global DEM avarage between three global DEM sources: (1) http://www.viewfinderpanoramas.org/dem3.html, (2) SRTMGL3 - http://e4ftl01.cr.usgs.gov/SRTM/SRTMGL3.003/2000.02.11/ and (3) GMTED2010 - https://lta.cr.usgs.gov/GMTED2010
## then derive standard DEM parameters

load("../Equi7/equi7t3.rda")

## Prepare the GMTED2010 (https://lta.cr.usgs.gov/GMTED2010):
for(j in 1:length(equi7t3)){
  if(!file.exists(paste0("GMTED2010_", names(equi7t3)[j], "_250m.sdat"))){ ## GMTED
    te <- extent(raster(paste0("DEM_", names(equi7t3)[j], "_250m.sdat")))
    gdalwarp("H:\\srtm15_plus\\GMTED2010_250m.tif", dstfile=paste0("GMTED2010_", names(equi7t3)[j], "_250m.sdat"), tr=c(250,250), r="near", t_srs=proj4string(equi7t3[[j]]), te=as.vector(te)[c(1,3,2,4)], srcnodata=-32768, dstnodata=-32767, ot="Int16", of="SAGA") ##
  }
}

## For Europe we could also use [http://www.eea.europa.eu/data-and-maps/data/eu-dem]
#te <- extent(raster(paste0("DEM_", names(equi7t3)[4], "_250m.sdat")))
#gdalinfo("eudem_dem_4258_europe.tif")
#gdalwarp("eudem_dem_4258_europe.tif", dstfile="eudem_dem_4258_europe.sdat", tr=c(250,250), r="average", t_srs=proj4string(equi7t3[["EU"]]), te=as.vector(te)[c(1,3,2,4)], dstnodata=-32767, ot="Int16", of="SAGA")

## Derive MEAN from 3 DEMs
for(j in 1:length(equi7t3)){
  if(!file.exists(paste0("MDEM_", names(equi7t3)[j], "_250m.sdat"))){
    rsaga.geoprocessor(lib="statistics_grid", module=4, param=list(GRIDS=paste0(paste0("GMTED2010_", names(equi7t3)[j], "_250m.sgrd"), "; ", paste0("DEM_", names(equi7t3)[j], "_250m.sgrd"), "; ", paste0("SRTMGL3_", names(equi7t3)[j], "_250m.sgrd")), MEAN=paste0("MDEM_", names(equi7t3)[j], "_250m.sgrd"), STDDEV=paste0("stdev_", names(equi7t3)[j], "_250m.sdat")), check.module.exists =FALSE, show.output.on.console = FALSE, warn=FALSE)
  }
}
## filter missing values / artifacts (if stdev > 50)?

## Export to GeoTiffs:
for(j in 1:length(equi7t3)){
   if(!file.exists( paste0("MDEM_", names(equi7t3)[j], "_250m.tif"))){
     if(any(names(equi7t3)[j] %in% c("AF","EU"))){
     gdal_translate(paste0("MDEM_", names(equi7t3)[j], "_250m.sdat"), paste0("MDEM_", names(equi7t3)[j], "_250m.tif"), a_nodata=-32767, ot="Int16", s_srs=proj4string(equi7t3[[j]]))
     } else {
       gdal_translate(paste0("MDEM_", names(equi7t3)[j], "_250m_f.sdat"), paste0("MDEM_", names(equi7t3)[j], "_250m.tif"), a_nodata=-32767, ot="Int16", s_srs=proj4string(equi7t3[[j]]))
     }
   }
}
## TH: Some tiles I had to clean up manually (digitize polygons, remove values and fill in the gaps using CLOSE GAPS using splines in SAGA GIS)

## Derive 6 standard DEM parameters:
source("DEM_parameters_parallel.R")
for(j in c(1,3,4,5,6,7)){
  saga_TA(inputFile=paste0("DEM_", names(equi7t3)[j], "_250m.tif.gz"), smaskFile=paste0("SMK_", names(equi7t3)[j], "_250m.tif.gz"))
}
