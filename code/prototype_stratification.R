rm(list=ls())
gc()

library(fields)
library(lubridate)
library(NISTunits)
library(raster)
library(ncdf4)
library(rgdal)


### bathymetry
setwd("~/Desktop/professional/biblioteca/data")
bathy <- nc_open('etopo1.nc')
topo <- ncvar_get(bathy, 'Band1')
topo_lat <- ncvar_get(bathy, 'lat')
topo_lon <- ncvar_get(bathy, 'lon')
nc_close(bathy)


################## geographic scope
lonbox_e <- -79 ### Florida Bay
lonbox_w <- -86 ### mouth of Mississippi River
latbox_n <- 38 ### northern coast
latbox_s <- 24 ### remove the Keys

ind_lat <- which(topo_lat<latbox_n & topo_lat>latbox_s)
ind_lon <- which(topo_lon<lonbox_e & topo_lon>lonbox_w)

topo_lat <- topo_lat[ind_lat]
topo_lon <- topo_lon[ind_lon]
topo <- topo[ind_lon,ind_lat]

### load map
# setwd("C:/Users/brendan.turley/Desktop/FL_habs/ne_10m_admin_0_countries")
# setwd("~/Desktop/professional/biblioteca/data/shapefiles/ne_10m_admin_0_countries")
# world <- readOGR('ne_10m_admin_0_countries.shp')
setwd("~/Desktop/professional/biblioteca/data/shapefiles/gshhg-shp-2.3.7/GSHHS_shp/h/")
world <- readOGR('GSHHS_h_L1.shp')
world <- crop(world, extent(-86, -79, 24.5, 28))
# setwd("~/Desktop/professional/biblioteca/data/shapefiles/Florida_Shoreline__1_to_40%2C000_Scale_-shp")
# FL <- readOGR('Florida_Shoreline__1_to_40%2C000_Scale_.shp')

### colorpalettes
### breaks and colors
temp_col <- colorRampPalette(c('gray20','purple','darkorange','gold'))
sal_col <- colorRampPalette(c('midnightblue','dodgerblue4','seagreen3','khaki1'))
chl_col <- colorRampPalette(c('honeydew2','darkseagreen3','forestgreen','darkslategrey'))
ox.col1 <- colorRampPalette(c(1,'firebrick4','red'))
ox.col2 <- colorRampPalette(c('darkgoldenrod4','goldenrod2','gold'))
# ox.col3 <- colorRampPalette(c('dodgerblue4','deepskyblue2','cadetblue1'))
ox.col3 <- colorRampPalette(c(1,'dodgerblue4','deepskyblue2','cadetblue1','azure'))
ex_col <- colorRampPalette(c('gray20','dodgerblue4','indianred3','gold1'))
strat_n_col <- colorRampPalette(c('purple4','purple2','orchid1','gray90'),interpolate='spline',bias=.9)
strat_p_col <- colorRampPalette(rev(c('darkgreen','green3','palegreen2','gray90')),interpolate='spline',bias=.9)



setwd('~/Desktop/professional/projects/Postdoc_FL/data/walton_smith/')
data <- read.csv('WS22141_SampleLog_Initial.csv')
### cruise name for file naming
cruise <- 'WS22141'
### only stations at depth
ind <- which(data$Depth!=0)
data2 <- data[ind,]
st_rm <- c('2','3','6.5','MR','9','9.5','10','12','16','18','21/LK','EK MID','EK OFF','WS','KW1','KW2')
data3 <- data2[!is.element(data2$Station,st_rm),]
data3$Date.GMT <- mdy(data3$Date.GMT)
data3 <- data3[order(data3$Station),]
st_kp <- data3$Station

### only stations at surface
ind <- which(data$Depth==0)
data2 <- data[ind,]
data4 <- data2[is.element(data2$Station,st_kp),]
data4$Date.GMT <- mdy(data4$Date.GMT)
data4 <- data4[order(data4$Station),]

### temperature stratification
temp_ind <- grep('temperature',names(data3),ignore.case = T)
data3$tempF <- NISTdegCtOdegF(data3[,temp_ind])
temp_ind <- grep('temperature',names(data4),ignore.case = T)
data4$tempF <- NISTdegCtOdegF(data4[,temp_ind])

temp_str <- data4$tempF - data3$tempF

bubblePlot(data3$Longitude.Decimal,data3$Latitude.Decimal,temp_str,asp=1)

### salinity stratification
sal_ind_bot <- grep('salinity',names(data3),ignore.case = T)
sal_ind_surf <- grep('salinity',names(data4),ignore.case = T)

sal_str <- data4[,sal_ind_surf] - data3[,sal_ind_bot]

bubblePlot(data3$Longitude.Decimal,data3$Latitude.Decimal,sal_str,asp=1)


### subset bathymetry
adj1 <- mean(diff(topo_lon))*5
adj2 <- mean(diff(topo_lat))*5
ind_lat <- which(topo_lat<(max(data3$Latitude.Decimal,na.rm=T)+adj1) & topo_lat>(min(data3$Latitude.Decimal,na.rm=T)-adj1))
ind_lon <- which(topo_lon<(max(data3$Longitude.Decimal,na.rm=T)+adj2) & topo_lon>(min(data3$Longitude.Decimal,na.rm=T)-adj2))
### bathymetry for spatial covariate model
topo_lat2 <- topo_lat[ind_lat]
topo_lon2 <- topo_lon[ind_lon]
topo2 <- topo[ind_lon,ind_lat]
topo2[which(topo2>0)] <- NA

### locations for krig
ctd.loc <- cbind(lon=data3$Longitude.Decimal,lat=data3$Latitude.Decimal)

loc.grid <- list(lon=topo_lon2,
                 lat=topo_lat2)
Z <- list(x=loc.grid$lon,y=loc.grid$lat,z=topo2)

xlims <- range(loc.grid$lon)
ylims <- range(loc.grid$lat)


my.krig2 <- spatialProcess(ctd.loc,temp_str,Z=-data3$Depth)
temp_kriged2 <- predictSurface(my.krig2, loc.grid, ZGrid=Z,extrap=T)
temp_se2 <- predictSurfaceSE(my.krig2, loc.grid, ZGrid=Z,extrap=T)

imagePlot(temp_kriged2,asp=1)
imagePlot(temp_se2,asp=1)

my.krig2 <- spatialProcess(ctd.loc,sal_str,Z=-data3$Depth)
sal_kriged2 <- predictSurface(my.krig2, loc.grid, ZGrid=Z,extrap=T)
sal_se2 <- predictSurfaceSE(my.krig2, loc.grid, ZGrid=Z,extrap=T)

imagePlot(sal_kriged2,asp=1)
imagePlot(sal_se2,asp=1)


### stratification
# sstrat_breaks <- pretty(sal_str[which(sal_str<=quantile(sal_str,quant,na.rm=T))],n=20)
sstrat_breaks <- pretty(sal_str[which(sal_str<=quantile(sal_str,.99,na.rm=T) & sal_str>=quantile(sal_str,.01,na.rm=T))],n=30)
if(any(sstrat_breaks==0) & sstrat_breaks[1]!=0 & sstrat_breaks[length(sstrat_breaks)]!=0){
  # limit <- mean(abs(range(sal_str,na.rm=T)))
  # sal_str[which(sal_str<(-limit))] <- -limit
  # sal_str[which(sal_str>limit)] <- limit
  # sstrat_breaks <- pretty(sal_str,n=30)
  sal_str[which(sal_str>sstrat_breaks[length(sstrat_breaks)])] <- sstrat_breaks[length(sstrat_breaks)]
  sal_str[which(sal_str<sstrat_breaks[1])] <- sstrat_breaks[1]
  sstrat_cols <- c(strat_n_col(length(which(sstrat_breaks<0))),
                   strat_p_col(length(which(sstrat_breaks>0))))
}else{
  sstrat_cols <- rev(strat_p_col(length(sstrat_breaks)-1))  
  sal_str[which(sal_str>sstrat_breaks[length(sstrat_breaks)])] <- sstrat_breaks[length(sstrat_breaks)]
  sal_str[which(sal_str<sstrat_breaks[1])] <- sstrat_breaks[1]
}
# tstrat_breaks <- pretty(temp_str[which(temp_str<=quantile(temp_str,quant,na.rm=T))],n=20)
tstrat_breaks <- pretty(temp_str[which(temp_str<=quantile(temp_str,.99,na.rm=T) & temp_str>=quantile(temp_str,.01,na.rm=T))],n=30)
if(any(tstrat_breaks==0) & tstrat_breaks[1]!=0 & tstrat_breaks[length(tstrat_breaks)]!=0){
  # limit <- mean(abs(range(temp_str,na.rm=T)))
  # temp_str[which(temp_str<(-limit))] <- -limit
  # temp_str[which(temp_str>limit)] <- limit
  temp_str[which(temp_str>tstrat_breaks[length(tstrat_breaks)])] <- tstrat_breaks[length(tstrat_breaks)]
  temp_str[which(temp_str<tstrat_breaks[1])] <- tstrat_breaks[1]
  tstrat_breaks <- pretty(temp_str,n=30)
  tstrat_cols <- c(strat_n_col(length(which(tstrat_breaks<0))),
                   strat_p_col(length(which(tstrat_breaks>0))))
}else{
  tstrat_cols <- strat_n_col(length(tstrat_breaks)-1)
  temp_str[which(temp_str>tstrat_breaks[length(tstrat_breaks)])] <- tstrat_breaks[length(tstrat_breaks)]
  temp_str[which(temp_str<tstrat_breaks[1])] <- tstrat_breaks[1]
}


imagePlot(temp_kriged2,
          breaks=tstrat_breaks,col=tstrat_cols,
          asp=1)
points(data3$Longitude.Decimal,data3$Latitude.Decimal,pch=16,col='gray50',cex=.5)

imagePlot(sal_kriged2,
          breaks=sstrat_breaks,col=sstrat_cols,
          asp=1)
points(data3$Longitude.Decimal,data3$Latitude.Decimal,pch=16,col='gray50',cex=.5)


### locations for krig
adj <- .1
resolution <- .01
ctd.loc <- cbind(lon=data3$Longitude.Decimal,lat=data3$Latitude.Decimal)
loc.grid <- list(lon=seq(min(data3$Longitude.Decimal,na.rm=T)-adj, max(data3$Longitude.Decimal,na.rm=T)+adj,resolution),
                 lat=seq(min(data3$Latitude.Decimal,na.rm=T)-adj, max(data3$Latitude.Decimal,na.rm=T)+adj,resolution))

xlims <- range(loc.grid$lon)
ylims <- range(loc.grid$lat)

my.krig <- spatialProcess(ctd.loc,temp_str)
temp_kriged <- predictSurface(my.krig, loc.grid, extrap=T)
temp_se <- predictSurfaceSE(my.krig, loc.grid, extrap=T)

imagePlot(temp_kriged,asp=1)
imagePlot(temp_se,asp=1)

my.krig <- spatialProcess(ctd.loc,sal_str)
sal_kriged <- predictSurface(my.krig, loc.grid, extrap=T)
sal_se <- predictSurfaceSE(my.krig, loc.grid, extrap=T)

imagePlot(sal_kriged,asp=1)
imagePlot(sal_se,asp=1)

imagePlot(temp_se2$z-temp_se$z,asp=1)
imagePlot(sal_se2$z-sal_se$z,asp=1)
