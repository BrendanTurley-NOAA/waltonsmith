rm(list=ls())
gc()

library(fields)
library(lubridate)
library(NISTunits)
library(raster)
library(ncdf4)
library(rgdal)

source('~/Desktop/professional/biblioteca/scripts/color_breaks.R')

data_plot <- function(longitude, latitude, data, color_fxn, n_breaks=15, title='', cex=2, xlab='',ylab=''){
  breaks <- pretty(data,n=n_breaks)
  cols <- color_fxn(length(breaks)-1)
  
  plot(longitude,latitude,pch=21,asp=1,las=1,cex=cex,
       bg=cols[as.numeric(cut(data,breaks))],
       col=cols[as.numeric(cut(data,breaks))],
       xlab=xlab,ylab=ylab)
  imagePlot(zlim=range(breaks,na.rm=T),breaks=breaks,col=cols,legend.only=TRUE,legend.width = 1.5)
  plot(world,add=T,col='gray80')
  contour(topo_lon,
          topo_lat,
          topo,
          add=T,levels=c(-200,-100,-50,-25,-10),col='gray70')
  mtext(title)
}

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
world <- crop(world, extent(-86, -79, 24.5, 30))
# setwd("~/Desktop/professional/biblioteca/data/shapefiles/Florida_Shoreline__1_to_40%2C000_Scale_-shp")
# FL <- readOGR('Florida_Shoreline__1_to_40%2C000_Scale_.shp')
### load shollow masks
setwd('~/Desktop/professional/projects/Postdoc_FL/data/')
masks <- read.csv('TB_SB_CH_masks.csv')

### colorpalettes
### breaks and colors
temp_col <- colorRampPalette(c('gray20','purple3','darkorange','gold'))
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
setwd('~/Downloads')
list.files()
data <- read.csv('WS23011_SampleLog_initial.csv')
### cruise name for file naming
cruise <- 'WS23011'
### check lat/lons
lons_cal <- -(data$Longitude.Deg+data$Longitude.Min/60)
lats_cal <- (data$Latitude.Deg+data$Latitude.Min/60)

plot(data$Longitude.Decimal,data$Latitude.Decimal,asp=1)


### surface
ind <- which(data$Depth==0)
data1 <- data[ind,]
# st_rm <- c('2','3','6.5','MR','9','9.5','10','12','16','18','21/LK','EK MID','EK OFF','WS','KW1','KW2')
st_rm <- c('2','3','6.5','MR','9','9.5','10','12','16','18','21/LK','EK MID','EK OFF','WS')
data1 <- data1[!is.element(data1$Station,st_rm),]
data1$Date.GMT <- mdy(data1$Date.GMT)

### check for repeated stations
ind1 <- which(duplicated(data1$Station))
# for(i in data1$Station[ind1]){
#   tmp <- data1[which(data1$Station==i),]
#   cat(paste('Station',i),'\n')
#   cat(which(data1$Station==i),'\n')
#   cat(tmp$Depth,'\n')
#   cat(paste(tmp$Date.GMT,tmp$Time..GMT.),'\n')
#   cat('\n\n')
# }


### only stations at depth
ind <- which(data$Depth!=0)
data2 <- data[ind,]
# st_rm <- c('2','3','6.5','MR','9','9.5','10','12','16','18','21/LK','EK MID','EK OFF','WS','KW1','KW2')
st_rm <- c('2','3','6.5','MR','9','9.5','10','12','16','18','21/LK','EK MID','EK OFF','WS')
data3 <- data2[!is.element(data2$Station,st_rm),]
data3$Date.GMT <- mdy(data3$Date.GMT)

### check for repeated stations
ind1 <- which(duplicated(data3$Station))
ind2 <- which(duplicated(data3$Station,fromLast=T))
as.vector(t(cbind(ind1,ind2)))
data3[as.vector(t(cbind(ind1,ind2))),]
remove <- rep(NA,length(ind1))
for(i in 1:length(ind1)){
  remove[i] <- ifelse(data3$Depth[ind1[i]]>data3$Depth[ind2[i]],ind2[i],ind1[i])
}
### remove shallower samples
if(!is.na(remove)){
  data3 <- data3[-remove,]  
}


data1 <- data1[is.element(data1$Station,data3$Station),]

dt_dz <- data3$Temperature.CTD.data - data1$Temperature.CTD.data
ds_dz <- data3$Salinity.CTD.data - data1$Salinity.CTD.data

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


### ----------------- Temperature krig -----------------
temp_ind <- grep('temperature',names(data3),ignore.case = T)
# data3$tempF <- NISTdegCtOdegF(data3[,temp_ind])

### no covariate
my.krig <- spatialProcess(ctd.loc,data3[,temp_ind],REML=T)
rmse1 <- sqrt(sum(resid.Krig(my.krig)^2)/length(resid.Krig(my.krig)))
### depth as covariate
my.krig2 <- spatialProcess(ctd.loc,data3[,temp_ind],Z=-data3$Depth,REML=T)
rmse2 <- sqrt(sum(resid.Krig(my.krig2)^2)/length(resid.Krig(my.krig2)))

### which is smallest?
cbind(rmse1,rmse2)

# temp_kriged2 <- predictSurface(my.krig, loc.grid, extrap=T)
# temp_se2 <- predictSurfaceSE(my.krig, loc.grid, extrap=T)
### depth as covariate
temp_kriged2 <- predictSurface(my.krig2, loc.grid, ZGrid=Z,extrap=T)
temp_se2 <- predictSurfaceSE(my.krig2, loc.grid, ZGrid=Z,extrap=T)

if(max(temp_kriged2$z,na.rm=T)>max(data3[,temp_ind],na.rm=T)){
  temp_kriged2$z[which(temp_kriged2$z>max(data3[,temp_ind],na.rm=T))] <- max(data3[,temp_ind],na.rm=T)
}
if(min(temp_kriged2$z,na.rm=T)<min(data3[,temp_ind],na.rm=T)){
  temp_kriged2$z[which(temp_kriged2$z<min(data3[,temp_ind],na.rm=T))] <- min(data3[,temp_ind],na.rm=T)
}

temp_breaks <- pretty(data3[,temp_ind],n=20)
temp_cols <- temp_col(length(temp_breaks)-1)

### convert to temp F
tempF <- NISTdegCtOdegF(temp_kriged2$z)
tempF_breaks <- pretty(NISTdegCtOdegF(data3[,temp_ind]),n=20)
tempF_cols <- temp_col(length(tempF_breaks)-1)

### ----------------- Salinity krig -----------------
sal_ind <- grep('salinity',names(data3),ignore.case = T)

### no covariate
my.krig <- spatialProcess(ctd.loc,data3[,sal_ind],REML=T)
rmse1 <- sqrt(sum(resid.Krig(my.krig)^2)/length(resid.Krig(my.krig)))
### depth as covariate
my.krig2 <- spatialProcess(ctd.loc,data3[,sal_ind],Z=-data3$Depth,REML=T)
rmse2 <- sqrt(sum(resid.Krig(my.krig2)^2)/length(resid.Krig(my.krig2)))

### which is smallest?
cbind(rmse1,rmse2)

sal_kriged2 <- predictSurface(my.krig, loc.grid, extrap=T)
sal_se2 <- predictSurfaceSE(my.krig, loc.grid, extrap=T)
### depth as covariate
# sal_kriged2 <- predictSurface(my.krig2, loc.grid, ZGrid=Z,extrap=T)
# sal_se2 <- predictSurfaceSE(my.krig2, loc.grid, ZGrid=Z,extrap=T)

if(max(sal_kriged2$z,na.rm=T)>max(data3[,sal_ind],na.rm=T)){
  sal_kriged2$z[which(sal_kriged2$z>max(data3[,sal_ind],na.rm=T))] <- max(data3[,sal_ind],na.rm=T)
}
if(min(sal_kriged2$z,na.rm=T)<min(data3[,sal_ind],na.rm=T)){
  sal_kriged2$z[which(sal_kriged2$z<min(data3[,sal_ind],na.rm=T))] <- min(data3[,sal_ind],na.rm=T)
}

sal_breaks <- pretty(data3[,sal_ind],n=20)
sal_cols <- sal_col(length(sal_breaks)-1)


### ----------------- DO krig -----------------
oxy_ind <- grep('oxygen',names(data3),ignore.case = T)
if(length(oxy_ind>1)){
  oxy_ind <- oxy_ind[1]
}
if(class(data3[,oxy_ind])!='numeric'){
  data3[,oxy_ind] <- as.numeric(data3[,oxy_ind])
}
### no covariate
my.krig <- spatialProcess(ctd.loc,data3[,oxy_ind],REML=T)
rmse1 <- sqrt(sum(resid.Krig(my.krig)^2)/length(resid.Krig(my.krig)))
### depth as covariate
my.krig2 <- spatialProcess(ctd.loc,data3[,oxy_ind],Z=-data3$Depth,REML=T)
rmse2 <- sqrt(sum(resid.Krig(my.krig2)^2)/length(resid.Krig(my.krig2)))
### temperatureas covariate
my.krig3 <- spatialProcess(ctd.loc,data3[,oxy_ind],Z=data3[,temp_ind],REML=T)
rmse3 <- sqrt(sum(resid.Krig(my.krig3)^2)/length(resid.Krig(my.krig3)))

### which is smallest?
cbind(rmse1,rmse2,rmse3)

# do_kriged2 <- predictSurface(my.krig, loc.grid, extrap=T)
# do_se2 <- predictSurfaceSE(my.krig, loc.grid, extrap=T)
### depth as covariate
# do_kriged2 <- predictSurface(my.krig2, loc.grid, ZGrid=Z,extrap=T)
# do_se2 <- predictSurfaceSE(my.krig2, loc.grid, ZGrid=Z,extrap=T)
### temperatureas covariate
do_kriged2 <- predictSurface(my.krig3, loc.grid, ZGrid=temp_kriged2,extrap=T)
do_se2 <- predictSurfaceSE(my.krig3, loc.grid, ZGrid=Z,extrap=T)

if(max(do_kriged2$z,na.rm=T)>max(data3[,oxy_ind])){
  do_kriged2$z[which(do_kriged2$z>max(data3[,oxy_ind]))] <- max(data3[,oxy_ind])
}
if(min(do_kriged2$z,na.rm=T)<min(data3[,oxy_ind])){
  do_kriged2$z[which(do_kriged2$z<min(data3[,oxy_ind]))] <- min(data3[,oxy_ind])
}

breaks <- pretty(data3[,oxy_ind],n=20)
cols <- c(ox.col1(length(breaks[breaks<2])),
          ox.col2(length(breaks[breaks>=2 & breaks<3.5])),
          ox.col3(length(breaks[breaks>=3.5])-1))


### ----------------- dt_dz krig -----------------
my.krig2 <- spatialProcess(ctd.loc,dt_dz)
dtdz_kriged2 <- predictSurface(my.krig2, loc.grid, extrap=T)
dtdz_se2 <- predictSurfaceSE(my.krig2, loc.grid, extrap=T)

# dtdz_brks <- pretty(dtdz_kriged2$z,n=20)
dtdz_brks <- pretty(dt_dz,n=20)
dtdz_cols <- c(strat_n_col(length(which(dtdz_brks<0))),
               strat_p_col(length(which(dtdz_brks>0))))

### ----------------- ds_dz krig -----------------
my.krig2 <- spatialProcess(ctd.loc,ds_dz)
dsdz_kriged2 <- predictSurface(my.krig2, loc.grid, extrap=T)
dsdz_se2 <- predictSurfaceSE(my.krig2, loc.grid, extrap=T)

# dsdz_brks <- pretty(dsdz_kriged2$z,n=20)
dtdz_brks <- pretty(ds_dz,n=20)
dsdz_cols <- c(strat_n_col(length(which(dsdz_brks<0))),
               strat_p_col(length(which(dsdz_brks>0))))


### ----------------- plots -----------------
setwd('~/Documents/R/Github/waltonsmith/figures')

png(paste0(cruise,'_bottom.png'), height = 11, width = 4, units = 'in', res=300)
par(mfrow=c(3,1),mar=c(4.5,4,2,1),oma=c(4,1,4,1))
imagePlot(temp_kriged2$x,
          temp_kriged2$y,
          tempF,
          col=tempF_cols,breaks=tempF_breaks,asp=1,
          xlab='',ylab='',las=1,
          xlim=xlims,ylim=ylims,
          nlevel=length(tempF_cols),legend.width=.7,legend.mar=3)
polygon(masks$longitude[c(1:8,1)],masks$latitude[c(1:8,1)],col='white',border='white')
polygon(masks$longitude[c(9:14,9)],masks$latitude[c(9:14,9)],col='white',border='white')
# contour(temp_kriged2$x,
#         temp_kriged2$y,
#         temp_kriged2$z,
#         levels=temp_breaks,add=T)
image(temp_se2,add=T,breaks=quantile(temp_se2$z,c(.5,1),na.rm=T),col='white')
image(topo_lon,topo_lat,topo,breaks=c(-1,100),col='white',add=T)
plot(world,col='gray70',add=T)
contour(topo_lon,topo_lat,topo,add=T,levels=c(-100,-50,-25,-10),col='gray40')
points(data3$Longitude.Decimal,data3$Latitude.Decimal,pch=16,col='gray50',cex=.5)
# mtext(expression(paste('Longitude (',degree,'W)')),1,line=3,cex=.75)
mtext(expression(paste('Latitude (',degree,'N)')),2,line=3,cex=.75)
mtext(expression(paste('Bottom Temperature (',degree,'F)')),adj=1,cex=.75)

imagePlot(sal_kriged2$x,
          sal_kriged2$y,
          sal_kriged2$z,
          col=sal_cols,breaks=sal_breaks,asp=1,
          xlab='',ylab='',las=1,
          xlim=xlims,ylim=ylims,
          nlevel=length(sal_cols),legend.width=.7,legend.mar=3)
polygon(masks$longitude[c(1:8,1)],masks$latitude[c(1:8,1)],col='white',border='white')
polygon(masks$longitude[c(9:14,9)],masks$latitude[c(9:14,9)],col='white',border='white')
# contour(sal_kriged2$x,
#         sal_kriged2$y,
#         sal_kriged2$z,
#         levels=sal_breaks,add=T)
# image(sal_se2,add=T,breaks=quantile(sal_se2$z,c(.5,1),na.rm=T),col='white')
image(temp_se2,add=T,breaks=quantile(temp_se2$z,c(.5,1),na.rm=T),col='white')
image(topo_lon,topo_lat,topo,breaks=c(-1,100),col='white',add=T)
plot(world,col='gray70',add=T)
contour(topo_lon,topo_lat,topo,add=T,levels=c(-100,-50,-25,-10),col='gray40')
points(data3$Longitude.Decimal,data3$Latitude.Decimal,pch=16,col='gray50',cex=.5)
# mtext(expression(paste('Longitude (',degree,'W)')),1,line=3,cex=.75)
mtext(expression(paste('Latitude (',degree,'N)')),2,line=3,cex=.75)
mtext('Salinity (PSU)',adj=1,cex=.75)

imagePlot(do_kriged2$x,
          do_kriged2$y,
          do_kriged2$z,
          col=cols,breaks=breaks,asp=1,
          xlab='',ylab='',las=1,
          xlim=xlims,ylim=ylims,
          nlevel=length(cols),legend.width=.7,legend.mar=3)
polygon(masks$longitude[c(1:8,1)],masks$latitude[c(1:8,1)],col='white',border='white')
polygon(masks$longitude[c(9:14,9)],masks$latitude[c(9:14,9)],col='white',border='white')
# contour(do_kriged2$x,
#         do_kriged2$y,
#         do_kriged2$z,
#         levels=breaks,add=T)
# image(do_se2,add=T,breaks=quantile(do_se2$z,c(.45,1),na.rm=T),col='white')
image(temp_se2,add=T,breaks=quantile(temp_se2$z,c(.5,1),na.rm=T),col='white')
image(topo_lon,topo_lat,topo,breaks=c(-1,100),col='white',add=T)
plot(world,col='gray70',add=T)
contour(topo_lon,topo_lat,topo,add=T,levels=c(-100,-50,-25,-10),col='gray40')
points(data3$Longitude.Decimal,data3$Latitude.Decimal,pch=16,col='gray50',cex=.5)
mtext(expression(paste('Longitude (',degree,'W)')),1,line=3,cex=.75)
mtext(expression(paste('Latitude (',degree,'N)')),2,line=3,cex=.75)
mtext(expression(paste('Bottom DO (mg l'^-1,')')),adj=1,cex=.75)

mtext('Walton Smith Bulletin',
      outer=T,line=1,side=3,font=2,at=.05,adj=0,cex=1.25)
mtext(paste('Collected:',
            paste(
              paste(month.abb[month(data3$Date.GMT[1])],
                    day(data3$Date.GMT[1])),
              paste(month.abb[month(data3$Date.GMT[nrow(data3)])],
                    day(data3$Date.GMT[nrow(data3)])),
              sep='-')),
      outer=T,line=-.1,side=3,at=.05,adj=0)
mtext(paste('Note: Data are early release and subject to further QA/QC, \nplease contact brendan.turley@noaa.gov with concerns \nProcessed: ',as.Date(Sys.time())),
      outer=T,line=2,side=1,col='red',font=2,at=.01,adj=0,cex=.75)
dev.off()


png(paste0(cruise,'_bottom_tracks.png'), height = 10, width = 9, units = 'in', res=300)
par(mfrow=c(2,2),mar=c(4.5,5,2,3),oma=c(4,1,4,2))
imagePlot(temp_kriged2$x,
          temp_kriged2$y,
          temp_kriged2$z,
          col=temp_cols,breaks=temp_breaks,asp=1,
          xlab='',ylab='',las=1,
          xlim=xlims,ylim=ylims,
          nlevel=length(temp_cols),legend.width=1.2,legend.mar=3)
polygon(masks$longitude[c(1:8,1)],masks$latitude[c(1:8,1)],col='white',border='white')
polygon(masks$longitude[c(9:14,9)],masks$latitude[c(9:14,9)],col='white',border='white')
# contour(temp_kriged2$x,
#         temp_kriged2$y,
#         temp_kriged2$z,
#         levels=temp_breaks,add=T)
image(temp_se2,add=T,breaks=quantile(temp_se2$z,c(.5,1),na.rm=T),col='white')
image(topo_lon,topo_lat,topo,breaks=c(-1,100),col='white',add=T)
plot(world,col='gray70',add=T)
contour(topo_lon,topo_lat,topo,add=T,levels=c(-100,-50,-25,-10),col='gray40')
points(data3$Longitude.Decimal,data3$Latitude.Decimal,pch=16,col='gray50',cex=.5)
# mtext(expression(paste('Longitude (',degree,'W)')),1,line=3,cex=.75)
mtext(expression(paste('Latitude (',degree,'N)')),2,line=3,cex=.75)
mtext(expression(paste('Bottom Temperature (',degree,'C)')),adj=1,cex=.75)

imagePlot(sal_kriged2$x,
          sal_kriged2$y,
          sal_kriged2$z,
          col=sal_cols,breaks=sal_breaks,asp=1,
          xlab='',ylab='',las=1,
          xlim=xlims,ylim=ylims,
          nlevel=length(sal_cols),legend.width=1.2,legend.mar=3)
polygon(masks$longitude[c(1:8,1)],masks$latitude[c(1:8,1)],col='white',border='white')
polygon(masks$longitude[c(9:14,9)],masks$latitude[c(9:14,9)],col='white',border='white')
# contour(sal_kriged2$x,
#         sal_kriged2$y,
#         sal_kriged2$z,
#         levels=sal_breaks,add=T)
# image(sal_se2,add=T,breaks=quantile(sal_se2$z,c(.5,1),na.rm=T),col='white')
image(temp_se2,add=T,breaks=quantile(temp_se2$z,c(.5,1),na.rm=T),col='white')
image(topo_lon,topo_lat,topo,breaks=c(-1,100),col='white',add=T)
plot(world,col='gray70',add=T)
contour(topo_lon,topo_lat,topo,add=T,levels=c(-100,-50,-25,-10),col='gray40')
points(data3$Longitude.Decimal,data3$Latitude.Decimal,pch=16,col='gray50',cex=.5)
# mtext(expression(paste('Longitude (',degree,'W)')),1,line=3,cex=.75)
# mtext(expression(paste('Latitude (',degree,'N)')),2,line=3,cex=.75)
mtext('Salinity (PSU)',adj=1,cex=.75)

imagePlot(do_kriged2$x,
          do_kriged2$y,
          do_kriged2$z,
          col=cols,breaks=breaks,asp=1,
          xlab='',ylab='',las=1,
          xlim=xlims,ylim=ylims,
          nlevel=length(cols),legend.width=1.2,legend.mar=3)
polygon(masks$longitude[c(1:8,1)],masks$latitude[c(1:8,1)],col='white',border='white')
polygon(masks$longitude[c(9:14,9)],masks$latitude[c(9:14,9)],col='white',border='white')
# contour(do_kriged2$x,
#         do_kriged2$y,
#         do_kriged2$z,
#         levels=breaks,add=T)
# image(do_se2,add=T,breaks=quantile(do_se2$z,c(.45,1),na.rm=T),col='white')
image(temp_se2,add=T,breaks=quantile(temp_se2$z,c(.5,1),na.rm=T),col='white')
image(topo_lon,topo_lat,topo,breaks=c(-1,100),col='white',add=T)
plot(world,col='gray70',add=T)
contour(topo_lon,topo_lat,topo,add=T,levels=c(-100,-50,-25,-10),col='gray40')
points(data3$Longitude.Decimal,data3$Latitude.Decimal,pch=16,col='gray50',cex=.5)
mtext(expression(paste('Longitude (',degree,'W)')),1,line=3,cex=.75)
mtext(expression(paste('Latitude (',degree,'N)')),2,line=3,cex=.75)
mtext(expression(paste('Bottom DO (mg l'^-1,')')),adj=1,cex=.75)

image(do_kriged2$x,
          do_kriged2$y,
          do_kriged2$z,
          col=cols,breaks=breaks,asp=1,
          xlab='',ylab='',las=1,
          xlim=xlims,ylim=ylims,
          nlevel=length(cols),legend.width=.7,legend.mar=3)
image(topo_lon,topo_lat,topo,breaks=c(-200,100),col='white',add=T)
plot(world,col='gray70',add=T)
contour(topo_lon,topo_lat,topo,add=T,levels=c(-100,-50,-25,-10),col='gray40')
points(data3$Longitude.Decimal,data3$Latitude.Decimal,pch=16,col=plasma(nrow(data3)),typ='o')
text(data3$Longitude.Decimal[1],data3$Latitude.Decimal[1],'Start',pos=2,font=2)
text(data3$Longitude.Decimal[nrow(data3)],data3$Latitude.Decimal[nrow(data3)],'End',pos=1,font=2)
mtext(expression(paste('Longitude (',degree,'W)')),1,line=3,cex=.75)
# mtext(expression(paste('Latitude (',degree,'N)')),2,line=3,cex=.75)
mtext('Cruise track',adj=1,cex=.75)

mtext('Walton Smith Bulletin',
      outer=T,line=1,side=3,font=2,at=.05,adj=0,cex=1.25)
mtext(paste('Collected:',
            paste(
              paste(month.abb[month(data3$Date.GMT[1])],
                    day(data3$Date.GMT[1])),
              paste(month.abb[month(data3$Date.GMT[nrow(data3)])],
                    day(data3$Date.GMT[nrow(data3)])),
              sep='-')),
      outer=T,line=-.1,side=3,at=.05,adj=0)
mtext(paste('Note: Data are early release and subject to further QA/QC, \nplease contact brendan.turley@noaa.gov with concerns \nProcessed: ',as.Date(Sys.time())),
      outer=T,line=2,side=1,col='red',font=2,at=.01,adj=0,cex=.75)
dev.off()


imagePlot(dtdz_kriged2$x,
          dtdz_kriged2$y,
          dtdz_kriged2$z,
          col=dtdz_cols,breaks=dtdz_brks,asp=1,
          xlab='',ylab='',las=1,
          xlim=xlims,ylim=ylims,
          nlevel=length(dtdz_cols),legend.width=.7,legend.mar=3)
polygon(masks$longitude[c(1:8,1)],masks$latitude[c(1:8,1)],col='white',border='white')
polygon(masks$longitude[c(9:14,9)],masks$latitude[c(9:14,9)],col='white',border='white')
image(temp_se2,add=T,breaks=quantile(temp_se2$z,c(.4,1),na.rm=T),col='white')
image(topo_lon,topo_lat,topo,breaks=c(-2,100),col='white',add=T)
plot(world,col='gray70',add=T)
contour(topo_lon,topo_lat,topo,add=T,levels=c(-100,-50,-25,-10),col='gray40')
points(data3$Longitude.Decimal,data3$Latitude.Decimal,pch=16,col='gray50',cex=.5)

imagePlot(dsdz_kriged2$x,
          dsdz_kriged2$y,
          dsdz_kriged2$z,
          col=dsdz_cols,breaks=dsdz_brks,asp=1,
          xlab='',ylab='',las=1,
          xlim=xlims,ylim=ylims,
          nlevel=length(dsdz_cols),legend.width=.7,legend.mar=3)
polygon(masks$longitude[c(1:8,1)],masks$latitude[c(1:8,1)],col='white',border='white')
polygon(masks$longitude[c(9:14,9)],masks$latitude[c(9:14,9)],col='white',border='white')
image(temp_se2,add=T,breaks=quantile(temp_se2$z,c(.4,1),na.rm=T),col='white')
image(topo_lon,topo_lat,topo,breaks=c(-2,100),col='white',add=T)
plot(world,col='gray70',add=T)
contour(topo_lon,topo_lat,topo,add=T,levels=c(-100,-50,-25,-10),col='gray40')
points(data3$Longitude.Decimal,data3$Latitude.Decimal,pch=16,col='gray50',cex=.5)


setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
png('new_plot.png',width=9,height=12,units='in',pointsize=12,res=300)
par(mfrow=c(2,2),mar=c(4,4,3,4),oma=c(5,0,3,1))
data_plot(data3$Longitude.Decimal,data3$Latitude.Decimal,data3$Temperature.CTD.data,temp_col,title=expression(paste('Bottom temperature (',degree,'F)')))
data_plot(data3$Longitude.Decimal,data3$Latitude.Decimal,data3$Salinity.CTD.data,sal_col,title='Bottom salinity (psu)')

plot(data3$Longitude.Decimal,data3$Latitude.Decimal,pch=21,asp=1,cex=2,las=1,
     bg=cols[as.numeric(cut(data3$Oxygen.mg.l..CTD.data,breaks))],
     col=cols[as.numeric(cut(data3$Oxygen.mg.l..CTD.data,breaks))],
     xlab='',ylab='')
imagePlot(zlim=range(breaks,na.rm=T),breaks=breaks,col=cols,legend.only=TRUE,legend.width = 1.5)
mtext(expression(paste('Bottom dissolved oxygen (mg l'^-1,')')))
plot(world,add=T,col='gray80')
contour(topo_lon,
        topo_lat,
        topo,
        add=T,levels=c(-200,-100,-50,-25,-10),col='gray70')
dev.off()