
library(akima)
library(fields)
library(lubridate)
library(ncdf4)
library(NISTunits)
library(oce)
library(raster)
library(rgdal)
library(scales)

source('~/Desktop/professional/biblioteca/scripts/color_breaks.R')

setwd("~/Desktop/professional/biblioteca/data")
bathy <- nc_open('etopo1.nc')
topo <- ncvar_get(bathy, 'Band1')
topo_lat <- ncvar_get(bathy, 'lat')
topo_lon <- ncvar_get(bathy, 'lon')
nc_close(bathy)

setwd("~/Desktop/professional/biblioteca/data/shapefiles/ne_10m_admin_0_countries")
world <- readOGR('ne_10m_admin_0_countries.shp')
world <- crop(world, extent(-100, -74, 24.5, 40))


### colors
d_col <- colorRampPalette(rev(c('purple4','purple2','orchid1','gray90')))
m_col <- colorRampPalette(rev(c('blue4','dodgerblue2','deepskyblue1','gray90')))
t_col <- colorRampPalette(c(1,'purple','darkorange','gold'))
s_col <- colorRampPalette(c('purple4','dodgerblue4','seagreen3','khaki1'))
c_col <- colorRampPalette(c('honeydew2','darkseagreen3','darkgreen'))
ox.col1 <- colorRampPalette(c(1,'firebrick4','red'))
ox.col2 <- colorRampPalette(c('darkgoldenrod4','goldenrod2','gold'))
ox.col3 <- colorRampPalette(c('midnightblue','dodgerblue4','deepskyblue2','cadetblue1','azure'))
o_breaks <- seq(0,10,by=.25)
o_cols <- c(ox.col1(length(o_breaks[o_breaks<2])),
          ox.col2(length(o_breaks[o_breaks>=2 & o_breaks<3.5])),
          ox.col3(length(o_breaks[o_breaks>=3.5])-1))


setwd('~/Desktop/noaa_ctd')
files <- list.files()
files <- files[grep('.cnv',files)]
write.table(files,paste(as.Date(Sys.time()),'files.txt',sep='_'),row.names = F,col.names = F,quote=F)

lonbox_e <- -81.5 ### Florida Bay
lonbox_w <- -87 ### mouth of Mississippi River
latbox_n <- 30.5 ### northern coast
latbox_s <- 24.3 ### southern edge of Key West

bottom_dat <- data.frame(filename=rep(NA,length(files)),
                         date_utc=rep(NA,length(files)),
                         lon=rep(NA,length(files)),
                         lat=rep(NA,length(files)),
                         ship=rep(NA,length(files)),
                         cruise=rep(NA,length(files)),
                         station=rep(NA,length(files)),
                         bot_t=rep(NA,length(files)),
                         bot_s=rep(NA,length(files)),
                         bot_c=rep(NA,length(files)),
                         bot_do=rep(NA,length(files)),
                         r_depth=rep(NA,length(files)),
                         m_depth=rep(NA,length(files)),
                         dd_dz1=rep(NA,length(files)),
                         dd_dz2=rep(NA,length(files)),
                         mld=rep(NA,length(files)),
                         surf_s=rep(NA,length(files)),
                         surf_c=rep(NA,length(files)),
                         pycno=rep(NA,length(files)),
                         z_cmax=rep(NA,length(files)),
                         chl_max=rep(NA,length(files)),
                         chl_p50=rep(NA,length(files))
)
for(i in 1:length(files)){
  data <- read.ctd(files[i])
  ### filter by location
  lat <- data@metadata$latitude
  lon <- data@metadata$longitude
  # if(lat<latbox_n & lat>latbox_s &
  #    lon<lonbox_e & lon>lonbox_w){
    ind_bot <- which.max(data@data$depth)
    bottom_dat$filename[i] <- as.character(files[i])
    if(is.na(data@metadata$date)){
      bottom_dat$date_utc[i] <- as.character(data@metadata$startTime)
    } else {
      bottom_dat$date_utc[i] <- as.character(data@metadata$date)
    }
    bottom_dat$lat[i] <- data@metadata$latitude
    bottom_dat$lon[i] <- data@metadata$longitude
    bottom_dat$ship[i] <- as.character(data@metadata$ship)
    bottom_dat$cruise[i] <- data@metadata$cruise
    bottom_dat$station[i] <- data@metadata$station
    bottom_dat$bot_t[i] <- data@data$temperature[ind_bot]
    bottom_dat$bot_s[i] <- data@data$salinity[ind_bot]
    bottom_dat$bot_c[i] <- data@data$fluorescence[ind_bot]
    bottom_dat$bot_do[i] <- data@data$oxygen[ind_bot]
    bottom_dat$m_depth[i] <- data@data$depth[ind_bot]
    bottom_dat$r_depth[i] <- data@metadata$waterDepth
    bottom_dat$surf_s[i] <- data@data$salinity[1]
    bottom_dat$surf_c[i] <- data@data$fluorescence[1]
    
    if(length(data@data$depth)>=4){
      ### stratification
      y <- smooth.spline(data@data$depth,data@data$density,df=length(data@data$density)/3)
      dd_dz <-  diff(y$y)/diff(data@data$depth)
      bottom_dat$dd_dz1[i] <- max(dd_dz,na.rm=T)
      bottom_dat$dd_dz2[i] <- mean(dd_dz[which(data@data$depth<=50)],na.rm=T)
      bottom_dat$pycno[i] <- data@data$depth[which.max(dd_dz)]
    }
    
    ### mld
    mld_ind <- which(data@data$density>mean(data@data$density[1:3],na.rm=T)+.125)[1]
    # mld_ind <- which(data@data$density>data@data$density[1]+.125)[1]
    bottom_dat$mld[i] <- data@data$depth[mld_ind]
    
    ### chl max
    chl_max <- which.max(data@data$fluorescence)
    bottom_dat$z_cmax[i] <- data@data$depth[chl_max]
    bottom_dat$chl_max[i] <- data@data$fluorescence[chl_max]
    bottom_dat$chl_p50[i] <- median(data@data$fluorescence,na.rm=T)
    
  # }
}
bottom_dat <- bottom_dat[!is.na(bottom_dat$date_utc),]
bottom_dat$date_utc <- ymd_hms(bottom_dat$date_utc)
bottom_dat$mld2 <- bottom_dat$mld/bottom_dat$r_depth
bottom_dat$mld2[which(is.na(bottom_dat$mld2))] <- 1
bottom_dat$z_cmax2 <- bottom_dat$z_cmax/bottom_dat$r_depth
bottom_dat$z_cmax2[which(is.na(bottom_dat$z_cmax2))] <- 1

write.csv(bottom_dat,'NOAA_NMFS_bot_dat.csv',row.names = F)

plot(bottom_dat$lon,bottom_dat$lat,asp=1,cex=bottom_dat$r_depth/100)
plot(bottom_dat$lon,bottom_dat$lat,asp=1,cex=bottom_dat$m_depth/100)
plot(world,add=T)
hist(bottom_dat$r_depth-bottom_dat$m_depth)

### pick cruise
table(bottom_dat$cruise,month(bottom_dat$date_utc),year(bottom_dat$date_utc))
sort(unique(bottom_dat$cruise))
cruise <- sort(unique(bottom_dat$cruise))[4]
ind <- which(bottom_dat$cruise==cruise & 
               bottom_dat$lat<latbox_n & bottom_dat$lat>latbox_s &
               bottom_dat$lon<lonbox_e & bottom_dat$lon>lonbox_w)
bot_plot <- bottom_dat[ind,]
plot(bot_plot$lon,bot_plot$lat,asp=1,cex=bot_plot$r_depth/100)
### fix NA mlds to bottom
bot_plot$mld[which(is.na(bot_plot$mld))] <- bot_plot$m_depth[which(is.na(bot_plot$mld))]


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

par(mfrow=c(2,2),mar=c(5,4,2,3.5))
data_plot(bot_plot$lon,bot_plot$lat,bot_plot$chl_max,c_col,title=expression(paste('Chlorophyll max (mg m'^-3,')')))
data_plot(bot_plot$lon,bot_plot$lat,bot_plot$chl_p50,c_col,title=expression(paste('Chlorophyll median (mg m'^-3,')')))
data_plot(bot_plot$lon,bot_plot$lat,bot_plot$z_cmax,m_col,title='Chlorophyll max depth (m)')
data_plot(bot_plot$lon,bot_plot$lat,bot_plot$z_cmax2,m_col,title='Scaled Chlorophyll max depth')

data_plot(bot_plot$lon,bot_plot$lat,bot_plot$bot_c,c_col,title=expression(paste('Bottom chlorophyll (mg m'^-3,')')))
data_plot(bot_plot$lon,bot_plot$lat,bot_plot$surf_c,c_col,title=expression(paste('Surface chlorophyll (mg m'^-3,')')))

bot_plot$bot_t_F <- NISTdegCtOdegF(bot_plot$bot_t)

setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
png('new_plot.png',width=9,height=12,units='in',pointsize=12,res=300)
par(mfrow=c(2,2),mar=c(4,4,3,4),oma=c(5,0,3,1))
data_plot(bot_plot$lon,bot_plot$lat,bot_plot$bot_t_F,t_col,title=expression(paste('Bottom temperature (',degree,'F)')))
data_plot(bot_plot$lon,bot_plot$lat,bot_plot$bot_s,s_col,title='Bottom salinity (psu)')
data_plot(bot_plot$lon,bot_plot$lat,bot_plot$bot_c,c_col,title=expression(paste('Bottom chlorophyll (mg m'^-3,')')))

plot(bot_plot$lon,bot_plot$lat,pch=21,asp=1,cex=2,las=1,
     bg=o_cols[as.numeric(cut(bot_plot$bot_do,o_breaks))],
     col=o_cols[as.numeric(cut(bot_plot$bot_do,o_breaks))],
     xlab='',ylab='')
imagePlot(zlim=range(o_breaks,na.rm=T),breaks=o_breaks,col=o_cols,legend.only=TRUE,legend.width = 1.5)
mtext(expression(paste('Bottom dissolved oxygen (mg l'^-1,')')))
plot(world,add=T,col='gray80')
contour(topo_lon,
        topo_lat,
        topo,
        add=T,levels=c(-200,-100,-50,-25,-10),col='gray70')

mtext('NMFS-SEAMAP Fall Plankton Cruise',
# mtext('NMFS Bottom Longline Survey',
      outer=T,line=1,side=3,font=2,at=.05,adj=0,cex=1.25)
mtext(paste(paste('Collected: R/V',
            bot_plot$ship[1], sep=' '),
            paste(
              paste(month.abb[month(bot_plot$date_utc[1])],
                    day(bot_plot$date_utc[1])),
              paste(month.abb[month(bot_plot$date_utc[nrow(bot_plot)])],
                    day(bot_plot$date_utc[nrow(bot_plot)])),
              sep=' - '),
            year(bot_plot$date_utc[1]),sep=', '),
      outer=T,line=-.1,side=3,at=.05,adj=0)
mtext('Bottom contours: 10, 25, 50, 100, 200 meters',
      outer=T,line=0,side=1,at=.05,adj=0)
# mtext(paste('Note: Data are early release and subject to further QA/QC, \nplease contact brendan.turley@noaa.gov for comments / concerns \nUpdated: ',as.Date(Sys.time())),
      # outer=T,line=4,side=1,col='red',font=2,at=.05,adj=0)
mtext(paste('Note: Data are early release and subject to further QA/QC'),
      outer=T,line=1,side=1,col='red',font=2,at=.05,adj=0)
mtext(paste('please contact brendan.turley@noaa.gov for comments / concerns'),
      outer=T,line=2,side=1,col='red',font=2,at=.05,adj=0)
mtext(paste('Updated: ',as.Date(Sys.time())),
      outer=T,line=3,side=1,col='red',font=2,at=.05,adj=0)
dev.off()

par(mfrow=c(2,2),mar=c(5,4,3,5))
data_plot(bot_plot$lon,bot_plot$lat,bot_plot$pycno,m_col,title='Pycnocline depth (m)')
data_plot(bot_plot$lon,bot_plot$lat,bot_plot$mld,m_col,title='Mixed-layer depth (m)')
data_plot(bot_plot$lon,bot_plot$lat,bot_plot$z_cmax,m_col,title='Chlorophyll max depth (m)')

par(mfrow=c(2,2),mar=c(5,4,3,5))
data_plot(bot_plot$lon,bot_plot$lat,bot_plot$dd_dz1,d_col,title='Density stratification')
data_plot(bot_plot$lon,bot_plot$lat,bot_plot$dd_dz2,d_col,title='Density stratification')
data_plot(bot_plot$lon,bot_plot$lat,bot_plot$mld,m_col,title='Mixed-layer depth (m)')
data_plot(bot_plot$lon,bot_plot$lat,bot_plot$mld2,m_col,title='Scaled MLD')



### Kriging locations
locs <- cbind(bot_plot$lon,bot_plot$lat)
### Kriging
parms <- c(8:11,14:19)
names(bot_plot)[parms]
system.time(
  for(i in parms){
    parm <- bot_plot[,i]
    print(names(bot_plot)[i])
    test <- tryCatch(spatialProcess(locs,parm,
                                    Distance = "rdist.earth"),
                     error=function (x) x,
                     warning=function(w) w)
    if(any(class(test)=='warning')){
      test <- tryCatch(Tps(locs,parm),
                       error=function (x) x,
                       warning=function(w) w)
    }
    if(any(class(test)=='warning')){
      test <- Krig(locs,parm,
                   cov.function = "stationary.cov",
                   Covariance="Matern",
                   GCV=F)
    }
    
    kriged <- predictSurface(test,
                             list(lon=seq(min(bot_plot$lon,na.rm=T), max(bot_plot$lon,na.rm=T),.01),
                                  lat=seq(min(bot_plot$lat,na.rm=T), max(bot_plot$lat,na.rm=T),.01)))
    imagePlot(kriged)
    surf_se <- predictSurfaceSE(test,
                                grid.list = list(lon=seq(min(bot_plot$lon,na.rm=T), max(bot_plot$lon,na.rm=T),.01),
                                                 lat=seq(min(bot_plot$lat,na.rm=T), max(bot_plot$lat,na.rm=T),.01)))
    assign(paste(names(bot_plot)[i],'kriged',sep='_'),kriged)
    assign(paste(names(bot_plot)[i],'se',sep='_'),surf_se)
  }
)
mld2_kriged$z[which(mld2_kriged$z<0)] <- 0
mld2_kriged$z[which(mld2_kriged$z>1)] <- 1


par(mfrow=c(3,2),mar=c(5,5,3,2),oma=c(5,0,1,0))
imagePlot(bot_t_se,asp=1)
imagePlot(bot_s_se,asp=1)
imagePlot(bot_c_se,asp=1)
imagePlot(bot_do_se,asp=1)
imagePlot(dd_dz2_se,asp=1)
imagePlot(mld2_se,asp=1)

par(mfrow=c(2,2))
imagePlot(surf_c_kriged,asp=1)
points(bot_plot$lon,bot_plot$lat,pch=20,col='green')
imagePlot(bot_c_kriged,asp=1)
points(bot_plot$lon,bot_plot$lat,pch=20,col='green')
imagePlot(surf_s_kriged,asp=1)
points(bot_plot$lon,bot_plot$lat,pch=20,col='green')
imagePlot(bot_s_kriged,asp=1)
points(bot_plot$lon,bot_plot$lat,pch=20,col='green')

par(mfrow=c(1,2))
imagePlot(surf_c_kriged$x,
          surf_c_kriged$y,
  surf_c_kriged$z-bot_c_kriged$z,asp=1)
points(bot_plot$lon,bot_plot$lat,pch=20,col='green')
imagePlot(surf_s_kriged$x,
          surf_s_kriged$y,
          surf_s_kriged$z-bot_s_kriged$z,asp=1)
points(bot_plot$lon,bot_plot$lat,pch=20,col='green')

par(mfrow=c(1,2))
imagePlot(dd_dz1_kriged,asp=1)
points(bot_plot$lon,bot_plot$lat,pch=20,col='green')
imagePlot(dd_dz2_kriged,asp=1)
points(bot_plot$lon,bot_plot$lat,pch=20,col='green')

par(mfrow=c(1,2))
imagePlot(mld_kriged,asp=1)
points(bot_plot$lon,bot_plot$lat,pch=20,col='green')
imagePlot(mld2_kriged,asp=1)
points(bot_plot$lon,bot_plot$lat,pch=20,col='green')



t_breaks <- breaks(bot_t_kriged$z)
t_cols <- t_col(length(t_breaks)-1)
s_breaks <- breaks(bot_s_kriged$z,.2)
s_cols <- s_col(length(s_breaks)-1)
c_breaks <- breaks(round(bot_c_kriged$z,2),.25)
c_cols <- c_col(length(c_breaks)-1)
d1_breaks <- breaks(round(dd_dz1_kriged$z,2),.01,decimal = T)
d1_cols <- d_col(length(d1_breaks)-1)
d2_breaks <- breaks(round(dd_dz2_kriged$z,2),.001,decimal = T)
d2_cols <- d_col(length(d2_breaks)-1)
m1_breaks <- breaks(mld_kriged$z,2)
m1_cols <- m_col(length(m1_breaks)-1)
m_breaks <- seq(0,1,.1)
m_cols <- m_col(length(m_breaks)-1)


setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
png('new_plot.png',width=9,height=12,units='in',pointsize=12,res=300)
par(mfrow=c(3,2),mar=c(5,5,3,2),oma=c(5,0,1,0))
imagePlot(bot_t_kriged,asp=1,breaks=t_breaks,col=t_cols)
contour(bot_t_kriged,levels=t_breaks,add=T)
# image(bot_t_se,breaks=c(1.2,10),col=alpha('gray50',.5),add=T)
contour(topo_lon,
        topo_lat,
        topo,
        add=T,levels=c(-100),col='gray70')
plot(world,add=T,col='gray80')
points(bot_plot$lon,bot_plot$lat,pch=20,col='gray70')
mtext(expression(paste('Temperature (',degree,'C)')),adj=1)
# mtext(expression(paste('Longitude (',degree,'W)')),1,2.5)
mtext(expression(paste('Latitude (',degree,'N)')),2,2)

mtext(bot_plot$ship[1],adj=0,line=1)
mtext(paste(paste(month.abb[month(bot_plot$date_utc[1])],day(bot_plot$date_utc[1])),
            paste(month.abb[month(bot_plot$date_utc[nrow(bot_plot)])],day(bot_plot$date_utc[nrow(bot_plot)])),
            sep=' - '),
      adj=0)

imagePlot(bot_s_kriged,asp=1,breaks=s_breaks,col=s_cols)
contour(bot_s_kriged,levels=s_breaks,add=T)
# image(bot_t_se,breaks=c(1.2,10),col=alpha('gray50',.5),add=T)
contour(topo_lon,
        topo_lat,
        topo,
        add=T,levels=c(-100),col='gray70')
plot(world,add=T,col='gray80')
points(bot_plot$lon,bot_plot$lat,pch=20,col='gray70')
mtext('Salinity',adj=1)

imagePlot(bot_c_kriged,asp=1,breaks=c_breaks,col=c_cols)
contour(bot_c_kriged,levels=c_breaks,add=T)
# image(bot_t_se,breaks=c(1.2,10),col=alpha('gray50',.5),add=T)
contour(topo_lon,
        topo_lat,
        topo,
        add=T,levels=c(-100),col='gray70')
plot(world,add=T,col='gray80')
points(bot_plot$lon,bot_plot$lat,pch=20,col='gray70')
mtext(expression(paste('Chlorophyll (mg m'^-3,')')),adj=1)
# mtext(expression(paste('Longitude (',degree,'W)')),1,2.5)
mtext(expression(paste('Latitude (',degree,'N)')),2,2)

imagePlot(bot_do_kriged,asp=1,breaks=o_breaks,col=o_cols)
contour(bot_do_kriged,levels=o_breaks,add=T)
# image(bot_t_se,breaks=c(1.2,10),col=alpha('gray50',.5),add=T)
contour(topo_lon,
        topo_lat,
        topo,
        add=T,levels=c(-100),col='gray70')
plot(world,add=T,col='gray80')
points(bot_plot$lon,bot_plot$lat,pch=20,col='gray70')
mtext(expression(paste('Dissolved Oxygen (mg l'^-1,')')),adj=1)

imagePlot(dd_dz1_kriged,asp=1,breaks=d1_breaks,col=d1_cols)
contour(dd_dz1_kriged,levels=d1_breaks,add=T)
# imagePlot(dd_dz2_kriged,asp=1,breaks=d2_breaks,col=d2_cols)
# contour(dd_dz2_kriged,levels=d2_breaks,add=T)
# image(bot_t_se,breaks=c(1.2,10),col=alpha('gray50',.5),add=T)
contour(topo_lon,
        topo_lat,
        topo,
        add=T,levels=c(-100),col='gray70')
plot(world,add=T,col='gray80')
points(bot_plot$lon,bot_plot$lat,pch=20,col='gray70')
mtext(expression(paste('Stratification (d',rho,' dz'^-1,')')),adj=1)
mtext(expression(paste('Longitude (',degree,'W)')),1,2.5)
mtext(expression(paste('Latitude (',degree,'N)')),2,2)

imagePlot(mld2_kriged,asp=1,breaks=m_breaks,col=m_cols)
contour(mld2_kriged,levels=m_breaks,add=T)
# imagePlot(mld_kriged,asp=1,breaks=m1_breaks,col=m1_cols)
# contour(mld_kriged,levels=m1_breaks,add=T)
# image(bot_t_se,breaks=c(1.2,10),col=alpha('gray50',.5),add=T)
contour(topo_lon,
        topo_lat,
        topo,
        add=T,levels=c(-100),col='gray70')
plot(world,add=T,col='gray80')
points(bot_plot$lon,bot_plot$lat,pch=20,col='gray70')
mtext('Scaled MLD',adj=1)
mtext(expression(paste('Longitude (',degree,'W)')),1,2.5)
# mtext(expression(paste('Latitude (',degree,'N)')),2,2)

mtext('R/V Oregon II Bulletin',
      outer=T,line=1,side=3,font=2,at=.05,adj=0,cex=1.25)
mtext(paste('Collected:',
            paste(
              paste(month.abb[month(bot_plot$date_utc[1])],
                    day(bot_plot$date_utc[1])),
              paste(month.abb[month(bot_plot$date_utc[nrow(bot_plot)])],
                    day(bot_plot$date_utc[nrow(bot_plot)])),
              sep='-')),
      outer=T,line=-.1,side=3,at=.05,adj=0)
mtext(paste('Note: Data are early release and subject to further QA/QC, \nplease contact brendan.turley@noaa.gov for comments/concerns \nUpdated: ',as.Date(Sys.time())),
      outer=T,line=3,side=1,col='red',font=2,at=.05,adj=0)
dev.off()


### any unacceptable plots
bot_do_kriged <- interp(bot_plot$lon,bot_plot$lat,bot_plot$bot_do,
                        xo=seq(min(bot_plot$lon,na.rm=T), max(bot_plot$lon,na.rm=T),.01),
                        yo=seq(min(bot_plot$lat,na.rm=T), max(bot_plot$lat,na.rm=T),.01),linear = F)

dd_dz2_kriged <- interp(bot_plot$lon,bot_plot$lat,bot_plot$dd_dz2,
                        xo=seq(min(bot_plot$lon,na.rm=T), max(bot_plot$lon,na.rm=T),.01),
                        yo=seq(min(bot_plot$lat,na.rm=T), max(bot_plot$lat,na.rm=T),.01),linear = F)
dd_dz2_kriged$z[which(dd_dz2_kriged$z<0)] <- 0
d2_breaks <- breaks(round(dd_dz2_kriged$z,2),.01,decimal = T)
d2_cols <- d_col(length(d2_breaks)-1)

mld2_kriged <- interp(bot_plot$lon,bot_plot$lat,bot_plot$mld2,
                      xo=seq(min(bot_plot$lon,na.rm=T), max(bot_plot$lon,na.rm=T),.01),
                      yo=seq(min(bot_plot$lat,na.rm=T), max(bot_plot$lat,na.rm=T),.01),linear = F)
mld2_kriged$z[which(mld2_kriged$z<0)] <- 0
mld2_kriged$z[which(mld2_kriged$z>1)] <- 1


setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
png('new_plot.png',width=9,height=12,units='in',pointsize=12,res=300)
par(mfrow=c(3,2),mar=c(5,5,3,2),oma=c(5,0,1,0))
imagePlot(bot_t_kriged,asp=1,breaks=t_breaks,col=t_cols)
contour(bot_t_kriged,levels=t_breaks,add=T)
# image(bot_t_se,breaks=c(1.2,10),col=alpha('gray50',.5),add=T)
contour(topo_lon,
        topo_lat,
        topo,
        add=T,levels=c(-100),col='green')
plot(world,add=T,col='gray80')
points(bot_plot$lon,bot_plot$lat,pch=20,col='green')
mtext(expression(paste('Temperature (',degree,'C)')),adj=1)
# mtext(expression(paste('Longitude (',degree,'W)')),1,2.5)
mtext(expression(paste('Latitude (',degree,'N)')),2,2)

mtext(bot_plot$ship[1],adj=0,line=1)
mtext(paste(paste(month.abb[month(bot_plot$date_utc[1])],day(bot_plot$date_utc[1])),
            paste(month.abb[month(bot_plot$date_utc[nrow(bot_plot)])],day(bot_plot$date_utc[nrow(bot_plot)])),
            sep=' - '),
      adj=0)

imagePlot(bot_s_kriged,asp=1,breaks=s_breaks,col=s_cols)
contour(bot_s_kriged,levels=s_breaks,add=T)
# image(bot_t_se,breaks=c(1.2,10),col=alpha('gray50',.5),add=T)
contour(topo_lon,
        topo_lat,
        topo,
        add=T,levels=c(-100),col='purple')
plot(world,add=T,col='gray80')
points(bot_plot$lon,bot_plot$lat,pch=20,col='purple')
mtext('Salinity',adj=1)

imagePlot(bot_c_kriged,asp=1,breaks=c_breaks,col=c_cols)
contour(bot_c_kriged,levels=c_breaks,add=T)
# image(bot_t_se,breaks=c(1.2,10),col=alpha('gray50',.5),add=T)
contour(topo_lon,
        topo_lat,
        topo,
        add=T,levels=c(-100),col='purple')
plot(world,add=T,col='gray80')
points(bot_plot$lon,bot_plot$lat,pch=20,col='purple')
mtext(expression(paste('Chlorophyll (mg m'^-3,')')),adj=1)
# mtext(expression(paste('Longitude (',degree,'W)')),1,2.5)
mtext(expression(paste('Latitude (',degree,'N)')),2,2)

imagePlot(bot_do_kriged,asp=1,breaks=o_breaks,col=o_cols)
contour(bot_do_kriged,levels=o_breaks,add=T)
# image(bot_t_se,breaks=c(1.2,10),col=alpha('gray50',.5),add=T)
contour(topo_lon,
        topo_lat,
        topo,
        add=T,levels=c(-100),col='orange')
plot(world,add=T,col='gray80')
points(bot_plot$lon,bot_plot$lat,pch=20,col='orange')
mtext(expression(paste('Dissolved Oxygen (mg l'^-1,')')),adj=1)

# imagePlot(dd_dz1_kriged,asp=1,breaks=d1_breaks,col=d1_cols)
# contour(dd_dz1_kriged,levels=d1_breaks,add=T)
imagePlot(dd_dz2_kriged,asp=1,breaks=d2_breaks,col=d2_cols)
contour(dd_dz2_kriged,levels=d2_breaks,add=T)
# image(bot_t_se,breaks=c(1.2,10),col=alpha('gray50',.5),add=T)
contour(topo_lon,
        topo_lat,
        topo,
        add=T,levels=c(-100),col='green')
plot(world,add=T,col='gray80')
points(bot_plot$lon,bot_plot$lat,pch=20,col='green')
mtext(expression(paste('Stratification (d',rho,' dz'^-1,')')),adj=1)
mtext(expression(paste('Longitude (',degree,'W)')),1,2.5)
mtext(expression(paste('Latitude (',degree,'N)')),2,2)

imagePlot(mld2_kriged,asp=1,breaks=m_breaks,col=m_cols)
contour(mld2_kriged,levels=m_breaks,add=T)
# image(bot_t_se,breaks=c(1.2,10),col=alpha('gray50',.5),add=T)
contour(topo_lon,
        topo_lat,
        topo,
        add=T,levels=c(-100),col='orange')
plot(world,add=T,col='gray80')
points(bot_plot$lon,bot_plot$lat,pch=20,col='orange')
mtext('Scaled MLD',adj=1)
mtext(expression(paste('Longitude (',degree,'W)')),1,2.5)
# mtext(expression(paste('Latitude (',degree,'N)')),2,2)
mtext(paste('Note: Data are early release and subject to further QA/QC, \nplease contact brendan.turley@noaa.gov for comments/concerns \nUpdated: ',as.Date(Sys.time())),
      outer=T,line=3,side=1,col='red',font=2,at=.05,adj=0)
dev.off()


# my.krig <- spatialProcess(locs,bot_plot$bot_t,
#                           Distance = "rdist.earth")
# temp_kriged <- predictSurface(my.krig,
#                               grid.list = list(lon=seq(min(bot_plot$lon,na.rm=T), max(bot_plot$lon,na.rm=T),.01),
#                                                lat=seq(min(bot_plot$lat,na.rm=T), max(bot_plot$lat,na.rm=T),.01)))
# my.krig <- spatialProcess(locs,bot_plot$bot_s,
#                           Distance = "rdist.earth")
# sal_kriged <- predictSurface(my.krig,
#                             list(lon=seq(min(bot_plot$lon,na.rm=T), max(bot_plot$lon,na.rm=T),.01),
#                                  lat=seq(min(bot_plot$lat,na.rm=T), max(bot_plot$lat,na.rm=T),.01)))
# my.krig <- spatialProcess(locs,bot_plot$bot_c,
#                           Distance = "rdist.earth")
# chl_kriged <- predictSurface(my.krig,
#                             list(lon=seq(min(bot_plot$lon,na.rm=T), max(bot_plot$lon,na.rm=T),.01),
#                                  lat=seq(min(bot_plot$lat,na.rm=T), max(bot_plot$lat,na.rm=T),.01)))
# my.krig <- spatialProcess(locs,bot_plot$bot_do,
#                           Distance = "rdist.earth")
# do_kriged <- predictSurface(my.krig,
#                             list(lon=seq(min(bot_plot$lon,na.rm=T), max(bot_plot$lon,na.rm=T),.01),
#                                  lat=seq(min(bot_plot$lat,na.rm=T), max(bot_plot$lat,na.rm=T),.01)))
# 
# my.krig <- spatialProcess(locs,bot_plot$dd_dz1,
#                           Distance = "rdist.earth")
# dddz1_kriged <- predictSurface(my.krig,
#                             list(lon=seq(min(bot_plot$lon,na.rm=T), max(bot_plot$lon,na.rm=T),.01),
#                                  lat=seq(min(bot_plot$lat,na.rm=T), max(bot_plot$lat,na.rm=T),.01)))
# 
# my.krig <- spatialProcess(locs,bot_plot$dd_dz2,
#                           Distance = "rdist.earth")
# dddz2_kriged <- predictSurface(my.krig,
#                                grid.list = list(lon=seq(min(bot_plot$lon,na.rm=T), max(bot_plot$lon,na.rm=T),.01),
#                                                 lat=seq(min(bot_plot$lat,na.rm=T), max(bot_plot$lat,na.rm=T),.01)))
# 
# my.krig <- spatialProcess(locs,bot_plot$mld2,
#                           Distance = "rdist.earth")
# mld_kriged <- predictSurface(my.krig,
#                                list(lon=seq(min(bot_plot$lon,na.rm=T), max(bot_plot$lon,na.rm=T),.01),
#                                     lat=seq(min(bot_plot$lat,na.rm=T), max(bot_plot$lat,na.rm=T),.01)))
# mld_kriged$z[which(mld_kriged$z<0)] <- 0
# mld_kriged$z[which(mld_kriged$z>1)] <- 1
# 
# bot_plot$mld[which(is.na(bot_plot$mld))] <- bot_plot$m_depth[which(is.na(bot_plot$mld))]
# my.krig <- spatialProcess(locs,bot_plot$mld,
#                           Distance = "rdist.earth")
# mld1_kriged <- predictSurface(my.krig,
#                              list(lon=seq(min(bot_plot$lon,na.rm=T), max(bot_plot$lon,na.rm=T),.01),
#                                   lat=seq(min(bot_plot$lat,na.rm=T), max(bot_plot$lat,na.rm=T),.01)))

t_breaks <- breaks(temp_kriged$z)
t_cols <- t_col(length(t_breaks)-1)
s_breaks <- breaks(sal_kriged$z,.2)
s_cols <- s_col(length(s_breaks)-1)
c_breaks <- breaks(round(chl_kriged$z,2),.25)
c_cols <- c_col(length(c_breaks)-1)
d1_breaks <- breaks(round(dddz1_kriged$z,2),.02,decimal = T)
d1_cols <- d_col(length(d1_breaks)-1)
d2_breaks <- breaks(round(dddz2_kriged$z,2),.01,decimal = T)
d2_cols <- d_col(length(d2_breaks)-1)
m1_breaks <- breaks(mld1_kriged$z,2)
m1_cols <- m_col(length(m1_breaks)-1)
m_breaks <- seq(0,1,.1)
m_cols <- m_col(length(m_breaks)-1)


setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
png('new_plot.png',width=9,height=12,units='in',pointsize=12,res=300)
par(mfrow=c(3,2),mar=c(5,5,3,2),oma=c(5,0,1,0))
imagePlot(temp_kriged,asp=1,breaks=t_breaks,col=t_cols)
contour(temp_kriged,levels=t_breaks,add=T)
contour(topo_lon,
        topo_lat,
        topo,
        add=T,levels=c(-100),col='green')
plot(world,add=T,col='gray80')
points(bot_plot$lon,bot_plot$lat,pch=20,col='green')
mtext(expression(paste('Temperature (',degree,'C)')),adj=1)
# mtext(expression(paste('Longitude (',degree,'W)')),1,2.5)
mtext(expression(paste('Latitude (',degree,'N)')),2,2)

mtext(bot_plot$ship[1],adj=0,line=1)
mtext(paste(paste(month.abb[month(bot_plot$date_utc[1])],day(bot_plot$date_utc[1])),
      paste(month.abb[month(bot_plot$date_utc[nrow(bot_plot)])],day(bot_plot$date_utc[nrow(bot_plot)])),
      sep=' - '),
      adj=0)

imagePlot(sal_kriged,asp=1,breaks=s_breaks,col=s_cols)
contour(sal_kriged,levels=s_breaks,add=T)
contour(topo_lon,
        topo_lat,
        topo,
        add=T,levels=c(-100),col='purple')
plot(world,add=T,col='gray80')
points(bot_plot$lon,bot_plot$lat,pch=20,col='purple')
mtext('Salinity',adj=1)

imagePlot(chl_kriged,asp=1,breaks=c_breaks,col=c_cols)
contour(chl_kriged,levels=c_breaks,add=T)
contour(topo_lon,
        topo_lat,
        topo,
        add=T,levels=c(-100),col='purple')
plot(world,add=T,col='gray80')
points(bot_plot$lon,bot_plot$lat,pch=20,col='purple')
mtext(expression(paste('Chlorophyll (mg m'^-3,')')),adj=1)
# mtext(expression(paste('Longitude (',degree,'W)')),1,2.5)
mtext(expression(paste('Latitude (',degree,'N)')),2,2)

imagePlot(do_kriged,asp=1,breaks=o_breaks,col=o_cols)
contour(do_kriged,levels=o_breaks,add=T)
contour(topo_lon,
        topo_lat,
        topo,
        add=T,levels=c(-100),col='orange')
plot(world,add=T,col='gray80')
points(bot_plot$lon,bot_plot$lat,pch=20,col='orange')
mtext(expression(paste('Dissolved Oxygen (mg l'^-1,')')),adj=1)

# imagePlot(dddz1_kriged,asp=1,breaks=d1_breaks,col=d1_cols)
# contour(dddz1_kriged,levels=d1_breaks,add=T)
imagePlot(dddz2_kriged,asp=1,breaks=d2_breaks,col=d2_cols)
contour(dddz2_kriged,levels=d2_breaks,add=T)
contour(topo_lon,
        topo_lat,
        topo,
        add=T,levels=c(-100),col='green')
plot(world,add=T,col='gray80')
points(bot_plot$lon,bot_plot$lat,pch=20,col='green')
mtext(expression(paste('Stratification (d',rho,' dz'^-1,')')),adj=1)
mtext(expression(paste('Longitude (',degree,'W)')),1,2.5)
mtext(expression(paste('Latitude (',degree,'N)')),2,2)

imagePlot(mld_kriged,asp=1,breaks=m_breaks,col=m_cols)
contour(mld_kriged,levels=m_breaks,add=T)
contour(topo_lon,
        topo_lat,
        topo,
        add=T,levels=c(-100),col='orange')
plot(world,add=T,col='gray80')
points(bot_plot$lon,bot_plot$lat,pch=20,col='orange')
mtext('Scaled MLD',adj=1)
mtext(expression(paste('Longitude (',degree,'W)')),1,2.5)
# mtext(expression(paste('Latitude (',degree,'N)')),2,2)
mtext(paste('Note: Data are early release and subject to further QA/QC, \nplease contact brendan.turley@noaa.gov for comments/concerns \nUpdated: ',as.Date(Sys.time())),
      outer=T,line=3,side=1,col='red',font=2,at=.05,adj=0)
dev.off()


quilt.plot(bot_plot$lon,bot_plot$lat,bot_plot$surf_c,asp=1,breaks=c_breaks,col=c_cols)
quilt.plot(bot_plot$lon,bot_plot$lat,bot_plot$bot_t,asp=1,breaks=t_breaks,col=t_cols)
quilt.plot(bot_plot$lon,bot_plot$lat,bot_plot$bot_s,asp=1,breaks=s_breaks,col=s_cols)
quilt.plot(bot_plot$lon,bot_plot$lat,bot_plot$bot_c,asp=1,breaks=c_breaks,col=c_cols)
quilt.plot(bot_plot$lon,bot_plot$lat,bot_plot$bot_do,asp=1,breaks=o_breaks,col=o_cols)
quilt.plot(bot_plot$lon,bot_plot$lat,bot_plot$dd_dz1,asp=1,breaks=d1_breaks,col=d1_cols)
quilt.plot(bot_plot$lon,bot_plot$lat,bot_plot$dd_dz2,asp=1,breaks=d2_breaks,col=d2_cols)
quilt.plot(bot_plot$lon,bot_plot$lat,bot_plot$mld,asp=1,breaks=m1_breaks,col=m1_cols)
quilt.plot(bot_plot$lon,bot_plot$lat,bot_plot$mld2,asp=1,breaks=m_breaks,col=m_cols)

bubblePlot(bot_plot$lon,bot_plot$lat,bot_plot$bot_do,asp=1,breaks=o_breaks,col=o_cols)
bubblePlot(bot_plot$lon,bot_plot$lat,bot_plot$mld2,asp=1,breaks=m_breaks,col=m_cols)

temp <- quilt.plot(bot_plot$lon,bot_plot$lat,bot_plot$mld2,plot=F)
obj <- list(x=temp$x,y=temp$y,z=temp$z)
loc <- expand.grid(list(lon=seq(min(bot_plot$lon,na.rm=T), max(bot_plot$lon,na.rm=T),.01),
                 lat=seq(min(bot_plot$lat,na.rm=T), max(bot_plot$lat,na.rm=T),.01)))
z.new<- interp.surface(obj,loc)
image.plot(obj,asp=1,breaks=m_breaks,col=m_cols)
quilt.plot(loc,z.new)
image.plot(z.new)

grid.list <- list(x=seq(min(bot_plot$lon,na.rm=T), max(bot_plot$lon,na.rm=T),.01),
                        y=seq(min(bot_plot$lat,na.rm=T), max(bot_plot$lat,na.rm=T),.01))
z.new<- interp.surface.grid(obj,grid.list)

image(z.new$z)

as.image(bot_plot$mld2)
### for FCWC
library(NISTunits)
temp_kriged$z <- NISTdegCtOdegF(temp_kriged$z)
topo <- NISTmeterTOft(topo)
mld1_kriged$z <- NISTmeterTOft(mld1_kriged$z)

t_breaks <- breaks(temp_kriged$z,2)
if(max(temp_kriged$z,na.rm=T)>t_breaks[length(t_breaks)]){
  t_breaks <- c(t_breaks,t_breaks[length(t_breaks)]+2)
}
t_cols <- t_col(length(t_breaks)-1)
m1_breaks <- breaks(mld1_kriged$z,5)
if(max(mld1_kriged$z,na.rm=T)>m1_breaks[length(m1_breaks)]){
  m1_breaks <- c(m1_breaks,m1_breaks[length(m1_breaks)]+5)
}
m1_cols <- m_col(length(m1_breaks)-1)

setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
png('fcwc_plot.png',width=9,height=12,units='in',pointsize=12,res=300)
par(mfrow=c(3,2),mar=c(5,5,3,2),oma=c(5,0,1,0))
imagePlot(temp_kriged,asp=1,breaks=t_breaks,col=t_cols)
contour(temp_kriged,levels=t_breaks,add=T)
contour(topo_lon,
        topo_lat,
        topo,
        add=T,levels=c(-400,-300,-200,-100,-50),col='green')
plot(world,add=T,col='gray80')
points(bot_plot$lon,bot_plot$lat,pch=20,col='green')
mtext(expression(paste('Temperature (',degree,'F)')),adj=1)
# mtext(expression(paste('Longitude (',degree,'W)')),1,2.5)
mtext(expression(paste('Latitude (',degree,'N)')),2,2)

mtext(bot_plot$ship[1],adj=0,line=1)
mtext(paste(paste(month.abb[month(bot_plot$date_utc[1])],day(bot_plot$date_utc[1])),
            paste(month.abb[month(bot_plot$date_utc[nrow(bot_plot)])],day(bot_plot$date_utc[nrow(bot_plot)])),
            sep=' - '),
      adj=0)

imagePlot(sal_kriged,asp=1,breaks=s_breaks,col=s_cols)
contour(sal_kriged,levels=s_breaks,add=T)
contour(topo_lon,
        topo_lat,
        topo,
        add=T,levels=c(-400,-300,-200,-100,-50),col='purple')
plot(world,add=T,col='gray80')
points(bot_plot$lon,bot_plot$lat,pch=20,col='purple')
mtext('Salinity',adj=1)

imagePlot(chl_kriged,asp=1,breaks=c_breaks,col=c_cols)
contour(chl_kriged,levels=c_breaks,add=T)
contour(topo_lon,
        topo_lat,
        topo,
        add=T,levels=c(-400,-300,-200,-100,-50),col='purple')
plot(world,add=T,col='gray80')
points(bot_plot$lon,bot_plot$lat,pch=20,col='purple')
mtext(expression(paste('Chlorophyll (mg m'^-3,')')),adj=1)
# mtext(expression(paste('Longitude (',degree,'W)')),1,2.5)
mtext(expression(paste('Latitude (',degree,'N)')),2,2)

imagePlot(do_kriged,asp=1,breaks=o_breaks,col=o_cols)
contour(do_kriged,levels=o_breaks,add=T)
contour(topo_lon,
        topo_lat,
        topo,
        add=T,levels=c(-400,-300,-200,-100,-50),col='orange')
plot(world,add=T,col='gray80')
points(bot_plot$lon,bot_plot$lat,pch=20,col='orange')
mtext(expression(paste('Dissolved Oxygen (mg l'^-1,')')),adj=1)

imagePlot(dddz1_kriged,asp=1,breaks=d1_breaks,col=d1_cols)
contour(dddz1_kriged,levels=d1_breaks,add=T)
# imagePlot(dddz2_kriged,asp=1,breaks=d2_breaks,col=d2_cols)
# contour(dddz2_kriged,levels=d2_breaks,add=T)
contour(topo_lon,
        topo_lat,
        topo,
        add=T,levels=c(-400,-300,-200,-100,-50),col='green')
plot(world,add=T,col='gray80')
points(bot_plot$lon,bot_plot$lat,pch=20,col='green')
mtext(expression(paste('Stratification (d',rho,' dz'^-1,')')),adj=1)
mtext(expression(paste('Longitude (',degree,'W)')),1,2.5)
mtext(expression(paste('Latitude (',degree,'N)')),2,2)

imagePlot(mld1_kriged,asp=1,breaks=m1_breaks,col=m1_cols)
contour(mld1_kriged,levels=m1_breaks,add=T)
contour(topo_lon,
        topo_lat,
        topo,
        add=T,levels=c(-400,-300,-200,-100,-50),col='orange')
plot(world,add=T,col='gray80')
points(bot_plot$lon,bot_plot$lat,pch=20,col='orange')
mtext('MLD (ft)',adj=1)
mtext(expression(paste('Longitude (',degree,'W)')),1,2.5)
# mtext(expression(paste('Latitude (',degree,'N)')),2,2)
mtext(paste('Note: Data are early release and subject to further QA/QC, \nplease contact brendan.turley@noaa.gov for comments/concerns \nUpdated: ',as.Date(Sys.time())),
      outer=T,line=3,side=1,col='red',font=2,at=.05,adj=0)
dev.off()

setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
png('fcwc_plot_temp.png',width=9,height=9,units='in',pointsize=12,res=300)
par(mfrow=c(1,1),mar=c(5,5,3,2),oma=c(5,0,1,0))
imagePlot(temp_kriged,asp=1,breaks=t_breaks,col=t_cols)
contour(temp_kriged,levels=t_breaks,add=T)
contour(topo_lon,
        topo_lat,
        topo,
        add=T,levels=c(-400,-300,-200,-100,-50),col='green')
plot(world,add=T,col='gray80')
points(bot_plot$lon,bot_plot$lat,pch=20,col='green')
mtext(expression(paste('Temperature (',degree,'F)')),adj=1)
mtext(expression(paste('Longitude (',degree,'W)')),1,2.5)
mtext(expression(paste('Latitude (',degree,'N)')),2,2)

mtext(bot_plot$ship[1],adj=0,line=1)
mtext(paste(paste(month.abb[month(bot_plot$date_utc[1])],day(bot_plot$date_utc[1])),
            paste(month.abb[month(bot_plot$date_utc[nrow(bot_plot)])],day(bot_plot$date_utc[nrow(bot_plot)])),
            sep=' - '),
      adj=0)

mtext(paste('Note: Data are early release and subject to further QA/QC, \nplease contact brendan.turley@noaa.gov for comments/concerns \nUpdated: ',as.Date(Sys.time())),
      outer=T,line=3,side=1,col='red',font=2,at=.05,adj=0)
dev.off()

setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
png('fcwc_plot_temp2.png',width=8,height=9,units='in',pointsize=12,res=300)
par(mfrow=c(1,1),mar=c(5,5,3,2),oma=c(5,0,1,0))
# imagePlot(temp_kriged,asp=1,breaks=c(0,1),col='white')
contour(temp_kriged,levels=t_breaks,asp=1)
contour(topo_lon,
        topo_lat,
        topo,
        add=T,levels=c(-400,-300,-200,-100,-50),col='gray80')
plot(world,add=T,col='gray80')
points(bot_plot$lon,bot_plot$lat,pch=20,col='gray80')
mtext(expression(paste('Temperature (',degree,'F)')),adj=1)
mtext(expression(paste('Longitude (',degree,'W)')),1,2.5)
mtext(expression(paste('Latitude (',degree,'N)')),2,2)

mtext(bot_plot$ship[1],adj=0,line=1)
mtext(paste(paste(month.abb[month(bot_plot$date_utc[1])],day(bot_plot$date_utc[1])),
            paste(month.abb[month(bot_plot$date_utc[nrow(bot_plot)])],day(bot_plot$date_utc[nrow(bot_plot)])),
            sep=' - '),
      adj=0)

mtext(paste('Note: Data are early release and subject to further QA/QC, \nplease contact brendan.turley@noaa.gov for comments/concerns \nUpdated: ',as.Date(Sys.time())),
      outer=T,line=3,side=1,col='red',font=2,at=.05,adj=0)
dev.off()



library(raster)
library(sp)
r <-raster(
  t(mld_kriged$z[,ncol(mld_kriged$z):1]),
  xmn=range(mld_kriged$x)[1], xmx=range(mld_kriged$x)[2],
  ymn=range(mld_kriged$y)[1], ymx=range(mld_kriged$y)[2], 
  crs=CRS("+proj=longlat +datum=WGS84")
)
plot(r)
KML(r, file='mld.kml',col=m_cols)


coordinates(bot_plot) <- ~lon + lat
proj4string(bot_plot)<-CRS("+proj=longlat +datum=WGS84")
writeOGR(bot_plot, dsn="mld2.kml", layer= "mld", driver="KML")



