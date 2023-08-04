rm(list=ls())
gc()

library(fields)
library(lubridate)
library(ncdf4)
library(terra)

### bathymetric data
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
setwd("~/Desktop/professional/biblioteca/data/shapefiles/gshhg-shp-2.3.7/GSHHS_shp/h/")
world <- vect('GSHHS_h_L1.shp')
world <- crop(world, ext(-86, -79, 24.5, 28))


### load data
setwd('~/Desktop/professional/projects/Postdoc_FL/data/walton_smith/')
setwd('~/Downloads')
data <- read.csv('WS23203_SampleLog.csv')
### cruise name for file naming
cruise <- 'WS23203'
### only stations at depth
ind <- which(data$Depth!=0)
data2 <- data[ind,]
st_rm <- c('2','3','6.5','MR','9','9.5','10','12','16','18','21/LK','EK MID','EK OFF','WS','KW1','KW2')
data3 <- data2[!is.element(data2$Station,st_rm),]
# data3$Date.GMT <- mdy(data3$Date.GMT)
data3$date <- mdy_hms(paste(data3$Date.GMT,data3$Time..GMT.))


################## geographic scope for hycom
lonbox_e <- -80.5 ### Florida Bay
lonbox_e <- (lonbox_e + 360)
lonbox_w <- -86 ### mouth of Mississippi River
lonbox_w <- (lonbox_w + 360)
latbox_n <- 30.5 ### northern coast
latbox_s <- 24.5 ### remove the Keys

# https://tds.hycom.org/thredds/catalogs/GLBy0.08/expt_93.0.html
### now
url <- 'https://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0'
data <- nc_open(url)

time <- ncvar_get(data,'time')
time2 <- as.Date(time/24,origin='2000-01-01 00:00:00')
time3 <- as.POSIXct(time*3600,origin='2000-01-01',tz='GMT')
# time2 <- as.Date(time/24,origin='2022-02-25 12:00:00')
# time3 <- as.POSIXct(time*3600,origin='2022-02-25 12:00:00',tz='GMT')

lat <- ncvar_get(data,'lat')
# ind_lat <- which(lat>=latbox_s & lat<=latbox_n)
lon <- ncvar_get(data,'lon')
# ind_lon <- which(lon>=lonbox_w & lon<=lonbox_e)

z <- ncvar_get(data,'depth')


ind_time <- sapply(data3$date,function(x) which.min(abs(x-time3)))
ind_lon <- sapply(data3$Longitude.Decimal+360,function(x) which.min(abs(x-lon)))
ind_lat <- sapply(data3$Latitude.Decimal,function(x) which.min(abs(x-lat)))

hycom_dat <- data.frame(bot_temp=rep(NA,length(ind_time)),
                        bot_sal=rep(NA,length(ind_time)))
for(i in 1:length(ind_time)){
  hycom_dat[i,1] <- ncvar_get(data,'water_temp_bottom',
                              start=c(ind_lon[i],ind_lat[i],ind_time[i]),
                              count=c(1,1,1))
  
  hycom_dat[i,2] <- ncvar_get(data,'salinity_bottom',
                              start=c(ind_lon[i],ind_lat[i],ind_time[i]),
                              count=c(1,1,1))
}

t_rsq <- round(summary(lm(hycom_dat$bot_temp~data3$Temperature.CTD.data))$adj.r.squared,2)
s_rsq <- round(summary(lm(hycom_dat$bot_sal~data3$Salinity.CTD.data))$adj.r.squared,2)

setwd('~/Documents/R/Github/waltonsmith/figures')
png(paste0(cruise,'_bottom_bias2.png'), height = 8, width = 4, units = 'in', res=300)
par(mfrow=c(2,1))
plot(data3$Temperature.CTD.data,hycom_dat$bot_temp,
     xlab='Walton Smith',ylab='HYCOM',asp=1,las=1)
mtext('Bottom temperature bias',adj=0)
abline(0,1,lty=1,col=2)
mtext(bquote(paste(R^2, ' = ', .(t_rsq))),
      adj=1)

plot(data3$Salinity.CTD.data,hycom_dat$bot_sal,
     xlab='Walton Smith',ylab='HYCOM',asp=1,las=1)
mtext('Bottom salinity bias',adj=0)
abline(0,1,lty=1,col=2)
mtext(bquote(paste(R^2, ' = ', .(s_rsq))),
      adj=1)
dev.off()

quant <- .99
strat_n_col <- colorRampPalette(c('dodgerblue4','deepskyblue3','lightskyblue1','gray90'))
strat_p_col <- colorRampPalette(c('gray90','rosybrown1','tomato2','red4'))

resid_t <- hycom_dat$bot_temp-data3$Temperature.CTD.data
# tr_breaks <- pretty(resid_t[which(resid_t<=quantile(resid_t,quant,na.rm=T))],n=20)
tr_breaks <- pretty(resid_t,n=20)
tr <- cut(resid_t,tr_breaks)
if(any(tr_breaks==0)){
  tr_breaks <- pretty(resid_t,n=20)
  tr_cols <- c(strat_n_col(length(which(tr_breaks<0))),
                   strat_p_col(length(which(tr_breaks>0))))
}else{
  tr_cols <- (strat_p_col(length(tr_breaks)-1))  
}

resid_s <- hycom_dat$bot_sal-data3$Salinity.CTD.data
# sr_breaks <- pretty(resid_s[which(resid_s<=quantile(resid_s,quant,na.rm=T))],n=20)
sr_breaks <- pretty(resid_s,n=20)
sr <- cut(resid_s,sr_breaks)
if(any(sr_breaks==0)){
  sr_breaks <- pretty(resid_s,n=20)
  sr_cols <- c(strat_n_col(length(which(sr_breaks<0))),
                   strat_p_col(length(which(sr_breaks>0))))
}else{
  sr_cols <- (strat_p_col(length(sr_breaks)-1))  
}


hist(resid_t)
hist(resid_s)

setwd('~/Documents/R/Github/waltonsmith/figures')
png(paste0(cruise,'_bottom_bias.png'), height = 12, width = 6, units = 'in', res=300)
par(mfrow=c(2,1))
info<- setupLegend()
plot(data3$Longitude.Decimal,data3$Latitude.Decimal,
     bg=tr_cols[as.numeric(tr)],pch=21,asp=1,cex=1.5,
     xlab='Longitude',ylab='Latitude',las=1)
plot(world,col='gray70',add=T)
contour(topo_lon,topo_lat,topo,add=T,levels=c(-100,-50,-25,-10),col='gray40')
addLegend(info,col=tr_cols,zlim=range(resid_t,na.rm=T))
mtext('Bottom temperature bias (HYCOM - observations)',adj=1)

info<- setupLegend()
plot(data3$Longitude.Decimal,data3$Latitude.Decimal,
     bg=sr_cols[as.numeric(sr)],pch=21,asp=1,cex=1.5,
     xlab='Longitude',ylab='Latitude',las=1)
plot(world,col='gray70',add=T)
contour(topo_lon,topo_lat,topo,add=T,levels=c(-100,-50,-25,-10),col='gray40')
addLegend(info,col=sr_cols,zlim=range(resid_s,na.rm=T))
mtext('Bottom salinity bias (HYCOM - observations)',adj=1)
dev.off()
