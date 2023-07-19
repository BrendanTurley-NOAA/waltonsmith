library(fields)
library(lubridate)
library(cmocean)
library(rerddap)
library(sf)
library(terra)

erddap_extract_oscar <- function(data, info, time){
  data_temp <- data$data
  lon <- data$summary$dim$longitude$vals
  lat <- data$summary$dim$latitude$vals
  ind <- which(data_temp$time==unique(data_temp$time)[time])
  u <- matrix(data_temp[ind,4], ### modified to extract correct column
              length(lon),
              length(lat))
  v <- matrix(data_temp[ind,5], ### modified to extract correct column
              length(lon),
              length(lat))
  month <- month(data_temp$time[ind][1])
  row_ind <- which(info$alldata$NC_GLOBAL$attribute_name=='title')
  col_ind <- which(colnames(info$alldata$NC_GLOBAL)=='value')
  name <- info$alldata$NC_GLOBAL[row_ind,col_ind]
  name <- unlist(strsplit(name,split=','))
  return(list(u = u,
              v = v,
              lon = lon,
              lat = lat,
              month = month,
              name = name))
  # setClass('erddap',slots=c(data='matrix',lon='array',lat='array'))
  # return(new('erddap',data=new_data,lon=lon,lat=lat))
  # return(list(new_data,lon,lat))
}

setwd("~/Desktop/professional/biblioteca/data/shapefiles/gshhg-shp-2.3.7/GSHHS_shp/h/")
world <- vect('GSHHS_h_L1.shp')

setwd("~/Desktop/professional/biblioteca/data/shapefiles/fknms_py2")
fknms <- vect('fknms_py.shp')
crs(fknms) <- crs(world)
fknms <- st_as_sf(fknms)

setwd('~/Desktop/professional/projects/Postdoc_FL/data/walton_smith/underway_data')
fil <- list.files()

# cruise "WS17177_out.csv" did not stay out very long, will discard
fil <- fil[-which(fil=="WS17177_out.csv")]

i=1
out_all <- out_fknms <- data.frame(matrix(NA,length(fil),8))
for(i in 1:length(fil)){
  data <- read.csv(fil[i])
  data <- data[which(data$sog>2 & data$sog<15),]
  if(any(is.na(data$lon.dd))){
    data <- data[-which(is.na(data$lon.dd)),]
  }
  point.sf <- st_as_sf(data, coords = c("lon.dd","lat.dd"))
  st_crs(point.sf) <- st_crs(fknms)
  # data <- st_filter(point.sf, fknms)
  mat <- st_intersects(point.sf, fknms)
  ind_fknms <- apply(mat, 1, any)
  data_fknms <- data[ind_fknms,]
  
  if(any(names(data)=='temp.c3p')){
    ind <- which(names(data)=='temp.c3p')
  } else {
    ind <- which(names(data)=='temp.tsg')
  }
  
  temp <- data[,ind]
  # st_geometry(temp) <- NULL
  temp[which(temp<0)] <- NA
  
  temp2 <- data_fknms[,ind]
  temp2[which(temp2<0)] <- NA
  
  out_all[i,] <- c(substr(fil[i],1,7),data$time[1],data$time[nrow(data)],quantile(temp,na.rm=T))
  out_fknms[i,] <- c(substr(fil[i],1,7),data_fknms$time[1],data_fknms$time[nrow(data_fknms)],quantile(temp2,na.rm=T))
  
  par(mfrow=c(1,2))
  hist(temp)
  hist(temp2)
  
  # plot(data$lon.dd,data$lat.dd,
  #      pch=16,col=cmocean('curl')(nrow(data)),
  #      asp=1)
}

dats <- dmy_hms(out_all$X2)
dats[is.na(dats)] <- mdy_hm(out_all$X2[is.na(dats)])

table(month(dats))
table(year(dats),month(dats))

ind <- which(month(dats)>5 & month(dats)<9)

brks <- seq(28.5,33.5,.25)
hists <- data.frame(matrix(NA,length(ind),length(brks)-1))
tmps <- list()
png('WS_underway_summer_clim_maps.png',width = 9, height = 9, res=300, units='in')
par(mar=c(2,2,2,2),mfrow=c(2,2))
for(i in ind){
  setwd('~/Desktop/professional/projects/Postdoc_FL/data/walton_smith/underway_data')
  data <- read.csv(fil[i])
  data <- data[which(data$sog>2 & data$sog<15),]
  if(any(is.na(data$lon.dd))){
    data <- data[-which(is.na(data$lon.dd)),]
  }
  
  if(any(names(data)=='temp.c3p')){
    ind_t <- which(names(data)=='temp.c3p')
  } else {
    ind_t <- which(names(data)=='temp.tsg')
  }
  temp <- data[,ind_t]
  # st_geometry(temp) <- NULL
  temp[which(temp<0)] <- NA
  tmps[[which(i==ind)]] <- temp
  
  h <- hist(temp,breaks=brks,plot=F)
  hists[which(i==ind),] <- h$density
  
  tmp_cut <- cut(temp,brks)
  cols <- cmocean('thermal')(length(levels(tmp_cut)))
  
  plot(data$lon.dd,data$lat.dd,
       pch=16,col=cols[as.numeric(tmp_cut)],
       asp=1,xlab='',ylab='')
  plot(world,add=T)
  mtext(paste0(out_all[i,1],", ",year(dats)[i],'-',month.abb[month(dats)[i]]))
  
}
setwd('~/Desktop/professional/projects/Postdoc_FL/figures/')
dev.off()
tmps <- unlist(tmps)
h2 <- hist(tmps,breaks=brks)

cols2 <- cmocean('thermal')(length(brks)-1)
at <- seq(0,15,1.5)

setwd('~/Desktop/professional/projects/Postdoc_FL/figures/')
png('WS_underway_summer_clim.png',width = 8, height = 9, res=300, units='in')
par(mar=c(5,5,2,2),mfrow=c(1,1))
plot(0,0,xlim=c(28.5,33.5),ylim=c(-.25,(nrow(hists)+1)*1.4),yaxt='n',xlab='Underway Temperature (C)',ylab='')
axis(2,at[1:length(ind)],paste0(year(dats)[ind],'-',month.abb[month(dats)[ind]]),las=2)
for(i in 1:nrow(hists)){
  rect(h$mids-.125,rep(at[i],15),h$mids+.125,hists[i,]+at[i],col=cols2)
  arrows(as.numeric(out_all[ind[i],4]),rep(at[i]-.125,1),as.numeric(out_all[ind[i],8]),rep(at[i]-.125,1),
         col=1,code=3,angle=90,length=.05,lend=2,lwd=2)
  rect(as.numeric(out_all[ind[i],5]),rep(at[i],15),as.numeric(out_all[ind[i],7]),rep(at[i]-.25,15),col='gray60')
  segments(as.numeric(out_all[ind[i],6]),rep(at[i]-.03,15),as.numeric(out_all[ind[i],6]),rep(at[i]-.21,15),col='gold',lwd=2,lend=2)
  
}
rect(h$mids-.125,rep(at[length(ind)+1],15),h$mids+.125,h2$density+at[length(ind)+1],col=cols2)
arrows(min(tmps,na.rm=T),rep(at[length(ind)+1]-.125,1),max(tmps,na.rm=T),rep(at[length(ind)+1]-.125,1),
       col=1,code=3,angle=90,length=.05,lend=2,lwd=2)
rect(quantile(tmps,.25,na.rm=T),rep(at[length(ind)+1],15),quantile(tmps,.75,na.rm=T),rep(at[length(ind)+1]-.25,15),col='gray60')
segments(median(tmps,na.rm=T),rep(at[length(ind)+1]-.03,15),median(tmps,na.rm=T),rep(at[length(ind)+1]-.21,15),col='gold',lwd=2,lend=2)
axis(2,at[length(ind)+1],'Overall',las=2)
dev.off()



brks <- seq(29.75,33.5,.25)
hists <- data.frame(matrix(NA,length(ind),length(brks)-1))
tmps <- list()
png('WS_underway_summer_clim_maps_fk.png',width = 9, height = 9, res=300, units='in')
par(mar=c(2,2,2,2),mfrow=c(2,2))
for(i in ind){
  setwd('~/Desktop/professional/projects/Postdoc_FL/data/walton_smith/underway_data')
  data <- read.csv(fil[i])
  data <- data[which(data$sog>2 & data$sog<15),]
  if(any(is.na(data$lon.dd))){
    data <- data[-which(is.na(data$lon.dd)),]
  }
  point.sf <- st_as_sf(data, coords = c("lon.dd","lat.dd"))
  st_crs(point.sf) <- st_crs(fknms)
  # data <- st_filter(point.sf, fknms)
  mat <- st_intersects(point.sf, fknms)
  ind_fknms <- apply(mat, 1, any)
  data_fknms <- data[ind_fknms,]
  
  if(any(names(data)=='temp.c3p')){
    ind_t <- which(names(data)=='temp.c3p')
  } else {
    ind_t <- which(names(data)=='temp.tsg')
  }
  temp <- data_fknms[,ind_t]
  # st_geometry(temp) <- NULL
  temp[which(temp<0)] <- NA
  tmps[[which(i==ind)]] <- temp

  h <- hist(temp,breaks=brks,plot=F)
  hists[which(i==ind),] <- h$density
  
  tmp_cut <- cut(temp,brks)
  cols <- cmocean('thermal')(length(levels(tmp_cut)))
  
  plot(data$lon.dd,data$lat.dd,
       typ='l',asp=1,xlab='',ylab='',lwd=2,col='orange3',
       xlim=range(data_fknms$lon.dd),ylim=range(data_fknms$lat.dd))
  plot(st_geometry(fknms),add=T,border='forestgreen',lwd=2)
  plot(world,add=T,col='gray70')
  points(data_fknms$lon.dd,data_fknms$lat.dd,
       pch=16,col=cols[as.numeric(tmp_cut)])
  mtext(paste0(out_fknms[i,1],", ",year(dats)[i],'-',month.abb[month(dats)[i]]))
  
  # plot(data_fknms$lon.dd,data_fknms$lat.dd,
  #      pch=16,col=cols[as.numeric(tmp_cut)],
  #      asp=1,xlab='',ylab='')
  # plot(world,add=T)
  # plot(st_geometry(fknms),add=T,border=3)
  # mtext(paste0(out[i,1],", ",year(dats)[i],'-',month.abb[month(dats)[i]]))
  
}
setwd('~/Desktop/professional/projects/Postdoc_FL/figures/')
dev.off()
tmps <- unlist(tmps)
h2 <- hist(tmps,breaks=brks)

cols2 <- cmocean('thermal')(length(brks)-1)
at <- seq(0,15,1.5)

setwd('~/Desktop/professional/projects/Postdoc_FL/figures/')
png('WS_underway_summer_clim_fk.png',width = 8, height = 9, res=300, units='in')
par(mar=c(5,5,2,2),mfrow=c(1,1))
plot(0,0,xlim=c(29.75,33.5),ylim=c(-.25,(nrow(hists)+1)*1.4),yaxt='n',xlab='Underway Temperature (C)',ylab='')
axis(2,at[1:length(ind)],paste0(year(dats)[ind],'-',month.abb[month(dats)[ind]]),las=2)
for(i in 1:nrow(hists)){
  rect(h$mids-.125,rep(at[i],15),h$mids+.125,hists[i,]+at[i],col=cols2)
  arrows(as.numeric(out_fknms[ind[i],4]),rep(at[i]-.125,1),as.numeric(out_fknms[ind[i],8]),rep(at[i]-.125,1),
         col=1,code=3,angle=90,length=.05,lend=2,lwd=2)
  rect(as.numeric(out_fknms[ind[i],5]),rep(at[i],15),as.numeric(out_fknms[ind[i],7]),rep(at[i]-.25,15),col='gray60')
  segments(as.numeric(out_fknms[ind[i],6]),rep(at[i]-.03,15),as.numeric(out_fknms[ind[i],6]),rep(at[i]-.21,15),col='gold',lwd=2,lend=2)
  
}
rect(h$mids-.125,rep(at[length(ind)+1],15),h$mids+.125,h2$density+at[length(ind)+1],col=cols2)
arrows(min(tmps,na.rm=T),rep(at[length(ind)+1]-.125,1),max(tmps,na.rm=T),rep(at[length(ind)+1]-.125,1),
       col=1,code=3,angle=90,length=.05,lend=2,lwd=2)
rect(quantile(tmps,.25,na.rm=T),rep(at[length(ind)+1],15),quantile(tmps,.75,na.rm=T),rep(at[length(ind)+1]-.25,15),col='gray60')
segments(median(tmps,na.rm=T),rep(at[length(ind)+1]-.03,15),median(tmps,na.rm=T),rep(at[length(ind)+1]-.21,15),col='gold',lwd=2,lend=2)
axis(2,at[length(ind)+1],'Overall',las=2)
dev.off()


setwd('~/Desktop/professional/projects/Postdoc_FL/figures/')
png('WS_underway_summer_uv_map.png',width = 9, height = 9, res=300, units='in')
par(mar=c(2,2,2,2),mfrow=c(2,2))
for(i in ind){
  time <- as.character(c(dmy(substr(out[i,2],1,11))-14,dmy(substr(out[i,3],1,11))))
  
  currents <- info('miamicurrents') #SSH derived
  ### west Florida Shelf
  lonbox_e <- -79 ### Florida Bay
  lonbox_w <- -90 ### mouth of Mississippi River
  latbox_n <- 30 ### northern coast
  latbox_s <- 23 ### southern edge of Key West
  latitude = c(latbox_s, latbox_n)
  longitude = c(lonbox_w, lonbox_e)
  
  uv_grab <- griddap(currents, latitude=latitude, longitude=longitude, time=time, fields=c('u_current','v_current'))
  lon_o <- sort(unique(uv_grab$data$lon))
  lat_o <- sort(unique(uv_grab$data$lat))
  
  n <- length(unique(uv_grab$data$time))
  u_array <- array(NA,dim=c(n,length(lon_o),length(lat_o)))
  v_array <- array(NA,dim=c(n,length(lon_o),length(lat_o)))
  for(i in 1:(n)){
    uv <- erddap_extract_oscar(uv_grab,currents,i) #SSH
    u_array[i,,] <- uv$u
    v_array[i,,] <- uv$v
  }
  lon_lat <- expand.grid(lon_o,lat_o)
  names(lon_lat) <- c('lon','lat')
  extend <- 1/4
  
  u <- apply(u_array,c(2,3),mean,na.rm=T)
  v <- apply(v_array,c(2,3),mean,na.rm=T)
  uv <- sqrt(u^2 + v^2)
  cols <- cmocean('speed')(20)
  
  imagePlot(lon_o,lat_o,
            uv,asp=1,breaks = seq(0,2,.1),col=cols)
  arrows(lon_lat$lon,lon_lat$lat, #RTOFs and SSH
         lon_lat$lon+(as.vector(u)*extend),lon_lat$lat+(as.vector(v)*extend),
         length=.025,col='gray10')
  plot(world,add=T)
  mtext(paste(time[1],'-',time[2]))
  
}
dev.off()





plot(rep(at[2],5),out[ind[2],c(4:8)],cex=2)




i=1
data <- read.csv(fil[i])
### make date times; check format b/c can change between cruises
data$time <- dmy_hms(data$time)
### copy to plot tracklines later
underway <- data
plot(data$lon.dd,data$lat.dd,
     pch=16,col=cmocean('curl')(nrow(data)),
     asp=1)


underway <- underway[which(underway$sog>2 & underway$sog<15),]

plot(underway$time,underway$temp.tsg,typ='l')
hist(underway$temp.tsg)

n <- 10
underway$temp.tsg <- filter(underway$temp.tsg,rep(1/n,n),'convolution')

tmp_cut <- cut(underway$temp.tsg,pretty(underway$temp.tsg,n=10))
cols <- cmocean('thermal')(length(levels(tmp_cut)))
plot(data$lon.dd,data$lat.dd,
     pch=16,col=cols[as.numeric(tmp_cut)],
     asp=1)
