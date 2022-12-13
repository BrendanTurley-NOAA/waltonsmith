rm(list=ls())

library(fields)
library(lubridate)
library(rerddap)
library(rgdal)
library(ncdf4)

erddap_extract <- function(data, info, parameter){
  data_temp <- data$data
  ind_extract <- which(names(data_temp)==parameter)
  time_step <- unique(data_temp$time)
  lon <- data$summary$dim$longitude$vals
  lat <- data$summary$dim$latitude$vals
  
  new_data <- array(data_temp[,ind_extract], 
                    c(length(lon),
                      length(lat),
                      length(time_step)))
  
  row_ind <- which(info$alldata$NC_GLOBAL$attribute_name=='title')
  col_ind <- which(colnames(info$alldata$NC_GLOBAL)=='value')
  name <- info$alldata$NC_GLOBAL[row_ind,col_ind]
  name <- unlist(strsplit(name,split=','))
  return(list(data = new_data,
              lon = lon,
              lat = lat,
              time = time_step,
              name = name))
  # setClass('erddap',slots=c(data='matrix',lon='array',lat='array'))
  # return(new('erddap',data=new_data,lon=lon,lat=lat))
  # return(list(new_data,lon,lat))
}

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

setwd("~/Desktop/professional/biblioteca/data")
bathy <- nc_open('etopo1.nc')
topo <- ncvar_get(bathy, 'Band1')
topo_lat <- ncvar_get(bathy, 'lat')
topo_lon <- ncvar_get(bathy, 'lon')
nc_close(bathy)

### load map
# setwd("C:/Users/brendan.turley/Desktop/FL_habs/ne_10m_admin_0_countries")
setwd("~/Desktop/professional/biblioteca/data/shapefiles/ne_10m_admin_0_countries")
world <- readOGR('ne_10m_admin_0_countries.shp')

# SST anomaly
cols_ssta_neg <- colorRampPalette(c('dodgerblue4','deepskyblue3','lightskyblue1','gray95'))
cols_ssta_pos <- colorRampPalette(c('gray95','rosybrown1','tomato2','red4'))
# breaks_ssta <- seq(-round(max(abs(sst_a_grab$data$sstAnom),na.rm=T),digits=1),
#                    round(max(abs(sst_a_grab$data$sstAnom),na.rm=T),digits=1),
#                    by=.2)


# mur_sst_a <- info('jplMURSST41anommday')
mur_sst <- info('jplMURSST41')
currents <- info('miamicurrents') #SSH derived
# ssha <- info('nesdisSSH1day')

### whole Gulf of Mexico
# latitude = c(17, 31)
# longitude = c(-98, -80)
### west Florida Shelf
lonbox_e <- -80.6 ### Florida Bay
lonbox_w <- -88 ### mouth of Mississippi River
latbox_n <- 30.5 ### northern coast
latbox_s <- 24.3 ### southern edge of Key West
latitude = c(latbox_s, latbox_n)
longitude = c(lonbox_w, lonbox_e)

time = c(paste(2022,"-12-04",sep=''), paste(2022,"-12-08",sep=''))
sst_grab <- griddap(mur_sst, latitude=latitude, longitude=longitude, time=time, fields='analysed_sst')
sst_1 <- erddap_extract(sst_grab,mur_sst,'analysed_sst')

uv_grab <- griddap(currents, latitude=latitude, longitude=longitude, time=time, fields=c('u_current','v_current'))
lon_o <- sort(unique(uv_grab$data$lon))
lat_o <- sort(unique(uv_grab$data$lat))

n <- length(sst_1$time)
u_array <- array(NA,dim=c(n,length(lon_o),length(lat_o)))
v_array <- array(NA,dim=c(n,length(lon_o),length(lat_o)))
for(i in 1:(n)){
  uv <- erddap_extract_oscar(uv_grab,currents,i) #SSH
  u_array[i,,] <- uv$u
  v_array[i,,] <- uv$v
}
lon_lat <- expand.grid(lon_o,lat_o)
names(lon_lat) <- c('lon','lat')
extend <- 3/4


setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
pdf('WS_sst_ssh.pdf', height = 8, width = 10,useDingbats = T)
par(mfrow=c(2,2),mar=c(4,4,1.5,1))
for(i in 1:n){
  imagePlot(sst_1$lon,sst_1$lat,sst_1$data[,,i],
            breaks=pretty(range(sst_1$data,na.rm=T),n=50),
            col = plasma(length(pretty(range(sst_1$data,na.rm=T),n=50))-1),asp=1)
  arrows(lon_lat$lon,lon_lat$lat, #RTOFs and SSH
         lon_lat$lon+(u_array[i,,]*extend),lon_lat$lat+(v_array[i,,]*extend),
         length=.025,col='gray10')
  mtext(sst_1$time[i])
}
dev.off()

jday <- yday(sst_1$time)
slope <- p_val <- matrix(NA,dim(sst_1$data)[1],dim(sst_1$data)[2])
for(i in 1:dim(sst_1$data)[2]){
  for(j in 1:dim(sst_1$data)[2]){
    sst <- sst_1$data[i,j,]
    if(all(!is.na(sst))){
      tmp_res <- summary(lm(sst~jday))
      slope[i,j] <- tmp_res$coefficients[2]
      p_val[i,j] <- tmp_res$coefficients[8]
    }
  }
}
range(slope,na.rm=T)
brks <- seq(-.75,.75,.01)
lm_neg <- colorRampPalette(c('dodgerblue4','deepskyblue3','lightskyblue1','gray95'))
lm_pos <- colorRampPalette(c('gray95','rosybrown1','tomato2','red4'))
cols <- c(lm_neg(length(which(brks<0))),
          lm_pos(length(which(brks>0))))
# slope[which(p_val>.1)] <- NA
imagePlot(sst_1$lon,sst_1$lat,slope,breaks=brks,col=cols)
# contour(sst_1$lon,sst_1$lat,p_val,levels=c(1,.1),add=T,col='gray80')


mur_sst <- info('jplMURSST41mday')
mur_ssta <- info('jplMURSST41anommday')

lonbox_e <- -80.6 ### Florida Bay
lonbox_w <- -99 ### mouth of Mississippi River
latbox_n <- 30.5 ### northern coast
latbox_s <- 24.3 ### southern edge of Key West
latitude = c(latbox_s, latbox_n)
longitude = c(lonbox_w, lonbox_e)

time = c(paste(2022,"-12-04",sep=''), paste(2022,"-12-08",sep=''))
sst_grab <- griddap(mur_ssta, latitude=latitude, longitude=longitude, time=time, fields='sstAnom')
sst_1 <- erddap_extract(sst_grab,mur_ssta,'sstAnom')

par(mfrow=c(2,3))
for(i in 1:6){
  imagePlot(sst_1$lon,
            sst_1$lat,
            sst_1$data[,,i],
            asp=1,breaks=brks,col=cols)
  mtext(month.name[i+3])
}

