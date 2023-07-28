rm(list=ls())

library(fields)
library(lubridate)
library(rerddap)
library(sf)
library(terra)
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
world <- vect('ne_10m_admin_0_countries.shp')

setwd("~/Desktop/professional/biblioteca/data/shapefiles/fknms_py2")
fknms <- vect('fknms_py.shp')
crs(fknms) <- crs(world)
fknms <- st_as_sf(fknms)

# SST anomaly
cols_ssta_neg <- colorRampPalette(c('dodgerblue4','deepskyblue3','lightskyblue1','gray95'))
cols_ssta_pos <- colorRampPalette(c('gray95','rosybrown1','tomato2','red4'))
# breaks_ssta <- seq(-round(max(abs(sst_a_grab$data$sstAnom),na.rm=T),digits=1),
#                    round(max(abs(sst_a_grab$data$sstAnom),na.rm=T),digits=1),
#                    by=.2)


# mur_sst_a <- info('jplMURSST41anommday')
mur_sst <- info('jplMURSST41')
# currents <- info('miamicurrents') #SSH derived
# ssha <- info('nesdisSSH1day')

### whole Gulf of Mexico
# latitude = c(17, 31)
# longitude = c(-98, -80)
### west Florida Shelf
lonbox_e <- -80 ### Florida Bay
lonbox_w <- -83.2 ### mouth of Mississippi River
latbox_n <- 26 ### northern coast
latbox_s <- 24 ### southern edge of Key West
latitude = c(latbox_s, latbox_n)
longitude = c(lonbox_w, lonbox_e)

time = c(paste(2023,"-01-01",sep=''), paste(2023,"-01-01",sep=''))
sst_grab <- griddap(mur_sst, latitude=latitude, longitude=longitude, time=time, fields='analysed_sst')

point.sf <- st_as_sf(sst_grab$data, coords = c("longitude","latitude"))
st_crs(point.sf) <- st_crs(fknms)
# data <- st_filter(point.sf, fknms)
mat <- st_intersects(point.sf, fknms)
ind_fknms <- apply(mat, 1, any)

# data_fknms <- sst_grab$data[ind_fknms,]

sst_day <- data.frame(matrix(NA,366*length(2003:2023),2))
m <- 1
n <- 0
for(i in 2003:2023){
  print(i)
  end <- ifelse(i==2023,'-07-25','-12-31')
  time <- c(paste(i,"-01-01",sep=''),paste(i,end,sep='')) 
  sst_grab <- griddap(mur_sst, latitude=latitude, longitude=longitude, time=time, fields='analysed_sst')
  dats <- unique(sst_grab$data$time)
  test <- matrix(sst_grab$data$analysed_sst,length(ind_fknms),length(dats))
  test <- test[ind_fknms,]
  n <- n + length(dats)
  sst_day[m:n,1] <- as.character(dats)
  sst_day[m:n,2] <- apply(test,2,mean,na.rm=T)
  m <- n + 1
  # dats <- seq(as.Date(paste(i,"-01-01",sep='')),as.Date(paste(i,end,sep='')),by='day')
  # for(j in 1:length(dats)){
  #   print(j)
  #   # time = c(as.character(dats[j]),as.character(dats[j]))
  #   # sst_grab <- griddap(mur_sst, latitude=latitude, longitude=longitude, time=time, fields='analysed_sst')
  #   sst_tmp <- sst_grab$data$analysed_sst[which(sst_grab$data$time==dats[j])]
  #   sst_day[n,1] <- as.character(dats[j])
  #   sst_day[n,2] <- mean(sst_tmp[ind_fknms],na.rm=T)
  #   n <- n + 1
  # }
  # cache_delete_all(force=T)
}
cache_list()
cache_delete_all(force=T)


time = c(paste(2023,"-01-01",sep=''), paste(2023,"-01-02",sep=''))
sst_grab <- griddap(mur_sst, latitude=latitude, longitude=longitude, time=time, fields='analysed_sst')

mean(data_fknms$analysed_sst,na.rm=T)




sst_1 <- erddap_extract(sst_grab,mur_sst,'analysed_sst')
sst_30 <- sst_1$data[,,(205-30):205]

imagePlot(sst_1$lon,sst_1$lat,sst_1$data[,,31],
          breaks=pretty(range(sst_1$data[,,31],na.rm=T),n=50),
          col = plasma(length(pretty(range(sst_1$data[,,31],na.rm=T),n=50))-1),asp=1)
mtext(sst_1$time[205])

jday <- yday(sst_1$time[(205-30):205])
slope <- p_val <- matrix(NA,dim(sst_30)[1],dim(sst_30)[2])
for(i in 1:dim(sst_30)[1]){
  for(j in 1:dim(sst_30)[2]){
    sst <- sst_30[i,j,]
    if(all(!is.na(sst))){
      tmp_res <- summary(lm(sst~jday))
      slope[i,j] <- tmp_res$coefficients[2]
      p_val[i,j] <- tmp_res$coefficients[8]
    }
  }
}
range(slope,na.rm=T)
brks <- seq(-.1,.1,.01)
lm_neg <- colorRampPalette(c('dodgerblue4','deepskyblue3','lightskyblue1','gray95'))
lm_pos <- colorRampPalette(c('gray95','rosybrown1','tomato2','red4'))
cols <- c(lm_neg(length(which(brks<0))),
          lm_pos(length(which(brks>0))))
slope[which(p_val>.1)] <- NA
imagePlot(sst_1$lon,sst_1$lat,slope,breaks=brks,col=cols)
contour(sst_1$lon,sst_1$lat,p_val,levels=c(.05,1),add=T,col='gray80')


sst_ts <- apply(sst_1$data,3,mean,na.rm=T)
plot(ymd_hms(sst_1$time),sst_ts,typ='l')
