library(lubridate)
library(rerddap)

cache_delete_all()

query <- ed_search_adv(query='Walton Smith',url='https://gcoos5.geos.tamu.edu/erddap/',
                       minTime ="2000-01-01T00:00:00Z",
                       maxTime="2023-07-23T00:00:00Z",
                       page_size = 4000,
                       minLon = -88, maxLon = -80,
                       minLat = 24, maxLat = 29)

query_res <- as.data.frame(query$info)
ws_info <- base::strsplit(as.character(query_res[,1]),',')
ws_date <- ymd(lapply(ws_info,function(x)x[4]))
ws_lonlat <- lapply(ws_info,function(x)x[5])
test <- base::strsplit(unlist(as.character(ws_lonlat)),' ')
test2 <- gsub('[A-Z]', '', unlist(test))
lons <- as.numeric(test2[seq(3,length(test2),3)])
lats <- as.numeric(test2[seq(2,length(test2),3)])
plot(-lons,lats,asp=1)

table(month(ws_date),year(ws_date))

ind_sum <- which(month(ws_date)>5 & month(ws_date)<9)
tail(query_res[ind_sum,],n=50)

plot(-lons[ind_sum],lats[ind_sum],asp=1)

profile_out <- data.frame(matrix(NA,dim(query$info)[1]*50,9))
sum_out <- data.frame(matrix(NA,dim(query$info)[1],17))
m <- 1
n <- 0
y <- seq(1/24,101/24,1/24)
# i=1
pb <- txtProgressBar(min = 0, max = nrow(query$info), initial = 0, style=3) 
for(i in 1:nrow(query$info)){
  id <- query$info$dataset_id[i]
  url <- paste0('https://gcoos5.geos.tamu.edu/erddap/tabledap/',id,'.csv')
  data <- read.csv(url)
  data <- type.convert(data[-1,],as.is=T)
  # data <- tabledap(query$info$dataset_id[i],url='https://gcoos5.geos.tamu.edu/erddap/')
  # lap <- ymd_hms(data$time_elapsed)-ymd_hms(data$time_elapsed)[2]
  # data <- data[which(lap>30),] # 30 or 45 sec
  bottom <- which.max(data$depth)
  data <- data[1:bottom,]
  data_tmp <- matrix(NA,nrow(data),8)
  if(nrow(data)>100 & min(data$depth,na.rm=T)<2){
    
    depth_m <- data$depth
    # depth_m[which(as.numeric(data$depth_qc)>3)] <- NA
    
    dr <- data$descent_rate
    # dr[which(as.numeric(data$descent_rate_qc)>3)] <- NA
    
    nas <- approx(depth_m,xout=which(is.na(depth_m)))
    depth_m[nas$x] <- nas$y
    nas <- approx(dr,xout=which(is.na(dr)))
    dr[nas$x] <- nas$y
    # if(is.na(dr[1])){
      # dr[1:3] <- dr[4]
    # }
    # sm_dr <- smooth.spline(dr,spar=.9)
    # sm_z <- smooth.spline(depth_m,spar=.9)
    
    # par(mfrow=c(3,1))
    # hist(dr)
    # abline(v=0,lwd=2,col=2)
    # plot(dr,typ='l')
    # plot(smooth.spline(dr,spar=.9)$y,typ='l')
    # abline(h=0,lwd=2,col=2)
    # plot(depth_m,typ='l')
    # plot(depth_m,smooth.spline(dr,spar=.9)$y,typ='l')
    
    if(max(depth_m)<4){
      quant <- .95
      cutoff <- .75
    } else {
      quant <- .9
      cutoff <- .85
    }
    top <- ceiling(rev(which(depth_m<2))[1]*cutoff)
    
    # threshold <- (quantile(dr,quant,na.rm=T)-sd(dr,na.rm=T))
    # inds <- NA
    # for(j in 1:(length(depth_m)-100)){
    #   deps <- depth_m[j:(j+100)]
    #   # deps <- sm_z$y[j:(j+100)]
    #   coefs <- coefficients(lm(deps ~ y))
    #   if(coefs[2]>=threshold){
    #     inds <- j
    #     break
    #   }
    # }
    # # inds3 <- which(sm_dr$y>0)[1]
    # inds2 <- ceiling(rev(which(depth_m<2))[1]*cutoff)
    # inds <- ifelse(inds<48,inds,inds-48)
    # # tops <- c(inds,inds2,inds3)
    # tops <- c(inds,inds2)
    # if(is.na(inds)){
    #   tops[1] <- 100000
    # }
    # if(is.na(inds2)){
    #   tops[2] <- 100000
    # }
    # if(diff(tops)>500 & depth_m[inds2]<2){
    #   top <- inds2
    # } else {
    #   top <- tops[which.min(tops)]
    # }
    # 
    # plot(depth_m)
    # abline(v=tops,col=c(1,2,4))
    
    data_tmp <- data[top:bottom,]
    depth_m <- depth_m[top:bottom]
    
    temperature_c <- data_tmp$sea_water_temperature
    # temperature_c[which(as.numeric(data_tmp$sea_water_temperature_qc)>3)] <- NA
    
    salinity_psu <- data_tmp$sea_water_salinity
    # salinity_psu[which(as.numeric(data_tmp$sea_water_salinity_qc)>3)] <- NA
    
    sigma_t <- data_tmp$sea_water_sigma_t
    # sigma_t[which(as.numeric(data_tmp$sea_water_sigma_t_qc)>3)] <- NA
    
    chl_fl <- data_tmp$chlorophyll_fluorescence
    # chl_fl[which(as.numeric(data_tmp$chlorophyll_fluorescence_qc)>3)] <- NA
    
    do <- data_tmp$dissolved_oxygen
    # do[which(as.numeric(data_tmp$dissolved_oxygen_qc)>3)] <- NA
    
    ox_sat <- data_tmp$oxygen_saturation
    # ox_sat[which(as.numeric(data_tmp$oxygen_saturation_qc)>3)] <- NA
    
    data_qc <- data.frame(data_tmp$depth,data_tmp$depth_qc,data_tmp$sea_water_temperature_qc,data_tmp$sea_water_salinity_qc,data_tmp$sea_water_sigma_t_qc,data_tmp$chlorophyll_fluorescence_qc,data_tmp$dissolved_oxygen_qc,data_tmp$oxygen_saturation_qc)
    vars <- c('depth_m','temperature_c','salinity_psu','sigma_t','chl_fl','do','ox_sat')
    err_mes <- NULL
    if(any(data_qc[,2:8]>1,na.rm=T)){
      ind_qc <- which(apply(data_qc[,2:8],2,sum,na.rm=T)>nrow(data_qc))
      for(j in (ind_qc+1)){
        mes <- paste(paste(names(table(data_qc[,j])),round(table(data_qc[,j])/nrow(data_qc),2),sep = '=='),collapse = '; ')
        err <- paste(vars[j-1],mes,sep = ': ')
        err_mes <- paste(err_mes, err, sep=' | ')
      }
    }
    # data_qc <- data.frame(as.numeric(data_tmp$depth),apply(data_qc[,2:8],1,max,na.rm=T),vars[apply(data_qc[,2:8],1,which.max)])q
    data_qc <- data.frame('depth' = data_tmp$depth,
                          'max_qc' = apply(data_qc[,2:8],1,max,na.rm=T))
    
    data_tmp <- data.frame(cbind(depth_m,temperature_c,salinity_psu,sigma_t,chl_fl,do,ox_sat,data_qc$max_qc))
    colnames(data_tmp) <- c('depth_m','temperature_c','salinity_psu','sigma_t','chl_fl','do','ox_sat','flag')
    
    if(any(data_tmp$depth_m<.6,na.rm=T)){
      data_tmp <- data_tmp[-which(data_tmp$depth_m<.6),]
    }
    
    z_cut <- cut(data_tmp$depth_m,
                 seq(round(min(data_tmp$depth_m,na.rm=T))-.5,
                     round(max(data_tmp$depth_m,na.rm=T))+.5,1))
    
    data_agg <- data.frame(matrix(NA,length(levels(z_cut)),ncol(data_tmp)))
    data_agg[,1] <- round(aggregate(data_tmp[,1],by=list(z_cut),mean,na.rm=T)$x)
    data_agg[,8] <- round(aggregate(data_tmp[,8],by=list(z_cut),max,na.rm=T)$x)
    for(j in 2:(ncol(data_tmp)-1)){
      data_agg[,j] <- aggregate(data_tmp[,j],by=list(z_cut),mean,na.rm=T)$x
    }
    names(data_agg) <- c('depth_m','temp_c','sal_psu','sigma_t','chl_rfu','do_mll','oxy_sat','flag')
    lth <- nrow(data_agg)
    n <- n + lth
    
    profile_out[m:n,1] <- paste0('aoml',sprintf('%03d',i))
    profile_out[m:n,2:9] <- data_agg
    m <- n + 1
    
    sum_out[i,7] <- max(data_agg$depth_m,na.rm=T) #6
    sum_out[i,8] <- min(data_agg$depth_m,na.rm=T) #6
    sum_out[i,9:12] <- data_agg[1,c(2:4,6)] #10
    sum_out[i,13:16] <- data_agg[nrow(data_agg),c(2:4,6)] #14
    sum_out[i,17] <- ifelse(is.null(err_mes),NA,err_mes)
    
  }
  sum_out[i,1] <- paste0('aoml',sprintf('%03d',i))
  sum_out[i,2] <- query$info$dataset_id[i]
  sum_out[i,3] <- data$time[1]
  sum_out[i,4] <- data$longitude[1]
  sum_out[i,5] <- data$latitude[1]
  sum_out[i,6] <- max(data$depth,na.rm=T) 
  
  rm(data,data_qc,data_tmp,data_agg,depth_m,temperature_c,salinity_psu,sigma_t,chl_fl,do,ox_sat,dr,top,bottom,lth,z_cut,threshold,ind_qc,inds,inds2)
  setTxtProgressBar(pb,i)
}
names(profile_out) <- c('profile_ind','depth_m','temp_c','sal_psu','sigma_t','chl_rfu','do_mll','oxy_sat','flag')
names(sum_out) <- c('profile_ind','dataset_id','date_utc',
                    'lon_dd','lat_dd',
                    'z_max_report','z_max_bin','z_min_bin',
                    'surf_temp_c','surf_sal_psu','surf_sigma_t','surf_do_mll',
                    'bot_temp_c','bot_sal_psu','bot_sigma_t','bot_do_mll')
profile_out <- profile_out[which(!is.na(profile_out$profile_ind)),]
sum_out <- sum_out[which(!is.na(sum_out$profile_ind)),]


output <- output[which(!is.na(output$ship)),]
output2 <- output[which(output$longitude<(-82) & output$longitude>(-87)),]
plot(output2$longitude,output2$latitude,typ='b')
unique(output2$date_utc)

# write.csv(output2,'202108_WFS_aoml_ctd.csv')

query$alldata[[1]]
i=1
data <- tabledap(query$info$dataset_id[i],url='https://gcoos5.geos.tamu.edu/erddap/')
depth <- as.numeric(data$depth)
dr <- as.numeric(data$descent_rate)
dr <- as.numeric(data$descent_rate)

nas <- approx(depth,xout=which(is.na(depth)))
depth[nas$x] <- nas$y
nas <- approx(dr,xout=which(is.na(dr)))
dr[nas$x] <- nas$y
if(is.na(dr[1])){
  dr[1] <- dr[2]
}

beg <- 45*24
y <- seq(1/24,101/24,1/24)
threshold <- (max(dr,na.rm=T)-sd(dr,na.rm=T)*1.5)
inds <- NA
n <- 1
for(i in beg:(length(depth)-100)){
  deps <- depth[i:(i+100)]
  coefs <- coefficients(lm(deps ~ y))
  if(coefs[2]>=threshold){
    inds <- i
    n <- n + 1
    break
  }
}
plot(depth)
abline(v=inds-36)
plot(dr)
abline(h=threshold)
abline(v=inds-36)
plot(filter(dr,rep(1/101,101)))
abline(h=threshold)
abline(v=inds-36)

top <- inds-36
bottom <- which.max(depth)

depth <- depth[top:bottom]
depth <- depth[-which(depth<.5)]

z_cut <- cut(depth,seq(round(min(depth))-.5,round(max(depth))+.5,1))
aggregate(depth,by=list(z_cut),mean,na.rm=T)



z_range <- range(depth)
x <- aggregate(depth,by=list(cut(depth,seq(0,z_range[2],1))),length)
quantile(x$x,.1)

plot(depth)
plot(dr)
plot(ymd_hms(data$time_elapsed),as.numeric(data$depth))
plot(smooth.spline(depth,spar=.5))

depth <- as.numeric(data$depth)
dr <- as.numeric(data$descent_rate)
nas <- approx(depth,xout=which(is.na(depth)))
depth[nas$x] <- nas$y
nas <- approx(dr,xout=which(is.na(dr)))
dr[nas$x] <- nas$y
dr[1] <- dr[2]
z_smooth <- smooth.spline(depth,dr,spar=.5)
which((z_smooth$y)>0)[1]

z_smooth <- smooth.spline(depth,spar=.5)
dr_smooth <- smooth.spline(dr,spar=.5)
z_smooth <- smooth.spline(depth,dr,spar=.5)
which((z_smooth$y)>0)
which(rev(z_smooth$y)<0)
plot((rev(z_smooth$y)))
plot(diff(rev(z_smooth$y)))
plot(z_smooth$y,dr_smooth$y)

plot(z_smooth$y)
plot(rev(dr_smooth$y))

which(rev(dr_smooth$y)<0)

plot(depth)
abline(v=length(depth)-604)

t.test(depth[1:500],depth[250:750])
boxplot(depth[1:500],depth[250:750])

plot(depth,col=2,typ='l')
points(filter(depth,rep(1/100,100)))

plot(dr)


plot(cumsum(dr),xlim=c(0,500))
plot(depth)




https://gcoos5.geos.tamu.edu/erddap/tabledap/WS23061_WS23061_WS23061_Kelble_Stn_RP3.csv?profile%2Ctime%2Ctime_elapsed%2Clatitude%2Clongitude%2Csea_water_pressure%2Csea_water_pressure_qc%2Csea_water_temperature%2Csea_water_temperature_qc%2Csea_water_temperature_2%2Csea_water_temperature_2_qc%2Csea_water_electrical_conductivity%2Csea_water_electrical_conductivity_qc%2Csea_water_electrical_conductivity_2%2Csea_water_electrical_conductivity_2_qc%2Cchlorophyll_fluorescence%2Cchlorophyll_fluorescence_qc%2CCDOM%2CCDOM_qc%2Cchlorophyll_concentration%2Cchlorophyll_concentration_qc%2Cbeam_attenuation%2Cbeam_transmission%2Cbeam_attenuation_qc%2Cbeam_transmission_qc%2Cdissolved_oxygen%2Cdissolved_oxygen_qc%2Coxygen_saturation%2Coxygen_saturation_qc%2Coxygen_saturation_2%2Coxygen_saturation_2_qc%2Cphotosynthetically_available_radiation%2Cphotosynthetically_available_radiation_qc%2Csurface_photosynthetically_available_radiation%2Csurface_photosynthetically_available_radiation_qc%2Caltimeter%2Caltimeter_qc%2Cdepth%2Cdepth_qc%2Csea_water_salinity%2Csea_water_salinity_qc%2Csea_water_salinity_2%2Csea_water_salinity_2_qc%2Csea_water_sigma_t%2Csea_water_sigma_t_qc%2Csea_water_sigma_t_2%2Csea_water_sigma_t_2_qc%2Cdescent_rate%2Cdescent_rate_qc%2Csound_velocity%2Csound_velocity_qc%2Cflag%2Cplatform%2Cinstrument%2Cinstrument1%2Cinstrument2%2Cinstrument3%2Cinstrument4%2Cinstrument5%2Cinstrument6%2Cinstrument7%2Cinstrument8%2Cinstrument9%2Cinstrument10%2Cinstrument11%2Cinstrument12%2Cinstrument13%2Ccrs&time%3E=2023-02-25T00%3A00%3A00Z&time%3C=2023-03-04T21%3A31%3A24Z

https://gcoos5.geos.tamu.edu/erddap/tabledap/WS23061_WS23061_WS23061_Kelble_Stn_RP3.csv?profile%2Ctime%2Clatitude%2Clongitude%2Csea_water_pressure%2Csea_water_pressure_qc%2Csea_water_temperature%2Csea_water_temperature_qc%2Csea_water_temperature_2%2Csea_water_temperature_2_qc%2Csea_water_electrical_conductivity%2Csea_water_electrical_conductivity_qc%2Csea_water_electrical_conductivity_2%2Csea_water_electrical_conductivity_2_qc%2Cchlorophyll_fluorescence%2Cchlorophyll_fluorescence_qc%2CCDOM%2CCDOM_qc%2Cchlorophyll_concentration%2Cchlorophyll_concentration_qc%2Cbeam_attenuation%2Cbeam_transmission%2Cbeam_attenuation_qc%2Cbeam_transmission_qc%2Cdissolved_oxygen%2Cdissolved_oxygen_qc%2Coxygen_saturation%2Coxygen_saturation_qc%2Coxygen_saturation_2%2Coxygen_saturation_2_qc%2Cphotosynthetically_available_radiation%2Cphotosynthetically_available_radiation_qc%2Csurface_photosynthetically_available_radiation%2Csurface_photosynthetically_available_radiation_qc%2Caltimeter%2Caltimeter_qc%2Cdepth%2Cdepth_qc%2Csea_water_salinity%2Csea_water_salinity_qc%2Csea_water_salinity_2%2Csea_water_salinity_2_qc%2Csea_water_sigma_t%2Csea_water_sigma_t_qc%2Csea_water_sigma_t_2%2Csea_water_sigma_t_2_qc%2Cdescent_rate%2Cdescent_rate_qc%2Csound_velocity%2Csound_velocity_qc%2Cflag%2Cplatform%2Cinstrument%2Cinstrument1%2Cinstrument2%2Cinstrument3%2Cinstrument4%2Cinstrument5%2Cinstrument6%2Cinstrument7%2Cinstrument8%2Cinstrument9%2Cinstrument10%2Cinstrument11%2Cinstrument12%2Cinstrument13%2Ccrs&time%3E=2023-02-25T00%3A00%3A00Z&time%3C=2023-03-04T21%3A31%3A24Z



url <- 'https://gcoos5.geos.tamu.edu/erddap/tabledap/WS23061_WS23061_WS23061_Kelble_Stn_RP3.csv'
data <- read.csv(url)

i=1
id <- query$info$dataset_id[i]
url <- paste0('https://gcoos5.geos.tamu.edu/erddap/tabledap/',id,'.csv')
system.time(data1 <- read.csv(url))
data1 <- type.convert(data1[-1,],as.is=T)
system.time(data2 <- tabledap(id,url='https://gcoos5.geos.tamu.edu/erddap/'))
