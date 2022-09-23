library(scales)
rescale(data@data$temperature,to = c(0, 1))

setwd('~/Desktop/noaa_ctd')
files <- list.files()
files <- files[grep('.cnv',files)]

setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
pdf('new_plot.pdf',width=9,height=12,pointsize=12)
par(mfrow=c(2,2))
setwd('~/Desktop/noaa_ctd')
for(i in 1:length(files)){
  data <- read.ctd(files[i])
  
  if(data@metadata$cruise=='2203' & data@metadata$latitude<30.5 & data@metadata$longitude>=-87){
    # plot(data@data$temperature,-data@data$depth,typ='o',col='red',pch=16)
    # plot(data@data$salinity,-data@data$depth,typ='o',col='purple',pch=16)
    # plot(data@data$fluorescence,-data@data$depth,typ='o',col='forestgreen',pch=16)
    # plot(data@data$oxygen,-data@data$depth,typ='o',col=4,pch=16)
    # abline(v=2,lty=2)
    
    tmp_brks <- pretty(data@data$temperature,min.n=5)
    sal_brks <- pretty(data@data$salinity,min.n=5)
    chl_brks <- pretty(data@data$fluorescence,min.n=5)
    oxy_brks <- pretty(data@data$oxygen,min.n=5)
    
    par(mar=c(6.5,5.5,8.5,2.5))
    plot(rescale(data@data$temperature,to = c(0, 1)),-data@data$depth,
         typ='l',lwd=2,col='red',xlab='',ylab='Depth (m)',xaxt='n',las=1)
    points(rescale(data@data$salinity,to = c(0, 1)),-data@data$depth,
           typ='l',lwd=2,col='purple')
    points(rescale(data@data$fluorescence,to = c(0, 1)),-data@data$depth,
           typ='l',lwd=2,col='green3')
    points(rescale(data@data$oxygen,to = c(0, 1)),-data@data$depth,
           typ='l',lwd=2,col=4)
    axis(1,seq(0,1,length.out=length(tmp_brks)),tmp_brks,line=0,col='red',col.axis = 'red',lwd=2)
    mtext('Temperature',1,line=0,col='red',cex=.75,adj=1,at=-.01)
    axis(1,seq(0,1,length.out=length(sal_brks)),sal_brks,line=3,col='purple',col.axis='purple',lwd=2)
    mtext('Salinity',1,line=3,col='purple',cex=.75,adj=1,at=-.01)
    axis(3,seq(0,1,length.out=length(chl_brks)),chl_brks,line=0,col='green3',col.axis='green3',lwd=2)
    mtext('Chlorophyll',3,line=0,col='green3',cex=.75,adj=1,at=-.01)
    axis(3,seq(0,1,length.out=length(oxy_brks)),oxy_brks,line=3,col=4,col.axis=4,lwd=2)
    mtext('Dissolved oxygen',3,line=3,col=4,cex=.75,adj=1,at=-.01)
    mtext(paste0('Lat ',round(data@metadata$latitude,3),'N',
                 ', Lon ',abs(round(data@metadata$longitude,3)),'W; ',
                 data@metadata$date,' UTC'),
          outer=F,adj=0,line=5.5,cex=.75)
    mtext(paste('R/V ',data@metadata$ship, '- Station:',data@metadata$station),
          outer=F,adj=0,line=6.5,cex=.75)
  }
}
dev.off()

lonbox_e <- -81.5 ### Florida Bay
lonbox_w <- -87 ### mouth of Mississippi River
latbox_n <- 30.5 ### northern coast
latbox_s <- 24.3 ### southern edge of Key West

setwd('~/Desktop/noaa_ctd')
data <- read.csv('NOAA_NMFS_bot_dat.csv')
ind <- which(data$cruise=='2203' & 
               data$lat<latbox_n & data$lat>latbox_s &
               data$lon<lonbox_e & data$lon>lonbox_w)
data <- data[ind,]

setwd("~/Desktop/professional/projects/Postdoc_FL/figures")
png('cruise_tracks.png',width=5,height=5,units='in',pointsize=12,res=300)
plot(data$lon,data$lat,typ='o',asp=1,col=4,bg='white',pch=21,cex=2,xlab='Longitude',ylab='Latitude',lwd=2)
plot(world,add=T,col='gray80')
contour(topo_lon,
        topo_lat,
        topo,
        add=T,levels=c(-200,-100,-50,-25,-10),col='gray70')
mtext('R/V GORDON GUNTER cruise track')
text(data$lon[1],data$lat[1],'start',pos=1,font=2)
text(data$lon[nrow(data)],data$lat[nrow(data)],'end',pos=4,font=2)
text(data$lon,data$lat,data$station,font=2,cex=.5)
dev.off()


files_plt <- data$filename[which(data$station>=85 & data$station<=90)]

tmp <- data[which(is.element(data$station,c(69,84,83,82,87,88))),]
tmp <- tmp[order(tmp$m_depth,decreasing = T),]
files_plt <- tmp$filename

setwd('~/Desktop/noaa_ctd')

n <- 0
plot(1,1,xlim=c(0,length(files_plt)),ylim=c(-200,0),typ='n',xlab='',ylab='Depth (m)')
for(i in 1:length(files_plt)){
  data <- read.ctd(files_plt[i])
    
  points(rescale(data@data$temperature,to = c(0, 1))+n,-data@data$depth,typ='l', col='red')
  points(rescale(data@data$salinity,to = c(0, 1))+n,-data@data$depth,typ='l',col='purple')
  # points(rescale(data@data$fluorescence,to = c(0, 1))+n,-data@data$depth,typ='l',col='green3')
  # points(rescale(data@data$oxygen,to = c(0, 1))+n,-data@data$depth,typ='l',col=4)
  n <- n + 1
  
}


n <- 0
plot(1,1,xlim=c(0,length(files_plt)),ylim=c(-200,0),typ='n',xlab='',ylab='Depth (m)')
for(i in 1:length(files_plt)){
  data <- read.ctd(files_plt[i])
  
  # points(rescale(data@data$temperature,to = c(0, 1))+n,-data@data$depth,typ='l', col='red')
  # points(rescale(data@data$salinity,to = c(0, 1))+n,-data@data$depth,typ='l',col='purple')
  points(rescale(data@data$fluorescence,to = c(0, 1))+n,-data@data$depth,typ='l',col='green3')
  points(rescale(data@data$oxygen,to = c(0, 1))+n,-data@data$depth,typ='l',col=4)
  n <- n + 1
  
}

m <- 1
n <- 0
data_d <- data.frame(matrix(NA,200*7,6))
for(i in 1:length(files_plt)){
  tmp <- read.ctd(files_plt[i])
  n <- n + length(tmp@data$depth)
  data_d[m:n,1] <- tmp@data$depth
  data_d[m:n,2] <- tmp@metadata$station
  data_d[m:n,3] <- tmp@data$temperature
  data_d[m:n,4] <- tmp@data$salinity
  data_d[m:n,5] <- tmp@data$fluorescence
  data_d[m:n,6] <- tmp@data$oxygen
  m <- n + 1
}
data_d <- data_d[1:n,]
names(data_d) <- c('depth','station','temp','sal','chl','do')
data_d$station <- as.factor(data_d$station)
data_d$station <- factor(data_d$station, levels=c('069','084','087','083','082','088'))
# levels(data_d$station) <- c('069','084','083','087','082','088')

t_breaks <- pretty(data_d$temp,n=10) 
t_cols <- t_col(length(t_breaks))

plot(as.numeric(data_d$station),-data_d$depth,pch=21,
     bg=t_cols[as.numeric(cut(data_d$temp,t_breaks))],
     col=t_cols[as.numeric(cut(data_d$temp,t_breaks))])

s_breaks <- pretty(data_d$sal,n=10) 
s_cols <- s_col(length(s_breaks))

plot(as.numeric(data_d$station),-data_d$depth,pch=21,
     bg=s_cols[as.numeric(cut(data_d$sal,s_breaks))],
     col=s_cols[as.numeric(cut(data_d$sal,s_breaks))])

c_breaks <- pretty(data_d$chl,n=10) 
c_cols <- c_col(length(c_breaks))

plot(as.numeric(data_d$station),-data_d$depth,pch=21,
     bg=c_cols[as.numeric(cut(data_d$chl,c_breaks))],
     col=c_cols[as.numeric(cut(data_d$chl,c_breaks))])


plot(as.numeric(data_d$station),-data_d$depth,pch=21,
     bg=o_cols[as.numeric(cut(data_d$do,o_breaks))],
     col=o_cols[as.numeric(cut(data_d$do,o_breaks))])
