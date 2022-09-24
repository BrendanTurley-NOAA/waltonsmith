library(akima)
library(fields)
library(oce)

source('~/Desktop/professional/projects/Postdoc_FL/scripts/FCWC/FCWC_process_htm.R')

### colors
d_col <- colorRampPalette(rev(c('purple4','purple2','orchid1','gray90')))
m_col <- colorRampPalette(rev(c('blue4','dodgerblue2','deepskyblue1','gray90')))
t_col <- colorRampPalette(c(1,'purple','darkorange','gold'))
s_col <- colorRampPalette(c('purple4','dodgerblue4','seagreen3','khaki1'))
c_col <- colorRampPalette(c('honeydew2','darkseagreen3','darkgreen'))
strat_n_col <- colorRampPalette(c('purple4','purple2','orchid1','gray90'))
strat_p_col <- colorRampPalette(rev(c('darkgreen','green3','palegreen2','gray90')))

ox.col1 <- colorRampPalette(c(1,'firebrick4','red'))
ox.col2 <- colorRampPalette(c('darkgoldenrod4','goldenrod2','gold'))
ox.col3 <- colorRampPalette(c('midnightblue','dodgerblue4','deepskyblue2','cadetblue1','azure'))
o_breaks <- seq(0,10,by=.25)
o_cols <- c(ox.col1(length(o_breaks[o_breaks<2])),
            ox.col2(length(o_breaks[o_breaks>=2 & o_breaks<3.5])),
            ox.col3(length(o_breaks[o_breaks>=3.5])-1))

###----------- profile plots -----------
setwd('~/Desktop/noaa_ctd')
files <- list.files()
files <- files[grep('.cnv',files)]

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

files_plt <- data$filename[which(data$station>=85 & data$station<=90)]

# tmp <- data[which(is.element(data$station,c(69,84,83,82,87,88))),]
# tmp <- tmp[order(tmp$m_depth,decreasing = T),]
# files_plt <- tmp$filename

setwd('~/Desktop/noaa_ctd')

m <- 1
n <- 0
data_d <- data.frame(matrix(NA,200*7,8))
for(i in 1:length(files_plt)){
  tmp <- read.ctd(files_plt[i])
  n <- n + length(tmp@data$depth)
  data_d[m:n,1] <- tmp@data$depth
  data_d[m:n,2] <- tmp@metadata$station
  data_d[m:n,3] <- tmp@metadata$longitude
  data_d[m:n,4] <- tmp@metadata$latitude
  data_d[m:n,5] <- tmp@data$temperature
  data_d[m:n,6] <- tmp@data$salinity
  data_d[m:n,7] <- tmp@data$fluorescence
  data_d[m:n,8] <- tmp@data$oxygen
  m <- n + 1
}
data_d <- data_d[1:n,]
names(data_d) <- c('depth','station','lon','lat','temp','sal','chl','do')

data_d$depth <- -data_d$depth
data_d <- data_d[order(data_d$lon),]
data_d$chl[which(data_d$chl<0)] <- 0

### find deepest samples to create a bottom for plotting
longitudes <- data_d$lon
depths <- data_d$depth

bottoms <- bottom_finder(data_d$lon,data_d$depth)

### create interpolated cross-sections for plotting
resolution <- 100
temp_in <- interp(data_d$lon,data_d$depth,data_d$temp,
                  yo=seq(min(data_d$depth,na.rm=T), max(data_d$depth,na.rm=T), length = resolution),
                  xo=seq(min(data_d$lon,na.rm=T), max(data_d$lon,na.rm=T), length = resolution*3),
                  duplicate='mean',linear=T)
sal_in <- interp(data_d$lon,data_d$depth,data_d$sal,
                 yo=seq(min(data_d$depth,na.rm=T), max(data_d$depth,na.rm=T), length = resolution),
                 xo=seq(min(data_d$lon,na.rm=T), max(data_d$lon,na.rm=T), length = resolution*3),
                 duplicate='mean',linear=T)
chl_in <- interp(data_d$lon,data_d$depth,log(data_d$chl+1,base=10),
                 yo=seq(min(data_d$depth,na.rm=T), max(data_d$depth,na.rm=T), length = resolution),
                 xo=seq(min(data_d$lon,na.rm=T), max(data_d$lon,na.rm=T), length = resolution*3),
                 duplicate='mean',linear=T)
do_in <- interp(data_d$lon,data_d$depth,data_d$do,
                yo=seq(min(data_d$depth,na.rm=T), max(data_d$depth,na.rm=T), length = resolution),
                xo=seq(min(data_d$lon,na.rm=T), max(data_d$lon,na.rm=T), length = resolution*3),
                duplicate='mean',linear=T)

t_breaks <- pretty(temp_in$z,n=20)
t_cols <- t_col(length(t_breaks)-1)
s_breaks <- pretty(sal_in$z,n=20)
s_cols <- s_col(length(s_breaks)-1)
c_breaks <- pretty(chl_in$z,n=20)
c_cols <- c_col(length(c_breaks)-1)

imagePlot(temp_in,breaks=t_breaks,col=t_cols)
polygon(c(bottoms[,2]-5,bottoms[,2],bottoms[,2]+5),c(bottoms[,1]-100,bottoms[,1],bottoms[,1]-100),col='wheat4')
points(data_d$lon,data_d$depth,pch=20,cex=.1,col='gray50')
imagePlot(sal_in,breaks=s_breaks,col=s_cols)
polygon(c(bottoms[,2]-5,bottoms[,2],bottoms[,2]+5),c(bottoms[,1]-100,bottoms[,1],bottoms[,1]-100),col='wheat4')
points(data_d$lon,data_d$depth,pch=20,cex=.1,col='gray50')
imagePlot(chl_in,breaks=c_breaks,col=c_cols)
polygon(c(bottoms[,2]-5,bottoms[,2],bottoms[,2]+5),c(bottoms[,1]-100,bottoms[,1],bottoms[,1]-100),col='wheat4')
points(data_d$lon,data_d$depth,pch=20,cex=.1,col='gray50')
imagePlot(do_in,breaks=o_breaks,col=o_cols)
polygon(c(bottoms[,2]-5,bottoms[,2],bottoms[,2]+5),c(bottoms[,1]-100,bottoms[,1],bottoms[,1]-100),col='wheat4')
points(data_d$lon,data_d$depth,pch=20,cex=.1,col='gray50')

o_inter <- interpBarnes(data_d$lon,data_d$depth,data_d$do,pregrid=F,
             yg=seq(min(data_d$depth,na.rm=T), max(data_d$depth,na.rm=T), length = resolution),
             xg=seq(min(data_d$lon,na.rm=T), max(data_d$lon,na.rm=T), length = resolution*3))
imagePlot(o_inter$xg,o_inter$yg,o_inter$zg,breaks=o_breaks,col=o_cols)
polygon(c(bottoms[,2]-5,bottoms[,2],bottoms[,2]+5),c(bottoms[,1]-100,bottoms[,1],bottoms[,1]-100),col='wheat4')
