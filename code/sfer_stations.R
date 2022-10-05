rm(list=ls())
gc()

library(fields)
library(lubridate)
library(NISTunits)
library(raster)
library(ncdf4)
library(rgdal)


setwd("~/Desktop/professional/biblioteca/data")
bathy <- nc_open('etopo1.nc')
topo <- ncvar_get(bathy, 'Band1')
topo_lat <- ncvar_get(bathy, 'lat')
topo_lon <- ncvar_get(bathy, 'lon')
nc_close(bathy)

topo <- NISTmeterTOft(topo)

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
setwd("~/Desktop/professional/biblioteca/data/shapefiles/Florida_Shoreline__1_to_40%2C000_Scale_-shp")
world <- readOGR('Florida_Shoreline__1_to_40%2C000_Scale_.shp')
world <- crop(world, extent(-86, -79, 24.5, 28))
# setwd("~/Desktop/professional/biblioteca/data/shapefiles/gshhg-shp-2.3.7/GSHHS_shp/h/")
# world <- readOGR('GSHHS_h_L1.shp')
# world <- crop(world, extent(-86, -79, 24.5, 28))

setwd('~/Desktop/professional/projects/Postdoc_FL/data')
cities <- read.csv('fl_cities.csv')

setwd('~/Desktop/professional/projects/Postdoc_FL/data/walton_smith')
data <- read.csv('WS22022_Sample_Log.csv')

plot(data$Longitude.Decimal,data$Latitude.Decimal)
stations <- sort(unique(data$Station))
ind <- c(1:67,87:89,95,106:108,118)
stations[ind]

orig <- is.element(data$Station,stations[ind])


setwd('~/Desktop/professional/projects/Postdoc_FL/figures')
png('WS_stations.png', height = 8, width = 8, units = 'in', res=300)
plot(data$Longitude.Decimal,data$Latitude.Decimal,
     asp=1,type='n',las=1,
     xlab='',ylab='')
plot(world,add=T,col='gray70')
# contour(topo_lon,topo_lat,topo,add=T,levels=c(-100,-50,-25,-10),col='gray40')
points(data$Longitude.Decimal[orig],data$Latitude.Decimal[orig],pch=21,bg='green3',cex=2)
points(data$Longitude.Decimal[!orig],data$Latitude.Decimal[!orig],pch=25,bg='orange2',cex=2)
legend('bottomleft',c('Pre-2018','New stations'),pch=c(21,25),pt.bg=c('green3','orange2'),bty='n',pt.cex=2)
mtext(expression(paste('Longitude (',degree,'W)')),1,line=3,cex=1)
mtext(expression(paste('Latitude (',degree,'N)')),2,line=3,cex=1)
dev.off()


cities <- cities[c(5:8,11:14),]
cities$pos <- c(rep(4,3),2,2,4,2,2)
cities$name[6] <- 'RSMAS'

setwd('~/Desktop/professional/projects/Postdoc_FL/figures')
png('WS_stations_4z.png', height = 20, width = 20, units = 'in', res=300)
plot(data$Longitude.Decimal,data$Latitude.Decimal,
     asp=1,las=1,
     xlab='',ylab='')
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "lightblue1")
plot(world,add=T,col='seashell2')
contour(topo_lon,topo_lat,topo,add=T,levels=c(-300,-200,-100,-50,-25,-10),col='gray40')
points(data$Longitude.Decimal,data$Latitude.Decimal,cex=2,pch=21,bg='plum1',lwd=1.5)
text(data$Longitude.Decimal,data$Latitude.Decimal,data$Station,cex=.5)
points(cities$longitude,cities$latitude,pch=21,bg='white',cex=2)
points(cities$longitude,cities$latitude,pch=16,cex=1)
text(cities$longitude,cities$latitude,cities$name,cex=2,pos=cities$pos,font=2)
# points(data$Longitude.Decimal[orig],data$Latitude.Decimal[orig],pch=21,bg='green3',cex=2)
# points(data$Longitude.Decimal[!orig],data$Latitude.Decimal[!orig],pch=25,bg='orange2',cex=2)
# legend('bottomleft',c('Pre-2018','New stations'),pch=c(21,25),pt.bg=c('green3','orange2'),bty='n',pt.cex=2)
mtext(expression(paste('Longitude (',degree,'W)')),1,line=3,cex=1)
mtext(expression(paste('Latitude (',degree,'N)')),2,line=3,cex=1)
dev.off()

