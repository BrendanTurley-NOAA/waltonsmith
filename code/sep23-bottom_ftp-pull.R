

ind <- which(year(bottom_dat$date_utc)==2023 & month(bottom_dat$date_utc)>=9)

dat_0923 <- bottom_dat[ind,]

plot(dat_0923$lon, dat_0923$lat, asp = 1)

plot(dat_0923$lon, dat_0923$lat, asp = 1,
     cex = dat_0923$bot_do)
hist(dat_0923$bot_do)

par(xpd=F)
plot(dat_0923$lon, dat_0923$lat,pch=21,asp=1,cex=2,las=1,
     bg=o_cols[as.numeric(cut(dat_0923$bot_do,o_breaks))],
     col=o_cols[as.numeric(cut(dat_0923$bot_do,o_breaks))],
     xlab='',ylab='')
lines(world)
plot(world,add=T,col='gray80')
contour(topo_lon,
        topo_lat,
        topo,
        add=T,levels=c(-200,-100,-50,-25,-10),col='gray70')

hist(dat_0923$pycno)
m_breaks <- pretty(dat_0923$pycno, n = 20)
m_cols <- m_col(length(m_breaks)-1)
plot(dat_0923$lon, dat_0923$lat,pch=21,asp=1,cex=2,las=1,
     bg=m_cols[as.numeric(cut(dat_0923$pycno,m_breaks))],
     col=m_cols[as.numeric(cut(dat_0923$pycno,m_breaks))],
     xlab='',ylab='')

hist(dat_0923$dd_dz1)
m_breaks <- pretty(dat_0923$dd_dz1, n = 20)
m_cols <- m_col(length(m_breaks)-1)
plot(dat_0923$lon, dat_0923$lat,pch=21,asp=1,cex=2,las=1,
     bg=m_cols[as.numeric(cut(dat_0923$dd_dz1,m_breaks))],
     col=m_cols[as.numeric(cut(dat_0923$dd_dz1,m_breaks))],
     xlab='',ylab='')

hist(dat_0923$dd_dz2)
m_breaks <- pretty(dat_0923$dd_dz2, n = 20)
m_cols <- m_col(length(m_breaks)-1)
plot(dat_0923$lon, dat_0923$lat,pch=21,asp=1,cex=2,las=1,
     bg=m_cols[as.numeric(cut(dat_0923$dd_dz2,m_breaks))],
     col=m_cols[as.numeric(cut(dat_0923$dd_dz2,m_breaks))],
     xlab='',ylab='')

hist(dat_0923$z_cmax2)
m_breaks <- pretty(dat_0923$z_cmax2, n = 20)
c_cols <- c_col(length(m_breaks)-1)
plot(dat_0923$lon, dat_0923$lat,pch=21,asp=1,cex=2,las=1,
     bg=c_cols[as.numeric(cut(dat_0923$z_cmax2,m_breaks))],
     col=c_cols[as.numeric(cut(dat_0923$z_cmax2,m_breaks))],
     xlab='',ylab='')

hist(dat_0923$mld2)
m_breaks <- pretty(dat_0923$mld2, n = 20)
m_cols <- m_col(length(m_breaks)-1)
plot(dat_0923$lon, dat_0923$lat,pch=21,asp=1,cex=2,las=1,
     bg=m_cols[as.numeric(cut(dat_0923$mld2,m_breaks))],
     col=m_cols[as.numeric(cut(dat_0923$mld2,m_breaks))],
     xlab='',ylab='')



library(curl)
url <- 'ftp url from notes'
h <- new_handle(dirlistonly=T)
con <- curl(url, 'r', h)
tbl <- read.table(con, stringsAsFactors = T, fill = T)
close(con)
head(tbl)
urls <- paste0(url, tbl[1,1])
fls <- basename(urls)
setwd("/Users/Brendan/Desktop/noaa_ctd")
curl_fetch_disk(urls, fls)
?curl_download
