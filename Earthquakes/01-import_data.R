# import earthquake data and plot example

library(data.table)
library(rnaturalearth)
library(rnaturalearthdata)
library(sp)
library(sf)
library(rgdal)
library(ggspatial)
library(maps)
library(ggplot2)

setwd('C:/Users/claudio/Documents/NR/ProjectPointProcess/ModelEvaluation/simstudy/Earthquakes/')


dt = fread('CaliforniaEQData.csv',sep = ' ')

setnames(dt,c('date','time','lat','lon','depth','mag','magt','nst','gap','clo','RMS','src','ID'))

dt[,date := as.Date(date)]

### attach utm coordinates ###

# most of southern Cal is zone 11, most of northern is 10:
crs_utm = sp::CRS('+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
crs_lonlat = sp::CRS("+proj=longlat +datum=WGS84")   

#create SpatialPoints object
sp_lonlat = data.frame(x = dt[,lon], y = dt[,lat])
names(sp_lonlat) = c('x','y') # required if vectors are provided
sp::coordinates(sp_lonlat) = c('x','y')
sp::proj4string(sp_lonlat) = crs_lonlat

sp_utm = sp::spTransform(sp_lonlat,crs_utm)
utmcs = as.data.table(coordinates(sp_utm))
setnames(utmcs,c('x','y'))

dt=  data.table(dt,utmcs)


dt_new = dt[,.(date,lat,lon,x,y,depth,mag)]

### plot ###

# get earth data as sf:

US = ne_countries(scale = 10,country = 'United States of America', returnclass = "sf")

states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))

lon_win = c(dt[,min(lon)] - 0.1,dt[,max(lon)] + 0.1)
lat_win = c(dt[,min(lat)] - 0.1,dt[,max(lat)] + 0.1)

pp = ggplot(data = US) + geom_sf(data = states,fill = 'antiquewhite') +
                         theme(panel.grid.major = element_line(color = gray(.5), 
                                                               linetype = 'dashed', 
                                                               size = 0.5), 
                               panel.background = element_rect(fill = 'aliceblue')) +
                         coord_sf(xlim = lon_win, ylim = lat_win, expand = FALSE)

pp = pp + geom_point(mapping = aes(x = lon,y = lat,size = mag),data = dt[year(date) > 2010],alpha = 0.5)

pp = pp + ggtitle('Earthquakes 2011 to today') + 
    xlab('Longitude') + ylab('Latitude') 


print(pp)

ggsave('EQ2011+.pdf')

save(dt_new,file = 'eqs1967+.RData')
