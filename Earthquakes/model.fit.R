
### Earthquakes from Entiat, Leavenworth, and Grand Coolee area
### 1969 through Oct 1986

### x has seven columns:
### time (in days from 01011970) x y depth mag yrmodyhrmn sec

library(data.table)
library(SAPP)
library(spatstat)
library(ETAS)

setwd('C:/Users/claudio/Documents/NR/ProjectPointProcess/ModelEvaluation/simstudy/Earthquakes/')

x <- read.table("cal.data")

eq_dt = as.data.table(x)

eq_dt[,date:= as.Date(paste0('19',V1,'-',V2,'-',V3),format = "%Y-%m-%d")]
eq_dt[,time:=paste0(V4,':',V5,':',V6)]


eq_dt[,paste0('V',1:6) := NULL]
setnames(eq_dt,paste0('V',7:9),c('lat','long','mag'))

eq_dt[,c('hour','min','sec') := NULL]

eq_dt[,long := -long]

setcolorder(eq_dt,c('date','time','long','lat','mag'))

##### plot ######

cal_map <- map_data("state",region = 'california')

# bounding box 

xlim <- c(eq_dt[,min(long)],eq_dt[,max(long)])
ylim =  c(eq_dt[,min(lat)],eq_dt[,max(lat)])

# with tolerance for plotting:

k <- 0.05

xlim_pl =  c(eq_dt[,min(long)-k],eq_dt[,max(long)+k])
ylim_pl =  c(eq_dt[,min(lat)-k],eq_dt[,max(lat)+k])




pp = ggplot(cal_map ) +a
geom_polygon(aes(x = long, y = lat, group = group),fill="white", colour = "black") + 
  geom_point(data = eq_dt,aes(x = long,y = lat,size = mag),shape = 1)  + coord_map(xlim = xlim_pl,ylim = ylim_pl) 

pp = pp + theme(
  panel.background = element_rect(fill = "lightblue",
                                  colour = "lightblue",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank())

pp
#####################


c0 = catalog(eq_dt,
             time.begin = '1962-01-01',study.start = '1963-01-01',
             lat.range = ylim,long.range = xlim)

fit0 = etas(c0)

par0 = round(fit0$param,2)

pred.years = 1968:1981

for(yy in pred.years)
{
  cc = catalog(eq_dt,
               time.begin = '1962-01-01',study.start = '1963-01-01',study.end = paste0(yy+1,'-01-01'),
               lat.range = ylim,long.range = xlim)
  
  fit = etas(cc,param0 = par0)
}


########## figure out cRS:

library(sp)
library(rgdal)

crs_utm = sp::CRS('+proj=utm +zone=10 +lon0=-121 +lat0=47.5 +units=m +no_defs')

coords = eq_dt[,.(x,y)]
coordinates(coords) = c('x','y')
proj4string(coords) = crs_utm
plot(coords)

# convert to lon lat:
crs_lonlat = sp::CRS("+proj=longlat +datum=WGS84")   
sp_utm = data.frame(lon = coords$x, lat =coords$y)
names(sp_utm) = c('lon','lat') 
sp::coordinates(sp_utm) = c('lon','lat')
sp::proj4string(sp_utm) = crs_utm

lonlatcoords = sp::spTransform(sp_utm, crs_lonlat)

######### check whether this is the right projection!!! ########

library(maps)

map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-60, 90), mar=c(0,0,0,0))
points(coords)

#########################
