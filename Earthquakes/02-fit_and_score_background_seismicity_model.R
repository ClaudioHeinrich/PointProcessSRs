# estimate model

rm(list = ls())

library(data.table)
library(rnaturalearth)
library(rnaturalearthdata)
library(sp)
library(sf)
library(rgdal)
library(ggspatial)
library(maps)
library(ggplot2)
library(spatstat)


setwd('C:/Users/claudio/Documents/NR/ProjectPointProcess/ModelEvaluation/simstudy/Earthquakes/')

load('eqs1967+.RData')

dt_new

# scale x and y:

dt_new[,x := x/1000][,y:=y/1000]

# 
# xr = dt_new[,range(x)]
# yr = dt_new[,range(y)]
# scale_factor = 1/(sqrt((xr[2]-xr[1])*(yr[2]-yr[1])))
# 
# transform_coords = function(x,y)
# {
#   return(data.table(x = x/(xr[2]-xr[1]), y = y/(yr)))
# }

source('../TreeData/functions/summaryScore.R')

window = owin(xrange = dt_new[,range(x)],yrange = dt_new[,range(y)])


get_scores = function(training_years = 1970:2000,
                      val_years = 2001,
                      sig = 8, # in km
                      a = 0.5,
                      ...)
{
  

  # 1. get earthquakes as ppp object
  
  dt_temp = dt_new[year(date) %in% training_years]
  
  utm_cor = dt_temp[,cor(x,y)]
  
  pp = ppp(dt_temp[,x],dt_temp[,y],window = window)
  
  # 2. create kernel smoothed image and density as convex combination:
  
  dd = density(pp,varcov = matrix(c(sig^2,sig^2*utm_cor,sig^2*utm_cor,sig^2),ncol = 2))
  #plot(dd)
  
  nu = pp$n / area(pp$window)
  
  dd_new = a*dd + (1-a) * nu
  
  # scale to number of years:
  
  dd_pred = length(val_years) * dd_new/length(training_years)
  
  sim = function(){return(rpoispp(lambda = dd_pred))}
  
  obs_dt = dt_new[year(date) %in% val_years]
  obs_ppp = ppp(obs_dt[,x],obs_dt[,y], xrange = dt_temp[,range(x)],yrange = dt_temp[,range(y)])
  
  scores = summaryScore(obs_ppp,predictive.model = sim,...)
  
  return(scores)
}


#### test run ####

validation_years = dt_new[,unique(year(date))]
as = seq(0.3,0.9,length.out = 7)
sigs = c(4,8,12,16)

res_dt = as.data.table(expand.grid(obs_year = validation_years,a = as,sig = sigs,est = c('K','intensity')))

## specify evaluation points for K function: ##

K_range = sqrt(area(window))/7.5 # a bit hacked: This ends somewhere not to far from spatstats recommendation for the configuration I tried...

K_eval_locs = seq(0,K_range,length.out = 10)


for(ii in 1:res_dt[,.N])
{
  print(paste0(ii,'/',res_dt[,.N]))
  obs_year = res_dt[ii,obs_year]
  sig = res_dt[ii,sig]
  est = res_dt[ii,est]
  if(est == 'K')
  {
    res_dt[ii,score:= get_scores(training_years = validation_years[validation_years != obs_year],
                                 val_years = obs_year,
                                 a = res_dt[ii,a],
                                 sig = sig,
                                 N_crps = 50,
                                 est = est,
                                 eval_points = K_eval_locs)]
  }
  if(est == 'intensity')
  {
    res_dt[ii,score:= get_scores(training_years = validation_years[validation_years != obs_year],
                                 val_years = obs_year,
                                 a = res_dt[ii,a],
                                 sig = sig,
                                 N_crps = 50,
                                 est = est
                                 )]
  }
  
  
}

save(res_dt,file = 'results.RData')
