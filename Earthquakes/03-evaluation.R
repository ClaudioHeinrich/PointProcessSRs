

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


load( 'results.RData')

res_dt[,a := as.factor(a)]

pp = ggplot(res_dt[est == 'K']) + geom_boxplot(aes(a,score)) + facet_grid(cols = vars(sig))
pp = pp +ggtitle('K-scores')
pp
ggsave('K_scores.pdf')


pp = ggplot(res_dt[est == 'intensity']) + geom_boxplot(aes(a,score)) + facet_grid(cols = vars(sig))
pp = pp +ggtitle('intensity scores')
pp
ggsave('intensity_scores.pdf')



### get mean ###

means = res_dt[,mean(score),.(est,sig,a)]

means[,sig := as.factor(sig)]
means[,a := as.numeric(as.character(a))]

setnames(means,'V1','ms')
pp = ggplot(means[est == 'K']) + geom_line(aes(x = a,y = ms,colour = sig,linetype = sig))
pp = pp + ggtitle('mean K score')
pp

ggsave('mean_K_score.pdf')

pp = ggplot(means[est == 'intensity']) + geom_line(aes(x = a,y = ms,colour = sig,linetype = sig))
pp = pp + ggtitle('mean intensity score')
pp

ggsave('mean_intensity_score.pdf')
