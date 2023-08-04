
#### compute the kernel estimator score and the K-function score ####

rm(list = ls())

library(spatstat)
library(RandomFields)
library(data.table)
library(ggplot2)
library(latex2exp)

set.seed(102030)

setwd('~/NR/ProjectPointProcess/PointProcessSRs/ScoringRules/')
load('setup_r1.RData')

#######################

N = 10000# samples to compute for CRPS

# K-score

eval_points = c(0.5,1) # points where the K estimator is to be evaluated. Needs to take very specific form due to the way K_hat is coded in spatstat
models= models[1:3]

K_hat = function(dat){return(as.function(Kinhom(dat))(eval_points))}


# get the length of the output of est:
l_est = length(est(sim(mod[1])))

#initialize data table
dt = data.table()

for(mm in mod){
  
  #simulate and fill:
  for(i in 1:N)
  {
    if(i%%100 == 0) cat(paste0('\r',mm,': ',round(100 *i/N,1),'%'))
    
    dat = sim(mm)
    dt = rbindlist(list(dt, data.table(model = mm, sample = i,r = eval_points,K = est(dat)))) 
  }
  
}

theoretical_values = data.table(r = eval_points,val = pi*eval_points^2)

dt[,rlab := paste0('r = ',r)]
theoretical_values[,rlab := paste0('r = ',r)]

theme_set(theme_bw(base_size = 16))

pp = ggplot(dt) + 
  geom_density(aes(color = model,fill = model,x = K),alpha = 0.4) + 
  scale_x_continuous(TeX('$\\widehat{K}(r)$')) + 
  facet_grid(cols = vars(rlab))

pp = pp + geom_vline(data= theoretical_values,aes(xintercept = val))
pp

plot_dir = c('../../ModelEvaluation/Figures/review1/')

ggsave(paste0(plot_dir,'K_hat_densities.jpg'),width = 10,height = 5)
