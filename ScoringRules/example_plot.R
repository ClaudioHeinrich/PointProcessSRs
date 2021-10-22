
### plot examples for all models ###


rm(list = ls())

library(spatstat)

library(data.table)

set.seed(102040)


setwd('~/NR/ProjectPointProcess/ModelEvaluation/simstudy/ScoringRules/')
load('setup.RData')

plot_dir = '../../Figures/'


par0 = list(mar = c(1,1,1,1),cex = 1.5)

par0$mar = c(1,1,1,1)

par0$cex = 1.5



for(mod in models)
{
  pdf(paste0(plot_dir,mod,'.pdf'))
  
  par(par0)
  
  pp = sim(mod)
  
  par('cex' = 3)
  
  plot(pp,main = mod,pch = 20)
  
  dev.off()
}


# The thomas model has a high variance in the number of points - even though the mean of this is 50
# to make it look not too different, rerole until you have between 45 and 55 points:

mod = 'ihT'

pdf(paste0(plot_dir,mod,'.pdf'))

  par(par0)
  par('pch' = 20)
  pp = sim(mod)
  
  if(abs(pp$n - n_exp) > 5)
  {
    while(abs(pp$n - n_exp) > 5)
    {
      pp = sim(mod)
    }
  }
  
  par('cex' = 3)
  
  plot(pp,main = mod,pch = 20)

dev.off()
