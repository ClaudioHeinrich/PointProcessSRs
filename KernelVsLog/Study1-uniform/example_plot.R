
### plot examples for all models ###


rm(list = ls())

library(spatstat)

library(data.table)


setwd('~/NR/ProjectPointProcess/ModelEvaluation/simstudy/ScoringRules/')
load('setup.RData')

plot_dir = '../../Figures/'


par0 = list(mar = c(1,1,1,1),cex = 1.5)

par0$mar = c(1,1,1,1)

par0$cex = 1.5


lambdas = list(0.5,0.4,0.6,function(x,y){return(0.05*x + 0.25)})


for(ind in 1:4)
{
  pdf(paste0(plot_dir,'P',ind,'.pdf'))
  
  par(par0)
  
  lambda = lambdas[[ind]]
  pp = rpoispp(lambda = lambda,
               win = owin(xrange = c(0,10), yrange = c(0,10)))
  
  par('cex' = 3)
  plot(pp,pch = 20,main = '')
  
  dev.off()
}

