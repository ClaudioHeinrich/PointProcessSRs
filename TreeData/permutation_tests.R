# This script runs permutation tests and plots the reslts.

rm(list = ls())

set.seed(102030)

library(spatstat)
library(data.table)
library(boot)
library(ggplot2)

setwd('~/NR/ProjectPointProcess/ModelEvaluation/simstudy/TreeData')

# setwd('~/R/PointProcess/Trees') # for working on server


#### directories ####

data_dir = './Data/'
function_dir = './functions/'
plot_dir = '../../Figures/'

# get data:

load(paste0(data_dir,'scores_multobsfitting.RData'))
load(paste0(data_dir,'trees.RData'))

source(paste0(function_dir,'PermutationTests.R'))


#####################

years = dt[,unique(year)]
models = dt[,unique(model)]
ests = 'K'


pt_dt = as.data.table(expand.grid(est = ests, model1 = models, model2 = models, year = years))

for(yy in years)
{
  for (mod1 in models)
  {
    for(mod2 in models)
    {
      print(c(yy,mod1,mod2))
      for(estim in ests)
      {
        a = dt[year == yy & model == mod1 & estimator == estim,scores]
        b = dt[year == yy & model == mod2 & estimator == estim,scores]
        
        pt = permutation_test_difference(a,b,N = 5e3)
        
        frac = (rank(c(pt$d_bar,pt$D))[1]-1)/(length(pt$D)) # get p_value
        
        pt_dt[year == yy & model1 == mod1 & model2 == mod2 & est == estim, p_val := frac]
        
      }
    }
  }
}

################
pt_dt[model1 == 'Cauchy']

pt_dt[model1 == 'LGCP']

####### plot results #############

library(ggplot2)

# specify font sizes

title_size = 20
labels_size = 16
ticks_size = 12

####################

for(yy in years)
{
  for(ee in ests)
  {
  p = ggplot(data = pt_dt[est == ee & year == yy & model1 != model2,],
             aes(x = model1,y = p_val,color = model2,shape = model2)) + geom_point(alpha = .75,size = 3)
  
  p = p+ geom_hline(yintercept = c(0.05,0.95),lty = 2)
  
    p = p + labs(title = paste0(ee,'-score, year ',yy),x = 'model 1',y = 'p-value')
  p = p + labs(color = 'model 2',shape = 'model 2') + theme(plot.title = element_text(size = title_size),
                                                            axis.title = element_text(size = labels_size),
                                                            axis.text = element_text(size = ticks_size))
  
  p = p + theme(legend.title = element_text(size = labels_size),
                legend.text = element_text(size = ticks_size))
  
  # plotting
  
  pdf(paste0(plot_dir,'perm_test_',ee,yy,'.pdf'))
    print(p)
  dev.off()
  
  }
}
