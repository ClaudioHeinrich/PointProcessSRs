
# This script runs permutation tests for the final paper version.

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

# remove the diagonal

#dt = dt[tr_win != val_win]


##### get bootstrap estimates #####

R = 500 # how often do you want to resample?

years = unique(dt[,year])
models = unique(dt[,model])
ests = unique(dt[,estimator])

dt_boot = data.table(expand.grid(year = years,
                                 model = models,
                                 est = ests,
                                 ind = 1:R))



boot_fct = function(data,indices)
{
  rv = NULL
  for(i in 1:ncol(data))
  {
    rv[i] = mean(unlist(lapply(X = data[indices,..i],FUN = as.numeric)))
  }
  return(rv)
}

for(yy in years)
{
  for(mm in models)
  {
    print(paste0('year = ',yy,', model = ',mm))
    for(ee in ests)
    {
      data = dt[year == yy & model == mm & estimator == ee,.(scores)]
      
      boot_samples = boot(data, statistic = boot_fct, R)
      
      dt_boot[year == yy & model == mm & est == ee,mean_score:= boot_samples$t]
      
    }
  }
}



### plotting ###

## finetune naming: ##

dt_boot[model == 'VarGamma',model := 'vG']
dt_boot[model == 'MatClust',model := 'Mat']
dt_boot[model == 'Thomas',model := 'Th']
dt_boot[model == 'Cauchy',model := 'Cau']


# specify font sizes

library(ggplot2)

theme_set(theme_bw(base_size = 28))



for(yy in years)
{
  
  for(ee in ests)
  {
    
    sc_p = ggplot(data = dt_boot[est == ee & year == yy],
                  aes(x = model, y = mean_score, color = model)) + geom_boxplot(show.legend = FALSE)
    
    # add custom title
    
    sc_p = sc_p + labs(title = paste0(yy),
                       x = '',
                       y = ifelse(yy == 1978,yes = 'mean K-function score',no = '')) 
    
    
    pdf(paste0(plot_dir,'mean_',ee,'_score_y',yy,'_new.pdf'))
    print(sc_p)
    dev.off()
  }
}

