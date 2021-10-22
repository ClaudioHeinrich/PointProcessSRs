
# This script fits several models to the data and compute scores.
# It generates the file scores.RData.
# This takes some time, so it's best to run this on a server in parallel, see below.
# For each year, the trees in a single plot are used to fit the predictive model and this model is evaluated against all other plots of the same year.
# This is repeated for each plot, i.e. each plot is used once to fit the predictive model for all other plots.


rm(list = ls())

set.seed(102030)

library(spatstat)
library(spatstat.utils)
library(data.table)
library(parallel)

#setwd('~/NR/ProjectPointProcess/ModelEvaluation/simstudy/TreeData')

setwd('~/R/PointProcess/Trees') # for working on server

data_dir = './Data'
dir.create(data_dir,showWarnings = FALSE)
plot_dir = './Figures'


#mc_cores = 1 # for running it on Windows



# mc_cores = 8

#### directories ####

data_dir = './Data/'
function_dir = './functions/'




##############################
### get tree data in shape ###
##############################

load(paste0(data_dir,"trees.Rdata"))

win = owin(xrange = c(0,6),yrange = c(0,6))

### need number of points in each plot - take only those that
### fall within [0,6]x[0,6] and are of species 1
n.x.78 <- NULL
x.78 <- list()
y.78 <- list()
for(i in 1:8) {
  ind <- which(trees[[i]]$species==1 & trees[[i]]$x>0 & trees[[i]]$x<6 &
               trees[[i]]$y>0 & trees[[i]]$y<6)
  x.78[[i]] <- trees[[i]]$x[ind]
  y.78[[i]] <- trees[[i]]$y[ind]
  n.x.78 <- c(n.x.78,length(x.78[[i]]))
}

n.x.90 <- NULL
x.90 <- list()
y.90 <- list()
for(i in 1:8) {
  ind <- which(trees[[i]]$species==1 & trees[[i]]$x>0 & trees[[i]]$x<6 &
               trees[[i]]$y>0 & trees[[i]]$y<6 & !is.na(trees[[i]]$D90))
  x.90[[i]] <- trees[[i]]$x[ind]
  y.90[[i]] <- trees[[i]]$y[ind]
  n.x.90 <- c(n.x.90,length(x.90[[i]]))
}

n.x.09 <- NULL
x.09 <- list()
y.09 <- list()
for(i in 1:8) {
  ind <- which(trees[[i]]$species==1 & trees[[i]]$x>0 & trees[[i]]$x<6 &
               trees[[i]]$y>0 & trees[[i]]$y<6 & !is.na(trees[[i]]$D09))
  x.09[[i]] <- trees[[i]]$x[ind]
  y.09[[i]] <- trees[[i]]$y[ind]
  n.x.09 <- c(n.x.09,length(x.09[[i]]))
}


# # plot 
# 
# par(mfrow=c(4,3),mex=0.25,mar=c(1,1,1,1)+0.01)
# for(i in 1:4) {
#   plot(x.78[[i]],y.78[[i]],type="p",pch=16,cex=0.5,
#        xlab="",ylab="",xlim=c(0,6),ylim=c(0,6),xaxs="i",
#        yaxs="i",axes=FALSE)
#   box()
#   plot(x.90[[i]],y.90[[i]],type="p",pch=16,cex=0.5,
#        xlab="",ylab="",xlim=c(0,6),ylim=c(0,6),xaxs="i",
#        yaxs="i",axes=FALSE)
#   box()
#   plot(x.09[[i]],y.09[[i]],type="p",pch=16,cex=0.5,
#        xlab="",ylab="",xlim=c(0,6),ylim=c(0,6),xaxs="i",
#        yaxs="i",axes=FALSE)
#   box()
# 
# }
# 
# 
# 
# years_char = c('78','90','09')





################################################################################
### fit model to 4 of the eight observations and predict the 4 other with it ###
################################################################################


estimators = c('intensity','K') # estimators you want to try

models = c('Pois','LGCP','Thomas','MatClust','VarGamma','Cauchy') # models you want to try

years = c(1978,1990,2009)

N_crps = 100 # MC samples generated from each predictive distribution

### get functions for score computation and multiple model fitting ###

source(paste0(function_dir,'summaryScore.R'))


source('mkppm.R')


#### transform tree data to list of ppps ####

lpp = list()

years_char = c('78','90','09')

for(y_ind in 1:3)
{
  x_data = get(paste0('x.',years_char[y_ind]))
  y_data = get(paste0('y.',years_char[y_ind]))
  
  l_temp = list()
  for(i in 1:8)
  {
    pp = ppp(x_data[[i]],y_data[[i]],xrange=c(0,6),yrange=c(0,6))
    l_temp[[i]] = pp
  }
  lpp[[y_ind]] = l_temp
}


names(lpp) = years


###################################################




get_scores = function(mod,est = 'K', yy = 1978, ntrain = 4, N_crps = 100,lpp)
{
  # get all subsets of ntrain indices out of 8
  scores = c()
    
  combs = combinat::combn(8,ntrain)
  
  for(col_ind in 1:dim(combs)[2])
  {
    print(col_ind)
    train_ind = combs[,col_ind] 
    val_ind = (1:8)[-train_ind]
    
    
    
    y_str = as.character(yy)
    
    Tr = lpp[[y_str]][train_ind]
    
    lambda_hat = 0
    for(i in 1:ntrain)
    {
      lambda_hat = lambda_hat + Tr[[i]]$n / (ntrain*36)
    }
    
    if(mod == 'Pois')
    {
      sim = function(){return(rpoispp(lambda = lambda_hat,win = Tr[[1]]$window))}
    }else{
    
    K_hat = mean_Kest(Tr)
    
    a = mkppm(X = Tr[[1]], Kfv = K_hat,lambda = lambda_hat,clusters = mod)
    
    # function to simulate predictive model:
    sim = function(){return(simulate(a,nsim = 1)[[1]])}
    }
    
    for(obs_ind in val_ind)
    {
    Obs = lpp[[y_str]][obs_ind]
    
    sc = summaryScore(Obs,
                      predictive.model = sim,
                      est = est,
                      N_crps = N_crps,
                      verbose = FALSE)
    scores = c(scores,sc)
    }
  }
  
  ret_val = data.table(year = yy,model = mod,estimator = est,scores = scores)
  
  return(ret_val)
}

# run in parallel:

dt = data.table()

for(est in estimators)
{
  for(yy in  years)
  {
  print(paste0('estimator = ',est,', year = ',yy))
  dt_temp = mclapply(X = models,
                     FUN = function(mod){return(get_scores(mod = mod,
                                                           est = est,
                                                           yy = yy,
                                                           ntrain = 4,
                                                           N_crps = N_crps,
                                                           lpp = lpp))},
                     mc.cores = length(models), mc.silent = FALSE)
  
  dt_temp = rbindlist(dt_temp)
  
  dt = rbindlist(list(dt,dt_temp))
  }
}


# save:

save(dt,file = paste0(data_dir,'scores_multobsfitting.RData'))

