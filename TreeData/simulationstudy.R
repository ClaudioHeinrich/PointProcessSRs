
# This script fits several models to the data and compute scores.
# It generates the file scores.RData.
# This takes some time, so it's best to run this on a server in parallel, see below.
# For each year, the trees in a single plot are used to fit the predictive model and this model is evaluated against all other plots of the same year.
# This is repeated for each plot, i.e. each plot is used once to fit the predictive model for all other plots.


rm(list = ls())

set.seed(102030)

library(spatstat)
library(data.table)
library(parallel)

setwd('~/NR/ProjectPointProcess/ModelEvaluation/simstudy/TreeData')

data_dir = './Data'
dir.create(data_dir,showWarnings = FALSE)
plot_dir = './Figures'


mc_cores = 1 # for running it on Windows


# setwd('~/R/PointProcess/Trees') # for working on server
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
### fit model to 1 of the eight observations and predict the 7 other with it ###
################################################################################


estimators = c('intensity','K','F','G','pcf') # estimators you want to try

models = c('Pois','LGCP','Thomas','MatClust','VarGamma','Cauchy') # models you want to try

years = c(1978,1990,2009)

N_crps = 100 # MC samples generated from each predictive distribution

### get function for score computation ###

source(paste0(function_dir,'summaryScore.R'))


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

#### parallelize along different observations: ####


parallel_score = function (train,mod)
{
  dt_temp = as.data.table(expand.grid(year = years,tr_win = train,val_win = 1:8,est = estimators))
  for(estimator in estimators)
  {
    for(yy in years)
    {
      print(c(estimator,yy,train))
      
        y_str = as.character(yy)
        
        Tr = lpp[[y_str]][[train]]
        
        # fit model:
        
        if(mod == 'Pois')
        {
          a = ppm(Tr ~ 1, interaction = NULL)
        }
        if(mod != 'Pois')
        {
          a = kppm(Tr ~ 1, clusters =  mod)   
        }
        
        # function to simulate predictive model:
        sim = function(){return(simulate(a,nsim = 1)[[1]])}
        
        Obs = lpp[[y_str]]
        
        sc = summaryScore(Obs,
                          predictive.model = sim,
                          est = estimator,
                          N_crps = N_crps,
                          verbose = FALSE)
            
        dt_temp[year == yy & tr_win == train & est == estimator,score := sc]
    }
  }
  return(dt_temp)
}

# run in parallel:

dt = data.table()

for(mod in models)
{
  print(mod)
  dt_temp = mclapply(X = 1:8,
                     FUN = function(train){return(parallel_score(train,mod))},
                     mc.cores = mc_cores, mc.silent = FALSE)
  
  dt_temp = rbindlist(dt_temp)
  
  dt_temp[,model := mod]
  
  dt = rbindlist(list(dt,dt_temp))
}


# save:

save(dt,file = paste0(data_dir,'scores.RData'))

