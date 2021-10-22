

rm(list = ls())

library(data.table)
library(spatstat)

setwd('~/NR/ProjectPointprocess/ModelEvaluation/simstudy/LTest')

set.seed(102030)

#specify different models and consider each as true and as predictive distribution and perform the L-test conditional on having the same number of earthquakes as in the observation

# consider for this test processes instead of binning 

#### specify window ###

limx = c(0,10)
limy = c(0,10)

x_range = max(limx) - min(limx)
y_range = max(limy) - min(limy)

win = owin(xrange = limx,yrange = limy)

vol = x_range * y_range

### set parameters for models ###


model_names = c('qiP','liP1','liP2','hP')

nexps = rep(100,4)

avec = c(0.1,0.2,0.3,0.4)

#specify intensities: 
qiPlambda = function (x,y) return((x/5-1)*avec[1] + 1)
qiPlambda_int = pracma::integral2(qiPlambda,
                                   xmin = limx[1],xmax = limx[2],
                                   ymin = limy[1],ymax = limy[2])$Q


liP1lambda = function (x,y) return((x/5-1)*avec[2] + 1)
liP1lambda_int = pracma::integral2(liP1lambda,
                                 xmin = limx[1],xmax = limx[2],
                                 ymin = limy[1],ymax = limy[2])$Q


liP2lambda = function (x,y) return((x/5-1)*avec[3] + 1)
liP2lambda_int = pracma::integral2(liP2lambda,
                                   xmin = limx[1],xmax = limx[2],
                                   ymin = limy[1],ymax = limy[2])$Q

hPlambda = function (x,y) return((x/5-1)*avec[4] + 1)
hPlambda_int = pracma::integral2(hPlambda,
                                   xmin = limx[1],xmax = limx[2],
                                   ymin = limy[1],ymax = limy[2])$Q





lambdas = list(function(x,y) return(nexps[1] *qiPlambda(x,y)/qiPlambda_int),
                  function(x,y) return(nexps[2] *liP1lambda(x,y)/liP1lambda_int),
                  function(x,y) return(nexps[3] *liP2lambda(x,y)/liP2lambda_int),
                  function(x,y) return(nexps[4] *hPlambda(x,y)/hPlambda_int))


# create multi-level list

models = list()
for(i in 1:length(model_names) )
{
  mod_i = list(nexp = nexps[i],lambda = lambdas[[i]])
  models = c(models,list( mod_i))
}

names(models) = model_names


##### functions #####


# function to simulate point patterns with these parameters

simPP = function(i)
{
  dat = rpoispp(lambda = models[[i]]$lambda,win = win)
  return(dat)
}

simPPn = function(i,n)
{
  dat = rpoint(n = n, f = models[[i]]$lambda,win = win)
  return(dat)
}

# function to get loglikelihood for ppp under model i

lll = function(ppp,i)
{
  l = models[[i]]$lambda
  nexp = models[[i]]$nexp
  if(is.numeric(l)){
    ll = vol - nexp  + ppp$n * log(l)
  }else{
    lambdas = l(ppp$x,ppp$y)
    ll = vol - nexp  + sum(log(lambdas))
  }
  return(ll)
}



# compute quantile score for observation y and predictive model i

gamma = function(y,i,MC_sam = 10)
{
  n = y$n
  
  LLH = NULL
  for(j in 1:MC_sam)
  {
    X = simPPn(i,n)
    LLH[j] = lll(X,i)
  }
  LLH_obs = lll(y,i) 
  
  return( mean(LLH < LLH_obs))
}




######################

N = 500

alpha = 0.05

dt = as.data.table(expand.grid(model = 1:4,observ = 1:4))

for(mod in 1:length(models))
{
  for(obs in 1:length(models))
  {
    print(c(mod,obs))
    
    gammas = NULL
    for(j in 1:N)
    {
      print(j)
      y = simPP(i = obs)
      gammas[j] = gamma(y,i = mod,MC_sam = N)
    }
    dt[model == mod & observ == obs,gamma := mean(gammas)]
    dt[model == mod & observ == obs, fails := mean(gammas <= alpha)]
    dt[model == mod & observ == obs, fails_ts := mean(gammas <= alpha/2 | gammas >= 1-alpha/2)]
  }
}


save(dt,N,file = 'SMtest.RData')


###################################

