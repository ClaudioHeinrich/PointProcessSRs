
### weird example for the log score ###

# We consider 2 Poisson models here, for both of them the intensity is a homogeneous centered Gaussian kernel,
# but they slightly differ in number of points. Even though we consider the mean log-likelihood over N samples (usually N = 100),
# the log-likelihood frequently prefers the wrong model.
# You can change the seed a couple of times and rerun everything to see this.


rm(list = ls())

library(spatstat)
library(ggplot2)
library(data.table)
library(pracma)

set.seed(1)

# set some variables:

N = 100 # the mean over how many logscores is computed

limx = c(-5,5) # range of considered spatial area
limy = c(-5,5)

win = owin(xrange = limx,yrange = limy)


### set parameters for models ###

n_mod = 2 # two models are considered.

# model 1 has as intensity the Gaussian kernel with these specifications:
sig1 = 1 
mu1 = c(0,0)
cor1 = 0
nexp1 = 100 # the Gaussian kernel is multiplied by this to get to the intensity, 
            # so this is only approximately the expected number of points (up to 4 digits after comma).
            # It is not exactly the number of points because the Gaussian kernel is not exactly 0 outside the
            # considered window.

#model 2 has slightly more points:
sig2 = 1
mu2 = c(0,0)
cor2 = 0
nexp2 = 110


########### functions for simulating the models and getting loglikelihood: #########

gaussian_kernel = function(x,y,mu=mu1,sig = sig1,cor = cor1)
{
  if(cor == 0)
  {
    exponent = -1/2 * ((x-mu[1])^2/sig^2 + (y-mu[2])^2/sig^2)
    ret_val = 1/(2*pi*sig^2) * exp(exponent) 
  } else { # the following lines are not used in this example script...
    mat = matrix(c(sig^2,cor*sig^2,
                   cor*sig^2,sig^2),
                 ncol = 2)
    mat_inv = solve(mat)
    
    # writing out the matrixmultiplication allows us to pass vectors in x and y:
    sum1_exponent = -1/2 * mat_inv[1,1]*(x - mu[1])^2
    sum2_exponent = - mat_inv[1,2] * (x - mu[1]) * (y - mu[2])
    sum3_exponent = -1/2 * mat_inv[2,2]*(y - mu[2])^2
    exponent = sum1_exponent + sum2_exponent + sum3_exponent
    ret_val = as.numeric((2*pi)^(-1) * det(mat)^(-1/2) * exp(exponent))
  }
  return(ret_val)
}


# get intensities:
lambdas = list()

lambdas[[1]] = function(x,y) return(nexp1 * gaussian_kernel(x,y))
lambdas[[2]] = function(x,y) return(nexp2 * gaussian_kernel(x,y))

# create list of models:

models = list()
for(i in 1:n_mod)
{
  # get the exact number of expected points by integrating the intensity over the window:
  nexp_exact =  pracma::integral2(lambdas[[i]],
                                  xmin = win$xrange[1],
                                  xmax = win$xrange[2],
                                  ymin = win$yrange[1],
                                  ymax = win$yrange[2])$Q
  
  mod_i = list(nexp = nexp_exact,lambda = lambdas[[i]])
  models = c(models,list( mod_i))
}

# function to simulate point patterns from these models:

simPP = function(i)
{
  dat = rpoispp(lambda = models[[i]]$lambda,win = win)
  return(dat)
}

#### uncomment this if you want to plot examples: ####

# for(i in 1:10)
# {
#   pp1 = simPP(1)
#   pp2 = simPP(2)
#   
#   par(mfrow = c(1,2))
#   plot(pp1,main = paste0('sample of mod 1, ',pp1$n,' points'))
#   plot(pp2,main = paste0('sample of mod 2, ',pp2$n,' points'))
# }



#################################

# function to get loglikelihood for point pattern ppp under model i.
# Likelihood-function is considered with respect to unit rate poisson process

lll = function(ppp,i)
{
  vol = volume(win) # not really necessary but also not harming things, I had it in so far so I left it in
  
  l = models[[i]]$lambda
  nexp = models[[i]]$nexp
  if(is.numeric(l)){# this case is not used in this example script...
    ll = vol - nexp  + ppp$n * log(l)
  }else{
    lambdas = l(ppp$x,ppp$y)
    ll = vol - nexp  + sum(log(lambdas))
  }
  return(ll)
}

####################################

# now, consider all four combinations of mod 1 and mod 2 as observation and prediction model 
# and compute N log scores

log_scores = as.data.table(expand.grid(obs_mod = 1:n_mod,
                                       pred_mod = 1:n_mod,
                                       index = 1:N))

for(pred in 1:n_mod)
{
  for(obs in 1:n_mod)
  {
    log_sc = NULL
    for(n in 1:N)
    {
      ppi = simPP(obs) # observation drawn from observation model
      log_sc[n] = -lll(ppi,pred) # negative llh under prediction model
    }
    log_scores[obs_mod == obs & pred_mod == pred,log_score := log_sc]
  }
}

# plot results:

log_scores[,c('obs_mod','pred_mod') := lapply(.SD,as.factor),.SDcols = c('obs_mod','pred_mod')]

pp1 = ggplot(log_scores[obs_mod == 1]) + 
  geom_boxplot(aes(group = pred_mod,y = log_score,fill = pred_mod)) +
  ggtitle('True model 1')
pp1

pp2 = ggplot(log_scores[obs_mod == 2]) + 
  geom_boxplot(aes(group = pred_mod,y = log_score,fill = pred_mod)) +
  ggtitle('True model 2')
pp2


# print mean scores:
print_dt = log_scores[,mean(log_score),by = .(obs_mod,pred_mod)]
setkey(print_dt,obs_mod)
print_dt
