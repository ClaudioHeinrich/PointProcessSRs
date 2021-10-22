
### run permutation tests for score differences ###


### setup ###

rm(list = ls())

library(spatstat)

library(data.table)

library(ggplot2)


setwd('~/NR/ProjectPointProcess/ModelEvaluation/simstudy/KernelVsLog')

set.seed(102030)

load('setup.RData')


# load scores:

load(file = paste0('scores.RData'))


############ function for permutation tests ##################

#' Run a permutation test of the pairwise difference between two vectors of numbers
#' @param a Vector, the scores from one method
#' @param b Vector, the scores from some other method
#' @param N Integer, the size of the permutation distribution
#' @return A list with the mean of the difference and the permutation distribution of that difference and the p-value
#' @examples
#' N = 1e2
#' trend  = 1:N
#' a = trend + .01 + rnorm(N, .001)
#' b = trend - .01 + rnorm(N, .001)
#' l = permutation_test_difference(a,b)
#' q = sum(l$D <= l$d_bar) / length(l$D)
#' @author Alex,Claudio
#' 
#' @export
#' 
permutation_test_difference = function(a,
                                       b,
                                       N = 5e3){
  n = length(a)
  d = a - b
  d_bar = mean(d)
  D = NULL
  for(i in 1:N){
    swap = rbinom(n,1,0.5)
    w_swap = which(swap == 1)
    d_i = d
    d_i[w_swap] = -d_i[w_swap]
    D[i] = mean(d_i)
  }
  
  p_val = sum(d_bar > sort(D))/N + sum(d_bar == sort(D))/(2*N)
  
  return(list(d_bar = d_bar, D = D,p_val = p_val))
}


#################################

# run permutation tests: consider means of N_pt scores (sampled from within dt) and see how often out of R times a permutation test at level alpha rejects. This divided by R approximates the power of the test#


obs_vec = 1:4

dt = dt[N_sam == 1 & obs %in% obs_vec]

N_tot = dt[type == 'log' & pr == 1 & obs == obs_vec[1],.N] # total number of available score evaluation

N_vec_pt = seq(5,50,5) # vector containing all numbers of N_pt for which this is done

R = 500 # number of permutation tests run for each model comparison

ests = c('kernel','log') # considered estimators

R_pt = 500 # number of resamples considered in each permutation test

########

# initialize data table:

p_vals = as.data.table(expand.grid(est = ests,pr_mod = 1:4,N_pt = N_vec_pt,ind = 1:R,obs_mod = obs_vec))

for(mod_obs in  obs_vec) # models that resembles observation
{
  print(paste0('observation model ',mod_obs,':'))
  # run stuff:
  
  for(NN_pt in N_vec_pt )
  {
    print(NN_pt)
    for(tt in ests)
    {
      print(paste0(tt,'...'))
      for(mod_pr in  1:4)
      {
        dt_obs = dt[pr == mod_obs & obs == mod_obs & type == tt]
        dt_pred = dt[pr == mod_pr & obs == mod_obs & type == tt]
        for(ii in 1:R)
        {
          if(ii %% 100 == 0) {print(paste0(ii + R*(mod_pr-1),'/',4*R))}
          ran_inds = sample(N_tot,NN_pt,replace = FALSE)
          
          true_sc = dt_obs[ran_inds,score]
          pred_sc = dt_pred[ran_inds,score]
          
          pp_val = permutation_test_difference(a=true_sc, b = pred_sc,N = R_pt)$p_val
        
          p_vals[est == tt & pr_mod == mod_pr & ind == ii & N_pt == NN_pt & obs_mod == mod_obs, p_val:= pp_val]
        }
      }
    }
  }
}
  
### save ###

save.image('pvalspermtest.RData')

