
#### compute the kernel estimator score and the K-function score ####

rm(list = ls())

library(spatstat)
library(RandomFields)
library(data.table)

set.seed(102030)

setwd('~/NR/ProjectPointProcess/PointProcessSRs/ScoringRules/')
load('setup_r1.RData')

#######################

N = 100 # samples to compute for CRPS


# K-score

eval_points = seq(0,2.5,by = 0.05) # points where the K estimator is to be evaluated. Needs to take very specific form due to the way K_hat is coded in spatstat

K_hat = function(dat){return(as.function(Kinhom(dat))(eval_points))}

NEst_Kscore = computeNEst(mod = models,est = K_hat,N = N)

K_scores = get_crps(dt = NEst_Kscore,models = models)

# get permutation scores:

aem_Ksc = get_crps(dt = NEst_Kscore,models = models,ret_abserrmat = TRUE)

K_sc_eb = get_crps_errbars(aem_Ksc,models)

perm_sc_K = list()

for(i in 1:length(models))
{
  print(i)
  mod = models[i]
  tmn = colnames(K_sc_eb)[(i-1)*length(models) + i]
  true_scs = as.vector(K_sc_eb[,..tmn])
  for(j in 1:length(models))
  {
    pmn = colnames(K_sc_eb)[(i-1)*length(models) + j]
    pred_scs = as.vector(K_sc_eb[,..pmn])
    
    perm_sc_K = c(perm_sc_K,permutation_test_difference(unlist(true_scs),unlist(pred_scs))$p_val)
    names(perm_sc_K)[(i-1)*length(models) + j] = paste0(tmn,' vs ',pmn)
  }
}

perm_sc_K = data.table(matrix(unlist(perm_sc_K),nrow = length(models)))
setnames(perm_sc_K,paste0('obs_',models))


##### Intensity score ####

lambda_hat = function(dat){return(as.matrix(density(dat)))}

NEst_lambdascore = computeNEst(mod = models,est = lambda_hat,N = N,silence = FALSE)

lambda_scores = get_crps(dt = NEst_lambdascore,models = models)


# get permutation scores:

aem_lsc = get_crps(dt = NEst_lambdascore,models = models,ret_abserrmat = TRUE)

l_sc_eb = get_crps_errbars(aem_lsc,models)

perm_sc_l = list()

for(i in 1:length(models))
{
  print(i)
  mod = models[i]
  tmn = colnames(l_sc_eb)[(i-1)*length(models) + i]
  true_scs = as.vector(l_sc_eb[,..tmn])
  for(j in 1:length(models))
  {
    pmn = colnames(l_sc_eb)[(i-1)*length(models) + j]
    pred_scs = as.vector(l_sc_eb[,..pmn])
    
    perm_sc_l = c(perm_sc_l,permutation_test_difference(unlist(true_scs),unlist(pred_scs))$p_val)
    names(perm_sc_l)[(i-1)*length(models) + j] = paste0(tmn,' vs ',pmn)
  }
}

perm_sc_l = data.table(matrix(unlist(perm_sc_l),nrow = length(models)))
setnames(perm_sc_l,paste0('obs_',models))

##### Intensity score without edge correction ####

lnec_hat = function(dat){return(as.matrix(density(dat,edge = FALSE)))}

NEst_lnecscore = computeNEst(mod = models,est = lnec_hat,N = N,silence = FALSE)

lnec_scores = get_crps(dt = NEst_lnecscore,models = models)


# get permutation scores:

aem_lsc = get_crps(dt = NEst_lnecscore,models = models,ret_abserrmat = TRUE)

l_sc_eb = get_crps_errbars(aem_lsc,models)

perm_sc_lnec = list()

for(i in 1:length(models))
{
  print(i)
  mod = models[i]
  tmn = colnames(l_sc_eb)[(i-1)*length(models) + i]
  true_scs = as.vector(l_sc_eb[,..tmn])
  for(j in 1:length(models))
  {
    pmn = colnames(l_sc_eb)[(i-1)*length(models) + j]
    pred_scs = as.vector(l_sc_eb[,..pmn])
    
    perm_sc_lnec = c(perm_sc_lnec,permutation_test_difference(unlist(true_scs),unlist(pred_scs))$p_val)
    names(perm_sc_lnec)[(i-1)*length(models) + j] = paste0(tmn,' vs ',pmn)
  }
}

perm_sc_lnec = data.table(matrix(unlist(perm_sc_lnec),nrow = length(models)))
setnames(perm_sc_lnec,paste0('obs_',models))



# generate data table containing the data for the plotting

score_dt = as.data.table(expand.grid(type = c('lambda','K','lnec'), obs = models, pr = models))

for(t in c('lambda','K','lnec'))
{
  for(o in models)
  {
    for(p in models)
    {
      o_ind = which(models == o)
      p_ind = which(models == p)
      
      if(t == 'lambda') sc_val = lambda_scores[o_ind,get(p)]
      if(t == 'K') sc_val = K_scores[o_ind,get(p)]
      if(t == 'lnec') sc_val = lnec_scores[o_ind,get(p)]
      
      score_dt[type == t & obs == o & pr == p, score := sc_val]
    }
  }
}



# save 

save(score_dt,K_scores,NEst_Kscore,perm_sc_K,
    lambda_scores,NEst_lambdascore,perm_sc_l,
    lnec_scores,NEst_lnecscore,perm_sc_lnec,
     file = 'scores_r1.RData')




