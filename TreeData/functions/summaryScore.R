


##### function for computing scores #####

# has an issue currently, in that it norms the W-space to 1, just multiply the score by the area at the end


#' @param obs A ppp or list of ppps containing the observation.
#' @param predictive.model A function taking no argument and returning a simulation of the predictive model.
#' @param est The estimator for the summary statistic the score is based on. Takes 'intensity', 'K', 'F','G','pcf' or any estimator specified as a functiontaking a ppp and returning a vector of values, 
#'            the estimated summary statistic values for different r. 
#' @param N_crps How many samples of the predictive model should be generated for the MC-spproximation of the CRPS?
#' @param eval_points points where the K/F/G/pcf estimator is to be evaluated. Needs to take very specific form due to the way K_hat is coded in spatstat.


summaryScore = function(obs,
                        predictive.model,
                        est = 'intensity',
                        N_crps = max(length(obs),100), 
                        verbose = FALSE,
                        eval_points = seq(0,2.5,by = 0.05),
                        volume = range(eval_points)[2] - range(eval_points)[1])
  
{

  
  # get est into the correct form:
  if(!is.function(est))
  {
    if(est %in% c('intensity','Intensity','lambda','int','Int')) 
    {
      est = function(pp) return(as.vector(as.matrix(density(pp))))
      volume = area(obs[[1]]) ## assumes that all observations live on the same window
      
    }else if(est %in% c('RipleyK','K','Kfun','KFun','Kfunction','KFunction')) 
    {
      est = function(pp){return(as.function(Kinhom(pp))(eval_points))}
      
    }else if(est %in% c('F')) 
    {
      est = function(pp){return(as.function(Fest(pp))(eval_points))}
      
    }else if(est %in% c('G')) 
    {
      est = function(pp){return(as.function(Fest(pp))(eval_points))}
    }else if(est %in% c('pcf','pair correlation function')) 
    {
      est = function(pp){return(as.function(pcf(pp))(eval_points))}
    }
  }
  
  
  # convert obs into a list of ppps:
  if(is.ppp(obs)) obs = list(obs)
  
  if(verbose) print('applying estimator for observations')
  
  # get data table with estimator values for each observation:
  obs_dt = list()
  for(i in 1:length(obs))
  {
    if(verbose & (i%%10 == 0)) print(paste0(i,'/',length(obs)))
    obs_dt[[i]] = est(obs[[i]])
  }
  
  obs_dt = as.data.table(obs_dt)
  setnames(obs_dt,paste0('obs',1:length(obs)))
  
  # get N_crps simulations and store their estimators in a data table
  
  if(verbose) print('getting estimator for simulations')
  
  pr_dt = list()
  for(i in 1:N_crps)
  {
    if(verbose & (i%%10 == 0)) print(paste0(i,'/',N_crps))
    
    pr_dt[[i]] = est(predictive.model())
  }
  
  pr_dt = as.data.table(pr_dt)
  setnames(pr_dt,paste0('pr',1:N_crps))
  
  
  # get mean differences between observation estimates and simulated estimates
  
  if(verbose) print('getting absolute error between obs and simulations')
  
  ae_obs_pr = matrix(0,ncol = length(obs),nrow = N_crps)
  
  for(obs_in in 1:length(obs))
  {
    if(verbose & (obs_in%%10 == 0)) print(paste0(obs_in,'/',length(obs)))
    for(pr_in in 1:N_crps)
    {
      ae_obs_pr[pr_in,obs_in] = mean(abs(obs_dt[[obs_in]] - pr_dt[[pr_in]]),na.rm = TRUE)  
    }
  }
  
  # get mean differences between observation estimates and simulated estimates
  
  if(verbose) print('getting absolute error between combinations of simulations')
  
  # compute upper triagonal matrix
  ae_pr_pr = matrix(0,ncol = N_crps,nrow = N_crps)
  for(pr_in1 in 1:(N_crps-1))
  {
    if(verbose & (pr_in1%%10 == 0)) print(paste0(pr_in1,'/',N_crps))
    
    for(pr_in2 in (pr_in1+1):N_crps)
    {
      ae_pr_pr[pr_in1,pr_in2] = mean(abs(pr_dt[[pr_in1]] - pr_dt[[pr_in2]]),na.rm = TRUE)  
    }
  }
  
  # get crps: 
  
  mae_pr = sum(ae_pr_pr)/(N_crps*(N_crps-1)/2)
  
  mae_obs = colMeans(ae_obs_pr)
  
  crps = volume*(mae_obs - mae_pr/2)
  
  return(crps)
}
