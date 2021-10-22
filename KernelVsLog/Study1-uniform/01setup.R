
### setup simulation study ###


rm(list = ls())

library(spatstat)

library(data.table)


setwd('~/NR/ProjectPointProcess/ModelEvaluation/simstudy/KernelVsLog/Study1-uniform/')

set.seed(102030)

# how many samples should we use for the mean score?

N_vec = c(10,100)

#### specify window ###

limx = c(0,10)
limy = c(0,10)

x_range = max(limx) - min(limx)
y_range = max(limy) - min(limy)

win = owin(xrange = limx,yrange = limy)

vol = x_range * y_range

### set parameters for models ###


model_names = c('hP','hP-','hP+','ihP')

nexps = c(50,40,60,50)

#inhomogeneous lambda: lambda = ax + b, needs to be vectorized function
ihlambda = function (x,y) return(1*x +5)
ihlambda_int = pracma::integral2(ihlambda,
                                 xmin = limx[1],xmax = limx[2],
                                 ymin = limy[1],ymax = limy[2])$Q




lambdas = as.list(nexps[1:3] / vol )
lambdas[[4]] = function(x,y) return(nexps[4] *ihlambda(x,y)/ihlambda_int)

# create multi-level list

models = list()
for(i in 1:length(model_names) )
{
  mod_i = list(nexp = nexps[i],lambda = lambdas[[i]])
  models = c(models,list( mod_i))
}

names(models) = model_names

################ functions ########################

# function to simulate point patterns with these parameters

simPP = function(i)
{
  dat = rpoispp(lambda = models[[i]]$lambda,win = win)
  return(dat)
}

# function to get loglikelihood for point pattern ppp under model i

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



computeNEst = function(pplist,mod,est,silence = TRUE)
{
  #mod is an index vector for the models, e.g. 1:6, pplist is a list of N*length(mod) point patterns
  N = length(pplist)/length(mod)
  # compute the estimator for N simulations of model mod.
  
  # get the length of the output of est:
  l_est = length(est(pplist[[1]]))
  
  #initialize data table
  dt = data.table()
  
  counter = 0
  for(mm in mod){
    
    if(! silence) print(mm)
    
    # initialize matrix:
    dt_temp = matrix(0,nrow = l_est,ncol = N)
    
    #simulate and fill:
    for(i in 1:N)
    {
      if(!silence & (i %% 10 == 0)) print(paste0(i,'/',N))
      
      dat = pplist[[counter*N + i]]
      dt_temp[,i] = est(dat)
    }
    
    dt_temp = as.data.table(dt_temp)
    setnames(dt_temp,c(paste0(mm,'_',1:N)))
    
    dt =c(dt,dt_temp)
    
    dt = as.data.table(dt)
    
    counter = counter + 1
  }
  return(dt)
}





get_crps = function(dt,mod,ret_abserrmat = FALSE)
{
  # the i,j-th entry of the absolute error matrix is the mean absolute difference between columns i and j
  # dt is the output from computeNEst and contains N realizations of the used estimator for each considered model.
  # For each model combination, N crps values are derived: Say the observation model is m_o and the forecast model is m_f, then
  # for each estimator est_i for i in m_o we compute 1/N sum_{j in m_f}|est_i - est_j| - 1/2 (1/(N(N-1))) sum_{j,k in m_f} |est_j - est_k|.
  # Note that, by this calculation, both the number of derived crps-values and the number of Monte-Carlo samples for the predictive distributions are N.
  # This got me confused once before...
  
  # get the integrated absolute error for all column combinations in dt:
  nn = ncol(dt)
  abs_err_matrix = matrix(0,nrow = nn,ncol = nn)
  
  # fill upper triangle matrix, thereafter symmetrize
  for(i in 1:(nn-1))
  {
    if(i %% 100 == 0) print(paste0(i,'/',nn-1))
    for(j in (i+1):nn)
    {
      abs_err_matrix[i,j] = mean(abs(dt[[i]]-dt[[j]]))
    }
  }
  
  # symmetrize:
  abs_err_matrix = abs_err_matrix + t(abs_err_matrix)
  
  
  crps = data.table()
  
  for(dat_mod  in mod)
  {
    for(fc_mod  in mod)
    {
      
      # get relevant sections of the matrix containing the mean absolute errors
      N = nn/length(mod)
      
      dt_fc_mat = abs_err_matrix[(dat_mod-1) * N + 1:N, (fc_mod-1) * N + 1:N]
      fc_fc_mat= abs_err_matrix[(fc_mod-1) * N + 1:N,(fc_mod-1) * N + 1:N]
      
      crps = data.table(crps,
                        setnames(data.table((N/(N-1))*rowMeans(dt_fc_mat) - 1/2 * mean(fc_fc_mat[fc_fc_mat>0])),
                                 paste0('ks_mod',fc_mod,'_obs',dat_mod)))
      
    }
  }
  
  crps = crps[,lapply(X = .SD,as.numeric),.SDcols = colnames(crps)]
  
  
   return(crps)
}


# save everything

save.image(file = 'setup.RData')



  