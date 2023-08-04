
### setup simulation study ###


rm(list = ls())

library(spatstat)

library(data.table)

#install.packages("../RandomFieldsUtils_1.2.5.tar.gz",repos = NULL, type = "source")
#install.packages("../RandomFields_3.3.14.tar.gz",repos = NULL, type = "source")


setwd('~/NR/ProjectPointProcess/PointProcessSRs/ScoringRules/')

set.seed(102030)

#### specify window ###

limx = c(0,10)
limy = c(0,10)

x_range = max(limx) - min(limx)
y_range = max(limy) - min(limy)

win = owin(xrange = limx,yrange = limy)

### set parameters ###

n_exp = 50 # number of expected points for all models, except hP+

models = c('hP','hP+','ihP','Str','ihT','LGCP')

# intensity for nonhomogeneous process:

lambda0 = function(x,y){return(sqrt(x^2+y^2))}

int_val = pracma::integral2(lambda0,
                            xmin = limx[1],xmax = limx[2],
                            ymin = limy[1],ymax = limy[2])

lambda =  function(x,y){return(n_exp * sqrt(x^2+y^2)/int_val$Q)}

# for Strauss process: approximate beta such that the number of expected points is roughly n_exp

gamma = 0.5
R = 1

N = 200

exp_n_Strauss = function(beta) # returns the mean squared error between n_exp and the number of points of a Strauss model with parameter beta 
{
  ret_val = 0
  for(i in 1:N){ret_val = ret_val + (1/N) * (rStrauss(beta = beta, gamma = gamma, R = R, W = win)$n - n_exp)^2}
  return(ret_val)
}

beta_vec = seq(0.8,1.4,by = 0.05)

exp_n = c()
for (b in beta_vec)
{
  print(b)
  exp_n = c(exp_n,exp_n_Strauss(b))
}

plot(beta_vec,exp_n)

beta = beta_vec[which.min(exp_n)] 


# for Thomas process 

m = 2 # mean number of points per cluster, the intensity of the parent process is scaled by 1/m
sc = 1 # scale parameter

# LGCP: homogeneous with 50 points in window:

# RandomFields required but deprecated. I filed a spatstat bug report.
#packageurl <- "https://cran.r-project.org/src/contrib/Archive/RandomFields/RandomFields_3.3.13.tar.gz"
#install.packages(packageurl, repos=NULL, type="source")
#install.packages("https://cran.r-project.org/src/contrib/Archive/RandomFieldsUtils/RandomFieldsUtils_1.1.0.tar.gz", repos = NULL, type = "source")


# ok, the RandomFields package seems to be a piece of work. We want exponential covariance,
# with variance v and mean mu.
# I think, the following should has 60 expected points in the window (area 100):

test = rLGCP(mu = -1/2 + log(0.6),v = 1,win = win)

########################


##### functions #####


# function to simulate point patterns with these parameters

sim = function(type)
{
  if(type == 'hP') dat = runifpoint(n = rpois(1,lambda = n_exp), win = win)
  if(type == 'hP+') dat = runifpoint(n = rpois(1,lambda = 1.2 * n_exp), win = win)
  if(type == 'ihP') dat = rpoispp(lambda = lambda,win = win)
  if(type == 'Str') dat = rStrauss(beta = beta, gamma = gamma, R = R, W = win)
  if(type == 'ihT')  dat = rThomas(kappa = function(x,y){return( (1/m)*lambda(x,y))},scale = sc, mu = m, win = win)  
  if(type == 'LGCP') dat = rLGCP(mu = -1/2 + log(0.6),v = 1,win = win)
  return(dat)
}


#' Auxiliary function, returns for several models multiple realizations of the estimator
#' 
#' @param mod Vector of models
#' @param est {The considered estimator. Needs to be a function that takes a ppp as input and returns a vector vec, the sum over this vector is then considered to be the integral over mathcal R
#'            i.e. each entry of the vector resembles \wh T(dat,r) for a fixed r. The length of vec may not depend on the input data.}
#' @param N The number of realizations of the estimator for each model
#' 
#' @return data table with length(mod) x N columns, each one labelled paste0(mod,i) where i = 1:N. Each column has n rows where n is the output length of est

computeNEst = function(mod,est,N,silence = TRUE)
{
  # get the length of the output of est:
  l_est = length(est(sim(mod[1])))
  
  #initialize data table
  dt = data.table()
  
  for(mm in mod){
    
    if(! silence) print(mm)
    
    # initialize matrix:
    dt_temp = matrix(0,nrow = l_est,ncol = N)
    
    #simulate and fill:
    for(i in 1:N)
    {
      if(!silence & (i %% 10 == 0)) print(paste0(i,'/',N))
      
      dat = sim(mm)
      dt_temp[,i] = est(dat)
    }
    
    dt_temp = as.data.table(dt_temp)
    setnames(dt_temp,c(paste0(mm,1:N)))
    
    dt =c(dt,dt_temp)
    
    dt = as.data.table(dt)
  }
  return(dt)
}


# computes the CRPS from the output of computeNEst  

get_crps = function(dt,models,ret_abserrmat = FALSE)
{
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
  
  
  crps = matrix(0,nrow = length(models),ncol = length(models))
  
  for(dat_mod_ind  in 1:length(models))
  {
    for(fc_mod_ind  in 1:length(models))
    {
      data_mod = models[dat_mod_ind]
      fc_mod = models[fc_mod_ind]
      
      # get relevant sections of the matrix containing the mean absolute errors
      N = nn/length(models)
      
      dt_fc_mat = abs_err_matrix[(dat_mod_ind-1) * N + 1:N, (fc_mod_ind-1) * N + 1:N]
      fc_fc_mat= abs_err_matrix[(fc_mod_ind-1) * N + 1:N,(fc_mod_ind-1) * N + 1:N]
      
      crps[dat_mod_ind,fc_mod_ind] = mean(dt_fc_mat[dt_fc_mat>0]) - 1/2 * mean(fc_fc_mat[fc_fc_mat>0])
      
    }
  }
  
  CRPS = data.table(crps)
  setnames(CRPS,models)
  
  if(ret_abserrmat) return(abs_err_matrix) else return(CRPS)
}




get_crps_errbars = function(aem,models)
{
  n = length(models)
  N = dim(aem)[1]/n
  
  mod_2_vec = paste0('obs_',rep(models,each = n),'_fc_',rep(models,n))
  
  ret_dt = data.table()
  
  for(dat_mod_ind  in 1:length(models))
  {
    for(fc_mod_ind  in 1:length(models))
    {
      data_mod = models[dat_mod_ind]
      fc_mod = models[fc_mod_ind]
      
      # get relevant sections of the matrix containing the mean absolute errors
      dt_fc_mat = aem[(dat_mod_ind-1) * N + 1:N, (fc_mod_ind-1) * N + 1:N]
      fc_fc_mat= aem[(fc_mod_ind-1) * N + 1:N,(fc_mod_ind-1) * N + 1:N]
      
      
      crps1 = apply(X = dt_fc_mat,MARGIN = 2,FUN = mean)
      crps2 = mean(fc_fc_mat[fc_fc_mat>0])
      
      crps = crps1 - 1/2 * crps2
      
      name = paste0('obs_',models[dat_mod_ind],'_fc_',models[fc_mod_ind])
      dt_temp = data.table(crps)
      setnames(dt_temp,name)
      
      ret_dt = data.table(ret_dt,dt_temp)
      
    }
  }
  return(ret_dt) 
}


#### permutation tests ####

#' Run a permutation test of the pairwise difference between two vectors of numbers
#' @param a Vector, the scores from one method
#' @param b Vector, the scores from some other method
#' @param N Integer, the size of the permutation distribution
#' @return A list with the mean of the difference and the permutation distribution of that difference
#' @examples
#' N = 1e2
#' trend  = 1:N
#' a = trend + .01 + rnorm(N, .001)
#' b = trend - .01 + rnorm(N, .001)
#' l = permutation_test_difference(a,b)
#' q = sum(l$D <= l$d_bar) / length(l$D)
#' @author Alex
#' 
permutation_test_difference = function(a,
                                       b,
                                       pval = TRUE,
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
  
  p_val = NULL
  if(pval)
  {
    p_val  = rank(x = c(d_bar,D))[1]/(N+1)
    
  }
  
  return(list(d_bar = d_bar, D = D,p_val = p_val))
}



# save everything 

save.image(file = 'setup_r1.RData')
