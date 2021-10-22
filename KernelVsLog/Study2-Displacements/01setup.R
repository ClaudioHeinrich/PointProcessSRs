
### setup simulation study ###


rm(list = ls())

library(spatstat)

library(data.table)

timer1 = Sys.time()

setwd('~/NR/ProjectPointProcess/ModelEvaluation/simstudy/KernelVsLog/Study2-Displacements//')

set.seed(1234)

# how many samples should we use for the mean score?

N_vec = c(10,100)

#### specify window ###

limx = c(-5,5)
limy = c(-5,5)

x_range = max(limx) - min(limx)
y_range = max(limy) - min(limy)

win = owin(xrange = limx,yrange = limy)

vol = x_range * y_range

### set parameters for models ###

n_mod = 6

#model_names = c('hP','hP-','hP+','ihP')

nexps = 100 # this is in fact the multiplicative constant used on the Gaussian kernel, not exactly the number of points. It would be exactly the number of points if the
           # Gaussian kernels determining the models would be exactly 0 out of the window, which is approximately the case.

# the true model has as intensity the Gaussian kernel with these specifications:
sig1 = 1 
mu1 = c(0,0)
cor1 = 0

mu2 = c(0.1,0) # model 2 has the kernel mean shifted by mu1
sig3 = 0.9 # sd of the kernel for model 3
sig4 = 1.1 # sd of the kernel for model 4
cor5 = 0.1 # correlation of model 5 

nexps6 = 105

gaussian_kernel = function(x,y,mu=mu1,sig = sig1,cor = cor1)
{
  if(cor == 0)
  {
  exponent = -1/2 * ((x-mu[1])^2/sig^2 + (y-mu[2])^2/sig^2)
  ret_val = 1/(2*pi*sig^2) * exp(exponent) 
  } else {
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


# get intensities
lambdas = list()

lambdas[[1]] = function(x,y) return(nexps * gaussian_kernel(x,y))
lambdas[[2]] = function(x,y) return(nexps * gaussian_kernel(x,y,mu = mu2))
lambdas[[3]] = function(x,y) return(nexps * gaussian_kernel(x,y,sig = sig3))
lambdas[[4]] = function(x,y) return(nexps * gaussian_kernel(x,y,sig = sig4))
lambdas[[5]] = function(x,y) return(nexps * gaussian_kernel(x,y,cor = cor5))
lambdas[[6]] = function(x,y) return(nexps6 * gaussian_kernel(x,y))

# create multi-level list

models = list()
nexpvec = c(rep(nexps,n_mod-1),nexps6)

for(i in 1:n_mod)
{
  mod_i = list(nexp = nexpvec[i],lambda = lambdas[[i]])
  models = c(models,list( mod_i))
}

### nexp is not exactly the number of expected points since the Gaussian kernels are not exactly 0 
# outside the window. We should correct for that:


for(ii in 1:length(models))
{
  FUN = models[[ii]]$lambda
  nexp_precise = pracma::integral2(FUN,
                                   xmin = win$xrange[1],
                                   xmax = win$xrange[2],
                                   ymin = win$yrange[1],
                                   ymax = win$yrange[2])
  models[[ii]] = c(models[[ii]],list(nexp_precise = nexp_precise$Q))
}

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
  nexp = models[[i]]$nexp_precise
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
      
      
      # this row takes out the diagonal of dt_fc_mat, to make things out-of-sample when dt = fc, and comparable if not.
      dt_fc_mat = t(matrix(t(dt_fc_mat)[row(dt_fc_mat) != col(dt_fc_mat)],nrow = N-1))
      
      
      crps = data.table(crps,
                        setnames(data.table(rowMeans(dt_fc_mat) - 1/2 * mean(fc_fc_mat[fc_fc_mat>0])),
                                 paste0('ks_mod',fc_mod,'_obs',dat_mod)))
      
    }
  }
  
  crps = crps[,lapply(X = .SD,as.numeric),.SDcols = colnames(crps)]
  
  
   return(crps)
}


timer1 = Sys.time() - timer1 

# save everything

save.image(file = 'setup.RData')


