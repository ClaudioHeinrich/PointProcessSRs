
### setup simulation study ###


rm(list = ls())

library(spatstat)

library(data.table)
library(ggplot2)



setwd('~/NR/ProjectPointProcess/PointProcessSRs/FEstimatorStudy/')

# for LGCPs:
#install.packages("../RandomFieldsUtils_1.2.5.tar.gz",repos = NULL, type = "source")
#install.packages("../RandomFields_3.3.14.tar.gz",repos = NULL, type = "source")


set.seed(102030)

#### specify window ###

limx = c(0,10)
limy = c(0,10)

x_range = max(limx) - min(limx)
y_range = max(limy) - min(limy)

win = owin(xrange = limx,yrange = limy)

### set parameters ###

n_exp = 50 # number of expected points for all models, except hP+

models = c('Str1','Str2','hP','LGCP1','LGCP2')


# for Strauss process: approximate beta such that the number of expected points is roughly n_exp

gamma1 = 0.4
gamma2 = 0.1 # This takes freaking forever to simulate, probably rejection sampling. We should take smaller range of interaction
R = 0.75


get_beta = FALSE # I've done this once and this doesn't change and takes quite some time...
generate_ppp = FALSE


if(get_beta)
{
  N = 100
  
  exp_n_Strauss = function(beta,gamma) # returns the mean squared error between n_exp and the number of points of a Strauss model with parameter beta 
  {
    ns = c()
    for(i in 1:N){ns = c(ns, rStrauss(beta = beta, gamma = gamma, R = R, W = win)$n)}
    return(ns)
  }
  
  beta_vec = seq(0.75,1.5,by = 0.05)
  
  n_dt = data.table()
  for(beta in beta_vec)
  {
    print(beta)
    ns1 = exp_n_Strauss(beta,gamma1)
    ns2 = exp_n_Strauss(beta,gamma2)
    n_dt = rbindlist(list(n_dt,data.table(beta = beta,n1 = ns1,n2 = ns2)))
  }
  
  means = n_dt[,lapply(.SD,mean),.SDcols = c('n1','n2'), by = beta]
 
  means = melt(means,'beta',c('n1','n2'))
  ggplot(means) + geom_line(aes(x = beta,color = variable,y = value)) + geom_hline(yintercept = n_exp)
  
  means[,error:= (value - n_exp)^2]
  betas = means[,.SD[which.min(error)], by = variable][,beta]
  
  beta1 = betas[1]
  beta2 = betas[2]
} else {
  beta1 = 0.9
  beta2 = 1.35
}

# The simulation of these

# LGCP: homogeneous with 50 points in window:

# RandomFields required but deprecated. I filed a spatstat bug report.
#packageurl <- "https://cran.r-project.org/src/contrib/Archive/RandomFields/RandomFields_3.3.13.tar.gz"
#install.packages(packageurl, repos=NULL, type="source")
#install.packages("https://cran.r-project.org/src/contrib/Archive/RandomFieldsUtils/RandomFieldsUtils_1.1.0.tar.gz", repos = NULL, type = "source")


# ok, the RandomFields package seems to be a piece of work. We want exponential covariance,
# with variance v and mean mu.
# The following has nexp expected points in the window (area 100):


sc1 = 0.5
sc2 = 2

test = rLGCP(mu = -1/2 + log(n_exp / 100),v = 1,scale = sc1,win = win)

plot(test)

test = rLGCP(mu = -1/2 + log(n_exp / 100),v = 1,scale = sc2,win = win)

plot(test)


# testing whether we got the mean correct:
#
# ns1 = 0
# ns2 = 0
# for(i in 1:500)
# {
#   cat(paste0('\r',round(100*i/500),'%'))
#   n1 = rLGCP(mu = -1/2 + log(n_exp / 100),v = 1,scale = sc1,win = win)$n
#   n2 = rLGCP(mu = -1/2 + log(n_exp / 100),v = 1,scale = sc2,win = win)$n
#   ns1 = c(ns1,n1)
#   ns2 = c(ns2,n2)
# }




########################
mods = c('Str1','Str2','hP','LGCP1','LGCP2')

if(generate_ppp)
{
  # how many point patterns:
  N = 10000
  
  #model names:
  
  # model simulation routines:
  sims = list(function() rStrauss(beta = beta1, gamma = gamma1, R = R, W = win),
              function() rStrauss(beta = beta2, gamma = gamma2, R = R, W = win),
              function() runifpoint(n = rpois(1,lambda = n_exp), win = win),
              function() rLGCP(mu = -1/2 + log(n_exp / 100),v = 1,scale = sc1,win = win),
              function() rLGCP(mu = -1/2 + log(n_exp / 100),v = 1,scale = sc2,win = win))
  
  ppp_list = list()           
  for(j in seq_along(mods))
  {
    print(mods[j])
    
    templist = list()
    for(i in 1:N)
    {
      templist = c(templist,list(sims[[j]]()))
    }
  
    ppp_list = c(ppp_list,list(templist))
  }
  
  names(ppp_list) = mods
  
  save(ppp_list,file = 'ppp_list.rds')
}else load(file = 'ppp_list.rds')

#################################################################
#################################################################
#################################################################

# Calculate the F-estimator for the saved list of point patterns:

N = 1000 # only use the first N point patterns (per model) in ppp_list
eval_points = seq(0,2.5,0.05)
F_hat = function(dat){return(as.function(Fest(dat))(eval_points))}

F_vals = data.table()
for(mod in mods)
{
  temp_list = ppp_list[[mod]]
  
  l_est = length(eval_points)
  # initialize matrix:
  dt_temp = matrix(0,nrow = l_est,ncol = N)
  
  #simulate and fill:
  for(i in 1:N)
  {
    dat = temp_list[[i]]
    dt_temp[,i] = F_hat(dat)
  }
  
  dt_temp = as.data.table(dt_temp)
  setnames(dt_temp,c(paste0(mod,1:N)))
  
  F_vals = c(F_vals,dt_temp)
  
}

# NAs possible, and they mean 1:
F_vals[is.na(F_vals)] = 1

# computes the CRPS from the output of computeNEst  

dt = as.data.table(F_vals)

weight_cutoffs = seq(0.25,1.5,0.25)
K = length(weight_cutoffs)
co_indices = match(weight_cutoffs,eval_points)

# get_crps = function(dt,models,ret_abserrmat = FALSE)
# {
  # get the integrated absolute error for all column combinations in dt:
  nn = ncol(dt)
  abs_errs = array(0,dim = c(nn,nn,K))
  
  # fill upper triangle matrix, thereafter symmetrize
  for(i in 1:(nn-1))
  {
    if(i %% 100 == 0) cat(paste0('\r',round(100*i/(nn-1)),'%'))
    for(j in (i+1):nn)
    {
      abs_errs[i,j,] = cumsum(abs(dt[[i]]-dt[[j]]))[co_indices]
    }
  }
  
  # symmetrize:
  
  for(k in 1:K){
    abs_errs[,,k] = abs_errs[,,k]  + t(abs_errs[,,k])
  }
  
  
  # function to get mean crps over a range of observations (specified by y_inds) using the simulated point patterns specified by x_inds,
  # over a range of weight functions (specified by R_inds) and a range of models specified by mods.
  #crps(y_inds,x_inds,R_inds = 1:K,mods = mods,first_term = TRUE,second_term = TRUE)
  
  
  # produces a list of length mods, each containing a vector of length K:
  crps_second_term = function(xinds,R_inds = 1:K,models = mods)
  {
    whichmods = match(models,mods)
    ret_list = list()
    for(wm in whichmods)
    {
      xinds_temp = xinds + N*(wm-1)
      temp = apply(abs_errs[xinds_temp,xinds_temp,],c(3),FUN = mean)
      # correct for the diagonal being included in the mean:
      temp = length(xinds_temp)/(length(xinds_temp)-1) * temp
      ret_list[[c("Str1", "Str2",  "hP",    "LGCP1", "LGCP2")[wm]]] = temp
    }
    return(ret_list)
  }
  
  
  ns = seq(25,500,25)
  rr = 100
  
  dt_2 = data.table()
  for(nnn in ns)
  {
    print(nnn)
    for(rrr in 1:rr)
    {
      xinds = sample(N,nnn,replace = TRUE)
      temp = as.data.table(crps_second_term(xinds))
      temp[,R:=weight_cutoffs][,r:=rrr][,n := nnn]
      dt_2 = rbindlist(list(dt_2,temp))
    }
  }
  
  
  
  crps_first_term = function(yinds,xinds,R_inds = 1:K,models = mods)
  {
    whichmods = match(models,mods)
    ret_dt = data.table()
    for(om in whichmods)
    {
      yinds_temp = yinds + N*(om-1)
      temp_dt = data.table()
      for(pm in whichmods)
      {
        xinds_temp = xinds + N*(pm-1)
        temp_yx = apply(abs_errs[yinds_temp,xinds_temp,],c(3),FUN = mean)
        temp_dt = rbindlist(list(temp_dt,
                                 data.table(value = temp_yx,
                                            pred_mod = c("Str1", "Str2",  "hP",    "LGCP1", "LGCP2")[pm],
                                            R=weight_cutoffs)))
      }
      temp_dt[,obs_mod := c("Str1", "Str2",  "hP",    "LGCP1", "LGCP2")[om]]
      ret_dt = rbindlist(list(ret_dt,temp_dt))
    }
    return(ret_dt)
  }
  
  
nxs = seq(25,500,25)
nys = c(10,50,100)
rr = 100

dt_2 = data.table() # second term 
dt_1 = data.table() # first term
for(nx in nxs)
{
  print(nx)
  for(rrr in 1:rr)
  {
    
    
    xinds = sample(N,nx,replace = TRUE)
    temp2 = as.data.table(crps_second_term(xinds))
    temp2[,R:=weight_cutoffs][,r:=rrr][,nx := nx]
    dt_2 = rbindlist(list(dt_2,temp2))
    
    for(ny in nys)
    {
    
      yinds = sample(N,ny,replace = TRUE)
      temp1 = crps_first_term(yinds,xinds)
      temp1[,r:=rrr][,nx := nx][,ny := ny]
      dt_1 = rbindlist(list(dt_1,temp1))
    }
  }
}

dt_2 = melt(dt_2,c('R','r','nx'))

setnames(dt_2,c('variable','value'),c('pred_mod','crps2'))
setnames(dt_1,'value','crps1')
dt_2[,pred_mod := as.character(pred_mod)]

dt = merge(dt_1,dt_2,by = c('R','r','nx','pred_mod'))


save(dt,dt_2,N,file = 'study_results.RData')
