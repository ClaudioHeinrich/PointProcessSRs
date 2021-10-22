

rm(list = ls())

library(data.table)

set.seed(102030)

setwd('~/NR/ProjectPointprocess/ModelEvaluation/simstudy/LTest')

# set parameters for models

n = 100 # number of bins

# intensity increases linearly from mu - a to mu + a
mu_vec = c(0.6,0.55,0.5,0.45)
a_vec = c(0,0.1,0.2,0.3) 

n_mod = length(mu_vec)

# parameters for the test

alpha = 0.05 # size of test
MC_sam = 1000 # for approximating gamma
N = MC_sam # repetitions of the test
randomized = FALSE # should we consider a randomized test to achieve maximal power at size alpha?


########## functions ###########

gamma_Poi = function(lam_vec,y, MC_sam = 1000)
{
  n = length(lam_vec)
  
  X = rpois(n*MC_sam,lambda = lam_vec)
  
  LLH = dpois(X,lambda = lam_vec,log = TRUE)
  
  LLH = colSums(matrix(LLH,nrow = n))
  
  LLH_obs = sum(dpois(y,lambda = lam_vec,log = TRUE))
  
  ret_val = mean(LLH < LLH_obs)
  
  return(ret_val)
}

# for resolving ties at random

gamma_Poi_ran = function(lam_vec,y, MC_sam = 1000)
{
  n = length(lam_vec)
  
  X = rpois(n*MC_sam,lambda = lam_vec)
  
  LLH = dpois(X,lambda = lam_vec,log = TRUE)
  
  LLH = colSums(matrix(LLH,nrow = n))
  
  LLH_obs = sum(dpois(y,lambda = lam_vec,log = TRUE))
  
  ties = which(abs(LLH - LLH_obs) < 1e-8)
  
  if(identical(ties,integer(0)))
  {
    ret_val = mean(LLH < LLH_obs)
  } else{
    ret_val = mean(LLH[-ties] < LLH_obs) + rbinom(n=1, size = length(ties), prob = 1/2) / length(ties)
  }
  
  return(ret_val)
}



######################

# get results:

dt = as.data.table(expand.grid(model = 1:4,observ = 1:4))

for(mod in 1:n_mod)
{
  mod_lam = seq(mu_vec[mod] - a_vec[mod],mu_vec[mod] + a_vec[mod], length.out = n)
  
  for(obs in 1:n_mod)
  {
    print(c(mod,obs))
    
    obs_lam = seq(mu_vec[obs] - a_vec[obs],mu_vec[obs] + a_vec[obs], length.out = n)
    
    gammas = NULL
    for(j in 1:N)
    {
      print(j)
      y = rpois(n = n,lambda = obs_lam)
      if(!randomized) gammas[j] = gamma_Poi(lam_vec = mod_lam, y = y, MC_sam = MC_sam)
      if(randomized) gammas[j] = gamma_Poi_ran(lam_vec = mod_lam, y = y, MC_sam = MC_sam)
    }
    dt[model == mod & observ == obs,gamma := mean(gammas)]
    dt[model == mod & observ == obs, fails := mean(gammas <= alpha)]
    dt[model == mod & observ == obs, fails_ts := mean(gammas <= alpha/2 | gammas >= 1-alpha/2)]
  }
}


### save ###

save.image(file = 'Ltest.RData')


# 
# fails1s = matrix(0,nrow = n_mod,ncol = n_mod)
# fails2s = matrix(0,nrow = n_mod,ncol = n_mod)
# gammas = matrix(0,nrow = n_mod,ncol = n_mod)
# 
# for(obs_mod in  1:n_mod)
# {
#   
#   obs_lam = seq(mu_vec[obs_mod] - a_vec[obs_mod],mu_vec[obs_mod] + a_vec[obs_mod], length.out = n)
#   
#   for(fc_mod in 1:n_mod)
#   {
#     print(c(obs_mod,fc_mod))
#     
#     fc_lam = seq(mu_vec[fc_mod] - a_vec[fc_mod],mu_vec[fc_mod] + a_vec[fc_mod], length.out = n)
#     
#     
#     fails_1s = 0
#     gamma_mean = 0
#     fails_2s = 0
#     for(i in 1:test_sam)
#     {
#       obs = rpois(n = n,lambda = obs_lam)
#       if(!randomized) gamma = gamma_Poi(lam_vec = fc_lam, y = obs, MC_sam = MC_sam)
#       if(randomized) gamma = gamma_Poi_ran(lam_vec = fc_lam, y = obs, MC_sam = MC_sam)
#       gamma_mean = gamma_mean + gamma/test_sam
#       if(gamma < alpha ) fails_1s = fails_1s + 1/test_sam
#       if(gamma < alpha/2 | gamma > 1-alpha/2 ) fails_2s = fails_2s + 1/test_sam
#     }
#     gammas[obs_mod,fc_mod] = gamma_mean
#     fails1s[obs_mod,fc_mod] = fails_1s
#     fails2s[obs_mod,fc_mod] = fails_2s
#   }
# }
# 
# fails1s
# fails2s
# gammas

# save

############
# plot gammas
# 
# plot_dir = '~/NR/ProjectPointProcess/ModelEvaluation/Figures/'
# 
# 
# 
# pdf(paste0(plot_dir,'LTestGammas.pdf'))
#   
#   col_vec = c('blue','red','black','green')
#   
#   library(latex2exp)
#   
#   x_lab = c(TeX('$F_1$'),TeX('$F_2$'),TeX('$F_3$'),TeX('$F_4$'))
#   
#   par(mar = c(4,2,2,2))
#   
#   plot(gammas[1,],ylim = range(gammas),pch = 15,type = 'b',col = col_vec[1],cex = 1.2,xaxt = 'n',xlab = 'pred. dist.',ylab = TeX('$\\overline{\\gamma}$'))
#   axis(side = 1,at = 1:4,labels = x_lab)
#   
#   points(x = 1,y = gammas[1,1],pch = 1,cex = 3,col = col_vec[1])
#   
#   for(i in 2:4)
#   {
#     lines(gammas[i,],pch = 14+i,type = 'b',col = col_vec[i],cex = 1.2)
#   
#     points(x = i,y = gammas[i,i],pch = 1,cex = 3,col = col_vec[i])
#   }
#   
#   legend("topright", legend = x_lab, pch=15:19, lty = 1, cex=1.5, bty='n',
#          col = col_vec,title = 'true dist.')
# 
# dev.off()
# 
# ### plot probabilities of failure ###
# 
# 
# pdf(paste0(plot_dir,'rp1s.pdf'))
#   
#   col_vec = c('blue','red','black','green')
#   
#   x_lab = c(TeX('$F_1$'),TeX('$F_2$'),TeX('$F_3$'),TeX('$F_4$'))
#   
#   par(mar = c(4,2,2,2))
#   
#   plot(fails1s[1,],ylim = range(fails1s),pch = 15,type = 'b',col = col_vec[1],cex = 1.2,xaxt = 'n',xlab = 'pred. dist.')
#   axis(side = 1,at = 1:4,labels = x_lab)
#   
#   points(x = 1,y = fails1s[1,1],pch = 1,cex = 3,col = col_vec[1])
#   
#   for(i in 2:4)
#   {
#     lines(fails1s[i,],pch = 14+i,type = 'b',col = col_vec[i],cex = 1.2)
#     
#     points(x = i,y = fails1s[i,i],pch = 1,cex = 3,col = col_vec[i])
#   }
#   
#   legend("topleft", legend = x_lab, pch=15:19, lty = 1, cex=1.5, bty='n',
#          col = col_vec,title = 'true dist.')
# 
# 
# dev.off()
# #
# 
# pdf(paste0(plot_dir,'rp2s.pdf'))
#   
#   col_vec = c('blue','red','black','green')
#   
#   x_lab = c(TeX('$F_1$'),TeX('$F_2$'),TeX('$F_3$'),TeX('$F_4$'))
#   
#   par(mar = c(4,2,2,2))
#   
#   plot(fails2s[1,],ylim = range(fails2s),pch = 15,type = 'b',col = col_vec[1],cex = 1.2,xaxt = 'n',xlab = 'pred. dist.')
#   axis(side = 1,at = 1:4,labels = x_lab)
#   
#   points(x = 1,y = fails2s[1,1],pch = 1,cex = 3,col = col_vec[1])
#   
#   for(i in 2:4)
#   {
#     lines(fails2s[i,],pch = 14+i,type = 'b',col = col_vec[i],cex = 1.2)
#     
#     points(x = i,y = fails2s[i,i],pch = 1,cex = 3,col = col_vec[i])
#   }
#   
#   legend("topleft", legend = x_lab, pch=15:19, lty = 1, cex=1.5, bty='n',
#          col = col_vec,title = 'true dist.')
#   
# dev.off()
# 

################### create example plots of models  #########

library(spatstat)


#function for gluing ppps:

glueppp = function(lppp)
{
  n = length(lppp)/5 # a ppp is a list of 5 elements
  
  #glue windows
  
  win_inds = ((1:n)-1)*5 +1
  x_range = NULL
  y_range = NULL
  
  for(i in win_inds)
  {
    x_range = c(x_range,lppp[[i]]$xrange)
    y_range = c(y_range,lppp[[i]]$yrange)
  }
  
  new_win = owin(xrange = range(x_range),yrange = range(y_range))
  
  # glue x and y vectors
  
  new_x = NULL
  x_inds = ((1:n)-1)*5 +3
  for(i in x_inds) new_x = c(new_x,lppp[[i]])
  
  new_y = NULL
  y_inds = ((1:n)-1)*5 +4
  for(i in y_inds) new_y = c(new_y,lppp[[i]])
  
  return(ppp(x = new_x,y = new_y,window = new_win))
}



#### specify models ####




for(i in 1:n_mod)
{
  mu = mu_vec[i]
  a = a_vec[i]
  lam_vec = seq(mu-a,mu+a,length.out = n)
  
  lppp = list()
  
  for(j in 1:n)
  {
    xcor = j %% 10
    xcor[xcor == 0] = 10
    ycor = (j-xcor)/10 + 1
    
    pp = rpoispp(lambda = lam_vec[j],win = owin(xrange = c(xcor-1,xcor),yrange = c(ycor-1,ycor)))
    lppp = c(lppp,pp)
  }
  pdf(paste0(plot_dir,'model',i,'Ltest.pdf'))
    par('mar' = c(1,1,1,1))
    par('cex' = 3)
    plot(glueppp(lppp),main = '',pch = 20)
  dev.off()
}



