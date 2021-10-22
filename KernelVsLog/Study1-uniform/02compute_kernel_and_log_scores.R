
### setup simulation study ###


rm(list = ls())

library(spatstat)
library(ggplot2)
library(data.table)


setwd('~/NR/ProjectPointProcess/ModelEvaluation/simstudy/KernelVsLog/Study2-Displacements/')

set.seed(1234)

load('setup.RData')

###############################
 



N = 3*max(N_vec) # get more scores to draw from when bootstrapping the mean

nm = length(models)

print(paste('simulate ',N,' point patterns from each predictive distribution...'))

for(i in 1:nm)
{
  pplist = list()
  for(j in 1:N)
  {
    pplist[[j]] = simPP(i)
  }
  assign(paste0('pp',i),pplist)
}

# log scores

log_scores = data.table()

for(pred_ind in 1:nm)
{
  for(obs_ind in 1:nm)
  {
    obs_pps = get(paste0('pp', obs_ind ))
    logs = unlist(lapply(X = obs_pps,FUN = function(x){-lll(x,pred_ind)}))
    temp_dt = data.table(obs = obs_ind,pred_ind,NN = 1:N,score = logs)
    
    log_scores = rbindlist(list(log_scores,temp_dt))
  }
}


##### Intensity score ####

# bandwidth selection. Standard is sigma = (1/8) * shortside(as.rectangle(w)) = 1.25

bws = c(0.1,0.25,0.5,1,1.5)

int_scores = list()

for(ii in 1:length(bws))
{
  bw = bws[ii]
  
  lambda_hat = function(dat){return(as.matrix(density.ppp(dat,sigma = bw)))}
  
  pplist = list()
  for(j in 1:nm)
  {
    pp = get(paste0('pp',j))
    pplist = c(pplist,pp)
  }
  
  print(paste0('bandwidth = ',bw,', computing estimators...'))
  NEst_lambdascore = computeNEst(pplist = pplist,mod = 1:nm,est = lambda_hat,silence = TRUE)
  
  print(paste0('bandwidth = ',bw,', computing crps...'))
  lambda_scores = get_crps(dt = NEst_lambdascore,mod = 1:nm,ret_abserrmat = TRUE)
  setnames(lambda_scores,substring(names(lambda_scores),4))
  int_scores[[ii]] = list(lambda_scores = lambda_scores,bandwidth = bw)
}


log_scores = dcast(log_scores,NN ~ obs + pred_ind)

new_names = unlist(lapply(paste0('obs',1:nm),function(x){paste0('mod',1:nm,'_',x)}))
setnames(log_scores,c('NN',new_names))
log_scores[,NN := NULL]

setcolorder(log_scores,colnames(int_scores[[1]]$lambda_scores))

###############

# the boot function sucks for this, so do-it-yourself it is:

dt = data.table()

R = 500

bt_dt = data.table()

for(N in N_vec)
{
  print(N)
  
  all_inds = c()
  
  counter = 0  
 
  while(counter < R)
  {
    if(counter %% 100 == 0) print(counter)
    
    inds = sample(3*max(N_vec),size = N)
    all_inds = c(all_inds,inds)
    
    dt_temp_log = log_scores[inds,lapply(.SD,mean)]
    dt_temp_log[,type := 'log'][,bw := 0][,NN := N]
    
    bt_dt = rbindlist(list(bt_dt,dt_temp_log))
    
    for(ii in 1:length(bws))
    {
      dt_temp_int = int_scores[[ii]]$lambda_scores[inds,lapply(.SD,mean)]
      dt_temp_int[,type := 'kernel'][,bw := bws[ii]][,NN := N]
      
      bt_dt = rbindlist(list(bt_dt,dt_temp_int))
    }
    
    counter = counter+1
  }
}     
  
  
setkey(bt_dt,NN,type,bw)   

bt_dt = melt(bt_dt,id.vars = c('type','bw','NN'),variable.name = 'mod_com') 


get_obs_mod = function(x)
{
  substring(x,9)
}

get_pr_mod = function(x)
{
  substring(x,4,4)
}

bt_dt[,obs := get_obs_mod(mod_com)][,pr:=get_pr_mod(mod_com)]
bt_dt[,mod_com := NULL]

setnames(bt_dt,'value','score')

dt = copy(bt_dt)


save(dt,file = 'scores.RData')

save(bws,log_scores,int_scores,file = 'raw_scores.RData')


