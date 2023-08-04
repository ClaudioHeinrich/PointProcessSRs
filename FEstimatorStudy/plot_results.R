

rm(list = ls())

library(spatstat)

library(data.table)
library(ggplot2)
library(latex2exp)
library(patchwork)


setwd('~/NR/ProjectPointProcess/PointProcessSRs/FEstimatorStudy/')

# if you don't have it anymore:
# for LGCPs:
#install.packages("../RandomFieldsUtils_1.2.5.tar.gz",repos = NULL, type = "source")
#install.packages("../RandomFields_3.3.14.tar.gz",repos = NULL, type = "source")


plot_dir = '../../ModelEvaluation/Figures/review1/FStudy/'

dir.create(plot_dir)

set.seed(102030)


load('study_results.RData')
load(file = 'ppp_list.rds')

# correct integrals: in abs_errs we used cumsum, now we want to go to means:

eval_points = seq(0,2.5,0.05)
dt[,crps1:= crps1/match(R,eval_points)]
dt[,crps2:= crps2/match(R,eval_points)]

### plot examples of the point patterns ###

indices = c(1,1,1,1,1)
models = c('Str1','Str2','hP','LGCP1','LGCP2')
model_names = c('Strauss 1','Strauss 2','Poisson','LGCP 1', 'LGCP 2') # label such that '2' is always further away from Poisson



for(i in seq_along(indices))
{
  ppp = ppp_list[[models[i]]][[indices[i]]]
  png(file=  paste0(plot_dir,'pp_',models[i],'.png'),width = 960,height = 960)
    par(cex.main = 5,cex = 2)
    plot(ppp,pch = 19,main = model_names[i])
  dev.off()
}

########### study results #############

dt[,crps := crps1 - 1/2 * crps2]


# comparison with correct model:
dt_new = dt[obs_mod == pred_mod]
dt_new[,pred_mod:=NULL]
dt_new = dt_new[,.(crps1 = mean(crps1),crps2 = mean(crps2)),by = .(R,obs_mod)]
setnames(dt_new,c('crps1'),c('crps1_truemod'))
dt_new[,crps_truemod := crps1_truemod - 1/2 * crps2]
dt_new[,crps2:= NULL]

dt = merge(dt,dt_new,by = c('R','obs_mod'))
dt[,crps_diff:= crps - crps_truemod]
dt[,crps1_diff:= crps1 - crps1_truemod]
dt[,crps2_diff:= crps2 - mean(crps2)]


### first plot ###
obs1 = 5
obs2 = 20
nxs = c(25,seq(100,500,100))
Rs = c(0.5,1.5)

Obs_Mod = 'hP'

plot_dt1 = dt[R %in% Rs & ny %in% c(obs1,obs2) & obs_mod == Obs_Mod][,.(R,r,nx,pred_mod,ny,crps1,crps1_diff)]
plot_dt2 = dt[R %in% Rs & ny == obs1 & obs_mod == Obs_Mod][,.(R,r,nx,pred_mod,ny,crps2,crps2_diff)]# obs_mod and ny don't matter
# for consistent naming:
plot_dt2[,crps1 := 1/2*crps2]
plot_dt2[,crps1_diff := 1/2*crps2_diff]
plot_dt2[,c('crps2','crps2_diff') := NULL]
plot_dt2[,ny:= 1000] # dummy value

plot_dt = rbindlist(list(plot_dt1,plot_dt2))
plot_dt = plot_dt[nx %in% nxs]
plot_dt[,crps1 := crps1/R] # divide by length of integral.
plot_dt[,pred_mod := model_names[match(pred_mod,models)]]
plot_dt[,pred_mod := factor(pred_mod,levels = c('Strauss 2','Strauss 1','Poisson','LGCP 1','LGCP 2'))]

plot_dt[,nx := as.factor(nx)]


# custom labels for facets using latex:
rlabels = TeX(c(paste0('$R = ',Rs[1],'$'),paste0('$R = ',Rs[2],'$')),output = 'character')
ylabels = TeX(c(paste0('$E_1(G,F),\ N = ',obs1,'$'),paste0('$E_1(G,F),\ N = ',obs2,'$'),'$E_2(F)$'),output = 'character')
yvals = c(obs1,obs2,1000)
plot_dt[,rlab := rlabels[match(R,Rs)]]
plot_dt[,ylab := ylabels[match(ny,yvals)]]
# keep correct order:
plot_dt[,ylab := factor(ylab,levels = ylabels)]

theme_set(theme_bw(base_size = 18))

pp = ggplot(plot_dt) + 
  geom_boxplot(aes(x = nx,y = crps1,color = pred_mod,fill = pred_mod),alpha = 0.5) + 
  facet_grid(cols = vars(ylab),rows = vars(rlab),labeller = label_parsed) +
  scale_y_continuous('') + 
  scale_color_discrete('pred. dist.') + scale_fill_discrete('pred. dist.')+
  scale_x_discrete('n')
pp

ggsave(paste0(plot_dir,'score_variability.png'),width = 15,dpi = 'retina')


########################

### plot percentage of correct classifications ###

Rs = c(0.25,0.5,1,1.5)
nxs = c(25,100,500)
Obs_Mod = 'hP'
models_new = c('Str2','Str1','hP','LGCP1','LGCP2')
model_names_new = model_names[match(models,models_new)]

dt_pc = dt[,.(percentage_correct = 100*sum(crps_diff >0)/.N),by = .(R,nx,ny,obs_mod,pred_mod)]

plot_dt = dt_pc[obs_mod == Obs_Mod & nx %in% nxs & R %in% Rs & obs_mod != pred_mod]

xlabels = TeX(paste0('$n = ',nxs,'$'),output = 'character')
plot_dt[,xlab := xlabels[match(nx,nxs)]]
plot_dt[,xlab := factor(xlab,levels = xlabels)]


model_labels=  TeX(model_names_new,output = 'character')
plot_dt[,pred_mod_lab := model_labels[match(pred_mod,models_new)]]
plot_dt[,pred_mod_lab := factor(pred_mod_lab,levels = model_labels[-which(Obs_Mod == models_new)])]

# # We now want to use the parse-labeller for the rows and the value-labeller for the columns:
# # Creating custom labeller functions
# row_labeller <- function(...){
#   return(label_parsed(...))
# }
# 
# col_labeller <- function(...){
#   return(label_value(...))
# }
# 
# # Creating a combined labeller function
# combined_labeller <- function(variable, value) {
#   if (variable == 'cyl')
#     return(row_labeller(variable,value))
#   else
#     return(col_labeller(variable,value))
# }

plot_dt[,R := as.factor(R)]

theme_set(theme_bw(base_size = 18))

pp = ggplot(plot_dt[ny <= 25]) + 
  geom_line(aes(x = ny, y = percentage_correct, color = R,linetype = R)) +
  facet_grid(cols = vars(pred_mod_lab),rows = vars(xlab),labeller = label_parsed) +
  scale_y_continuous('percentage correctly classified') + 
  scale_x_continuous('N')
pp
  

ggsave(paste0(plot_dir,'percentage_correct.png'),width = 15,dpi = 'retina')
# plot_dt[,crps2_centralized:= crps2 - mean(crps2),by = .(R,pred_mod)]


# discrimination ability: #

# dt[,crps := crps1 - 1/2*crps2]
# dt[,pred_mod := factor(pred_mod,levels = models)]
# dt[,obs_mod := factor(obs_mod,levels = models)]
# 
# 
# dt[obs_mod == pred_mod, crps_diff := NA]
# 
# plot_dt = dt[ ny == 50 & nx == 500]
# plot_dt[,crps:= crps/R]
# plot_dt[,R := as.factor(R)]
# 
# 
# pp = ggplot(plot_dt) + geom_boxplot(aes(x = R,y = crps_diff,color = pred_mod,fill = pred_mod)) +facet_grid(cols = vars(obs_mod))
# pp


# percentage correct #

################### 

# calculate F-function estimates for all models


N = 100 # only use the first N point patterns (per model) in ppp_list

mods = c('Str1','Str2','hP','LGCP1','LGCP2')


eval_points = seq(0,2,0.01)
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

dt_F = data.table(x = eval_points)
for(i in seq_along(mods))
{
  indices = N*(i-1) + 1:N
  tempmat = matrix(unname(unlist(F_vals[indices])),ncol = N)
  tempmat[is.na(tempmat)] = 1
  dt_F[,(mods[i]) := rowMeans(tempmat)]
}

dt_F = melt(dt_F,'x')
dt_F[,variable := model_names[match(variable,mods)]]
model_names_new = model_names[c(2:1,3:5)]
dt_F[,variable := factor(variable, levels = model_names_new)]

theme_set(theme_bw(base_size = 18))
pp = ggplot(dt_F) + 
  geom_line(aes(x = x,y = value,color = variable,linetype = variable)) +
  scale_x_continuous('r') + 
  scale_y_continuous(TeX('${{F}}_r$')) + 
  scale_color_discrete('Model') + scale_linetype_discrete('Model')
pp

ggsave(paste0(plot_dir,'F_fcts.pdf'))

# combine plot of F-fcts and point patterns


indices = c(1,1,1,1,1)
models = c('Str1','Str2','hP','LGCP1','LGCP2')
model_names = c('Strauss 1','Strauss 2','Poisson','LGCP 1', 'LGCP 2') # label such that '2' is always further away from Poisson



theme_set(theme_bw(base_size = 25))
pp = ggplot(dt_F) + 
  geom_line(aes(x = x,y = value,color = variable,linetype = variable)) +
  scale_x_continuous('r') + 
  scale_y_continuous(TeX('')) + #scale_y_continuous(TeX('${{F}}_r$')) + 
  scale_color_discrete('Model') + scale_linetype_discrete('Model') 
pp


indices = rep(1,5)

for(i in seq_along(indices))
{
  i_new = i
  ppp = ppp_list[[models[i_new]]][[indices[i_new]]]
  dt_temp = data.table(x = ppp$x,y = ppp$y)
  assign(paste0('pp',i),value = ggplot(dt_temp) + 
           geom_point(aes(x=x,y = y)) + 
           ggtitle(model_names[i]) + 
           scale_x_continuous('') + scale_y_continuous('') + 
           theme(panel.grid = element_blank(),
                 axis.ticks = element_blank(),
                 axis.text = element_blank()))
}

library(patchwork)

p = pp2 + pp1 + pp3 + pp4 + pp5 + pp + plot_layout(ncol = 6)
p
ggsave(paste0(plot_dir,'pp_combined.pdf'),height = 5,width = 25)
