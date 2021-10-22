
######## plot results of permutation test analysis #########

library(ggplot2)

setwd('~/NR/ProjectPointProcess/ModelEvaluation/simstudy/KernelVsLog')

set.seed(102030)

#############################

load('pvalspermtest.RData') # get results

##############################

# specify font sizes for plots

title_size = 20
labels_size = 16
ticks_size = 12


## analysis ##

alpha = 0.05 # choose size of the test

pt_ana = as.data.table(expand.grid(est = ests,N_pt = N_vec_pt,pr_mod = 1:4,obs_mod = c(1,4)))


for(oo in obs_vec)
{
  for(NN_pt in N_vec_pt )
  {
    print(NN_pt)
    for(tt in ests)
    {
      print(paste0(tt,'...'))
      for(mod_pr in  1:4)
      {
        frac = p_vals[est == tt & pr_mod == mod_pr & N_pt == NN_pt & obs_mod == oo & p_val <= alpha,.N]/R
        pt_ana[est == tt & pr_mod == mod_pr & N_pt == NN_pt & obs_mod == oo, power:=frac]
      }
    }
  }
}

### plot 1: observed model 1 ###

for(obs in obs_vec)
{
  
  data = pt_ana[N_pt <= 50 & obs_mod == obs & pr_mod != obs_mod]
  
  pp = ggplot(data = data, aes(x = N_pt,y = power,color = factor(pr_mod))) + 
    geom_point(mapping = aes(shape = est)) +
    geom_line(mapping = aes( linetype = est)) 
  
  
  
  
  # get title and labels correct
  
  pp = pp + labs(title = substitute('True distribution F'[oo],list(oo = obs)),
                 x = 'n',
                 shape = 'score',
                 linetype = 'score',
                 color = 'pred. model')
  
  # change labels 
  
  pp = pp + scale_shape_discrete(labels = c('intensity','log')) + 
    scale_linetype_discrete(labels = c('intensity','log'))
  
  
  pp = pp +  theme(plot.title = element_text(size = title_size),
          axis.title = element_text(size = labels_size),
          axis.text = element_text(size = ticks_size),
          legend.title = element_text(size = labels_size),
          legend.text = element_text(size = ticks_size))
  
  
  # generate plot
  
  pdf(paste0('../../Figures/perm_test_log_vs_int',obs,'.pdf'))
    print(pp)
  dev.off()

}
