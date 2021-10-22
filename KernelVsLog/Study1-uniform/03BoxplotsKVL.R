
rm(list = ls())

library(spatstat)

library(data.table)

library(ggplot2)

setwd('~/NR/ProjectPointProcess/ModelEvaluation/simstudy/KernelVsLog/Study1-uniform/')

set.seed(102030)

load('setup.RData')

# specify font sizes

theme_set(theme_bw(base_size = 22))

# load data

load(file = paste0('scores.RData'))

pdf_dir = '../../../Figures/study2/uniform/'
dir.create(pdf_dir,showWarnings = F)

# For every score apply an affine linear function to scale it to [0,1] so that they can be plotted on the same range.

dt[,obs := as.factor(obs) ][,pr := as.factor(pr)][,bw := as.factor(bw) ][,NN := as.factor(NN) ]


for(oo in unique(dt[,obs]))
{
  for(bw0 in unique(dt[,bw]))
  {
    for(N in unique(dt[,NN]))
    {
      rr = dt[obs == oo & bw == bw0 & NN == N,range(score)]
      fac = 
      addcon = -rr[1]
      
      dt[obs == oo & bw == bw0 & NN == N, fac := 1/(rr[2] - rr[1])]
      dt[obs == oo & bw == bw0 & NN == N, addcon := - rr[1]]
      dt[obs == oo & bw == bw0 & NN == N, scaled_score := (score + addcon) * fac]
    }
  }
}

# plotting

for(oo in 1:4) # observation model
{
  
  for(N in N_vec)
  {
    dt_temp = dt[obs == oo & NN == N]
    
    bws = as.numeric(as.character(unique(dt_temp[type != 'log',bw])))
    
    score_labels = as.expression(c('log', bquote("int., "~sigma == .(bws[1])),bquote("int., "~sigma == .(bws[2])),bquote("int., "~sigma == .(bws[3]))))
    
    bq = function(x){as.expression(bquote(x))}
    
    model_labels = as.expression(c(bquote(F[1]),bquote(F[2]),bquote(F[3]),bquote(F[4])))
    
    pp = ggplot(dt_temp) + 
      geom_boxplot(aes(x = bw,y = scaled_score,fill = as.factor(pr)),position=position_dodge(0.5),width = 0.33) +
      scale_x_discrete(labels = score_labels,name = '') +
      scale_y_discrete(labels = '',name = 'mean score',expand = c(0.1,0.1)) + 
      scale_fill_discrete(name = 'pred. model', labels = model_labels) + 
      ggtitle(paste0('N = ',N))
    
    # adjust so that they look nice when plotted next to each other:remove legend from first and adjust plot width:
    if(N == 10) 
    {pp = pp + theme(legend.position = 'none')
      ww = 8.5 # width of pdf
    }
    if(N == 100) ww = 10
    
    pdf(paste0(pdf_dir,'score_poisson_obsmod',oo,'_N',N,'.pdf'),width = ww, height = 7)
      print(pp)
    dev.off()
  }
}


# #############################################
# ########### logarithmic score ###############
# #############################################
# 
# #### generate boxplot for log-score without bootstrap: ####
# 
# p_log = ggplot(dt[type == 'log' & obs == 1 & N_sam == 1, ], aes( factor(pr),  score)) + geom_boxplot()
# 
# p_log = p_log + labs(title = 'log-score',
#                      x = 'pred. model',
#                      y = expression('score'))  + 
#   theme(plot.title = element_text(size = title_size),
#         axis.title = element_text(size = labels_size),
#         axis.text = element_text(size = ticks_size))
# 
# # check:
# 
# print(p_log)
# 
# #generate plot:
# 
# pdf(paste0('../../Figures/log_score_poisson.pdf'))
#   print(p_log)
# dev.off()
# 
# 
# ###### generate bootstrap plots: ##########
# 
# for(N in N_vec)
# {
#   #### generate boxplot for mean log-score with bootstrap: ####
#   
#   p_log_m = ggplot(dt[type == 'log' & obs == 1 & N_sam == N, ], aes( factor(pr),  score)) + geom_boxplot()
#   
#   # custom title and labels:
#   
#   p_log_m = p_log_m + labs(title = paste0('mean log-score, ',N,' obs.'),
#                        x = 'pred. model',
#                        y = expression('mean score'))  + 
#     theme(plot.title = element_text(size = title_size),
#           axis.title = element_text(size = labels_size),
#           axis.text = element_text(size = ticks_size))
#   
#   
#   # check: 
#   
#   print(p_log_m)
#   
#   # generate plot
#   
#   pdf(paste0('../../Figures/mean_log_score_poisson_N',N,'.pdf'))
#     print(p_log_m)
#   dev.off()
#   
# }
# 
# 
# 
# #############################################
# ########### intensity score ###############
# #############################################
# 
# #### generate boxplot for kernel score without bootstrap: ####
# 
# p_kernel = ggplot(dt[type == 'kernel' & obs == 1 & N_sam == 1], aes( factor(pr),  score)) +  geom_boxplot()
# 
# # custom title and labels:
# 
# p_kernel = p_kernel + labs(title = 'intensity score',
#                            x = 'pred. model',
#                            y = expression('score'))  + 
#   theme(plot.title = element_text(size = title_size),
#         axis.title = element_text(size = labels_size),
#         axis.text = element_text(size = ticks_size))
# 
# 
# # check: 
# 
# print(p_kernel)
# 
# # generate plot:
# 
# pdf(paste0('../../Figures/kernel_score_poisson.pdf'))
#   print(p_kernel)
# dev.off()
# 
# #### generate boxplots for mean kernel score with bootstrap: ####
# 
# for(N in N_vec)
# {
#   
#   p_kernel_m = ggplot(dt[type == 'kernel' & obs == 1 & N_sam == N ], aes( factor(pr),  score)) + geom_boxplot()
#   
#   # custom title and labels:
#   
#   p_kernel_m = p_kernel_m + labs(title = paste0('mean intensity score, ',N,' obs.'),
#                            x = 'pred. model',
#                            y = expression('mean score'))  + 
#     theme(plot.title = element_text(size = title_size),
#           axis.title = element_text(size = labels_size),
#           axis.text = element_text(size = ticks_size))
#   
#   
#   # check: 
#   
#   print(p_kernel_m)
#   
#   # generate plots:
#   
#   pdf(paste0('../../Figures/mean_kernel_score_poisson_N',N,'.pdf'))
#     print(p_kernel_m)
#   dev.off()
#  
# }
# 
