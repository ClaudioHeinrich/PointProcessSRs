
# generates the boxplots and scatterplots for the paper

rm(list = ls())

library(spatstat)

library(data.table)

library(ggplot2)

setwd('~/NR/ProjectPointProcess/ModelEvaluation/simstudy/KernelVsLog/Study2-Displacements/')

set.seed(123)

load('setup.RData')

# specify font sizes

theme_set(theme_bw(base_size = 22))

# load data

load(file = paste0('scores.RData'))

pdf_dir = paste0('../../../Figures/study2/displacements/nexp',nexps,'_new/')
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
      
      dt[obs == oo & bw == bw0 & NN == N, fac := 1/(rr[2] - rr[1])]
      dt[obs == oo & bw == bw0 & NN == N, addcon := - rr[1]]
      dt[obs == oo & bw == bw0 & NN == N, scaled_score := (score + addcon) * fac]
    }
  }
}

# plotting

# for(oo in 1:6) # observation model
# {

oo = 1
  
  for(N in N_vec)
  {
    dt_temp = dt[obs == oo & NN == N]
    
    bws = as.numeric(as.character(unique(dt_temp[type != 'log',bw])))
    
    score_labels = as.expression(c('log', 
                                   bquote(~sigma == .(bws[1])),
                                   bquote(~sigma == .(bws[2])),
                                   bquote(~sigma == .(bws[3])),
                                   bquote(~sigma == .(bws[4])),
                                   bquote(~sigma == .(bws[5]))))
    
    bq = function(x){as.expression(bquote(x))}
    
    model_labels = as.expression(c(paste0('N(0,1), N = ',nexps),
                                   bquote(mu == .(mu2)),
                                   bquote(eta == .(sig3)),
                                   bquote(eta == .(sig4)),
                                   bquote(rho == .(cor5)),
                                   paste0('N = ',nexps6)
                                   ))
    
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


save.image(file = paste0(pdf_dir,'image.RData'))

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

######## generate scatterplots #########

dt_sub = dt[NN == 10 & obs == 1]

pred_mod = 5
dt1 = dt_sub[pr == 1]
dt2 = dt_sub[pr == pred_mod]

setkey(dt1,type,bw)
setkeyv(dt2,key(dt1))

setnames(dt1,'scaled_score','xscore')

dt1 = data.table(dt1,dt2[,.(scaled_score)])
setnames(dt1,'scaled_score','yscore')


N = 10

bws = c(0.1,0.4,1.6) #  as.numeric(as.character(unique(dt_sub[type != 'log',bw])))

dt1 = dt1[bw %in% c(0,bws)]

score_labels = as.expression(c('log', 
                                 bquote(sigma == .(bws[1])),
                                 bquote(sigma == .(bws[2])),
                                 bquote(sigma == .(bws[3]))))
  
setkey(dt1,bw)

dt1$bw = factor(dt1$bw,labels = score_labels)

theme_set(theme_bw(base_size = 36))

pp = ggplot(dt1) + geom_point(aes(x = xscore,y = yscore)) + facet_grid(cols = vars(bw), labeller = label_parsed)

pp = pp + theme(axis.text = element_blank(),axis.ticks = element_blank())

pp = pp + xlab('mean score under correct prediction') + ylab(bquote('mean score,  '~rho == 0.1))

pdf(paste0(pdf_dir,'scatter_plot_model5_N',N,'.pdf'),width = 28, height = 7)
  print(pp)
dev.off()
  
