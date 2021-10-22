
rm(list = ls())

library(data.table)
library(spatstat)
library(ggplot2)
library(latex2exp)

setwd('~/NR/ProjectPointprocess/ModelEvaluation/simstudy/LTest')

set.seed(102030)

# generate whether results of L-test should be plotted or of S-test/M-test

type = 'Ltest' # takes 'Ltest' or SMtest

load(file = paste0(type,'.RData'))

# specify font sizes

title_size = 20
labels_size = 16
ticks_size = 12

# expected quantile score plot:

gam_p = ggplot(data = dt,aes(x = model,
                             y = gamma,
                             color = factor(observ),
                             shape = factor(observ))) + geom_line() + geom_point(size = 2.5) 

# add custom title

gam_p = gam_p + labs(title = 'mean quantile score',
                     x = 'predictive model',
                     y = expression(paste('mean(',gamma,')')),
                     color = 'true model',
                     shape = 'true model') + 
  theme(plot.title = element_text(size = title_size),
        axis.title = element_text(size = labels_size),
        axis.text = element_text(size = ticks_size),
        legend.title = element_text(size = labels_size),
        legend.text = element_text(size = ticks_size))


# encircle true model
gam_p = gam_p + geom_point(mapping = aes(x = model,y = gamma,color = factor(observ)), 
             data = dt[model == observ], shape = 1,size = 6,show.legend = FALSE)



print(gam_p)

# generate plots

pdf(paste0('../../Figures/',type,'gamma.pdf'))
  print(gam_p)
dev.off()

#something along the following lines should work to replace the tick labels by F_1,...,F_4, but it doesn't

# pred_mod = c(expression(F[1]),expression(F[2]),expression(F[3]),expression(F[4]))
# 
# gam_p  = gam_p + scale_x_discrete(labels = pred_mod)

#### rejection probability plot: ####

rej_p = ggplot(data = dt,aes(x = model,
                             y = fails,
                             color = factor(observ),
                             shape = factor(observ))) + geom_line() + geom_point(size = 2.5)


# add custom title/labels

rej_p = rej_p + labs(title = 'rejection probability',
                     x = 'predictive model',
                     y = '',
                     color = 'true model',
                     shape = 'true model') + 
  theme(plot.title = element_text(size = title_size),
        axis.title = element_text(size = labels_size),
        axis.text = element_text(size = ticks_size),
        legend.title = element_text(size = labels_size),
        legend.text = element_text(size = ticks_size))


#encircle true model

rej_p = rej_p + geom_point(mapping = aes(x = model,y = fails,color = factor(observ)), 
                           data = dt[model == observ], shape = 1,size = 6,show.legend = FALSE)

# add abline

rej_p = rej_p + geom_hline(yintercept = 0.05, linetype = 2, alpha = 0.75, show.legend = FALSE)

# check result

print(rej_p)

# generate plot

pdf(paste0('../../Figures/',type,'rejp.pdf'))
  print(rej_p)
dev.off()


######################################

#### rejection probability plot for two-sided test: ####

if(type == 'Ltest') title_str = 'rejection prob., 1st simulation study'
if(type == 'SMtest') title_str = 'rejection prob., 2nd simulation study'

rejts_p = ggplot(data = dt,aes(x = model,
                               y = fails_ts,
                               color = factor(observ),
                               shape = factor(observ))) + geom_line() + geom_point(size = 2.5)


# add custom title/labels

rejts_p = rejts_p + labs(title = title_str,
                         x = 'predictive model',
                         y = '',
                         color = 'true model',
                         shape = 'true model') + 
  theme(plot.title = element_text(size = title_size),
        axis.title = element_text(size = labels_size),
        axis.text = element_text(size = ticks_size),
        legend.title = element_text(size = labels_size),
        legend.text = element_text(size = ticks_size))


#encircle true model

rejts_p = rejts_p + geom_point(mapping = aes(x = model,y = fails_ts,color = factor(observ)), 
                               data = dt[model == observ], shape = 1,size = 6,show.legend = FALSE)

# add abline

rejts_p = rejts_p + geom_hline(yintercept = 0.05, linetype = 2, alpha = 0.75, show.legend = FALSE)

# check result

print(rejts_p)

# generate plot

pdf(paste0('../../Figures/',type,'rejpts.pdf'))
  print(rejts_p)
dev.off()