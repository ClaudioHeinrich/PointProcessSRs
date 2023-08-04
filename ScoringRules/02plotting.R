

# make fancy plots

rm(list = ls())

library(data.table)
library(spatstat)
library(ggplot2)
library(latex2exp)


setwd('~/NR/ProjectPointProcess/PointProcessSRs/ScoringRules/')

set.seed(102030)


load(file = 'setup_r1.RData')
load(file = 'scores_r1.RData')

# base font size:

bs = 24

### plots ###

# kernel estimator score:

l_p = ggplot(data = score_dt[type == 'lambda'],aes(x = obs, y = score,color = pr,shape = pr,fill = pr)) + 
  geom_point(size = 5,alpha = 0.7)


# manually adjust scales:

# old configuration:
# shapes = c(16,17,15,23,25)
# cols = c('black','goldenrod','forestgreen','red','purple')

# new configuration, trying to make 'similar processes' look similar
shapes = c(23,15,17,16,25,7)
cols = c('forestgreen','blue2','purple','goldenrod','firebrick3','darkblue')
perm = 1:6 # experiment by permuting the symbol and color assignment

l_p = l_p + scale_shape_manual(values = shapes[perm]) + 
  scale_color_manual(values = cols[perm]) + 
  scale_fill_manual(values = cols[perm])


l_p = l_p +labs(title = 'mean intensity score',
                x = 'true dist.',
                y = 'mean score',
                color = 'pred. dist.',
                shape = 'pred. dist.',
                fill = 'pred. dist.')

# mark correct model

l_p = l_p + geom_point(data = score_dt[type == 'lambda' & obs == pr],
                       mapping = aes(x = obs,y = score),
                       shape = 3,size = 2.5,colour = 'black',show.legend = FALSE,
                       position = position_nudge(x = -0.2,y = 0))

# make background white and adjust text sizes

l_p = l_p + theme_bw(base_size = bs)

# check

print(l_p)

# generate plot
ggsave(l_p,file = '../../ModelEvaluation/Figures/review1/lambda_scores.pdf',height = 7)




# K-function score:

K_p = ggplot(data = score_dt[type == 'K'],aes(x = obs, y = score,color = pr,shape = pr,fill = pr)) + 
  geom_point(size = 5,alpha = 0.7)

K_p = K_p + scale_shape_manual(values = shapes[perm]) + 
  scale_color_manual(values = cols[perm]) + 
  scale_fill_manual(values = cols[perm])



K_p = K_p +labs(title = 'mean K-function score',
                x = 'true dist.',
                y = 'mean score',
                color = 'pred. dist.',
                shape = 'pred. dist.',
                fill = 'pred. dist.') 
 

# encircle correct model

K_p = K_p + geom_point(data = score_dt[type == 'K' & obs == pr],
                       mapping = aes(x = obs,y = score),
                       shape = 3,size = 2.5,colour = 'black',show.legend = FALSE,
                       position = position_nudge(x = -0.2,y = 0))

# make background white

K_p = K_p + theme_bw(base_size = bs) 
# check

print(K_p)

# generate plot
ggsave(K_p,file = '../../ModelEvaluation/Figures/review1/K_scores.pdf',height = 7)


#Print out permutation tests: