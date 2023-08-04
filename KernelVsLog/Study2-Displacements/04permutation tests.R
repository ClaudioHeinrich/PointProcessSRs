
### run permutation tests for score differences ###


### setup ###

rm(list = ls())

library(spatstat)

library(data.table)

library(ggplot2)


setwd('~/NR/ProjectPointProcess/PointProcessSRs/KernelVsLog/Study2-Displacements/')

set.seed(102030)

load('setup.RData')
pdf_dir = paste0('../../../ModelEvaluation/Figures/review1/study2/displacements/nexp',nexps,'_new/')
dir.create(pdf_dir,showWarnings = F)



# load scores:

load(file = paste0('raw_scores.RData'))


############ function for permutation tests ##################

#' Run a permutation test of the pairwise difference between two vectors of numbers
#' @param a Vector, the scores from one method
#' @param b Vector, the scores from some other method
#' @param N Integer, the size of the permutation distribution
#' @return A list with the mean of the difference and the permutation distribution of that difference and the p-value
#' @examples
#' N = 1e2
#' trend  = 1:N
#' a = trend + .01 + rnorm(N, .001)
#' b = trend - .01 + rnorm(N, .001)
#' l = permutation_test_difference(a,b)
#' q = sum(l$D <= l$d_bar) / length(l$D)
#' @author Alex,Claudio
#' 
#' @export
#' 
permutation_test_difference = function(a,
                                       b,
                                       N = 5e3){
  n = length(a)
  d = a - b
  d_bar = mean(d)
  D = NULL
  for(i in 1:N){
    swap = rbinom(n,1,0.5)
    w_swap = which(swap == 1)
    d_i = d
    d_i[w_swap] = -d_i[w_swap]
    D[i] = mean(d_i)
  }
  
  p_val = sum(d_bar > D)/N + sum(d_bar == D)/(2*N)
  
  return(list(d_bar = d_bar, D = D,p_val = p_val))
}


#################################

# run permutation tests: consider means of N_pt scores (sampled from within dt) and see how often, out of R times, a permutation test at level alpha rejects. 
# This divided by R approximates the power of the test


dt = melt(log_scores,variable.name = 'mod_com') 


get_obs_mod = function(x)
{
  substring(x,9)
}

get_pr_mod = function(x)
{
  substring(x,4,4)
}


dt[,obs := get_obs_mod(mod_com)][,pr:=get_pr_mod(mod_com)]

dt[,mod_com := NULL][,type := 'log'][,bw := 0]

for(ii in 1:length(int_scores))
{
  dt_temp = int_scores[[ii]]$lambda_scores
  dt_temp = melt(dt_temp,variable.name = 'mod_com') 
  dt_temp[,obs := get_obs_mod(mod_com)][,pr:=get_pr_mod(mod_com)]
  
  dt_temp[,mod_com := NULL][,type := 'int'][,bw := int_scores[[ii]]$bandwidth]
  dt = rbindlist(list(dt,dt_temp),use.names = T)
  
}


oo = 1 # observation model
N_pt = 100 # how many permutation tests are run
N_sample = 100 # how many draws used in each permutation test

p_vals = as.data.table(expand.grid(pr_mod = 1:n_mod,N_obs = 1:15,ind = 1:N_pt,obs_mod = oo,bw = unique(dt[,bw])))

p_vals[bw == 0, type := 'log']
p_vals[bw !=0, type := 'int']

p_vals = p_vals[pr_mod != obs_mod]


# get random draws:
N_max = log_scores[,.N]

all_inds = c()
counter = 0
while(counter < N_pt)
{
  inds = sample(N_max,size = max(N_vec)) 
  all_inds = c(all_inds,inds)
  counter = counter + 1
}

all_inds = matrix(all_inds,nrow = max(N_vec))

get_pval = function(row_ind)
{
  sample_ind = p_vals[row_ind,ind]
  nobs = p_vals[row_ind,N_obs]
  pm = p_vals[row_ind,pr_mod]
  om = p_vals[row_ind,obs_mod]
  sc = p_vals[row_ind,type]
  bw0 = p_vals[row_ind,bw]
  
  inds = as.vector(all_inds[,sample_ind])[1:nobs]
  
  scores_pm = dt[pr == pm & obs == om & type == sc & bw == bw0][inds,value]
  scores_om = dt[pr == om & obs == om & type == sc & bw == bw0][inds,value]
  
  pval = permutation_test_difference(a = scores_om, b = scores_pm, N = N_sample)$p_val
  p_vals[row_ind, p_val := pval]
}


for(i in 1:p_vals[,.N])
{
  if(i %% 500 == 0)print(paste0(i,'/',p_vals[,.N]))
  get_pval(i)
}


p_val_means = p_vals[,mean(p_val),.(pr_mod,N_obs,obs_mod,bw,type)]

bws = unique(p_val_means[,bw])


p_val_means[,bw := as.factor(bw)]

score_labels = as.expression(c('log',
                               bquote(sigma == .(bws[2])),
                               bquote(sigma == .(bws[3])),
                               bquote(sigma == .(bws[4])),
                               bquote(sigma == .(bws[5])),
                               bquote(sigma == .(bws[6]))))

# it's rather complicated to get expressions into the facet labeling, but let's go

p_val_means[,pr_mod := as.factor(pr_mod)]

model_labels = as.expression(c(bquote(mu[x] == .(mu2[1])),
                               bquote(eta == .(sig3)),
                               bquote(eta == .(sig4)),
                               bquote(rho == .(cor5)),
                               bquote(N == .(nexps6))))


model_labels = as.expression(c(bquote(mu[x] == .(mu2[1])),
                               bquote(eta == .(sig3)),
                               bquote(eta == .(sig4)),
                               bquote(rho == .(cor5)),
                               bquote(N == .(nexps6))))

model_labels_new = TeX(c(paste0('$F_2:\\ \\mu_x = ',mu2[1],'$'),
                         paste0('$F_3:\\ \\eta = ',sig3,'$'),
                         paste0('$F_4:\\ \\eta = ',sig4,'$'),
                         paste0('$F_5:\\ \\rho = ',cor5,'$'),
                         paste0('$F_6:\\ N = ',nexps6,'$')))

p_val_means$pr_mod = factor(p_val_means$pr_mod,
                            levels = 2:6,
                            labels = model_labels_new)



# specify font sizes

bws = 0.1 * 2^(0:4)

theme_set(theme_bw(base_size = 32))

pdf(paste0(pdf_dir,'pvals.pdf'),width = 25)
  
  pp = ggplot(p_val_means[bw %in% c(0,bws)]) + 
    geom_line(aes(x = N_obs,y = V1,color = bw, linetype = bw),linewidth = 1) + 
    geom_hline(yintercept = 0.05,linetype = 'dashed',alpha = 0.75) + 
    facet_grid(cols = vars(pr_mod),labeller = label_parsed) +
    xlab('Number of observed point patterns') + ylab('mean p-value') +
    scale_color_discrete(name = 'score',labels = score_labels) + 
    scale_linetype(name = 'score',labels = score_labels) +
    theme(legend.text.align = 0)
  print(pp)

dev.off()


save.image('pvalspermtest.RData')

