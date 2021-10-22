
# some data exploration, e.g. we generate plots of random samples to be compared with the true data

rm(list = ls())

library(spatstat)
library(data.table)
library(parallel)


set.seed(102030)

setwd('~/NR/ProjectPointProcess/ModelEvaluation/simstudy/TreeData')

# setwd('~/R/PointProcess/Trees') # for working on server


#### directories ####

data_dir = './Data/'
function_dir = './functions/'
plot_dir = '../../Figures/'

# get data:

load(paste0(data_dir,'scores.RData'))
load(paste0(data_dir,'trees.RData'))

source(paste0(function_dir,'summaryScore.R'))

##########################

models = c('Pois','LGCP','MatClust','Thomas','VarGamma','Cauchy')

##########################

###########################
### transform tree data ###
###########################


win = owin(xrange = c(0,6),yrange = c(0,6))

### need number of points in each plot - take only those that
### fall within [0,6]x[0,6] and are of species 1
n.x.78 <- NULL
x.78 <- list()
y.78 <- list()
for(i in 1:8) {
  ind <- which(trees[[i]]$species==1 & trees[[i]]$x>0 & trees[[i]]$x<6 &
                 trees[[i]]$y>0 & trees[[i]]$y<6)
  x.78[[i]] <- trees[[i]]$x[ind]
  y.78[[i]] <- trees[[i]]$y[ind]
  n.x.78 <- c(n.x.78,length(x.78[[i]]))
}

n.x.90 <- NULL
x.90 <- list()
y.90 <- list()
for(i in 1:8) {
  ind <- which(trees[[i]]$species==1 & trees[[i]]$x>0 & trees[[i]]$x<6 &
                 trees[[i]]$y>0 & trees[[i]]$y<6 & !is.na(trees[[i]]$D90))
  x.90[[i]] <- trees[[i]]$x[ind]
  y.90[[i]] <- trees[[i]]$y[ind]
  n.x.90 <- c(n.x.90,length(x.90[[i]]))
}

n.x.09 <- NULL
x.09 <- list()
y.09 <- list()
for(i in 1:8) {
  ind <- which(trees[[i]]$species==1 & trees[[i]]$x>0 & trees[[i]]$x<6 &
                 trees[[i]]$y>0 & trees[[i]]$y<6 & !is.na(trees[[i]]$D09))
  x.09[[i]] <- trees[[i]]$x[ind]
  y.09[[i]] <- trees[[i]]$y[ind]
  n.x.09 <- c(n.x.09,length(x.09[[i]]))
}


# plot 

par(mfrow=c(4,3),mex=0.25,mar=c(1,1,1,1)+0.01)
for(i in 1:8) {
  plot(x.78[[i]],y.78[[i]],type="p",pch=16,cex=0.5,
       xlab="",ylab="",xlim=c(0,6),ylim=c(0,6),xaxs="i",
       yaxs="i",axes=FALSE)
  box()
  plot(x.90[[i]],y.90[[i]],type="p",pch=16,cex=0.5,
       xlab="",ylab="",xlim=c(0,6),ylim=c(0,6),xaxs="i",
       yaxs="i",axes=FALSE)
  box()
  plot(x.09[[i]],y.09[[i]],type="p",pch=16,cex=0.5,
       xlab="",ylab="",xlim=c(0,6),ylim=c(0,6),xaxs="i",
       yaxs="i",axes=FALSE)
  box()
  
}

##############

years = c(1978,1990,2009)
years_char = c('78','90','09')


# transform data to list of ppps


lpp = list()

for(y_ind in 1:3)
{
  x_data = get(paste0('x.',years_char[y_ind]))
  y_data = get(paste0('y.',years_char[y_ind]))
  
  l_temp = list()
  for(i in 1:8)
  {
    pp = ppp(x_data[[i]],y_data[[i]],xrange=c(0,6),yrange=c(0,6))
    l_temp[[i]] = pp
  }
  lpp[[y_ind]] = l_temp
}


names(lpp) = years



#######################################
### get parameters of fitted models ###
#######################################


pl = 1:8

pr_mod = list()

i = 1
for(yy in years)
{
  for(plpl in pl)
  {
    print(c(yy,plpl))
    for(modmod in models)
    {
      y_str = as.character(yy)
      
      Tr = lpp[[y_str]][[plpl]]
      
      if(modmod == 'Pois')
      {
        a = ppm(Tr ~ 1, interaction = NULL) 
      } 
      if((modmod != 'Pois') & (modmod != 'LGCP2'))
      {
        a = kppm(Tr,clusters = modmod) 
      }
      if(modmod == 'LGCP2')
      {
        a = kppm(Tr,clusters = 'LGCP') 
      }
      pr_mod[[i]] = list(year = yy,pl = plpl,mod = a,mod_name = modmod)
      i = i+1
    }
  }
}
  

####################################

# generate simulations from all models, this takes some time!

N = 50 # how many do you want?

for(index in 1:length(pr_mod))
{
  print(paste0(index,'/',length(pr_mod)))
    
  a = pr_mod[[index]]$mod
  
  b = simulate(a,nsim = N)
  
  ns = NULL
  for(i in 1:N)
  {
    ns[i] = b[[i]]$n
  }
  
  pr_mod[[index]] = c(pr_mod[[index]],ns = list(ns))
}

###################

### analyse the distribution of the number of observed points ###

plot_data = data.table(expand.grid(i = 1:N, pl = 1:8,year = years,model = models ))

for(index in 1:length(pr_mod))
{
  print(index)
  yy = pr_mod[[index]]$year
  plpl = pr_mod[[index]]$pl
  modmod = pr_mod[[index]]$mod_name
  ns = pr_mod[[index]]$ns
  
  plot_data[model == modmod & year == yy & pl == plpl, n := ns]
  
  # get true number of observations:
  if(modmod == 'Pois')
  {
    plot_data[model == modmod & year == yy & pl == plpl,n_true := pr_mod[[index]]$mod$Q$data$n]
  }
}


# boxplots for number of points

library(ggplot2)

for(yy in years)
{
  for(ii in 1:8)
  {
    p = ggplot(data = plot_data[pl == ii & year == yy],mapping = aes(x = factor(model),y = n,color = factor(model))) + geom_boxplot()
    p = p + labs(title = c(ii,yy))
    
    # add horizontal line at true value:
    p = p + geom_abline(slope = 0,intercept = plot_data[pl == ii & year == yy & model == 'Pois' & i == 1,n_true])

  print(p)
  }
}



dt_test =  plot_data[,mean(n),by = .(model,year,pl)]

dt_test[,n_true := plot_data[model == 'Pois' & i == 1,n_true]]


############# look at simulations from the different models ##################

y = 1978

n = 8

n_mod = length(models)

for(obs in 1:8)
{

pdf(paste0(plot_dir,'pointsamples_obs',obs,'.pdf'),width = 56,height = 7*(n_mod +1))

  par(mfrow = c(5,8),mex=0.25,mar=c(1,1,1,1)+0.01)
  
  # plot true data:
  
  for(i in 1:n) {
    plot(x.78[[i]],y.78[[i]],type="p",pch=16,cex=1.5,
         xlab="",ylab="",xlim=c(0,6),ylim=c(0,6),xaxs="i",
         yaxs="i",axes=FALSE)
    box()
  }
  
  indices = ((obs-1) * n_mod + 1) : (obs * n_mod)
  
  for(index in indices)
  {
    print(pr_mod[[index]]$mod_name)
    for(i in 1:n)
    {
    a = pr_mod[[index]]$mod
    
    b = simulate(a,nsim = 1)[[1]]
    plot(b,pch=16,cex=1.5,
         xlab="",ylab="",xlim=c(0,6),ylim=c(0,6),xaxs="i",
         yaxs="i",axes=FALSE,main =  pr_mod[[index]]$mod_name)
    }
  }

dev.off()
} 


#### for understanding results of intensity score, plot intensities ####

for(obs in 1:8)
{
  
  pdf(paste0(plot_dir,'intensities_obs',obs,'.pdf'),width = 56,height = 7*(n_mod +1))
  
  par(mfrow = c(5,8),mex=0.25,mar=c(1,1,1,1)+0.01)
  
  
  
  for(i in 1:n) {
    ppi = ppp(x.78[[i]],y.78[[i]],window = owin(xrange = c(0,6),yrange = c(0,6)))
    
    ddi = density(ppi)
    
    plot(ddi,zlim = c(0,25),axes=FALSE,main = '')
    
  }
  
  indices = ((obs-1)*n_mod +1) : (obs*n_mod)
  
  for(index in indices)
  {
    for(i in 1:n)
    {
      a = pr_mod[[index]]$mod
      
      b = simulate(a,nsim = 1)[[1]]
      
      c = density(b)
      
      plot(c,zlim = c(0,25),axes=FALSE,main = pr_mod[[index]]$mod_name)
    }
  }
  
  dev.off()
} 


