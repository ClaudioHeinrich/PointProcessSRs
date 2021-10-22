
# compute the AIC for the predictive models

rm(list = ls())

set.seed(102030)

library(spatstat)
library(data.table)
library(parallel)

setwd('~/NR/ProjectPointProcess/ModelEvaluation/simstudy/TreeData')

data_dir = './Data'
dir.create(data_dir,showWarnings = FALSE)
plot_dir = './Figures'

# setwd('~/R/PointProcess/Trees') # for working on server


#### directories ####

data_dir = './Data/'
function_dir = './functions/'

models = c('Pois','LGCP','Thomas','MatClust','VarGamma','Cauchy')


###############################################################
### Plot data from 1978, 1990, and 2009 in one large figure ###
###############################################################

load(paste0(data_dir,"trees.Rdata"))

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


# # plot 
# 
# par(mfrow=c(4,3),mex=0.25,mar=c(1,1,1,1)+0.01)
# for(i in 1:4) {
#   plot(x.78[[i]],y.78[[i]],type="p",pch=16,cex=0.5,
#        xlab="",ylab="",xlim=c(0,6),ylim=c(0,6),xaxs="i",
#        yaxs="i",axes=FALSE)
#   box()
#   plot(x.90[[i]],y.90[[i]],type="p",pch=16,cex=0.5,
#        xlab="",ylab="",xlim=c(0,6),ylim=c(0,6),xaxs="i",
#        yaxs="i",axes=FALSE)
#   box()
#   plot(x.09[[i]],y.09[[i]],type="p",pch=16,cex=0.5,
#        xlab="",ylab="",xlim=c(0,6),ylim=c(0,6),xaxs="i",
#        yaxs="i",axes=FALSE)
#   box()
# 
# }


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

# get hyperframe

lpp_temp = list()
for(y_ind in 1:3)
{
  x_data = get(paste0('x.',years_char[y_ind]))
  y_data = get(paste0('y.',years_char[y_ind]))
  
  for(i in 1:8)
  {
    pp = list(ppp(x_data[[i]],y_data[[i]],xrange=c(0,6),yrange=c(0,6)))
    lpp_temp = c(lpp_temp,pp)
  }
}


lpph1 = as.hyperframe(expand.grid(obs = 1:8,year = years))
lpph2 = hyperframe(pp = lpp_temp)

lpph = cbind(lpph1,lpph2 )


################
### get AICs ###
################



dt_aic = as.data.table(expand.grid(model = models,tr_win = 1:8,year = years))

for(train in 1:8)
{
  for(mod in models)
  {
    print(c(train,mod))
    for(yy in years)
    {
      y_str = as.character(yy)
      
      Tr = lpp[[y_str]][[train]]
      
      if(mod == 'Pois')
      {
        a = ppm(Tr ~ 1, interaction = NULL)
      }
      if(mod != 'Pois')
      {
        a = kppm(Tr ~ 1, clusters =  mod, method = 'palm')   
      }
      
      dt_aic[year == yy & tr_win == train & model == mod, aic:= AIC(a)]
      
    }
  }
}



dt_aic[,mean(aic),by = .(year,model)]

# check out aics for each year:

for(yy in years)
{
  print(paste0('year = ',yy,':'))
  print(dt_aic[year == yy,mean(aic),by = .(year,model)])
} 


save(dt_aic,file = paste0(data_dir,'aic.RData'))

# yy = 2009
# 
# data = dt[year == yy & est == 'int']
# 
# p = ggplot(data = data,aes(x = Model,y = score,color = Model)) + geom_boxplot()
# 
# print(p)
