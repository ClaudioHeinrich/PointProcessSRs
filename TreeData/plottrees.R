# This script plots trees.

rm(list = ls())

set.seed(102030)

library(spatstat)
library(data.table)
library(boot)
library(ggplot2)

setwd('~/NR/ProjectPointProcess/ModelEvaluation/simstudy/TreeData')



# setwd('~/R/PointProcess/Trees') # for working on server


#### directories ####

data_dir = './Data/'
function_dir = './functions/'
plot_dir = '../../Figures/'



###############################################################
### Plot data from 1978, 1990, and 2009 in one large figure ###
###############################################################

load(paste0(data_dir,"trees.Rdata"))

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

# plot 2 example plots:

pdf(file="../../Figures/Treestwo.pdf",height=4,width=6)
  par(mfrow=c(2,3),mex=0.25,mar=c(1,1,1,1)+0.01)
  for(i in c(1,3)) 
    {
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
dev.off()

