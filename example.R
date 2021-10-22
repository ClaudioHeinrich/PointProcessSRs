



rm(list = ls())

library(spatstat)

library(data.table)


# area:

limx = c(0,10)
limy = c(0,10)

x_range = max(limx) - min(limx)
y_range = max(limy) - min(limy)

win = owin(xrange = limx,yrange = limy)


n_exp = 50

pp1 = runifpoint(n = rpois(1,lambda = n_exp), win = win)


# intensity for nonhomogeneous process:

lambda0 = function(x,y){return(sqrt(x^2+y^2))}

int_val = pracma::integral2(lambda0,
                            xmin = limx[1],xmax = limx[2],
                            ymin = limy[1],ymax = limy[2])

lambda =  function(x,y){return(n_exp * sqrt(x^2+y^2)/int_val$Q)}


pp2 = rpoispp(lambda = lambda,win = win)


# Strauss process: find approximately for which beta n_exp is the expected number of points 

N = 50

exp_n_Strauss = function(beta) # returns the mean squared error between n_exp and the number of points of a Strauss model with parameter beta 
{
  ret_val = 0
  for(i in 1:N){ret_val = ret_val + (1/N) * (rStrauss(beta = beta, gamma = 0.2, R = 0.75, W = win)$n - n_exp)^2}
  return(ret_val)
}

beta_vec = seq(0.5,2,by = 0.1)

exp_n = c()
for (b in beta_vec)
{
  print(b)
  exp_n = c(exp_n,exp_n_Strauss(b))
}

plot(beta_vec,exp_n)

beta = beta_vec[which.min(exp_n)] 


pp3 = rStrauss(beta = beta, gamma = 0.2, R = 0.75, W = win)


# draw from each point pattern N times, compute a kernel density estimator and compare to the estimation from the predictive distribution

models = c('homPois','inhomPois','Strauss','inhomNS')

# plot examples:

pp1 = runifpoint(n = rpois(1,lambda = n_exp), win = win)
pp2 = dat = rpoispp(lambda = lambda,win = win)
pp3 =  dat = rStrauss(beta = beta, gamma = 0.2, R = 0.75, W = win)

nclust <- function(x0, y0, radius, n) {
  return(runifdisc(n, radius, centre=c(x0, y0)))
}

m = 4

pp4 = rThomas(kappa = function(x,y){return( (1/m)*lambda(x,y))},scale = 1, mu = m,win = win)  


pp5 = runifpoint(n = rpois(1,lambda = 1.2*n_exp), win = win)





cex_val = 1.5
par(mar = c(1,1,1,1))

setwd('~/NR/ProjectPointProcess/Figures')

pdf('homPois.pdf')
par(mar = c(1,1,1,1))
par(cex = cex_val)

plot(pp1,main = 'hP')

dev.off()

pdf('inhomPois.pdf')
par(mar = c(1,1,1,1))
par(cex = cex_val)

plot(pp2,main = 'ihP')

dev.off()


pdf('Strauss.pdf')
par(mar = c(1,1,1,1))
par(cex = cex_val)

plot(pp3,main = 'Str')
par(mar = c(1,1,1,1))
dev.off()

pdf('ihThomas.pdf')
par(mar = c(1,1,1,1))
par(cex = cex_val)

plot(pp4,main = 'ihTh')

dev.off()


pdf('homPoisPlus.pdf')
par(mar = c(1,1,1,1))
par(cex = cex_val)

plot(pp5,main = 'hP+')

dev.off()
