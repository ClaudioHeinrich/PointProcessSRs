# this script tempers a bit with internal spatstat functions in order to make it possible to pass an fv object to kppm, rather than a ppp.
# this allows it to fit a cluster model to multiple observations using minimum contrast, by taking the mean of the K-function over all observations.
# The functions defined here are much less robust than the original spatstat functions, but work on the tree simulation study.

###################################################

mean_Kest = function(lppp,rmax = NULL,...)
{
  # Takes a list lppp of ppps and returns the mean Kest as fv, mimicing to be a K-function such that it can be passed to clusterfit. rmax and ... are passed on to Kest.
  
  np = length(lppp)
  
  
  # Attention! The r values contained in the fv returned by Kest applied to all observations need to be identical! 
  # This is e.g. satisfied when all observations have the same support window owin, check Kest and rmax.rule
  if(is.null(rmax)){
    win1  = lppp[[1]]$window
    for(i in 2:np)
    {
      if(!identical(lppp[[i]]$window,win1)) stop('When the point patterns have different support you need to specify rmax')
    }
  }
  
  
  # get K function estimate = mean of K function etimates for all windows.
  K_hat = 0 
  for(i in 1:np)
  {
    
    k_new = Kest(lppp[[i]],rmax = rmax,...)
    
    K_hat = K_hat + k_new/np
  }
  
  # reset r
  
  K_hat$r = k_new$r
  
  K_hat = as.fv(K_hat)
  
  # set fvnames of K_hat, important for correct plotting behaviour and for further use of K_hat
  for(a in c('.x','.y','.'))
  {
    fvnames(K_hat,a)<-fvnames(k_new,a)
  }
  # set all atributes 
  attributes(K_hat) = attributes(k_new)
  
  return(K_hat)
}

######################################################

mkppmMinCon = function (X, Xname, Kfv, po, lambda = NULL, clusters,rmax = NULL) 
{
  # alteration of the spatstat internal function kppmMinCon. This version takes a K-function estimate (function value object) which can originate from multiple observations.
  # Works only for homogeneous models.
  
  stationary <- TRUE
  
  mcfit <- clusterfit(Kfv, clusters, lambda = lambda, dataname = Xname)
  
  fitinfo <- attr(mcfit, "info")
  attr(mcfit, "info") <- NULL
  
  Fit <- list(method = "mincon", statistic = 'K', Stat = fitinfo$Stat, 
              StatFun = fitinfo$StatFun, StatName = fitinfo$StatName, 
              FitFun = fitinfo$FitFun, statargs = list(), mcfit = mcfit)
  out <- list(Xname = Xname, X = X, stationary = TRUE, 
              K = Kfv, 
              clusters = clusters, modelname = fitinfo$modelname, 
              isPCP = fitinfo$isPCP, po =po, lambda = lambda, 
              mu = mcfit$mu, par = mcfit$par, par.canon = mcfit$par.canon, 
              clustpar = mcfit$clustpar, clustargs = mcfit$clustargs, 
              modelpar = mcfit$modelpar, covmodel = mcfit$covmodel, 
              Fit = Fit)
  
  return(out)
}

######################################################


mkppm = function (Kfv,lambda = NULL, X,
                  trend = ~1, 
                  clusters = c("Thomas", "MatClust", "Cauchy", "VarGamma", "LGCP"), 
                  rmax = max(Kfv$r)) 
{
  #fits a homogeneous cluster point process to a K-function using the method of minimum contrast. The K-function can be estimated from multiple observations.
  # the variable X is supposed to contain one of the observations and is used as data sample etc. for the output which ought to be of type 'kppm'
  cl <- match.call()
  callstring <- paste(short.deparse(sys.call()), collapse = "")
  Xname <- short.deparse(substitute(X))
  clusters <- match.arg(clusters)
  
  if(trend != ~1) stop('At the current time, only trend = ~1 is supported.')
  
  statistic <- 'K'
  
  ClusterArgs <- list(method = 'minkon', 
                      improve.type = 'none',
                      improve.args = list(),
                      weightfun = NULL, 
                      control = list(),
                      algorithm = 'Nelder-Mead',
                      statistic = statistic, 
                      statargs = list(), 
                      rmax = rmax)
  
  po <- ppm(Q = X, trend = trend, covariates = NULL, 
            forcefit = TRUE, rename.intercept = FALSE, covfunargs = NULL, 
            use.gam = FALSE, nd = NULL, eps = NULL)
  
  out <- mkppmMinCon(X = X, Xname = Xname,po = po, Kfv = Kfv,lambda = lambda, clusters = clusters, rmax = rmax)
  
  out <- append(out, list(ClusterArgs = ClusterArgs, call = cl, 
                          callframe = parent.frame(), callstring = callstring))
  
  class(out) = 'kppm'
  
  return(out)
}

