library(progress)

# Only return the minimized MISE and the bandwidth associated
NW.minimizeMISE = function(Kernel = Knorm, h.test = seq(0.01, 0.15, 0.001), M = 100, ni = 100){
  set.seed(122020)
  best = Inf
  best_h = 0
  
  pb = progress_bar$new(
    format = "[NW MISE] Best: :bmise (:bh) | Current: :cmise (:ch) [:bar] :current/:total |:percent | :elapsed", 
    total = length(h.test)
  )
  Sys.sleep(0.1)
  pb$tick(0)

  for (h in h.test){
    
    ISE = matrix(NA, M, 1) 
    
    for (i in 1:M){
      boot.X = runif(ni)
      boot.Y = Y(boot.X)
      boot.ISE = integrate(
        function(xi) (sapply(xi, function(x) NW.regEst(x, boot.X, boot.Y, h, Kernel)) - m(xi))^2,
        lower=0, upper=1)$value
      ISE[i,] = boot.ISE
    }
    
    MISE = round(mean(ISE), 6)
    
    if (MISE < best){
      best = MISE
      best_h = h
    }
    pb$tick(tokens=list(bmise=best, bh=best_h, cmise=MISE, ch=h))
  }
  sprintf("[NW MISE]: %.6f with %.6f as bandwhidth", best, best_h)
  return (c(best, best_h))
}

# This one returns MSSE, bandwidth, Bias and also the mean of the regressions (for variance)
NW.minimizeMSSE = function(Kernel = Knorm, h.test = seq(0.01, 0.15, 0.001), M = 100, ni = 100){
  set.seed(122020)
  boot.x = seq(0, 1, length = ni)
  res = data.frame(MSSE=Inf, h=NA, Bias=matrix(NA, ni, 1), MReg=matrix(NA, ni, 1))
  
  pb = progress_bar$new(
    format = "[NW MSSE] Best: :bmsse (:bh) | Current: :cmsse (:ch) [:bar] :current/:total |:percent | :elapsed", 
    total = length(h.test)
  )
  Sys.sleep(0.1)
  pb$tick(0)
  for (h in h.test){
    SE = matrix(NA, M, ni)
    Reg = matrix(NA, M, ni)
    
    for (i in 1:M){
      boot.X = runif(ni)
      boot.Y = Y(boot.X)
      boot.reg = sapply(boot.x, function(xi) NW.regEst(xi, boot.X, boot.Y, h, Kernel))
      Reg[i,] = boot.reg
      SE[i,] = (boot.reg - m(boot.x))^2
    }

    MReg = colMeans(Reg)
    Bias = MReg - m(boot.x)
    MSE = colMeans(SE)
    MSSE = round(mean(MSE),6)
    
    if (MSSE < mean(res$MSSE)){
      res$MSSE = MSSE
      res$h = h
      res$Bias = Bias
      res$MReg = MReg
    }
    pb$tick(tokens=list(bmsse=res$MSSE, bh=res$h, cmsse=MSSE, ch=h))
  }
  return (res)
}

# -----
# Haven't found a way to compute MISE with LM (yet) 
# But as MSSE ~= MISE, this is fine.
# This one returns MSSE, bandwidth, Bias, the mean of the M regressions
# And also the min and max values taken
LM.minimizeMSSE = function(poly.test = seq(1, 20, 1), M = 200){
  set.seed(122020)
  res = data.frame(MSSE=Inf, poly=NA, MSE=matrix(Inf, n, 1), Bias=matrix(Inf, n, 1), MReg=matrix(NA, n, 1), 
                   low = matrix(NA, n, 1), up = matrix(NA, n, 1))
  
  pb = progress_bar$new(
    format = "[LM MSSE] Best: :bmsse (:bp) | Current: :cmsse (:cp) | Bias: :bias [:bar] :current/:total |:percent | :elapsed", 
    total = length(poly.test)
  )
  pb$tick(0)
  
  for (p in poly.test){
    Reg = matrix(NA, M, n)
    SE = matrix(NA, M, n)
    
    for (i in 1:M){
      boot.X = sort(X(n))
      boot.Y = Y(boot.X)
      boot.reg = lm(boot.Y ~ poly(boot.X, degree = p, raw=T))
      boot.pred = predict(boot.reg, newdata=poly(CP.x, degree=p), type="response")
      Reg[i,] = boot.pred
      SE[i,] = (boot.pred - m(CP.x))**2
    }
    
    # Add "confidence interval"
    # (this is simply the min / max values taken not real CI)
    CI = matrix(NA, n, 2)
    
    for (i in 1:n){
      CI[i,1] = min(Reg[,i])
      CI[i,2] = max(Reg[,i])
    }
  
    MReg = colMeans(Reg)
    Bias = MReg - m(CP.x)
    
    MSE = colMeans(SE)
    MSSE = round(mean(MSE),6)
    
    if (mean(MSE) < mean(res$MSE)){
      res$MSE = MSE
      res$MSSE = MSSE
      res$poly = p
      res$Bias = Bias
      res$MReg = MReg
      res$low = CI[,1]
      res$up = CI[,2]
    }

    pb$tick(tokens=list(bias=sum(abs(Bias)), bmsse=res$MSSE, bp=res$poly, cmsse=MSSE, cp=p))
  }
  return (res)
}
