library(progress)

# Only return the minimized MISE and the bandwidth associated
NW.minimizeMISE = function(Kernel = Knorm, h.test = seq(0.01, 0.15, 0.001), M = 100, ni = 100){
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
# Haven't found a way to compute MISE with LM
# But as MSSE ~= MISE, this is fine.
# This one returns MSSE, bandwidth, Bias and also the mean of the regressions (for variance)
LM.minimizeMSSE = function(poly.test = seq(1, 20, 1), M = 200, ni = 100, raw=F){
  
  boot.X = runif(ni)
  res = data.frame(MSSE=Inf, poly=NA, Bias=matrix(NA, ni, 1), MReg=matrix(NA, ni, 1))
  
  pb = progress_bar$new(
    format = "[LM MSSE] Best: :bmsse (:bp) | Current: :cmsse (:cp) [:bar] :current/:total |:percent | :elapsed", 
    total = length(poly.test)
  )
  Sys.sleep(0.1)
  pb$tick(0)
  for (p in poly.test){
    SE = matrix(NA, M, ni)
    Reg = matrix(NA, M, ni)
    
    for (i in 1:M){
      boot.Y = Y(boot.X)
      boot.reg = lm(boot.Y ~ poly(boot.X, degree = p, raw=raw))$fitted
      SE[i,] = (boot.reg - m(boot.X))^2
      Reg[i,] = boot.reg
    }
    
    MReg = colMeans(Reg)
    Bias = MReg - m(boot.X)
    MSE = colMeans(SE)
    MSSE = round(mean(MSE),6)
    
    if (MSSE < mean(res$MSSE)){
      res$MSSE = MSSE
      res$poly = p
      res$Bias = Bias
      res$MReg = MReg
    }
    pb$tick(tokens=list(bmsse=res$MSSE, bp=res$poly, cmsse=MSSE, cp=p))
  }
  return (res)
}
