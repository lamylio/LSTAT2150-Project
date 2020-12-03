###
# Made by Lionel Lamy
# 
# Nadayara-Watson & Linear Model with poly
# Returns: minimized MISE/MSSE and bandwidth/poly associated
##

library(progress)

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


NW.minimizeMSSE = function(Kernel = Knorm, h.test = seq(0.01, 0.15, 0.001), M = 100, ni = 100){
  best = Inf
  best_h = 0
  
  pb = progress_bar$new(
    format = "[NW MSSE] Best: :bmsse (:bh) | Current: :cmsse (:ch) [:bar] :current/:total |:percent | :elapsed", 
    total = length(h.test)
  )
  Sys.sleep(0.1)
  pb$tick(0)
  for (h in h.test){
    SE = matrix(NA, M, ni)
    for (i in 1:M){
      boot.x = seq(0, 1, length = ni) 
      boot.X = runif(ni)
      boot.Y = Y(boot.X)
      SE[i,] = (sapply(boot.x, function(xi) NW.regEst(xi, boot.X, boot.Y, h, Kernel)) - m(boot.x))^2
    }

    MSE = colMeans(SE)
    MSSE = round(mean(MSE),6)
    
    if (MSSE < best){
      best = MSSE
      best_h = h
    }
    pb$tick(tokens=list(bmsse=best, bh=best_h, cmsse=MSSE, ch=h))
  }
  sprintf("[NW MSSE]: %.6f with %.6f as bandwhidth", best, best_h)
  return (c(best, best_h))
}

# -----
# Haven't found a way to compute MISE with LM
# But as MSSE ~= MISE, this is fine.

LM.minimizeMSSE = function(poly.test = seq(5, 20, 0.01), M = 200, ni = 100, raw=F){
  best = Inf
  best_poly = 0
  
  pb = progress_bar$new(
    format = "[LM MSSE] Best: :bmsse (:bp) | Current: :cmsse (:cp) [:bar] :current/:total |:percent | :elapsed", 
    total = length(poly.test)
  )
  Sys.sleep(0.1)
  pb$tick(0)
  for (p in poly.test){
    SE = matrix(NA, M, ni)
    for (i in 1:M){
      boot.X = runif(ni)
      boot.Y = Y(boot.X)
      SE[i,] = (lm(boot.Y ~ poly(boot.X, degree = p, raw=raw))$fitted - m(boot.X))^2
    }
    
    MSE = colMeans(SE)
    MSSE = round(mean(MSE),6)
    
    if (MSSE < best){
      best = MSSE
      best_poly = p
    }
    pb$tick(tokens=list(bmsse=best, bp=best_poly, cmsse=MSSE, cp=p))
  }
  sprintf("[LM MSSE]: %.6f with %f as poly", best, best_poly)
  return (c(best, best_poly))
}
