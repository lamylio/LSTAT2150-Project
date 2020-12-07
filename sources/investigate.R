# Compute the MSE by M Monte-carlo simulations

NW.investigation = function(M=5000, n=100){
  set.seed(122020)
  SE = matrix(NA, M, length(CP.to_check))
  
  for (i in 1:M){
    boot.X = runif(n)
    boot.Y = Y(boot.X)
    SE[i,] = (sapply(CP.to_check,function(xi) NW.regEst(xi, boot.X, boot.Y, NW.h, Knorm))- m(CP.to_check))^2
  }
  
  MSE = colMeans(SE)
  return(MSE)
}

# ---

LM.investigation = function(M=5000, n=100){
  set.seed(122020)
  SE = matrix(NA, M, length(CP.to_check))
  
  for (i in 1:M){
    boot.X = sort(runif(n))
    boot.Y = Y(boot.X)
    boot.reg = lm(boot.Y ~ poly(boot.X, degree = LM.poly, raw=ifelse(n<50, T, F)))
    #boot.pred = predict(boot.reg, newdata=poly(CP.x, degree=LM.poly))[c(1,CP.to_check*100)]
    boot.pred = fitted(boot.reg)[c(1,CP.to_check*n)] 
    SE[i,] = (boot.pred - m(CP.to_check))^2
  }
  
  MSE = colMeans(SE)
  return(MSE)
}

# ---
# Create the dataframe for each sample size

CP.investigation.toframe = function(investigation_method){
  set.seed(122020)
  res = NULL
  pb = progress_bar$new(
    format = "[Investigation] Sample: :sample | [:bar] :current/:total |:percent | :elapsed", 
    total = length(CP.investigation.samples)
  )
  pb$tick(0)
  Sys.sleep(0.1)
  for (s in CP.investigation.samples){
    pb$tick(tokens=list(sample=s))
    res = cbind(res, investigation_method(n=s))
  }
  return (res)
}
