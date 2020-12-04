# Compute the MSE by M Monte-carlo simulations

NW.investigation = function(M=5000, n=100){
  
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
  
  SE = matrix(NA, M, length(CP.to_check))
  
  for (i in 1:M){
    boot.X = runif(n)
    # We need to sort in order to be able to retrieve the right x's
    boot.Y = Y(boot.X)[order(boot.X)]
    boot.X = sort(boot.X)
    boot.reg = lm(boot.Y ~ poly(boot.X, degree = LM.poly))
    boot.pred = fitted(boot.reg)[c(1,CP.to_check*n)]  # get em (add 1 bc index starts at 1 not 0)
    SE[i,] = (boot.pred - m(CP.to_check))^2
  }
  
  MSE = colMeans(SE)
  return(MSE)
}

# ---
# Create the dataframe for each sample size

CP.investigation.toframe = function(investigation_method){
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


# ---
# Drawing the evolution of the MSE

CP.investigation.drawMSE = function(investigation_frame, title){
  
  middle = nrow(investigation_frame)-1
  jpeg(paste0("./plots/investigation/",title,".jpg"), quality = 100, width = 1080, height = 720)
  plot(CP.investigation.samples, 
    head(investigation_frame, 1), col = 1, lwd = 2, lty=2,
    type="l", xlab="Sample size", ylab="MSE", ylim = c(-0.0001, min(max(investigation_frame), 0.4))
  )
  
  for (i in 2:middle){
    lines(CP.investigation.samples, investigation_frame[i,], col = i, lwd = 2)
  }
  
  lines(CP.investigation.samples, tail(investigation_frame, 1), col = 1, lwd = 2, lty=3)
  abline(h = 0, lty = 2)
  legend("topright", 
    legend = CP.to_check, 
    col = c(1, 2:middle, 1),
    lty = c(2, rep(1,middle-1), 3),
    lwd = 2,
    cex = 1.5 #0.8
  )
  title(title)
}
