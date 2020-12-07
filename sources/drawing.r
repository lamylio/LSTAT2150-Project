# ---
# Draw the comparison of LM and NW

CP.drawComparison = function(BestNW, BestLM, n=100, save=F){
  set.seed(122020)
  if(save) {jpeg("./plots/comparison.jpg", quality = 100, width = 1080, height = 720)}
  plot(CP.X, CP.Y, pch="+", col = "grey", xlab="X", ylab="Y")
  lines(CP.x, BestNW, col = 2, lwd=ifelse(save, 2, 1))
  lines(CP.x, BestLM, col=4, lwd=ifelse(save, 2, 1), pch=19)
  lines(CP.x, m(CP.x), col = 1, lwd=ifelse(save, 2, 1))
  legend("topleft", 
         legend = c("True", "NW", "LM"),
         col = c(1, 2, 4),
         lty = c(1, 1, 1),
         #pch = c(NA, 19, 19),
         cex = ifelse(save, 1.5, 0.8)
  )
  title("NW v.s LM")
  if(save) {dev.off()}
}




# ---
# Draw the evolution of the MSE

CP.drawMSE = function(investigation_frame, title, save=F){
  set.seed(122020)
  middle = nrow(investigation_frame)-1
  if(save) {jpeg(paste0("./plots/investigation/",title,".jpg"), quality = 100, width = 1080, height = 720)}
  plot(CP.investigation.samples, 
       head(investigation_frame, 1), col = 1, lwd = ifelse(save, 2, 1), lty=2,
       type="l", xlab="Sample size", ylab="MSE", ylim = c(-0.0001, min(max(investigation_frame), 0.4))
  )
  
  for (i in 2:middle){
    lines(CP.investigation.samples, investigation_frame[i,], col = i, lwd = ifelse(save, 2, 1))
  }
  
  lines(CP.investigation.samples, tail(investigation_frame, 1), col = 1, lwd = ifelse(save, 2, 1), lty=3)
  abline(h = 0, lty = 2)
  legend("topright", 
         legend = CP.to_check, 
         col = c(1, 2:middle, 1),
         lty = c(2, rep(1,middle-1), 3),
         lwd = 2,
         cex = ifelse(save, 1.5, 0.8)
  )
  title(title)
  if(save) {dev.off()}
}