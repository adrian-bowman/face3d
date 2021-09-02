"resample.face3d" <- function(path, n = 50, method = "spline", threshold = 1e-5){

   chd <- arclengths(path)
   
   if (method == "spline") {
      xr    <- spline(chd, path[,1], method = "natural", n = n)$y
      yr    <- spline(chd, path[,2], method = "natural", n = n)$y
      zr    <- spline(chd, path[,3], method = "natural", n = n)$y
      pathr <- cbind(xr, yr, zr)
      niter <- 1
      DIST <- 0.01
      while (DIST > threshold) {
         niter   <- niter + 1
         chd     <- cumchord(pathr)
         xr      <- spline(chd, pathr[,1], method = "natural", n = n)$y
         yr      <- spline(chd, pathr[,2], method = "natural", n = n)$y
         zr      <- spline(chd, pathr[,3], method = "natural", n = n)$y
         pathr   <- cbind(xr, yr, zr)
         chdnew  <- cumchord(pathr)
         DIST    <- sum(abs(diff(chd) - diff(chdnew)))
      }
      #nrs <- c(nr.iter,diff(chdnew)[1])
      #names(nrs) <- c("nr of iter","chord dist")
      #print(nrs)
   }
   
   else if (method == "linear") {
      xr    <- approx(chd, path[,1], n = n)$y
      yr    <- spline(chd, path[,2], n = n)$y
      zr    <- spline(chd, path[,3], n = n)$y
      pathr <- cbind(xr, yr, zr)
   }
   
   return(pathr)
}

# "cumchord" <- function(X) cumsum(sqrt(apply((X-rbind(X[1,],X[-(dim(X)[1]),]))^2,1,sum)))
