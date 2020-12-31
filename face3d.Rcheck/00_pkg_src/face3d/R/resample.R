"resample.face3d" <- function(shape, n = 50, threshold = 1e-5){
   CHD <- cumchord(shape)
   xr  <- spline(CHD,shape[,1],method = "natural",n=n)$y
   yr  <- spline(CHD,shape[,2],method = "natural",n=n)$y
   zr  <- spline(CHD,shape[,3],method = "natural",n=n)$y
   shaper <- cbind(xr,yr,zr)
   nr.iter <- 1
   DIST <- 0.01
   while (DIST > threshold) {
      nr.iter <- nr.iter + 1
      CHD     <- cumchord(shaper)
      xr      <- spline(CHD,shaper[,1],method = "natural",n=n)$y
      yr      <- spline(CHD,shaper[,2],method = "natural",n=n)$y
      zr      <- spline(CHD,shaper[,3],method = "natural",n=n)$y
      shaper  <- cbind(xr,yr,zr)
      CHDnew  <- cumchord(shaper)
      DIST    <- sum(abs(diff(CHD)-diff(CHDnew)))
   }
   #nrs <- c(nr.iter,diff(CHDnew)[1])
   #names(nrs) <- c("nr of iter","chord dist")
   #print(nrs)
   return(shaper)
}

"cumchord" <- function(X) cumsum(sqrt(apply((X-rbind(X[1,],X[-(dim(X)[1]),]))^2,1,sum)))
