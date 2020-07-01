gpa.face3d <- function(x, scale = TRUE, model.mesh = FALSE, tol = 1e-5, monitor = FALSE) {
   
   n                     <- dim(x)[3]
   gpa                   <- procGPA(x, scale = scale, distances = FALSE, pcaoutput = FALSE)
   rnms                  <- rownames(x)
   rownames(gpa$rotated) <- rnms
   rownames(gpa$mshape)  <- rnms
   if (!model.mesh)
      return(invisible(list(rotated = gpa$rotated, mean = gpa$mshape)))
   
   x         <- gpa$rotated
   xp        <- x
   mshape    <- gpa$mshape
   wts       <- area.face3d(as.face3d(mshape, model.mesh = TRUE))$points
   size.x    <- mean(apply(x, 3, function(x) area.face3d(as.face3d(x, model.mesh = TRUE))$area))
   mn        <- apply(mshape,   2, function(x) sum(wts  * x) / sum(wts ))
   mshape    <- sweep(mshape,   2, mn,   "-")
   dst       <- mean(apply(xp, 3, function(x) mean(sweep((x - mshape)^2, 1, wts, "*"))))
   if (monitor) cat("iteration", 0, ":", dst, "\n")
   
   dst0 <- 1000
   dst  <- 100
   iter <- 0
   while (abs(dst - dst0) / dst0 > tol) {
      wts    <- area.face3d(as.face3d(mshape, model.mesh = TRUE))$points
      mn     <- apply(mshape,   2, function(x) sum(wts  * x) / sum(wts ))
      mshape <- sweep(mshape,   2, mn,   "-")
      for (i in 1:n)
         xp[ , , i] <- opa.face3d(x[ , , i], mshape, weights = wts, model.mesh = TRUE)
      size.xp <- mean(apply(xp, 3, function(x) area.face3d(as.face3d(x, model.mesh = TRUE))$area))
      xp <- sweep(xp, 3, size.x / size.xp, "*")
      
      mshape <- apply(xp, 1:2, mean)
      # clear3d()
      # for (i in 1:n) points3d(xp[ , , i], col = i)
      dst0   <- dst
      dst    <- mean(apply(xp, 3, function(x) mean(sweep((x - mshape)^2, 1, wts, "*"))))
      iter   <- iter + 1
      if (monitor) cat("iteration", iter, ":", dst, "\n")
   }

   return(invisible(list(rotated = xp, mean = mshape, weights = wts)))
   
}

# gpa.face3d <- function(x, scale = TRUE, model.mesh = FALSE, monitor = FALSE) {
#    
#    gpa                  <- procGPA(x, scale = scale, distances = FALSE, pcaoutput = FALSE)
#    rownames(gpa$mshape) <- rownames(x)
#    if (!model.mesh)
#       return(invisible(list(rotated = gpa$rotated, mean = gpa$mshape)))
#    
#    # clear3d()
#    # points3d(gpa$mshape)
#    
#    mshape0           <- gpa$mshape
#    rownames(mshape0) <- rownames(x)
#    dst               <- 1000
#    iter              <- 0
#    while (dst > 1e-3) {
#       wts               <- area.face3d(as.face3d(mshape0, model.mesh = TRUE))$points
#       mn                <- apply(x, 2:3, function(x) sum(wts  * x) / sum(wts ))
#       x1                <- sweep(x, 2:3, mn, "-")
#       gpa               <- procGPA(sweep(x1, 1, sqrt(wts), "*"), scale = scale, distances = FALSE, pcaoutput = FALSE)
#       mshape1           <- sweep(gpa$mshape, 1, sqrt(wts), "/")
#       rownames(mshape1) <- rownames(x)
#       mshape1           <- opa.face3d(mshape1, mshape0, scale = FALSE, model.mesh = TRUE)
#       # clear3d()
#       # points3d(gpa$mshape)
#       # points3d(mshape1, col = "red")
#       dst               <- mean(apply((mshape0 - mshape1)^2, 1, function(x) sqrt(sum(x))))
#       mshape0           <- mshape1
#       iter              <- iter + 1
#       if (monitor) cat("iteration", iter, ":", dst, "\n")
#    }
#    wts <- sqrt(area.face3d(as.face3d(mshape0, model.mesh = TRUE))$points)
#    
#    return(invisible(list(rotated = sweep(gpa$rotated, 1, wts, "/"), mean = mshape1, weights = wts)))
#    
# }
