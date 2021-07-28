gpa.face3d <- function(x, triangles, match.ids, scale = TRUE, tol = 1e-5, monitor = 0) {
   
   if (missing(match.ids))    match.ids <- 1:dim(x)[1]
   if (is.logical(match.ids)) match.ids <- which(match.ids)
   if (is.character(match.ids)) {
      if (is.null(rownames(x))) stop("'match.ids' specifies rownames but 'x' has none.")
      match.ids <- match(match.ids, rownames(x))
   }

   gpa <- procGPA(x[match.ids, , ], scale = scale, distances = FALSE, pcaoutput = FALSE)
   n   <- dim(x)[3]
   y   <- x
   if (length(match.ids) < dim(x)[1]) {
      for (i in 1:n) y[ , , i] <- opa.face3d(x[match.ids, , i], gpa$rotated[ , , i], x[ , , i], scale = scale)
   }
   else {
      y <- gpa$rotated
   }
   dimnames(y) <- dimnames(x)

   if (missing(triangles))
      return(invisible(list(aligned = y, mean = apply(y, 1:2, mean))))
   
   x         <- y
   xp        <- x
   mshape    <- gpa$mshape
   rownames(mshape) <- rownames(y)
   wts       <- areas(as.face3d(list(vertices = mshape, triangles = triangles)))$points
   size.x    <- mean(apply(x, 3, function(x)
                       areas(as.face3d(list(vertices = x, triangles = triangles)))$area))
   mn        <- apply(mshape,   2, function(x) sum(wts  * x) / sum(wts ))
   mshape    <- sweep(mshape,   2, mn,   "-")
   dst       <- mean(apply(xp, 3, function(x) mean(sweep((x - mshape)^2, 1, wts, "*"))))
   if (monitor > 0) cat("iteration", 0, ":", dst, "\n")
   
   dst0 <- 1000
   dst  <- 100
   iter <- 0
   while (abs(dst - dst0) / dst0 > tol) {
      wts    <- areas(as.face3d(list(vertices = mshape, triangles = triangles)))$points
      mn     <- apply(mshape,   2, function(x) sum(wts  * x) / sum(wts ))
      mshape <- sweep(mshape,   2, mn,   "-")
      for (i in 1:n)
         xp[ , , i] <- opa.face3d(x[ , , i], mshape, weights = wts, triangles = triangles)
      size.xp <- mean(apply(xp, 3,
                            function(x) areas(as.face3d(list(vertices = x, triangles = triangles)))$area))
      xp <- sweep(xp, 3, size.x / size.xp, "*")
      
      mshape <- apply(xp, 1:2, mean)
      # clear3d()
      # for (i in 1:n) points3d(xp[ , , i], col = i)
      dst0   <- dst
      dst    <- mean(apply(xp, 3, function(x) mean(sweep((x - mshape)^2, 1, wts, "*"))))
      iter   <- iter + 1
      if (monitor > 0) cat("  iteration", iter, ":", dst, "\n")
   }
   if (monitor > 0) cat("  completed.", "\n")
   
   return(invisible(list(aligned = xp, mean = mshape, weights = wts)))
   
}
