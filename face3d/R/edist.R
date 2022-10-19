"edist" <- function(x1, x2, type = "pairwise", index = FALSE) {
  
   if (missing(x2))   x2 <- rep(0, 3)
   if (is.vector(x1)) x1 <- t(as.matrix(x1))
   if (is.vector(x2)) x2 <- t(as.matrix(x2))
   
   if (is.matrix(x1) | is.matrix(x2)) {
      if (is.face3d(x1)) x1 <- x1$vertices
      if (is.face3d(x2)) x2 <- x2$vertices
      if (ncol(x1) != 3) stop("the dimension of x1 do not correspond to 3d co-ordinates.")
      if (ncol(x2) != 3) stop("the dimension of x2 do not correspond to 3d co-ordinates.")
      if (type == "pairwise") {
         dst <- fields::rdist(x1, x2)
         if (any(dim(dst) == 1)) dst <- c(dst)
      }
      else if (type == "min") {
         dst  <- fields::rdist(x1, x2)
         if (index) ind  <- apply(dst, 1, which.min)
         dst <- apply(dst, 1, min)
         if (index) dst <- list(distance = dst, index = ind)
      }
      else if (type == "minsum")
         dst <- sum(c(apply(x1, 1, function(x) min(fields::rdist(t(x), x2)))))
      else
         stop("invalid value of the type argument.")
      invisible(dst)
   }
}
