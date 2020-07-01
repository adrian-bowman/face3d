pca.face3d <- function(x, group, weights, alpha = 0) {

   group.missing <- missing(group)
   if (!group.missing && length(unique(group)) != 2)
      stop("if group is specified it must define only two groups.")
   if (is.list(x) && all(c("aligned", "mean", "weights") %in% names(x))) {
      psc     <- x$aligned
      weights <- x$weights
   }
   else {
      psc <- x
      if (missing(weights)) weights <- rep(1, nrow(x))
   }
   if (length(weights) != nrow(psc))
      stop("'weights' does not have the correct length.")
   mns <- apply(psc, 2:3, function(x) sum(x * weights) / sum(weights))
   if (max(abs(mns)) > .Machine$double.eps * 1000) warning("the aligned shapes are not centred.")
   
   pms         <- apply(psc, 1:2, mean)
   k           <- nrow(psc)
   d           <- ncol(psc)
   n           <- dim(psc)[3]
   X           <- t(apply(psc, 3, c))
   X           <- sweep(X, 2, apply(X, 2, mean))
   X           <- sweep(X, 2, rep(sqrt(weights), d), "*")

   # psc.vec     <- matrix(0, n, d * k)
   # for (i in 1:n) psc.vec[i, ] <- c(psc[ , , i])
   # pms.vec     <- c(pms)
   # psc.vec.cen <- psc.vec - rep(1,n) %*% t(pms.vec)
   # psc.vec.cen <- sweep(psc.vec.cen, 2, rep(weights, d), "*")
   psc.vec.cen <- X
   
   if (group.missing) {
      if (alpha == 0) {
         Bpow     <- diag(rep(1, 3 * k))
         Bpow.inv <- diag(rep(1, 3 * k))
      }
      else {
         Be          <- ginvtps.face3d(pms)
         Be          <- Be[1:k,1:k]
         BEpow       <- bepower3d.face3d(Be, alpha)
         BEpow.inv   <- bepower3d.face3d(Be, -alpha)
         Bpow        <- kronecker(diag(1,3), BEpow)
         Bpow.inv    <- kronecker(diag(1,3), BEpow.inv)
      }
      psc.alpha   <- t(Bpow %*% t(psc.vec.cen))
      SVD         <- eigen(var(psc.alpha))
      evals       <- abs(SVD$values)
      evecs       <- SVD$vectors
      pct         <- round(100 * evals / sum(evals), 2)
      cumpct      <- cumsum(pct)
      scores      <- t(t(evecs) %*% Bpow %*% t(psc.vec.cen))
      evecs       <- Bpow.inv %*% evecs
      pca         <- list()
      pca$evecs   <- evecs # rescales GPCA eigenvectors 
      pca$scores  <- scores # GPC scores
      pca$sd      <- apply(scores, 2, sd)
      pca$mean    <- pms # Procrustes Mean Shape
      pca$percent <- pct # GPC eigenvalues for each GPC rescaled to percentages
   }
   else {
      n          <- table(group)
      code       <- levels(factor(as.character(group)))
      Sigma      <- ((n[1] - 1) * cov(X[group == code[1], ]) +
                     (n[2] - 1) * cov(X[group == code[2], ])) / (sum(n) - 2)
      eig        <- eigen(Sigma)
      pca        <- list(evecs = eig$vectors, scores = X %*% eig$vectors,
                         mean = pms, percent = 100 * eig$values / sum(eig$values))
      pca$sd     <- apply(pca$scores, 2, sd)
      pca$group  <- group
   }
   
   pca$weights <- weights
   invisible(pca)
}
