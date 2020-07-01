subspaces.face3d <- function(x, weights) {
  
   if (is.list(x) && all(c("rotated", "weights") %in% names(x))) {
      psc     <- x$rotated
      weights <- x$weights
   }
   else {
      psc <- x
      if (missing(weights)) weights <- rep(1, k)
   }
   
   pms        <- apply(psc, 1:2, mean)
   k          <- dim(psc)[1]
   d          <- dim(psc)[2]
   n          <- dim(psc)[3]
   psc.aff    <- array(0, c(k, d, n))
   psc.nonaff <- array(0, c(k, d, n))
   rsq        <- matrix(0, n, d)
  
   for (i in 1:n) {
      for (j in 1:d) {
         mod.ij              <- lm(psc[ , j, i] ~ pms, weights = weights)
         psc.aff[ , j, i]    <- mod.ij$fitted.values
         psc.nonaff[ , j, i] <- pms[ , j] + mod.ij$residuals
         rsq[i, j]           <- summary(mod.ij)[["r.squared"]]
      }
   }
   rsq.stats <- c(mean(rsq), 1 - mean(rsq))
  
   return(invisible(list(affine = psc.aff, nonaffine = psc.nonaff, rsq = rsq)))
}

########################################################################################

"gpca.face3d" <- function(psc, alpha = 0) 
{
  pms         <- apply(psc, 1:2, mean)
  k           <- nrow(pms)
  d           <- ncol(psc)
  n           <- dim(psc)[3]
  psc.vec     <- matrix(0,n,d*k)
  for (i in 1:n) psc.vec[i,] <- c(psc[,,i])
  pms.vec     <- c(pms)
  psc.vec.cen <- psc.vec - rep(1,n) %*% t(pms.vec)
  Be          <- ginvtps.face3d(pms)
  Be          <- Be[1:k,1:k]
  BEpow       <- bepower3d.face3d(Be,alpha)
  BEpow.inv   <- bepower3d.face3d(Be,-alpha)
  {
    if (alpha == 0) {
      Bpow     <- diag(rep(1,3*k))
      Bpow.inv <- diag(rep(1,3*k))
    }
    else {
      Bpow     <- kronecker(diag(1,3),BEpow)
      Bpow.inv <- kronecker(diag(1,3),BEpow.inv)
    }
  }
  
  psc.alpha   <- t(Bpow %*% t(psc.vec.cen))
  SVD         <- eigen(var(psc.alpha))
  evals       <- abs(SVD$values)
  evecs       <- SVD$vectors
  pct         <- round(evals/sum(evals)*100,2)
  cumpct      <- cumsum(pct)
  scores      <- t(t(evecs) %*% Bpow %*% t(psc.vec.cen))
  evecs       <- Bpow.inv %*% evecs
  
  GPCA        <- list()
  GPCA$evecs  <- evecs # rescales GPCA eigenvectors 
  GPCA$scores <- scores # GPC scores
  GPCA$sd     <- apply(scores, 2, sd)
  GPCA$mean   <- pms # Procrustes Mean Shape
  GPCA$pct    <- pct # GPC eigenvalues for each GPC rescaled to percentages
  GPCA$cumpct <- cumpct # cumulative GPC eigenvalues for each GPC rescaled to percentages
  
  return(GPCA)
  
}

########################################################################################

"predict.gpca" <- function(GPCA, type = "sd", mag = 1, pc = 1)
{
  
  if (length(mag) != length(pc)) 
    stop("length of mag and pc must me the same.")
  
  # evecs  <- GPCA$evecs # GPC eigenvalues for each GPC rescaled to percentages
  scores <- GPCA$scores
  pms    <- GPCA$pms
  
  k      <- length(pc) 
  k1     <- nrow(GPCA$evecs)
  evecs  <- array(0,c(k1/3,3,k)) # reformatting loadings to an array
  lower  <- array(0,c(k1/3,3,k)) # reformatting mean minus magn times loadings to an array
  upper  <- array(0,c(k1/3,3,k)) # reformatting mean plus magn times loadings to an array
  
  for (i in 1:k) {
    evecs.i <- GPCA$evecs[,i]
    evecs.i <- matrix(evecs.i,k1/3,3)
    evecs[,,i] <- evecs.i/sqrt(sum(mod3d.face3d(evecs.i)))
    if (type == "sd") mag.i <- c(-mag[i]*sd(scores[,i]),mag[i]*sd(scores[,i]))
    if (type == "range") {
      range.i <- range(scores[,i])
      mag.i   <- c(mag[i]*range.i[1],mag[i]*range.i[2])
    }
    lower[,,i] <- pms + mag.i[1]*evecs[,,i]
    upper[,,i] <- pms + mag.i[2]*evecs[,,i]
  }
  
  GPCA        <- list()
  GPCA$evecs  <- evecs
  GPCA$lower  <- lower
  GPCA$upper  <- upper
  
  return(GPCA) 
  
}

################################################################################

"bepower3d.face3d" <- function(Be,alpha)
{
  k <- nrow(Be)
  if (alpha == 0) {
    BEpow <- diag(rep(1,k))
  }
  else {
    SVD   <- eigen(Be,symm = TRUE)
    evals <- SVD$values
    evecs <- SVD$vectors
    ID    <- which(evals > 0)
    k1    <- k - length(ID) 
    evals <- c(evals[ID]^(-alpha/2),rep(0,k1))
    BEpow <- evecs %*% diag(evals) %*% t(evecs)
  }
  return(BEpow)
  
}

################################################################################

"ginvtps.face3d" <- function(pms){
  
  require(MASS)
  
  k1    <- nrow(pms)
  Sdx   <- outer(pms[,1],pms[,1],"-")
  Sdy   <- outer(pms[,2],pms[,2],"-")
  Sdz   <- outer(pms[,3],pms[,3],"-")
  Sphi  <- sqrt(Sdx * Sdx + Sdy * Sdy + Sdz * Sdz)
  ONE   <- matrix(1, k1, 1)
  ZERO  <- matrix(0, 4, 4)
  Q     <- rbind(Sphi, t(ONE), t(pms))
  O     <- cbind(ONE, pms)
  O     <- as.matrix(O)
  U     <- rbind(O, ZERO)
  GAMMA <- cbind(Q, U)
  Gi    <- ginv(GAMMA)
  return(Gi)
  
}
################################################################################

"mod3d.face3d" <- function(evecs.mat){
  scale <- (evecs.mat[,1])^2 + (evecs.mat[,2])^2 + (evecs.mat[,3])^2
  return(scale)
}

################################################################################
################################################################################

