sm.map        <- function(x, y, weights, fixed, nseg = 25, deg = 3, pdiff = 2,
                          ig.a = 1e-3, ig.b = 1e-3, lambdas = exp(seq(-10, 10, length = 50)),
                          prior = function(lambda) 1,
                          ngrid = 401, eval.points = seq(min(x), max(x), length = ngrid)) {

  if (!requireNamespace("MASS", quietly = TRUE)) stop("the MASS package is required.")

  B           <- bbase(x, xl = min(eval.points), xr = max(eval.points), nseg = nseg, deg = deg)
  D           <- diff(diag(ncol(B)), diff = pdiff)
  DtD         <- crossprod(D)              #compute penalty matrix
  weights[1]     <- 10
  weights[length(weights)] <- 10
  weights     <- MASS::ginv(diag(weights))
  BtB         <- t(B) %*% weights %*% B    #compute b-spline basis
  bty         <- t(B) %*% weights %*% y
  P.eigen     <- eigen(BtB + DtD)          # Eigen composition s.t. P.eigen = RoLoRot
  if (any(P.eigen$values < sqrt(.Machine$double.eps) * max(P.eigen$values)))
    stop("Singularity detected. No well-defined estimate.")
  Mt          <- t(P.eigen$vectors) * (1 / sqrt(P.eigen$values))   #Mt = Lo^(-.5)Rot
  
  #BM = BRoLo^(-.5) = ULVt   so   B= ULVtLo^(-.5)R0t  so Xtinv = RoLo^(-.5)V
  Q.svd       <- svd(B %*% t(Mt), nu = ncol(B), nv = ncol(B))   
  d           <- c(pmin(Q.svd$d,1), rep(0, ncol(B) - length(Q.svd$d)))^2
  e           <- 1 - d
  Xtinv       <- t(t(P.eigen$vectors) * sqrt(1 / P.eigen$values)) %*% Q.svd$v
  sel         <- 1:length(Q.svd$d)
  log.det.XtX <- sum(log(P.eigen$values))
  rank.D      <- sum(e > sqrt(.Machine$double.eps))
  z           <- drop(t(Q.svd$u) %*% y)

 #see notes but need to solve for z above since can then solve for alpha

# Function to compute the log-posterior
  loglik <- function(lambda) {
  	# Get the coefficient
    coef         <- drop(Xtinv[,sel]%*%(z*sqrt(d[sel]) / (d[sel]+lambda*e[sel])))
    residuals    <- y-B%*%coef
    # Get the posterior determinant
    log.post.det <- -sum(log(d+lambda*e))-log.det.XtX
    # Compute the log-posterior
    0.5 * rank.D * log(lambda) + 0.5 * log.post.det - 
    (ig.a + length(y)/2) * log(2*ig.b + sum(y*residuals)) + log(prior(lambda))
  } 
  
  
# Find the best lambda and compute alpha coeff for best lambda
  logliks  <- sapply(lambdas, loglik)
  lambda   <- lambdas[which.max(logliks)]
  alpha    <- drop(Xtinv[ , sel] %*% (z * sqrt(d[sel]) / (d[sel] + lambda * e[sel])))
  fitted   <- drop( B %*% alpha)
    
    
# Adjusting for fixed points add into DTD penalty matrix??
   
   # beta  <- fitted
   # B1    <- BtB
   # fixed <- matrix(c(fixed), ncol = 2)
   # A     <- bbase(fixed[ , 1:1], xl = min(eval.points), xr = max(eval.points), nseg = nseg, deg = deg)
   # beta  <- beta + B1 %*% t(A) %*% solve(A %*% B1 %*% t(A)) %*% (fixed[ , 2] - A %*% beta)
   # est <- c(B %*% beta)

   # i     <-1
   # x     <- matrix(x, ncol = 1)
   # ndim <- 1
   # mat   <- 0 + lambda * DtD[[1]]
   # B1    <- solve(BtB + mat) 
   # beta  <- as.vector(B1 %*% bty)
    # if (any(fixed[,1] <= min(x) - 0.05 * diff(range(x[,i]))) | 
           # any(fixed[,1] >= max(x) + 0.05 * diff(range(x[,i]))))
          # stop("fixed points must be inside the range of the data.")
       # fixed <- matrix(c(fixed), ncol = ndim + 1)
       # A     <- bbase(fixed[ , 1:ndim], xl = min(x[,i]) - 0.05 * diff(range(x[,i])),
                              # xr = max(x[,i]) + 0.05 * diff(range(x[,i])), 
                              # nseg = nseg, deg = bdeg)
       # beta  <- beta + 
                 # B1 %*% t(A) %*% solve(A %*% B1 %*% t(A)) %*% (fixed[ , ndim + 1] - A %*% beta)
       # edf   <- NA
             
 
 
return(list(best.lambda = lambda, trial.lambdas = lambdas, logliks = logliks,
            alpha = alpha, fitted = fitted, eval.points = eval.points))
}





bbase<-function(x, xl = min(x), xr = max(x), nseg = 10, deg = 3) {
    dx <- (xr - xl) / nseg
    knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
    P <- outer(x, knots, tpower, deg)
    n <- dim(P)[2]
    D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
    B <- (-1) ^ (deg + 1) * P %*% t(D)
    B
}

tpower <- function(x, t, p){
    (x - t) ^ p * (x > t)
}
