


smquad.face3d <- function(x, pdist, y, lambda, df = 5, method = "df",
                        weights, eval.points, ngrid,
                        nseg, bdeg = 3, pord = 2, display = "lines", fixed) {
    n <- length(x)
    if (missing(ngrid))  ngrid  <- switch(1, 100, 20, 20)
    lambda.select <- missing(lambda)
    
    if (missing(nseg))   nseg   <- switch(1, 100, 17, 10)

    if (missing(weights)) weights <- rep(1, length(y))
    weights <- diag(weights)
    
    if (missing(eval.points)) eval.points <- seq(min(x), max(x), length = ngrid)
   
    # Compute B-spline basis
    B0 <- bbase(x, xl = min(x), xr = max(x), nseg = nseg, deg = bdeg)
    m  <- dim(B0)[2]
    B  <- cbind(B0, sweep(B0, 1, pdist, "*"), sweep(B0, 1, pdist^2, "*"))

    # Construct penalty matrices
    # Matrix D --expresses roughness- based on the number of splines- called P here, second order differences (pord=2)      
    P <- diff(diag(m), diff = pord)
    P <- t(P) %*% P 

    btb <- t(B) %*% weights %*% B
   
    # bty <- t(B) %*% weights %*% log(y)
    # taking the log(y) created NA's
    bty <- t(B) %*% weights %*% y
    df <- 2 * df
    
    # Identify lambda 
    if (lambda.select) {
       if (method == "df") {
          lambda.df <- function(lambda, btb, bty, P) {
             #mat  <- 0
             # for (i in 1:length(P))
             #    mat <- mat + lambda[i] * P[[i]]
             #mat  <- lambda * P
             
             mat <- matrix(0, nrow = 3 * m, ncol = 3 * m)
             mat[1:m, 1:m]         <- lambda * P
             mat[m + 1:m, m + 1:m] <- lambda * P
             mat[2 * m + 1:m, 2 * m + 1:m] <- lambda * P
             B1   <- solve(btb + mat) 
             # B2   <- solve(btp + mat )
             beta <- as.vector(B1 %*% bty)
             sum(diag(btb %*% B1))
          }
          lambda <- 1
          # print(lambda.df(lambda, btb, bty, P))
          # print(adf)
          while (lambda.df(lambda, btb, bty, P) <= df) lambda <- lambda / 10
          lower <- lambda
          lambda <- 1
          while (lambda.df(lambda, btb, bty, P) >= df) lambda <- lambda * 10
          upper <- lambda
          lambda.crit <- function(lambda, btb, bty, P, df)
             lambda.df(lambda, btb, bty, P) - df
          result <- uniroot(lambda.crit, interval = c(lower, upper), btb, bty, P, df)
          lambda <- result$root
       }
    }
    
    # Fit
    mat <- matrix(0, nrow = 3 * m, ncol = 3 * m)
    mat[1:m, 1:m]         <- lambda * P
    mat[m + 1:m, m + 1:m] <- lambda * P
    mat[2 * m + 1:m, 2 * m + 1:m] <- lambda * P
    B1   <- solve(btb + mat) 
    beta <- as.vector(B1 %*% bty)
    edf  <- sum(diag(btb %*% B1))        

     
   # Adjust for fixed points
    if (!missing(fixed) && all(!is.na(fixed))) {
       if (any(fixed[,1] < min(eval.points)) | 
           any(fixed[,1] > max(eval.points)))
          stop("fixed points must be inside the range of eval.points.")
       fixed <- matrix(c(fixed), ncol = 2)
       A.beg     <- bbase(fixed[ , 1], xl = min(x), xr = max(x), nseg = nseg, deg = bdeg)
       # mat <- matrix(0, nrow = 2 , ncol = 2 * m + 1)
       mat   <- matrix(0, nrow = 2 , ncol = 3 * m)
       mat[1:2, m + 1:m]       <- A.beg
       #mat[1:2, 2 * m + 1:m]       <- A.beg
       #mat[m + 1:2, m + 1:m] <- sweep(A.beg , 1 , pdist , "*")
       #mat[2m:1, 2m:1]        <- pdist^2
       A     <- mat
       beta  <- beta + 
                 B1 %*% t(A) %*% solve(A %*% B1 %*% t(A)) %*% (fixed[ , 2] - A %*% beta)
                 
       edf   <- NA
    }  

    b   <- bbase(eval.points, xl = min(x), xr = max(x), nseg = nseg, deg = bdeg)
    #fixed points fixed
    #set smoothpath npts length to 103
    # c2 is 103 in length and being divided from 40 in length because in the original cvture is a single number
    #cvture  <- 2 * beta[2 * m + 1]
    c2      <- 2 * beta[2 * m + 1:m]
    #cat("curvature:", cvture, "\n")
    mu  <- c(B %*% beta)
    #est <- -c(b %*% beta[m + 1:m]) / cvture
    est <- -c(b %*% beta[m + 1:m]) / c2
  
    # Plot data and fit
    if (display != "none") {
       plot(x, pdist, col = topo.colors(20)[cut(y, 20, labels = FALSE)])
       lines(eval.points, est, col = "red")
       }

    # Return list
    pp <- list(x = x, pdist = pdist, y = y, nseg = nseg,
               bdeg = bdeg, pord = pord, beta = beta, B = B, b = b,
               lambda = lambda, eval.points = eval.points, estimate = est,
               fitted = mu, btb = btb, P = P, bty = bty, df = edf, curvature = c(B0 %*% beta[2 * m + 1:m]))
    return(invisible(pp))
}

bbase <- function(x, xl = min(x), xr = max(x), nseg = 10, deg = 3) {
# Construct B-spline basis
    dx <- (xr - xl) / nseg
    knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
    P <- outer(x, knots, tpower, deg)
    n <- dim(P)[2]
    D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
    B <- (-1) ^ (deg + 1) * P %*% t(D)
    B
}

tpower <- function(x, t, p)
# Truncated p-th power function
    (x - t) ^ p * (x > t)
