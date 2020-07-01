sm.psp <- function(x, y, lambda, df, method = "df",
                         weights, eval.points, ngrid,
                         nseg, bdeg = 3, pord = 2, display = "none",
                         period = rep(NA, ncol(x)),
                         increasing = FALSE, decreasing = FALSE, kappa = lambda * 100,
                         fixed, negative = TRUE) {

    if (is.vector(x)) x <- matrix(x, ncol = 1)
    ndim <- ncol(x)
    n    <- nrow(x)  
    if (missing(ngrid)) ngrid <- switch(ndim, 100, 20, 20)
    if (missing(df))    df    <- switch(ndim, 5, 12, 20)
    missing.eval.points       <- missing(eval.points)
    lambda.select             <- missing(lambda)
    if (missing(nseg)) nseg   <- switch(ndim, 100, 17, 10)

#    if (ndim == 1 & (!missing(fixed)) & (increasing | decreasing))
#       stop("monotonic increasing estimation is not available with fixed points.")
    if (ndim == 1 & increasing & decreasing)
       stop("only one of increasing and decreasing can be set to TRUE.")
    if (missing(weights)) weights <- rep(1, length(y))
    weights <- diag(weights)
    
    xrange <- matrix(nrow = 2, ncol = ndim)
    for (i in 1:ndim) {
      if (is.na(period[i]))
         # xrange[ , i] <- c(min(x[ , i]) - 0.05 * diff(range(x[,i])),
         #                   min(x[ , i]) + 0.05 * diff(range(x[,i])))
         xrange[ , i] <- range(x[ , i])
      else
         xrange[ , i] <- c(0, period[i])
    }

    # Compute B-spline basis
    
    # tim <- proc.time()
    b <- list(length = ndim)
    m <- vector(length = ndim)
    for (i in 1:ndim) {
       b[[i]] <- bbase(x[,i], xl = xrange[1, i], xr = xrange[2, i], 
                       nseg = nseg, deg = bdeg)
       m[i]   <- dim(b[[i]])[2]
    }
    B <- b[[1]]
    if (ndim > 1)
       B <- t(apply(cbind(b[[1]], b[[2]]), 1, 
                            function(x) c(x[1:m[1]] %x% x[-(1:m[1])])))
    if (ndim == 3)
       B <- t(apply(cbind(B,  b[[3]]), 1, 
                function(x) c(x[1:(m[1]*m[2])] %x% x[-(1:(m[1]*m[2]))])))
    # cat("bases:", proc.time()[3] - tim[3], "\n")
    # tim <- proc.time()

    # Construct penalty matrices
    P <- list(length = ndim)
    for (i in 1:ndim) {
       P[[i]] <- diff(diag(m[1]), diff = pord)
       if (!is.na(period)) {
          z      <- c(1, rep(0, m[i] - 4), -1)
          P[[i]] <- rbind(P[[i]], c(z, 0, 0))
          P[[i]] <- rbind(P[[i]], c(0, z, 0))
          P[[i]] <- rbind(P[[i]], c(0, 0, z))
       }
       P[[i]] <- crossprod(P[[i]])
    }
    if (ndim >= 2) {
       P[[1]] <- P[[1]] %x% diag(m[2])
       P[[2]] <- diag(m[2]) %x% P[[2]]
    }
    if (ndim == 3) {
       P[[1]] <- P[[1]] %x% diag(m[3])
       P[[2]] <- P[[2]] %x% diag(m[3])
       P[[3]] <- diag(m[1]) %x% diag(m[2]) %x% P[[3]]
    }
    # cat("penalties:", proc.time()[3] - tim[3], "\n")
    # tim <- proc.time()
    
    btb  <- t(B) %*% weights %*% B
    bty  <- t(B) %*% weights %*% y
    
    # Identify lambda (common value across dimensions)
    
    if (lambda.select) {
       if (method == "df") {
          lambda.df <- function(lambda, btb, bty, P) {
             mat  <- 0
             for (i in 1:length(P))
                mat <- mat + lambda * P[[i]]
             B1   <- solve(btb + mat) 
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
       lambda <- rep(lambda, ndim)
    }
    
    # Fit
    
    btb   <- t(B) %*% weights %*% B
    bty   <- t(B) %*% weights %*% y
    # cat("matrices:", proc.time()[3] - tim[3], "\n")
    # tim <- proc.time()
    mat <- 0
   
    for (i in 1:length(P))
       mat <- mat + lambda[i] * P[[i]]
    B1  <- solve(btb + mat) 
    # cat("inversion:", proc.time()[3] - tim[3], "\n")
    # tim <- proc.time()
    E     <- B1 %*% t(B)
    beta  <- as.vector(E %*% y)
    # cat("solution:", proc.time()[3] - tim[3], "\n")
    # tim <- proc.time()
    # f     <- lsfit(rbind(B, P1, P2, P3), c(y, nix), intercept = FALSE)
    # h     <- hat(f$qr)[1:m]
    # beta  <- f$coef

    mu  <- c(B %*% beta)
    edf <- sum(diag(btb %*% B1))
    # cat("fitting:", proc.time()[3] - tim[3], "\n")
    # tim <- proc.time()
    
    # Adjust for monotonicity and for fixed points - one covariate only
    
    if (ndim == 1 & !is.na(increasing) & (increasing | decreasing | !negative)) {
       # if (!missing(fixed))
       #     stop("fixed values cannot be used with increasing/decreasing/non-negative constraints")
       D1    <- diff(diag(m[1]), diff = 1)
       delta <- 1
       while (delta > 1e-5) {
          mat1 <- mat
       	  if (increasing | decreasing) {
             if (increasing) v <- as.numeric(diff(beta) <= 0)
             if (decreasing) v <- as.numeric(diff(beta) >= 0)
             mat1 <- mat1 + kappa * t(D1) %*% diag(v) %*% D1
       	  }
          if (!negative) {
             # A    <- bbase(seq(0, 1, length = 20), xl = 0, xr = 1, nseg = nseg, deg = bdeg)
             # v    <- as.numeric(c(A %*% beta) < 0)
             # mat1 <- mat1 + kappa * (t(A) %*% diag(v) %*% A)
             v    <- as.numeric(beta < 0)
             mat1 <- mat1 + kappa * diag(v)
          }
          B1       <- solve(t(B) %*% weights %*% B + mat1)
          beta.old <- beta
          beta     <- as.vector(B1 %*% t(B) %*% weights %*% y)
          delta    <- sum((beta - beta.old)^2) / sum(beta.old^2)
       }
    }
   
    if (ndim == 1 & !missing(fixed) && all(!is.na(fixed))) {
       if (any(fixed[ , 1] <= xrange[1, i]) | 
           any(fixed[ , 1] >= xrange[2, i]))
          stop("fixed points must be inside the range of the data.")
       fixed <- matrix(c(fixed), ncol = ndim + 1)
       A     <- bbase(fixed[ , 1:ndim], xl = xrange[1, i], xr = xrange[2, i], 
                              nseg = nseg, deg = bdeg)
       beta  <- beta + 
                 B1 %*% t(A) %*% solve(A %*% B1 %*% t(A)) %*% (fixed[ , ndim + 1] - A %*% beta)
       edf   <- NA
    }
    
    # Cross-validation and dispersion
    # r     <- (y - mu ) / (1 - h)
    # cv    <- sqrt(sum(r ^2))
    # sigma <- sqrt(sum((y - mu) ^2) / (n - sum(h)))

    # Evaluate the estimate at eval.points
    # evp is the version of the avaluation points returned by the function
    if (missing.eval.points) {
       eval.points <- list(length = ndim)
       for (i in 1:ndim)
          eval.points[[i]] <- seq(min(x[,i]), max(x[,i]), length = ngrid)
       evp <- if (ndim == 1) eval.points[[1]] else eval.points
       eval.points <- as.matrix(expand.grid(eval.points))
    }
    else if (ndim == 1) {
       evp         <- eval.points
       eval.points <- matrix(eval.points, ncol = 1)
    }
    else {
    	 if (is.list(eval.points)) {
          evp         <- eval.points
          eval.points <- as.matrix(expand.grid(eval.points))
    	 }
       else
          evp <- eval.points
    }
    
    for (i in 1:ndim) {
       b[[i]] <- bbase(eval.points[,i], xl = xrange[1, i], xr = xrange[2, i], 
                                        nseg = nseg, deg = bdeg)
       m[i]   <- dim(b[[i]])[2]
    }
    B <- b[[1]]
    if (ndim > 1)
       B <- t(apply(cbind(b[[1]], b[[2]]), 1, 
                            function(x) c(x[1:m[1]] %x% x[-(1:m[1])])))
    if (ndim == 3)
       B <- t(apply(cbind(B,  b[[3]]), 1, 
                function(x) c(x[1:(m[1]*m[2])] %x% x[-(1:(m[1]*m[2]))])))
    est <- c(B %*% beta)
    S   <- B %*% E
    
    
    if (missing.eval.points & (ndim > 1)) est <- array(est, dim = rep(ngrid, ndim))
    if (!missing.eval.points) display <- "none"
    
    # cat("estimate:", proc.time()[3] - tim[3], "\n")
    # tim <- proc.time()

    # Plot data and fit
    
    if (display != "none") {
       if (ndim == 1) {
          plot(x, y, main = '', xlab = '', ylab = '')
          lines(eval.points[,1], est, col = 'blue')
          # if (se > 0 ) {
          #    Covb <- solve(btb + P[[1]])
          #    Covz <- sigma^2 * B %*% Covb %*% t(B)
          #    seb  <- se * sqrt(diag(Covz))
          #    lines(u, est + seb, lty = 2, col = 'red')
          #    lines(u, est - seb, lty = 2, col = 'red')
          # }
       }
       else if (ndim == 2) {
       	  if (display == "image")
       	     image(evp[[1]], evp[[2]], est)
       	  else
             persp(evp[[1]], evp[[2]], est, 
                   ticktype = "detailed", col = "green", d = 10, theta = 30)
       }
    }

    # Return list
    pp <- list(x = x, y = y,
               eval.points = evp, estimate = est, muhat = mu,
               df = edf, lambda = lambda,
               beta = beta, B = B, S = S,
               nseg = nseg, bdeg = bdeg, pord = pord)
    return(invisible(pp))
}

ps.matrices <- function(x, xrange, ndims, nseg, bdeg = 3, pord = 2, period = NA,
												decompose =  TRUE) {

	# Compute a set of basis functions and a penalty matrix associated with x.
	# An intercept term and the main effect of any interaction terms are removed.

	ndimx <- ncol(x)
	if (ndimx > 3) stop("terms with more than three dimensions cannot be used.")
	n    <- nrow(x)

	if (missing(nseg)) nseg <- rep(switch(ndimx, 100, 17, 7), ndimx)

	# Compute B-spline basis

	b <- list(length = ndimx)
	m <- vector(length = ndimx)
	for (i in 1:ndimx) {
		b[[i]] <- bbase(x[,i], xl = xrange[i , 1], xr = xrange[i, 2], nseg = nseg[i],
										deg = bdeg)
		m[i]   <- ncol(b[[i]])
	}

	B <- b[[1]]
	if (ndimx > 1)
		B <- t(apply(cbind(b[[1]], b[[2]]), 1,
								 function(x) c(x[1:m[1]] %x% x[-(1:m[1])])))
	if (ndimx == 3)
		B <- t(apply(cbind(B,  b[[3]]), 1,
								 function(x) c(x[1:(m[1]*m[2])] %x% x[-(1:(m[1]*m[2]))])))

	# Construct smoothness penalty matrices
	P <- list()
	for (i in 1:ndimx) {
		P[[i]] <- diff(diag(m[i]), diff = pord)
		if (!is.na(period[i])) {
			z      <- c(1, rep(0, m[i] - 4), -1)
			P[[i]] <- rbind(P[[i]], c(z, 0, 0))
			P[[i]] <- rbind(P[[i]], c(0, z, 0))
			P[[i]] <- rbind(P[[i]], c(0, 0, z))
		}
		P[[i]] <- crossprod(P[[i]])
	}
	if (ndimx >= 2) {
		P[[1]] <- P[[1]] %x% diag(m[2])
		P[[2]] <- diag(m[2]) %x% P[[2]]
	}
	if (ndimx == 3) {
		P[[1]] <- P[[1]] %x% diag(m[3])
		P[[2]] <- P[[2]] %x% diag(m[3])
		P[[3]] <- diag(m[1]) %x% diag(m[2]) %x% P[[3]]
	}
	pmat <- matrix(0, nrow = ncol(B), ncol = ncol(B))
	for (i in 1:ndimx)
		pmat <- pmat + P[[i]]

	#     Construct anova constraint penalty matrices
	if (length(ndims) == 1) {
		# Sum of coefficients constraint
		# cmat <- matrix(1, nrow = prod(m), ncol = prod(m))
		# Sum of estimated values constraint
		Bsum <- apply(B, 2, sum)
		cmat <- Bsum %o% Bsum
		# Corner point constraint (first coefficient is 0
		# cmat <- diag(c(1, rep(0, ncol(B) - 1)))
		pmat <- pmat + cmat
	}
	else if (length(ndims) == 2) {
		if (all(ndims == c(1, 1))) ind <- c(m[1], m[2])
		if (all(ndims == c(1, 2))) ind <- c(m[1], m[2] * m[3])
		if (all(ndims == c(2, 1))) ind <- c(m[1] * m[2], m[3])
		pmat <- pmat + matrix(1, nrow = ind[1], ncol = ind[1]) %x% diag(ind[2])
		pmat <- pmat + diag(ind[1]) %x% matrix(1, nrow = ind[2], ncol = ind[2])
	}
	else if (length(ndims) == 3) {
		pmat <- pmat + matrix(1, nrow = m[1], ncol = m[1]) %x% diag(m[2]) %x% diag(m[3])
		pmat <- pmat + diag(m[1]) %x% matrix(1, nrow = m[2], ncol = m[2]) %x% diag(m[3])
		pmat <- pmat + diag(m[1]) %x% diag(m[2]) %x% matrix(1, nrow = m[3], ncol = m[3])
	}

	result <- list(B = B, P = pmat, xrange = xrange, nseg = nseg, bdeg = bdeg, pord = pord)
	if (length(ndims) == 1) result$cmat <- cmat
	invisible(result)
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
