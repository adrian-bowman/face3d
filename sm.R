sm <- function(x, y, data, subset, weights, bdeg = 3, pord = 2, h, model,
                    random, ...) {

   weights.missing <- missing(weights)
   random.missing  <- missing(random)

   if (!missing(y)) {
      x.name <- deparse(substitute(x))
      y.name <- deparse(substitute(y))
      if (weights.missing) weights <- NA
      if (missing(model)) model <- "none"
      return(sm.regression(x, y, h, model = model, weights = weights, ...))
   }
   else if (class(x) != "formula") {
      x.name <- deparse(substitute(x))
      if (weights.missing) weights <- NA
      if (missing(model)) model <- "none"
      return(sm.density(x, h, model = model, weights = weights, xlab = x.name, ...))
   }
      
   opt <- sm.options(list(...))
   replace.na(opt, display,   "none")
   replace.na(opt, reference, "none")
   replace.na(opt, panel,     FALSE)
   pam.formula    <- x
   data.missing   <- missing(data)
   subset.missing <- missing(subset)
   
   if (!data.missing) {
      mf  <- match.call(expand.dots = FALSE)
      nms <- c("x", "data", "subset", "weights")
      if (!random.missing) nms <- c(nms, "random")
      m   <- match(nms, names(mf), 0L)
      mf  <- mf[c(1L, m)]
      mf$drop.unused.levels <- TRUE
      mf[[1L]] <- quote(stats::model.frame)
      data <- eval(mf, parent.frame())
   }
   if (!random.missing) random <- as.factor(as.character(data[ , "(random)"]))

   terms.obj        <- terms(pam.formula, specials = "s")
   vars.inf         <- if (data.missing) eval.parent(attr(terms.obj, "variables"))
                       else eval(attr(terms.obj, "variables"), data)
   term.labels      <- attr(terms.obj, "term.labels")
   s.ind            <- attr(terms.obj, "specials")$s
   response.ind     <- attr(terms.obj, "response")
   involved         <- attr(terms.obj, "factors")
   terms.linear     <- matrix(c(involved[s.ind, ]), ncol = length(term.labels))
   terms.linear     <- which(apply(terms.linear, 2, sum) == 0)
   nterms           <- length(term.labels)
   terms.smooth     <- which(!(1:nterms %in% terms.linear))
   rhs.linear       <- paste(term.labels[terms.linear], collapse = " + ")
   rhs.linear       <- if (nchar(rhs.linear) == 0) "1" else rhs.linear
   formula.linear   <- as.formula(paste("~", rhs.linear))
   names(vars.inf)  <- rownames(involved)
   bricks.type      <- sapply(vars.inf[-response.ind], mode)
   ind              <- (bricks.type == "numeric") & 
                       sapply(vars.inf[-response.ind], is.factor)
   bricks.type[ind] <- "factor"
   Xlinear          <- vars.inf[-response.ind][bricks.type != "list"]
   names(Xlinear)   <- names(bricks.type)[bricks.type != "list"]
   
   ylab             <- attr(terms.obj, "variables")
   ylab             <- strsplit(deparse(ylab), ",")[[1]][1]
   ylab             <- substr(ylab, 6, nchar(ylab))
   y                <- unlist(vars.inf[[response.ind]])
   X                <- list()
   xlabels          <- list()
   xlab             <- list()
   xdims            <- list()
   ndims            <- list()
   df               <- list()
   nseg             <- list()
   lambda           <- list()
   period           <- list()
   increasing       <- list()
   xrange           <- list()
   fixed            <- list()
   fac              <- list()
   xmissing         <- FALSE
   
   if (length(terms.smooth) < 1) stop("there must be at least one smooth term.")
   
   if (any(apply(involved, 2, sum) > 3))
      stop("four-way interactions not yet implemented.")

   for (i in 1:length(terms.smooth)) {
   	
      inv     <- which(involved[ , terms.smooth[i]] == 1)
      ilinear <- which(bricks.type[names(inv)] == "numeric")
      ifactor <- which(bricks.type[names(inv)] == "factor")
      if (length(ilinear) > 0)
         stop("interactions with linear terms are not yet implemented.")
      if (length(ifactor) > 1)
         stop("interactions with more than one factor are not yet implemented.")
      else if (length(ifactor) == 1) {
         fact     <- names(bricks.type)[ifactor]
         inv      <- inv[-match(fact, names(inv))]
         fac[[i]] <- Xlinear[[fact]]
      }
      else
         fac[[i]] <- NA

      nvars           <- length(inv)
      X[[i]]          <- matrix(nrow = length(y), ncol = 0)
      xlabels[[i]]    <- vector("character")
      xlab[[i]]       <- vector("character")
      xdims[[i]]      <- numeric()
      ndims[[i]]      <- numeric()
      df[[i]]         <- numeric()
      lambda[[i]]     <- numeric()
      period[[i]]     <- numeric()
      increasing[[i]] <- logical()
      nseg[[i]]       <- numeric()
      xrange[[i]]     <- matrix( , nrow = 0, ncol = 2)
      fixed[[i]]      <- matrix( , nrow = 0, ncol = 2)
      for (j in inv) {
         lambda[[i]]  <- c(lambda[[i]], vars.inf[[j]]$lambda)
         nseg[[i]]    <- c(nseg[[i]],   vars.inf[[j]]$nseg)
         xlabels[[i]] <- c(xlabels[[i]], vars.inf[[j]]$variables)
      	 newvar       <- if (data.missing) eval.parent(parse(text = vars.inf[[j]]$variables[1]))
      	                 else              eval(parse(text = vars.inf[[j]]$variables[1]), data)
      	 xdims[[i]]   <- if (is.matrix(newvar)) c(xdims[[i]], ncol(newvar)) else c(xdims[[i]], 1)
         if (length(vars.inf[[j]]$variables) > 1) {
            for (k in 2:length(vars.inf[[j]]$variables)) {
            	   newvar <- if (data.missing) cbind(newvar, eval.parent(parse(text = vars.inf[[j]]$variables[k])))
                         else              cbind(newvar, eval(parse(text = vars.inf[[j]]$variables[k]), data))
               xdims[[i]] <- c(xdims[[i]], ncol(newvar) - sum(xdims[[i]]))
            }
         }
         if (is.matrix(newvar)) {
            nms <- colnames(newvar)
            if (any(is.null(colnames(newvar)))) 
                            nms <- paste(vars.inf[[j]]$variables, 
                                         "[", 1:ncol(newvar), "]", sep = "")
         }
         else 
            nms <- vars.inf[[j]]$variables
         xlab[[i]]    <- c(xlab[[i]], nms)
         # xlab[[i]]    <- c(xlab[[i]], vars.inf[[j]]$variables)
         newvar       <- as.matrix(newvar)
         ndims.new    <- ncol(newvar)
         ndims[[i]]   <- c(ndims[[i]], ndims.new)
         prd          <- vars.inf[[j]]$period
         if (length(prd) == 1 && is.na(prd)) prd <- rep(NA, ndims.new)
         if (length(prd) != ndims.new)
            stop("period does not match the columns of x.")
         period[[i]]  <- c(period[[i]], prd)
         if (any(!is.na(prd))) {         
            for (k in 1:ndims.new)
      	       if (!is.na(prd[k])) newvar[ , k] <- newvar[ , k] %% prd[k]
         }
         incr            <- vars.inf[[j]]$increasing
         increasing[[i]] <- c(increasing[[i]], incr)
         xrng <- vars.inf[[j]]$xrange
         if ((ndims.new == 1) & (length(xrng) == 2))
            xrng <- matrix(xrng, nrow = 1)
         if (!is.matrix(xrng))
            xrng <- matrix(NA, nrow = ndims.new, ncol = 2)
         if (nrow(xrng) != ndims.new)
            stop("xrange does not match columns of x.")
         for (k in 1:ndims.new) {
            if (any(is.na(xrng[k, ]))) {
               if (!is.na(prd[k]))
                  xrng[k, ] <- c(0, prd[k])
               else
                  xrng[k, ] <- c(min(newvar[ , k], na.rm = TRUE), max(newvar[ , k], na.rm = TRUE))
            # xrange <- t(apply(xrange, 1, function(x) c(x[1] - 0.05 * diff(x), x[2] + 0.05 * diff(x))))
            }
         }
         xrange[[i]]  <- rbind(xrange[[i]], xrng)
         fixed[[i]]   <- rbind(fixed[[i]], vars.inf[[j]]$fixed)
         X[[i]]       <- cbind(X[[i]], newvar)
         df.new       <- vars.inf[[j]]$df
         if (is.na(df.new)) df.new <- switch(ndims.new, 6, 12, 18)
         df[[i]]      <- c(df[[i]], df.new)
      }
#      if (any(is.na(nseg[[i]])) | prod(nseg[[i]]) > 400)
      if (any(is.na(nseg[[i]])))
         nseg[[i]] <- rep(switch(sum(ndims[[i]]), 100, 17, 7), sum(ndims[[i]]))
      if (any(is.na(X[[i]]))) xmissing  <- TRUE
   }
   
   # Remove observations which have missing data.
   # This may not be necessary as model.frame above may do it already.
   ind.missing <- lapply(X, function(x) apply(x, 1, function(z) any(is.na(z))))
   ind.missing <- cbind(is.na(y), matrix(unlist(ind.missing), ncol = length(X)))
   ind.missing <- apply(ind.missing, 1, any)
   if (any(ind.missing)) {
      y <- y[!ind.missing]
      for (i in 1:length(X)) X[[i]] <- as.matrix(X[[i]][!ind.missing, ])
      cat("warning: missing data removed.\n")
   }

   if (opt$verbose > 1) tim <- proc.time()
   	
   P <- list(length = length(terms.smooth))
   if (data.missing) {
      B.linear <- model.matrix(formula.linear, parent.frame())
   	  if (nrow(B.linear) == 0) B.linear <- matrix(1, nrow = length(y), ncol = 1)
      }
   else
      B.linear <- model.matrix(formula.linear, data)

   # If the data argument has not been used, but subset has, subset the data.
   if (data.missing & !subset.missing) {
   	  y <- y[subset]
      for (i in 1:length(X)) X[[i]] <- as.matrix(X[[i]][subset, ])
      B.linear <- as.matrix(B.linear[subset, ])
   }
   
   B <- B.linear
   m <- ncol(B)

   for (i in 1:length(terms.smooth)) {
      mat    <- ps.matrices(X[[i]], xrange[[i]], ndims = ndims[[i]], 
                     nseg = nseg[[i]], period = period[[i]])
   	  if (all(is.na(fac[[i]]))) {
         B      <- cbind(B, mat$B)
         m      <- c(m, ncol(mat$B))
         P[[i]] <- mat$P
   	  }
      else {
         Btemp <- matrix(nrow = length(y), ncol = 0)
         for (j in levels(fac[[i]]))
            Btemp <- cbind(Btemp, mat$B * as.numeric(fac[[i]] == j))
         B      <- cbind(B, Btemp)
         m      <- c(m, ncol(Btemp))
         nlevs  <- length(levels(fac[[i]]))
         pdim   <- nlevs * ncol(mat$P)
         P[[i]] <- matrix(0, nrow = pdim, ncol = pdim)
         for (j in 1:nlevs) {
         	ind <- (j - 1) * ncol(mat$P) + (1:ncol(mat$P))
            P[[i]][ind, ind] <- mat$P
         }
         P[[i]] <- P[[i]] + matrix(1, ncol = nlevs, nrow = nlevs) %x% diag(ncol(mat$B))
      }
      xrange[[i]] <- mat$xrange
   }
      
   if (opt$verbose > 1) {
      cat("Timings:\nconstructing matrices", (proc.time() - tim)[1], "seconds\n")
      tim <- proc.time()
   }

   b.ind <- list(length = length(m))
   for (i in 1:length(terms.smooth))
      b.ind[[i]] <- (cumsum(m)[i] + 1):cumsum(m)[i + 1]
   
   if (weights.missing) {
      # btb   <- t(B)  %*% B
      btb   <- crossprod(B)
      # Does crossprod also work with a vector, below?
      bty   <- t(B)  %*% y
   }
   else if (is.vector(weights)) {
      btb   <- t(B * weights)  %*% B
      bty   <- t(B * weights)  %*% y
   }
   else if (is.matrix(weights)) {
      btb   <- t(B) %*% weights %*% B
      bty   <- t(B) %*% weights %*% y
   }
   else
      stop("the weights argument is inappropriate.")

   if (opt$verbose > 1) {
      cat("matrix products", (proc.time() - tim)[1], "seconds\n")
      tim <- proc.time()
   }

   # Select the smoothing parameters, if required
   for (i in 1:length(terms.smooth)) {
      if (any(is.na(df[[i]]))) df[[i]] <- switch(sum(ndims[[i]]), 6, 12, 18)
      # code doesn't currently handle more than one df for terms with more than one variable.
      df[[i]] <- sum(df[[i]])
      if (df[[i]] > prod(nseg[[i]] + 3))
         stop(paste("df is too large for the value of nseg in term", i))
      if (any(is.na(lambda[[i]])))
      		lambda[[i]] <- lambda.select(btb[b.ind[[i]], b.ind[[i]]], P[[i]], df[[i]])
   }

   if (opt$verbose > 1) {
      cat("selecting smoothing parameters", (proc.time() - tim)[1], "seconds\n")
      tim <- proc.time()
   }

   # Fit the model
   Pall  <- matrix(0, nrow = ncol(B), ncol = ncol(B))
   for (i in 1:length(terms.smooth))
      Pall[b.ind[[i]], b.ind[[i]]] <- lambda[[i]] * P[[i]]
   B1       <- solve(btb + Pall) 
   alpha    <- as.vector(B1 %*% bty)
   df.model <- sum(btb * t(B1))
   df.error <- length(y) - sum(btb * (2 * B1 - B1 %*% btb %*% B1))
   
   # Force the estimate to pass through fixed points (1 covariate only)
   if (length(terms.smooth) == 1 & ndims[[1]] == 1 & all(!is.na(fixed[[1]]))) {
   	fxd <- fixed[[1]]
   	if (any(fxd[,1] < xrange[[1]][1]) | 
   			any(fxd[,1] > xrange[[1]][2]))
   		stop("fixed points must be inside the range of the data.")
   	A     <- cbind(1, ps.matrices(as.matrix(fxd[ , 1]), xrange[[1]], ndims[[1]],
   																nseg[[1]])$B)
   	alpha <- alpha +  B1 %*% t(A) %*% solve(A %*% B1 %*% t(A)) %*% 
   		(fxd[ , 2] - A %*% alpha)
   }   
   
   if (!random.missing) {
   	  nrnd   <- nlevels(random)
   	  n      <- length(y)
   	  ind    <- cbind(1:n, as.numeric(random))
   	  utu    <- diag(table(as.numeric(random)))
   	  btu    <- t(apply(B, 2, function(x) tapply(x, random, sum)))
   	  uty    <- tapply(y, random, sum)
      sig    <- 1
      sigr   <- 1
      sigo   <- 2
      sigro  <- 2
      p1     <- alpha
      p2     <- rep(0, nrnd)
      m1     <- c(B1 %*% bty)
      M1     <- B1 %*% btu
      while (abs(sig - sigo)/sigo > 1e-6 | abs(sigr - sigro)/sigro > 1e-6) {
         invu  <- solve(utu + diag(rep(sig^2/sigr^2, nrnd)))
         m2    <- as.vector(invu %*% uty)
         M2    <- invu %*% t(btu)
         p1    <- solve(diag(sum(m)) - M1 %*% M2) %*% (m1 - M1 %*% m2)
         p2    <- solve(diag(nrnd)   - M2 %*% M1) %*% (m2 - M2 %*% m1)
         ress  <- sum(y^2) + c(t(p1) %*% btb %*% p1 + t(p2) %*% utu %*% p2 - 2 * t(p1) %*% bty +
                     2 * t(p1) %*% btu %*% p2 - 2 * t(p2) %*% uty)
         nu    <- sum(diag(invu)) / sigr^2
         sigo  <- sig
         sigro <- sigr
         sig   <- sqrt(ress / (n - nrnd + nu))
         sigr  <- sqrt(sum(p2^2) / (nrnd - nu))
         # print(c(sig, sigr))
      }
      alpha     <- p1
      sigma     <- sig
      sigma.r   <- sigr
      N1        <- solve(diag(sum(m)) - M1 %*% M2)
      U1        <- solve(utu + diag(nrnd) * sigma^2 / sigma.r^2)
      cov.alpha <- sigma^2 * btb -
      	           sigma^2 * btu %*% U1 %*% t(btu) +
      	           sigma.r^2 * btu %*% t(btu) -
      	           sigma.r^2 * btu %*% utu %*% U1 %*% t(btu) - 
      	           sigma^2 * btu %*% U1 %*% t(btu) +
      	           sigma^2 * btu %*% U1 %*% utu %*% U1 %*% t(btu) -
      	           sigma.r^2 * btu %*% U1 %*% utu %*% t(btu) +
      	           sigma.r^2 * btu %*% U1 %*% utu %*% utu %*% U1 %*% t(btu)
      cov.alpha <- N1 %*% B1 %*% cov.alpha %*% B1 %*% t(N1)
      rss       <- NULL
      R.squared <- NULL
      mu         <- c(B %*% alpha)
   }
   else {
     mu         <- c(B %*% alpha)
     sigma      <- sqrt(sum((y - mu)^2) / df.error)
     cov.alpha  <- B1 %*% btb %*% t(B1) * sigma^2
     rss        <- sum((y - mu)^2)
     tss        <- sum((y - mean(y))^2)
     R.squared  <- 100 * (tss - rss) / tss
   }
   
   # Increasing function: specific to nseg 17 and so 20 basis fns in each dim.
   if (random.missing && (length(ndims) == 1) && ndims[[1]] == 2 &&
   		 increasing[[1]] && all(nseg[[1]] == 17) && bdeg == 3) {
   	 D1    <- diff(diag(20)) %x% diag(20)
   	 D2    <- diag(20)       %x% diff(diag(20))
   	 delta <- 1
   	 while (delta > 1e-5) {
   	 	 v1        <- as.numeric(c(D1 %*% alpha[-1]) <= 0)
   	 	 v2        <- as.numeric(c(D2 %*% alpha[-1]) <= 0)
   	 	 mat1      <- 100 * lambda[[1]] * t(D1) %*% diag(v1) %*% D1
   	 	 mat2      <- 100 * lambda[[1]] * t(D2) %*% diag(v2) %*% D2
   	 	 mat       <- matrix(0, nrow = ncol(B), ncol = ncol(B))
   	 	 mat[2:401, 2:401] <- mat1 + mat2
   	 	 B1        <- solve(btb + Pall + mat)
   	 	 alpha.old <- alpha
   	 	 alpha     <- as.vector(B1 %*% bty)
   	 	 delta     <- sum((alpha - alpha.old)^2) / sum(alpha.old^2)
   	 }
   	 mu       <- c(B %*% alpha)
   	 df.model <- NULL
   	 df.error <- NULL
   }
   
	# P[[1]] <- P[[1]] %x% diag(m[2])
	# P[[2]] <- diag(m[2]) %x% P[[2]]

   if (opt$verbose > 1) {
      cat("fitting", (proc.time() - tim)[1], "seconds\n")
      tim <- proc.time()
   }
   
   if (opt$verbose > 1) {
      cat("summaries", (proc.time() - tim)[1], "seconds\n")
      tim <- proc.time()
   }

   # If there is only one term, include the mean
   # if (nterms == 1) b.ind[[1]] <- c(1, b.ind[[1]])

   result <- list(fitted = mu, alpha = alpha, m = m, B = B, 
                  bty = bty, btb = btb, B1 = B1, Pall = Pall, xlabels = xlabels,
                  linear.matrix = B.linear,
                  terms.linear = terms.linear, terms.smooth = terms.smooth,
                  xlab = xlab, ylab = ylab, term.labels = term.labels,
                  lambda = lambda, ndims = ndims, xdims = xdims,
                  y = y, X = X, fac = fac, Xlinear = Xlinear,
                  bricks.type = bricks.type,
                  sigma = sigma, cov.alpha = cov.alpha, b.ind = b.ind,
                  df = df, df.model = df.model, df.error = df.error, 
                  rss = rss, R.squared = R.squared, xrange = xrange,
                  nseg = nseg, bdeg = bdeg, period = period, 
   							  increasing = increasing, pam.formula = pam.formula,
                  formula.linear = formula.linear,
                  involved = involved, nterms = nterms)
   if (!weights.missing) result$weights <- weights
   if (!random.missing)  result$sigma.r <- sigma.r
   class(result) <- "pam"
   
   # if (nterms == 1 & ndims[[1]] <= 2) {
   # 	  if (missing(model)) model <- "none"
   # 	  colnames(result$X[[1]]) <- xlab[[1]]
   #    if (ndims[[1]] == 1) sm.regression(result$X[[1]], result$y, h = h, model = model,
   #                                          xlab = xlab[[1]], ylab = ylab, ...)
   #    else                 sm.regression(result$X[[1]], result$y, h = h, model = model, zlab = ylab, ...)
   #    return(invisible(result))
   # }
   # else
   if (opt$display != "none") plot(result, ...)
   
   # if (nterms == 1 & ndims[[1]] <= 2) {
      # if (opt$panel) {
         # replace.na(opt, df.max, switch(ndims[[1]], 20, 50, 100))
      	 # df.min  <- switch(ndims[[1]], 2, 4, 8) + 0.1
         # df.max  <- if (!opt$panel) df[[1]] else min(length(y) - 5, opt$df.max)
         # df.min  <- if (!opt$panel) df[[1]] else df.min
         # Pall    <- rbind(0, cbind(0, P[[1]]))
         # llambda    <- 0
         # llambda.df <- lambda.df(exp(max(llambda)), btb, Pall)
         # while (min(llambda.df) >= df.min) {
            # llambda    <- c(llambda, max(llambda) + log(10))
            # llambda.df <- c(llambda.df, lambda.df(exp(max(llambda)), btb, Pall))
         # }
         # while (max(llambda.df) <= df.max) {
            # llambda    <- c(llambda, min(llambda) - log(10))
            # llambda.df <- c(llambda.df, lambda.df(exp(min(llambda)), btb, Pall))
         # }
         # df.fun <- approxfun(llambda.df, llambda)
   
         # sm.pam.draw <- function(pam.panel) {
         	# plot(pam.panel$model, options = pam.panel$opt)
            # # title(pam.panel$df)
            # pam.panel
         # }
         # sm.pam.redraw <- function(pam.panel) {
            # # pam.panel$model$lambda <- lambda.select(pam.panel$model$btb, pam.panel$model$bty,
            # #                                         Pall, pam.panel$df)
            # pam.panel$model$lambda <- exp(pam.panel$df.fun(pam.panel$df))
            # B1 <- solve(pam.panel$model$btb + pam.panel$model$lambda * pam.panel$Pall)
            # pam.panel$model$alpha  <- as.vector(B1 %*% pam.panel$model$bty)
            # pam.panel$model$fitted <- c(pam.panel$model$B %*% pam.panel$model$alpha)
            # pam.panel$opt$se       <- pam.panel$se
            # pam.panel$opt$theta    <- pam.panel$theta
            # pam.panel$opt$phi      <- pam.panel$phi
            # rp.control.put(pam.panel$panelname, pam.panel)
            # rp.tkrreplot(pam.panel, plot)
            # pam.panel
         # }
         # opt1 <- opt
         # opt1$panel <- FALSE
         # pam.panel <- rp.control(model = result, opt = opt1, Pall = rbind(0, cbind(0, P[[1]])),
                                 # df = opt$df, df.fun = df.fun, theta = opt$theta, phi = opt$phi)
         # rp.tkrplot(pam.panel, plot, sm.pam.draw, hscale = opt$hscale, vscale = opt$vscale, pos = "right")
         # rp.slider(pam.panel, df, df.min, df.max, sm.pam.redraw, showvalue = TRUE)
         # rp.checkbox(pam.panel, se, sm.pam.redraw, title = "Standard errors")
         # if (ndims[[1]] == 2) {
            # rp.slider(pam.panel, theta, -180, 180, sm.pam.redraw, "persp angle 1")
            # rp.slider(pam.panel, phi,      0,  90, sm.pam.redraw, "persp angle 2")
         # }
      # }
      # else if (opt$display != "none")
         # plot(result, ...)
   # }
   
   invisible(result)
}

#----------------------------------------------------------------------------

lambda.df <- function(lambda, btb, P, df) {
	B1   <- solve(btb + lambda * P)
	sum(diag(btb %*% B1)) - df
}

lambda.select <- function(btb, P, df, method = "df") {
#    This currently uses the same lambda in all dimensions
  if (method == "df") {
    lambda <- c(1, 1)
    f      <- rep(lambda.df(lambda[1], btb, P, df), 2)
    mult   <- if (f[1] < 0) 0.1 else 10
    while (sign(f[1]) == sign(f[2])) {
      lambda <- c(lambda[2], lambda[2] * mult)
      f      <- c(f[2], lambda.df(lambda[2], btb, P, df))
    }
    if (diff(f) > 0) {
      lambda <- rev(lambda)
      f      <- rev(f)
    }
    lambda <- uniroot(lambda.df, interval = lambda, btb, P, df,
    									f.lower = f[1], f.upper = f[2])$root

       # 	  lambda <- 1
       # 	  f.new  <- 0
       # 	  while (f.new <= df) {
       # 	  	lower   <- lambda
       # 	  	f.lower <- f.new
       # 	  	lambda  <- lambda / 10
       # 	  	f.new   <- lambda.df(lambda, btb, P)
       # 	  }
       # 	  lambda <- 1
       # 	  f.new  <- df + 1
       # 	  while (f.new >= df) {
       # 	  	upper   <- lambda
       # 	  	f.upper <- f.new
       # 	  	lambda  <- lambda * 10
       # 	  	f.new   <- lambda.df(lambda, btb, P)
       # 	  }
       # print(c(lower, upper))
       # # 	  
       # # 	  lower <- lambda / 10
       # #  	lambda  <- 0.1
       # #  	f.upper <- df + 1
       # #  	while (f.upper >= df) {
       # #  		lambda  <- lambda * 10
       # #  		f.upper <- lambda.df(lambda, btb, P)
       # #  	}
       # # 	 upper <- lambda * 10
       # # 	 f.lower <- lambda.df(lambda, btb, P)
       # #    lower  <- lambda
       # #    lambda <- 1
       # #    f.lambda <- df - 1
       # #    while (lambda.df(lambda, btb, P) >= df) lambda <- lambda * 10
       # #    upper  <- lambda
       #    lambda.crit <- function(lambda, btb, P, df)
       #       lambda.df(lambda, btb, P) - df
       #    cat("entering uniroot ... \n")
       #    result <- uniroot(lambda.crit, interval = c(lower, upper), 
       #    									f.lower = f.lower, f.upper = f.upper, btb, P, df)
       #    # cat("result$root", result$root, "\n")
       #    lambda <- result$root
       #    cat("leaving uniroot ... \n")
    }
  lambda
}

#----------------------------------------------------------------------------

predict.pam <- function(model, newdata, se.fit = FALSE, verbose = 1, deriv = 0) {

   if (!is.list(newdata)) {
      newdata <- list(newdata)
      # names(newdata) <- vars.inf[[2]]$variables
      names(newdata) <- model$xlabels[[1]]
      print(newdata)
   }
   if (!all(unlist(model$xlabels) %in% names(newdata)))
      stop("some required variables are not present in the new data.")

   nnew <- if (is.matrix(newdata[[1]])) nrow(newdata[[1]])
           else                         length(newdata[[1]])
   X    <- list()

   for (i in 1:length(model$terms.smooth)) {
   	  ii     <- model$terms.smooth[i]
      inv    <- which(model$involved[ , ii] == 1) - 1
      nvars  <- length(inv)
      X[[i]] <- matrix(nrow = nnew, ncol = 0)
      for (j in 1:length(inv)) {
      	 newvar  <- eval(parse(text = model$xlabels[[i]][j]), newdata)
         X[[i]]  <- cbind(X[[i]], newvar)
      }
   }

   inrange <- rep(TRUE, nnew)
   for (i in 1:length(model$terms.smooth))
      for (j in 1:ncol(X[[i]]))
   	     inrange <- inrange & X[[i]][ , j] >= model$xrange[[i]][j, 1] & 
   	                          X[[i]][ , j] <= model$xrange[[i]][j, 2]
   outrange <- which(!inrange)
   if (length(outrange) > 0 & verbose > 0)
      warning("some evaluation points are out of range and have been omitted.")
   nnew <- length(which(inrange))
   
   B.linear <- model.matrix(model$formula.linear, newdata)
   if (nrow(B.linear) == 0) B.linear <- matrix(1, nrow = nnew, ncol = 1)
   
   B    <- matrix(c(B.linear[inrange, ]), ncol = ncol(B.linear))
   for (i in 1:length(model$terms.smooth)) {
      mat    <- ps.matrices(as.matrix(X[[i]][inrange, ]), xrange = model$xrange[[i]], 
                     ndims = model$ndims[[i]], nseg = model$nseg[[i]],
      							 period = model$period[[i]])
      B      <- cbind(B, mat$B)
   }
   
   fv          <- rep(NA, nnew)
   fv[inrange] <- c(B %*% model$alpha)

   if (model$nterms == 1 & model$ndims[[1]] == 1 & deriv > 0) {
     mat         <- ps.matrices(as.matrix(X[[1]][inrange, ]), xrange = model$xrange[[1]], 
                           ndims = model$ndims[[1]], nseg = model$nseg[[1]], bdeg = model$bdeg - deriv, 
                           period = model$period[[1]])
     h           <- (model$xrange[[1]][,2] - model$xrange[[1]][,1]) / mat$nseg
     D           <- diff(diag(length(model$alpha[-1])), differences = deriv)
     fv[inrange] <- model$alpha[1] + c(mat$B %*% D %*% model$alpha[-1]) / (h^deriv)
   }
   # if (model$nterms == 1 & model$ndims[[1]] == 1 & deriv > 0) {
   #   mat    <- ps.matrices(as.matrix(X[[1]][inrange, ]), xrange = model$xrange[[1]], 
   #                         ndims = model$ndims[[1]], nseg = model$nseg[[1]], bdeg = model$bdeg - deriv, 
   #                         period = model$period[[1]])
   #   alpha1 <- diff(model$alpha[-1], differences = deriv)
   #   h      <- model$xrange[[1]][,2] - model$xrange[[1]][,1]
   #   fv[inrange] <- c(mat$B %*% alpha1) / (h / mat$nseg)^deriv
   # }
   
   results     <- list(fit = fv)
   
   if (se.fit)
      results$se.fit = sqrt(diag(B %*% model$cov.alpha %*% t(B)))
   
   return(invisible(results))
}

#----------------------------------------------------------------------------

anova.pam <- function(model, terms = 1:length(model$b.ind), method = "QF", verbose = 1) {

   terms.obj <- terms(model$pam.formula, specials = "s")
   involved  <- attr(terms.obj, "factors")
   ord       <- attr(terms.obj, "order")
   ind       <- vector("logical", length = 0)
   for (i in terms) {
      inv    <- which(involved[ , i] > 0)
      ind <- c(ind, any(apply(involved, 2, function(x) all(x[inv] > 0)) & (ord > ord[i])))
   }
   terms <- terms[!ind]

   n         <- length(model$y)
   btb       <- model$btb
   B1        <- model$B1
   p         <- vector(length = 0)
   
   for (i in terms) {

      ind  <- model$b.ind[[i]]
      ind  <- ind[ind != 1]
      # weighted by covariance matrix
      eig  <- eigen(model$cov.alpha[ind, ind] / model$sigma^2)
      pos  <- which(abs(eig$values) > .Machine$double.eps)
      n0   <- length(eig$values) - length(pos)
      if (n0 > 0 & verbose > 1)
         cat(n0, "eigenvalues omitted out of", length(eig$values), "\n")
      inv  <- eig$vectors[, pos]  %*% diag(1/eig$values[pos]) %*% t(eig$vectors[, pos])
      # A1   <- B1[ , ind] %*% inv %*% B1[ind, ]
      # simple unweighted form
      A1   <- B1[ , ind] %*% B1[ind, ]
      Fobs <- c(model$y %*% model$B %*% A1 %*% t(model$B) %*% model$y) / model$rss
      # Fobs <- sum(model$alpha[ind]^2) / model$rss
      
      if (method == "QF") {
         A2   <- 2 * B1 - B1 %*% btb %*% B1
         A3   <- A1 + Fobs * A2
         A4   <- A3 %*% btb %*% A3 - 2 * Fobs * A3
         A5   <- A4 %*% btb %*% A3 + Fobs^2 * A3 - Fobs * A4
         k1   <-     sum(btb * A3) - n * Fobs
         k2   <- 2 * sum(btb * A4) + n * Fobs^2
         k3   <- 8 * sum(btb * A5) - n * Fobs^3
         aa   <- abs(k3 / (4 * k2))
         bb   <- (8 * k2^3) / k3^2
         cc   <- k1 - aa * bb
         p    <- c(p, 1 - pchisq(-cc / aa, bb))
      }

      else if (method == "F") {
         # P1   <- diag(nrow(P1)) - P1
         # P2   <- diag(nrow(P2)) - P2
         # P1   <- t(P1) %*% P1
         # P2   <- t(P2) %*% P2
         # df1  <- n - sum(diag(P1))
         # df2  <- n - sum(diag(P2))
         # rss1 <- c(y %*% P1 %*% y)
         # rss2 <- c(y %*% P1 %*% y)
         # Fobs <- (rss1 - rss2) * (n - df2) / (rss2 * (df2 - df1))
      	 df.m <- sum(btb * A1)
      	 # print(c(model$df[[i]], df.m))
      	 Fobs <- Fobs * model$df.error / df.m
         p    <- c(p, 1 - pf(Fobs, model$df.m, model$df.error))
         # print(c(Fobs, p))
      	 # print(c(Fobs, model$df.error, model$df.model, df.m))
         # return(invisible(list(p = pval, Fobs = Fobs, df1 = df1, df2 = df2)))
      }
      else 
         stop("method not recognised.")
         
   }
   
   names(p) <- model$term.labels[terms]
   pmat <- matrix(round(p, 3), ncol = 1, dimnames = list(names(p), "p-value"))
   if (verbose > 0) print(pmat)
   return(invisible(list(p = p)))

}

#----------------------------------------------------------------------------

summary.pam <- function(model, verbose = 1) {
	
   lcoefs   <- model$alpha[1:model$m[1]]
   se       <- sqrt(model$cov.alpha[cbind(1:model$m[1], 1:model$m[1])])
   d.linear <- data.frame(Estimate = lcoefs, se = se,
                     row.names = paste("   ", colnames(model$linear.matrix)))
   if (verbose > 0) {
      cat("Linear terms:\n")
      print(d.linear)
   }
   
   nterms <- length(model$xlabels)
   adf    <- rep(0,  nterms)
   nvars  <- rep(1,  nterms)
   xl     <- rep("", nterms)
   low    <- rep(NA, nterms)
   high   <- rep(NA, nterms)
   for (i in 1:nterms) {
      l            <- paste("s(", model$xlabels[[i]], ")", sep = "")
      if (length(l) == 2) l <- paste(l[1], l[2], sep = ":")
      if (length(l) == 3) l <- paste(l[1], l[2], l[3], sep = ":")
      xl[i]        <- l
      ind          <- model$b.ind[[i]]
      btb.i        <- model$btb * 0
      btb.i[, ind] <- model$btb[ , ind]
      adf[i]       <- round(sum(btb.i * t(model$B1)), 1)
      nvars[i]     <- length(model$xlabels[[i]])
      mu.i         <- model$B[ , ind] %*% model$alpha[ind]
      low[i]       <- signif(min(mu.i))
      high[i]      <- signif(max(mu.i))
   }
   adf             <- c(adf,   round(model$df.model, 1))
   nvars           <- c(nvars, length(unique(unlist(model$xlabels))))
   low             <- c(low,   min(model$fitted))
   high            <- c(high,  max(model$fitted))
   d.smooth        <- data.frame(adf, nvars, high - low,
                           row.names = paste("   ", 1:length(adf), c(xl, "model")))
   names(d.smooth) <- c("df", "nvars", "range")
   if (verbose > 0) {
      cat("\nSmooth terms:\n")
      print(d.smooth)
      cat("\nError s.d.: ",  signif(model$sigma))
      if ("sigma.r" %in% names(model))
         cat("\nR.e.  s.d.: ",  signif(model$sigma.r))
      if (!is.null(model$R.squared))
         cat("\nR-squared :  ", round(model$R.squared, 1), "%", sep = "")
   }
   
   invisible(list(linear = d.linear, smooth = d.smooth))
}
   
#----------------------------------------------------------------------------

plot.pam <- function(model, components = 1:length(model$xlabels), plotinfo, 
                   options = list(), ...) {

   opt <- if (length(options) > 0) sm.options(options) else sm.options(list(...))

   if (any(components > model$nterms)) stop("some components do not match model terms.")
   
   ngrid <- opt$ngrid
   replace.na(opt, reference,               "none")
   replace.na(opt, nlevels,                 20)
   replace.na(opt, panel,                   TRUE)
   replace.na(opt, panel.plot,              TRUE)
   replace.na(opt, display,                 "image")
   replace.na(opt, eqscplot,                FALSE)
   replace.na(opt, deriv,                   NA)
   replace.na(opt, deriv.order,             0)
   if (!("include.lower.terms" %in% names(list(...)))) {
     if (model$nterms == 1) opt$include.lower.terms <- TRUE
   }
   replace.na(opt, include.lower.terms,     FALSE)
   inc.mn <- opt$include.lower.terms & (is.na(opt$deriv) | (opt$deriv.order == 0))
   replace.na(opt, include.mean, inc.mn)
   replace.na(opt, include.lower.reference, FALSE)
   replace.na(opt, se,                      FALSE)
   replace.na(opt, lwd,                     2)
   if (opt$reference != "none") opt$se <- TRUE
   # col.pal <- if (requireNamespace("RColorBrewer", quietly = TRUE)) 
                 # rev(colorRampPalette(RColorBrewer::brewer.pal(9, "BuGn"))(100))
                 # else topo.colors(100)
   if (requireNamespace("colorspace", quietly = TRUE)) 
      col.pal <- if (opt$include.lower.terms) colorspace::heat_hcl(100) else colorspace::diverge_hcl(100)
   else
      col.pal <- topo.colors(100)
   replace.na(opt, col.palette, col.pal)
   replace.na(opt, order, 1:3)
   if (!("vscale" %in% names(list(...)))) opt$vscale <- opt$hscale
   
   if (!is.na(opt$deriv)) {
   	  if (length(components) > 1)
         stop("only one component can be specified when a derivative is requested.")
      if (!(opt$deriv %in% model$xlabels[[components]]))
      	 stop("the requested derivative is not in the specified model term.")
   	  if (opt$deriv.order == 0) opt$deriv <- NA
   }
   
   # To handle additive terms
   if (is.character(components)) {
      comps      <- as.numeric(unlist(strsplit(components, "+", fixed = TRUE)))
      components <- comps[1]
      additive   <- TRUE
   }
   else
      additive <- FALSE
   # Code not yet followed through
   
   display           <- opt$display
   se                <- opt$se
   ylim              <- opt$ylim
   partial.residuals <- opt$partial.residuals
   ngrid.missing     <- is.na(ngrid)
   
   if (missing(plotinfo)) {
      plotinfo <- list()
      for (i in components) {
   
         x     <- model$X[[i]]
         xdims <- model$xdims[[i]]
         
         if (is.vector(x)) x <- matrix(x, ncol = 1)
         ndim <- sum(model$ndims[[i]])
         if (ngrid.missing) ngrid <- switch(ndim, 100, 20, 15)
         u   <- list(length = ndim)
         for (j in 1:ndim)
            u[[j]] <- seq(model$xrange[[i]][j, 1], model$xrange[[i]][j, 2], length = ngrid)
         U   <- as.matrix(expand.grid(u))

         if (opt$include.lower.terms) {
         	  ind.terms <- which(sapply(model$xlabels, function(x) all(x %in% model$xlabels[[i]])))
         	  if (!is.na(opt$deriv)) {
         	  	ind.t     <- integer(0)
         	    for (k in ind.terms)
         	       if (opt$deriv %in% model$xlabels[[k]]) ind.t <- c(ind.t, k)
         	    ind.terms <- ind.t
         	  }
         }
         else
            ind.terms <- i

         deriv <- if (is.na(opt$deriv)) 0 else opt$deriv.order
         	
         B   <- matrix(nrow = nrow(U), ncol = 0)
         ind <- integer(0)
         
         for (j in ind.terms) {
         	  Uj <- matrix(nrow = nrow(U), ncol = 0)
         	  for (k in 1:length(model$xlabels[[j]])) {
         	     indk <- match(model$xlabels[[j]][k], model$xlabels[[i]])
         	     cs   <- cumsum(model$xdims[[i]])[indk]
         	     indk <- (cs - model$xdims[[i]][indk] + 1):cs
               Uj   <- cbind(Uj, U[ , indk])
            }
            # Uj  <- U[, match(model$xlabels[[j]], model$xlabels[[i]])]
            # if (!is.matrix(Uj)) Uj <- matrix(Uj, ncol = 1)
         	  ind.d <- which(model$xlabels[[j]] == opt$deriv)
         	  bdg   <- rep(model$bdeg, sum(model$ndims[[j]]))
         	  bdg[ind.d] <- bdg[ind.d] - deriv
         	  Bj <- ps.matrices(Uj, model$xrange[[j]], model$ndims[[j]], nseg = model$nseg[[j]],
         	  									bdeg = bdg, period = model$period[[j]], penalty = FALSE)$B
         	  if (!is.na(opt$deriv) & deriv > 0) {
         	  	 h  <- (model$xrange[[j]][ind.d, 2] - model$xrange[[j]][ind.d, 1]) / model$nseg[[j]][ind.d]
               I  <- diag(model$nseg[[j]][1] + model$bdeg)
         	     D1 <- diff(I, differences = deriv)
         	     D  <- 1
         	     for (k in 1:length(bdg)) {
         	     	  Mk <- if (k == ind.d) D1 else I
         	     	  D  <- D %x% Mk
         	     }
         	     Bj <- Bj %*% D / (h^deriv)
            }
         	  if (!opt$include.lower.reference & j == i) {
               Bi   <- Bj
               indi <- model$b.ind[[i]]
            }
            B   <- cbind(B, Bj)
            ind <- c(ind, model$b.ind[[j]])
         }
                       
         if (!all(is.na(model$fac[[i]]))) {
            nlevs <- length(levels(model$fac[[i]]))
            Btemp <- B
            if (nlevs > 1)
               for (j in 2:nlevs) B <- cbind(B, Btemp)
         }
       
         if (opt$include.mean) {
            B   <- cbind(1, B)
            ind <- c(1, ind)
         }
         est <- c(B %*% model$alpha[ind])
         if (ndim > 1) est <- array(est, dim = rep(ngrid, ndim))
      
         mui <- if (is.na(opt$deriv)) c(model$B[ , ind] %*% model$alpha[ind])
                else rep(NA, length(model$y))
         model.fitted <- model$fitted

         result <- list(u = u, est = est, ndim = ndim, component = i, mui = mui, 
                        x = x, xdims = xdims, xlab = model$xlab[[i]], ylab = model$ylab,
                        model.fitted = model.fitted)
         if (se) {
            if (opt$reference == "no effect" & !opt$include.lower.reference) {
               est <- c(Bi %*% model$alpha[indi])
               if (ndim > 1) est <- array(est, dim = rep(ngrid, ndim))
            	 B   <- Bi
            	 ind <- indi
            }
            result$st.error <- sqrt(diag(B %*% model$cov.alpha[ind, ind] %*% t(B)))
            if (opt$reference == "no effect") {
            	 dm <- if (ndim > 1) dim(est) else length(est)
               result$reference <- est / array(result$st.error, dm)
            }
         }
         plotinfo[[length(plotinfo) + 1]] <- result
      }
   }
   else {
      mui <- plotinfo[[1]]$mui
      est <- plotinfo[[1]]$est
      ngrid <- length(plotinfo[[1]]$u[[1]])
   }
   
   if (display == "none") return(plotinfo)
   
   # One smooth term and one factor: sm.ancova
   if (length(which(model$bricks.type == "list"))   == 1 &&
       length(which(model$bricks.type == "factor")) == 1 &&
       plotinfo[[1]]$ndim == 1) {
      Xl    <- model$Xlinear[[1]]
      lvls  <- levels(Xl)
      nlvls <- length(lvls)
      cfs   <- model$alpha[1:model$m[1]]
      est   <- matrix(plotinfo[[1]]$est, nrow = length(est), ncol = nlvls)
      for (i in 1:nlvls)
         est[ , i] <- est[ , i] + 
              cfs %*% model$linear.matrix[min(which(Xl == lvls[i])), ]
      if (length(model$term.labels) == 3) {
         for (i in 1:nlvls) {
            cols      <- ncol(B) / nlvls
            cols      <- (i - 1) * cols + 1:cols
            ind       <- model$b.ind[[2]][cols]
            est[ , i] <- est[ , i] + c(B[ , cols] %*% model$alpha[ind])
         }
      }
      if (any(is.na(ylim))) ylim <- range(model$y, est)
      replace.na(opt, col, 1 + 1:nlvls)
   	  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("the ggplot2 package is required.")
   	  xgrid <- plotinfo[[1]]$u[[1]]
   	  dfrm  <- data.frame(xgrid = rep(xgrid, nlvls), est = c(est), lvl = factor(rep(1:nlvls, each = length(xgrid))))
      dfrm1 <- data.frame(x = plotinfo[[1]]$x, y = model$y, lvl = model$Xlinear[[1]])
      plt   <- ggplot2::ggplot(dfrm, ggplot2::aes(x = xgrid, y = est, colour = lvl)) +
               ggplot2::geom_line() + ggplot2::geom_point(data = dfrm1, ggplot2::aes(x, y, colour = lvl)) +
               ggplot2::labs(colour = names(model$Xlinear)[1]) +
               ggplot2::xlab(model$xlabels[[1]]) + ggplot2::ylab(model$ylab)
      print(plt)
   	  # matplot(plotinfo[[1]]$u[[1]], est, type = "l", lty = 1, ylim = ylim,
   	                  # xlab = model$xlabels[[1]], ylab = "component")
   	  # points(x, model$y, col = opt$col[as.numeric(Xl)])
      # # for (i in 1:nlvls)
      # #    lines(plotinfo[[1]]$u[[1]], est + shift[i],
      # #              col = opt$col[i], lty = opt$lty, lwd = opt$lwd)
      # title(paste(c(model$term.labels[model$terms.linear], ":", lvls), collapse = " "))
      # title("put the levels in different colours?", line = 1, cex.main = 0.7)
   	  return(invisible(plotinfo))
   }
   
   if (!is.na(opt$deriv)) partial.residuals <- FALSE

   # 1-d terms
   ndim1 <- vector(length = 0)
   for (i in 1:length(plotinfo))
   	  if (plotinfo[[i]]$ndim == 1) ndim1 <- c(ndim1, i)

   if (length(ndim1) > 0) {
      est.all  <- numeric()
      u.all    <- numeric()
      x.all    <- numeric()
      xx.all   <- numeric()
      pres.all <- numeric()
      y.all    <- numeric()
      se1.all  <- numeric()
      se2.all  <- numeric()
      lbl.x    <- character()
      lbl.u    <- character()
      lbl      <- character()
      xrng     <- vector("list", length(ndim1))
      for (i in 1:length(ndim1)) {
      	 pinf      <- plotinfo[[ndim1[i]]]
      	 uu        <- pinf$u[[1]]
      	 x.add     <- uu
         if (partial.residuals)
            x.add  <- c(c(pinf$x), uu)
         y.add     <- pinf$est
         if (partial.residuals)
         	   y.add <- c(model$y - model.fitted + pinf$mui, y.add)
         u.all     <- c(u.all, uu)
         est.all   <- c(est.all, pinf$est)
         lblx.add  <- rep("line", length(uu))
         if (partial.residuals)
         	  lblx.add  <- c(rep("data", length(pinf$x)), lblx.add)
         if (se | opt$reference != "none") {
         	  se1.add  <-  2 * pinf$st.error
         	  se2.add  <- -2 * pinf$st.error
         	  if (opt$reference == "none") {
         	  	se1.add <- se1.add + pinf$est
         	  	se2.add <- se2.add + pinf$est
         	  }
         	  se1.all  <- c(se1.all,  se1.add)
         	  se2.all  <- c(se2.all,  se2.add)
         	  x.add    <- c(x.add, uu, rev(uu))
         	  y.add    <- c(y.add, se1.add, rev(se2.add))
         	  lblx.add <- c(lblx.add, rep("se", 2 * length(uu)))
         }
         x.all     <- c(x.all, x.add)
         y.all     <- c(y.all, y.add)
         lbl.x     <- c(lbl.x, lblx.add)
         lbl       <- c(lbl,   rep(model$xlabels[[components[ndim1][i]]], length(x.add)))
         xrng[[i]] <- model$xrange[[ndim1[i]]]
      }
   	  if (any(is.na(ylim)))        ylim <- range(y.all)
      if (se)                      ylim <- range(ylim, se1.all,  se2.all)
      if (!requireNamespace("ggplot2", quietly = TRUE)) stop("the ggplot2 package is required.")
   	  dfrm  <- data.frame(x.all, y.all, lbl.x, lbl)
   	  # dfrm1 <- data.frame(u.all, est.all, se1.all, se2.all)
   	  plt   <- ggplot2::ggplot(dfrm, ggplot2::aes(x.all, y.all))
   	  if (se | opt$reference != "none")
   	  	plt <- plt + ggplot2::geom_polygon(data = subset(dfrm, lbl.x == "se"), fill = "lightblue")
   	  if (partial.residuals) plt <- plt + ggplot2::geom_point(data = subset(dfrm, lbl.x == "data"))
   	  if (all(!is.na(opt$xlim))) plt <- plt + ggplot2::xlim(opt$xlim[1], opt$xlim[2])
   	  plt   <- plt + ggplot2::ylim(ylim[1], ylim[2])
   	  plt   <- plt + ggplot2::geom_line(data = subset(dfrm, lbl.x == "line"), col = "blue", size = opt$lwd)
   	  plt   <- plt + ggplot2::facet_wrap( ~ lbl, scales = "free_x")
   	  plt   <- plt + ggplot2::xlab("") + ggplot2::ylab("partial residuals")
   	  print(plt)
   }
   
   if (length(ndim1) < length(plotinfo)) {
      if (length(ndim1) > 0)
         ndim2 <- (1:length(plotinfo))[-ndim1]
      else
         ndim2 <- 1:length(plotinfo)
      for (i in ndim2) {
      	 env <- environment()
         with(plotinfo[[i]], {
         # del1  <- u[[1]][2] - u[[1]][1]
         # del2  <- u[[2]][2] - u[[2]][1]
         # ugrid <- as.matrix(expand.grid(u[[1]], u[[2]]))
         
         if (ndim == 2) {
            mask  <- sm.mask(x[ , opt$order[1:2]], cbind(u[[opt$order[1]]], u[[opt$order[2]]]),
                             mask.method = opt$mask.method)
    	      if (any(is.na(ylim))) ylim <- range(est * mask, na.rm = TRUE)
            if (!opt$include.lower.terms) ylim <- c(-1, 1) * max(abs(ylim))
            if (length(xlab) == 1)
               xlab <- if (!is.null(colnames(x))) colnames(x) else paste(xlab, 1:2, sep = "-")
            clr <- "green"
            if (exists("reference")) {
            	sdiff <- array(c(reference[-ngrid, -ngrid], reference[    -1, -ngrid],
            									 reference[-ngrid,     -1], reference[    -1,     -1]),
            								 dim = c(ngrid - 1, ngrid - 1, 4))
            	sdiff <- apply(sdiff, 1:2, function(x) 
            		if (length(which(is.na(x))) > 1) NA else mean(x, na.rm = TRUE))
            	sdiff <- matrix(c(sdiff), nrow = ngrid - 1, ncol = ngrid - 1)
            	se.breaks <- c(-3, -2, 2, 3)
            	col.pal   <- rev(rainbow(length(se.breaks) + 1, start = 0/6, end = 4/6))
            	se.breaks <- c(min(-3, sdiff, na.rm = TRUE) - 1, se.breaks, 
            								 max(3, sdiff, na.rm = TRUE) + 1)
            	clr <- col.pal[cut(c(sdiff), se.breaks, labels = FALSE)]
            	clr <- matrix(clr, ngrid - 1, ngrid - 1)
            }
            if (display %in% c("persp", "lines")) {
               if (any(is.na(opt$x1lim))) opt$x1lim <- range(u[[opt$order[1]]])
               if (any(is.na(opt$x2lim))) opt$x2lim <- range(u[[opt$order[2]]])
               persp(u[[opt$order[1]]], u[[opt$order[2]]], aperm(est * mask, opt$order[1:2]), 
                     xlab = xlab[opt$order[1]], ylab = xlab[opt$order[2]], zlab = ylab,
                     ticktype = "detailed", col = c(clr),
                     xlim = opt$x1lim, ylim = opt$x2lim, zlim = ylim,
                     d = 10, theta = opt$theta, phi = opt$phi)
            }
            else if ((display %in% "rgl") && (require(rgl) & require(rpanel))) {
            	ylim <- range(model$y, aperm(est * mask, opt$order[1:2]), na.rm = TRUE)
            	yy <- model$y - model.fitted + mui
            	if (exists("reference"))
            		 clr <- col.pal[cut(c(reference), se.breaks, labels = FALSE)]
            	opt$scaling <- rp.plot3d(x[, 1], yy, x[, 2],
            					xlab = xlab[opt$order[1]], ylab = ylab, zlab = xlab[opt$order[2]],
            					xlim = opt$x1lim, ylim = ylim, zlim = opt$x2lim,
            					size = 0.5, col = "black")
            	surf.ids <- sm.surface3d(cbind(u[[opt$order[1]]], u[[opt$order[2]]]), 
            													 aperm(est * mask, opt$order[1:2]), opt$scaling, 
            													 col = clr, col.mesh = opt$col.mesh, 
            													 alpha = 0.7, alpha.mesh = 1, lit = FALSE)
            }
            else {
               layout(matrix(c(1, 2), ncol = 2), widths = c(7, 1))
               par(mar = c(3, 3, 1, 0.5) + 0.1, mgp = c(1.5, 0.2, 0), tcl = -0.2)
               if (opt$eqscplot)
                  MASS::eqscplot(x[ , opt$order[1:2]], type = "n", xlab = xlab[1], ylab = xlab[2])
               else {
                  plot(x[ , opt$order[1]], x[ , opt$order[2]], type = "n", axes = FALSE,
                       xlab = xlab[1], ylab = xlab[2])
                  usr <- par("usr")
                  rect(usr[1], usr[3], usr[2], usr[4], col = grey(0.9), border = NA)
                  grid(col = "white", lty = 1)
                  axis(1, lwd = 0, lwd.ticks = 2, col = grey(0.6), col.ticks = grey(0.6),
                          col.axis = grey(0.6), cex.axis = 0.8)
                  axis(2, lwd = 0, lwd.ticks = 2, col = grey(0.6), col.ticks = grey(0.6),
                          col.axis = grey(0.6), cex.axis = 0.8)
               }
               # brks[is.infinite(brks) & (brks > 0)] <- max(y, model$y, na.rm = TRUE) + 1
               # brks[is.infinite(brks) & (brks < 0)] <- min(y, model$y, na.rm = TRUE) - 1

               dfrm <- data.frame(x = rep(u[[opt$order[1]]], length(u[[opt$order[2]]])),
   	                              y = rep(u[[opt$order[2]]], each = length(u[[opt$order[1]]])),
   	                              z = c(t(aperm(est * mask, opt$order[1:2]))))
               ind  <- !is.na(dfrm$x + dfrm$y + dfrm$z)
   	           dfrm <- subset(dfrm, ind)
    	       if (requireNamespace("akima", quietly = TRUE)) {
  	           grdx  <- seq(min(dfrm$x), max(dfrm$x), length = 200)
  	           grdy  <- seq(min(dfrm$y), max(dfrm$y), length = 200)
  	           intrp <- akima::interp(dfrm$x, dfrm$y, dfrm$z, grdx, grdy)
  	           dfrmx <- dfrm$x
  	           dfrmy <- dfrm$y
  	           dfrm  <- data.frame(x = rep(intrp$x, length(intrp$y)),
  	                               y = rep(intrp$y, each = length(intrp$x)),
  	                               z = c(intrp$z))
  	           if (exists("reference"))
  	             intrpr <- akima::interp(dfrmx, dfrmy, c(reference)[ind], grdx, grdy)
  	           dfrm <- list(x = grdx, y = grdy, z = intrp$z)
  	           if (exists("reference")) dfrm$r <- intrpr$z
    	       }
   	         else {
   	            dfrm   <- as.list(dfrm)
   	            dfrm$z <- matrix(dfrm$z, nrow = dim(model$y)[1])
   	         }
   	           
               image(intrp$x, intrp$y, t(intrp$z), zlim = ylim, col = opt$col.palette, add = TRUE,
                     xlab = xlab[opt$order[1]], ylab = xlab[opt$order[2]])
               if (exists("reference")) {
                 lvls <- pretty(c(2, max(c(2, dfrm$r), na.rm = TRUE)))
                 mmx <- max(c(2, dfrm$r), na.rm = TRUE)
                 if (mmx >= 2) {
                   lvls <- if (trunc(mmx) > 5) pretty(c(2, trunc(mmx))) else 2:trunc(mmx)
                   # contour(mx[ , 1], mx[ , 2], mr, add = TRUE, col = "blue", levels = lvls, lty = 1)
                   contour(dfrm$x, dfrm$y, dfrm$r, add = TRUE, col = "blue", levels = lvls, lty = 1)
                 }
                 lvls <- pretty(c(-2, min(c(-2, dfrm$r), na.rm = TRUE)))
                 mmn <- min(c(-2, dfrm$r), na.rm = TRUE)
                 if (mmn <= -2) {
                   lvls <- if (trunc(mmn) < -5) pretty(c(-2, trunc(mmn))) else (-2):trunc(mmn)
                   # contour(mx[ , 1], mx[ , 2], mr, add = TRUE, col = "blue", levels = lvls, lty = 2)
                   contour(dfrm$x, dfrm$y, dfrm$r, add = TRUE, col = "blue", levels = lvls, lty = 2)
                 }
               }

               if (is.function(opt$foreground.plot)) opt$foreground.plot()

               del  <- 0.04 * diff(ylim)
               brks <- seq(ylim[1] - del, ylim[2] + del, length = length(opt$col.palette) + 1)
               rp.colour.key(opt$col.palette, brks, par.mar = c(3, 0, 1, 2.5) + 0.1)
               mtext(ylab, side = 4, line = 1.1, font = 1)
               layout(1)

            }
         }
         
         else if (ndim == 3) {

         	  mask  <- sm.mask(x[ , opt$order[1:2]], cbind(u[[opt$order[1]]], u[[opt$order[2]]]),
                             mask.method = opt$mask.method)
            mdl   <- list(x = cbind(u[[opt$order[1]]], u[[opt$order[2]]]), z = u[[opt$order[3]]],
                          y = aperm(est, opt$order))
            mdl$y <- array(c(mdl$y) * c(mask), dim = dim(mdl$y))
            
            if (opt$reference == "no effect") mdl$reference <- array(c(reference) * c(mask), dim = dim(reference))
            fg.plot <- if (is.function(opt$foreground.plot)) opt$foreground.plot else NULL
            if (length(xlab) < 3) {
               xl <- character(0)
               for (j in 1:length(xdims)) {
                  if (xdims[j] > 1) {
                  	 cs  <- cumsum(xdims[1:j])
                  	 ind <- (cs - xdims[j] + 1):cs
                     xln <- if (!is.null(colnames(x[ , ind]))) colnames(x[ , ind])
                            else paste(xlab, 1:xdims[j], sep = "-")
                     xl <- c(xl, xln)
                  }
                  else
                     xl <- c(xl, xlab[j])
               }
               xlab <- xl
            }
            ylb <- ylab
            if (opt$deriv.order == 1) ylb <- paste(ylb, "(1st deriv.)")
            if (opt$deriv.order == 2) ylb <- paste(ylb, "(2nd deriv.)")
            rp.plot4d(x[ , opt$order[1:2]], x[ , opt$order[3]], mui, mdl,
                  ylab = ylb,
                  zlab = xlab[opt$order[3]], x1lab = xlab[opt$order[1]], x2lab = xlab[opt$order[2]],
                  col.palette = opt$col.palette, hscale = opt$hscale, vscale = opt$vscale,
                  foreground.plot = fg.plot, display = display)
            if (length(components) == 1 & ndim == 3) assign("mdl", mdl, envir = env)
         }
      })
      }
   }
   
   if (length(components) == 1 & ndim == 3) {
   	  plotinfo <- plotinfo[[1]]
      plotinfo$model <- mdl
   }
   invisible(plotinfo)
}

sm.pam.colour.chart <- function(panel) {
  par(mar = c(5, 1, 4, 2) + 0.1)
  rp.colour.chart(panel$col.palette, panel$ylim)
  panel
  }

rp.colour.chart <- function(cols, zlim)  {
   ngrid <- length(cols)
   plot(0:1, zlim, type = "n", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
   axis(4)
   xvec <- rep(0, ngrid)
   yvec <- seq(zlim[1], zlim[2], length = ngrid + 1)
   rect(xvec, yvec[-length(yvec)], xvec + 1, yvec[-1], col = cols, border = NA)
   box()
   }

#----------------------------------------------------------------------------
fitted.pam <- function(model) model$fitted

residuals.pam <- function(model) model$y - model$fitted

#----------------------------------------------------------------------------
sm.mask <- function(x, eval.points, mask.method = "hull") {
	
   ngrid       <- nrow(eval.points)
   grid.points <- cbind(rep(eval.points[, 1], ngrid), 
                           rep(eval.points[, 2], rep(ngrid, ngrid)))
                           
   if (mask.method == "hull") {
      hull.points <- as.matrix(x[order(x[, 1], x[, 2]), ])
      dh          <- diff(hull.points)
      hull.points <- hull.points[c(TRUE, !((dh[, 1] == 0) & (dh[, 2] == 0))), ]
      hull.points <- hull.points[chull(hull.points), ]
      nh          <- nrow(hull.points)
      gstep       <- matrix(rep(eval.points[2, ] - eval.points[1, ], nh),
                        ncol = 2, byrow = TRUE)
      hp.start    <- matrix(rep(eval.points[1, ], nh), ncol = 2, byrow = TRUE)
      hull.points <- hp.start + gstep * round((hull.points - hp.start)/gstep)
      hull.points <- hull.points[chull(hull.points), ]
      D           <- diff(rbind(hull.points, hull.points[1, ]))
      temp        <- D[, 1]
      D[, 1]      <- D[, 2]
      D[, 2]      <- (-temp)
      C           <- as.vector((hull.points * D) %*% rep(1, 2))
      C           <- matrix(rep(C, ngrid^2), nrow = ngrid^2, byrow = TRUE)
      D           <- t(D)
      wy          <- ((grid.points %*% D) >= C)
      wy          <- apply(wy, 1, all)
      wy[wy]      <- 1
      wy[!wy]     <- NA
      mask        <- matrix(wy, ncol = ngrid)
   }
   else if (mask.method == "near") {
   	  del1  <- eval.points[2, 1] - eval.points[1, 1]
   	  del2  <- eval.points[2, 2] - eval.points[1, 2]
      mask  <- apply(grid.points, 1,
                   function(z) any(((z[1] - x[,1])/del1)^2 + ((z[2] - x[,2])/del2)^2 < 4^2))
      mask  <- matrix(as.numeric(mask), ncol = ngrid)
      mask[mask == 0] <- NA
   }
   else
      mask <- matrix(1, ncol = ngrid, nrow = ngrid)
      
   return(invisible(mask))
}

#----------------------------------------------------------------------------
s <- function(..., lambda = NA, df = NA, period = NA, increasing = FALSE,
							     xrange = NA, nseg = NA,
                   fixed = c(NA, NA)) {
   vars.list <- as.list(substitute(list(...)))[-1]
   nvar <- length(vars.list)
   if (nvar > 3)
      stop("smooth terms can be constructed from only 1, 2 or 3 variables.")
   variables <- character(0)
   for (i in 1:nvar) variables <- c(variables, deparse(vars.list[[i]]))
   list(variables = variables, lambda = lambda, df = df, period = period,
        increasing = increasing, xrange = xrange, nseg = nseg, fixed = fixed)
}

#----------------------------------------------------------------------------

ps.matrices <- function(x, xrange, ndims, nseg, bdeg = 3, pord = 2, period = NA,
                            decompose =  TRUE, penalty = TRUE) {

    # Compute a set of basis functions and a penalty matrix associated with x.
    # An intercept term and the main effect of any interaction terms are removed.
    
    ndimx <- ncol(x)
    if (ndimx > 3) stop("terms with more than three dimensions cannot be used.")
    n    <- nrow(x)
    
    if (missing(nseg)) nseg <- rep(switch(ndimx, 100, 17, 7), ndimx)
    if (length(bdeg) < ndimx) bdeg <- rep(bdeg[1], ndimx)
    
    # Compute B-spline basis
    
    b <- list(length = ndimx)
    m <- vector(length = ndimx)
    for (i in 1:ndimx) {
       b[[i]] <- bbase(x[,i], xl = xrange[i , 1], xr = xrange[i, 2], nseg = nseg[i], 
                       deg = bdeg[i])
       m[i]   <- ncol(b[[i]])
    }

    B <- b[[1]]
    if (ndimx > 1)
       B <- t(apply(cbind(b[[1]], b[[2]]), 1,
                            function(x) c(x[1:m[1]] %x% x[-(1:m[1])])))
    if (ndimx == 3)
       B <- t(apply(cbind(B,  b[[3]]), 1, 
                function(x) c(x[1:(m[1]*m[2])] %x% x[-(1:(m[1]*m[2]))])))
    
    result <- list(B = B, xrange = xrange, nseg = nseg, bdeg = bdeg, pord = pord)
    if (!penalty) return(invisible(result))

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
        
    result$P <- pmat
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

