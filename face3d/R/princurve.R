princurve.face3d <- function(x, df = 12, maxiter = 30, periodic = FALSE,
                             monitor = FALSE, threshold = 1e-3) {
   
   if (!requireNamespace("princurve", quietly = TRUE))
      stop("the princurve package is required.")
   
   if (periodic) x <- rbind(x, x[1, ])
   n         <- nrow(x)
   p         <- ncol(x)
   dist.old  <- sum(diag(var(x)))
   pcrv      <- princurve::project_to_curve(x, x, stretch = 0)
   period    <- if (periodic) max(pcrv$lambda) else NA
   
   iter      <- 0
   s         <- matrix(NA, nrow = n, ncol = p)
   converged <- (abs((dist.old - pcrv$dist)/dist.old) <= threshold)
   while (!converged && iter < maxiter) {
      iter <- iter + 1
      for (j in 1:p)
         s[ , j] <- smooth.face3d(pcrv$lambda, x[ , j], df = df,
                           period = period,
                           eval.points = pcrv$lambda, display = "none")$estimate
      dist.old <- pcrv$dist
      pcrv     <- princurve::project_to_curve(x, s, stretch = 0)
      period   <- max(pcrv$lambda)
      cnvg     <- abs(dist.old - pcrv$dist) / dist.old
      if (monitor) cat(iter, cnvg, "\n")
      converged <- (cnvg <= threshold)
   }
   if (iter == maxiter) cat("maximum number of iterations reached.\n")
   
   estimate  <- pcrv$s[pcrv$ord, ]
   arclength <- pcrv$lambda[pcrv$ord]
   if (periodic) {
      estimate  <- rbind(estimate, pcrv$s[pcrv$ord[1], ])
      arclength <- c(arclength, period)
   }
   
   result <- list(fitted           = pcrv$s,
                  fitted.arclength = pcrv$lambda,
                  ord              = pcrv$ord,
                  dist             = pcrv$dist,
                  estimate         = estimate,
                  arclength        = arclength)
   if (periodic) result$period <- period
   
   invisible(result)
}
