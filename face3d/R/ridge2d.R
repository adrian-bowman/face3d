 drv <- function(alpha, model, deriv, lambda) {
	if (deriv > 2) stop("deriv must be 1 or 2.")
	Bs     <- ps.matrices(matrix(model$eval.points[[1]]),
												xrange = matrix(range(model$eval.points[[1]]), nrow = 1), ndims = 1,
												nseg = model$nseg, bdeg = model$bdeg)
	Bd     <- ps.matrices(matrix(alpha), xrange = matrix(range(model$eval.points[[2]]), nrow = 1), ndims = 1,
												nseg = model$nseg, bdeg = model$bdeg - deriv)
	P      <- diff(diag(length(alpha)), diff = 1)
	P      <- crossprod(P)
	beta   <- matrix(model$beta, ncol = 20)
   alpha1 <- diff(beta, differences = deriv)
   h      <- diff(range(model$eval.points[[2]]))
   derv   <- Bd$B %*% alpha1 / (h / Bd$nseg)^deriv
   derv   <- rowSums(derv * Bs$B) / length(alpha)
   result <- if (deriv == 1) derv - 2 * lambda * P %*% alpha else diag(derv) - 2 * lambda * P
   result
}

ridge2d.face3d <- function(s, d, crv, lambda = 1, df = 12, ngrid = 50, endpoints.fixed = rep(TRUE, 2), monitor = FALSE) {
   
   if (length(endpoints.fixed) == 1) endpoints.fixed <- rep(endpoints.fixed, 2)
	evp    <- list(s = seq(min(s), max(s), length = ngrid),
                  d = seq(min(d), max(d), length = 100))
	x      <- cbind(s, d)
	model  <- smooth.face3d(x, crv, df = df, eval.points = evp)
	if (monitor)
	   image(model$eval.points[[1]], model$eval.point[[2]],
            matrix(model$estimate, nrow = length(model$eval.points[[1]])),
            col = topo.colors(20), xlab = "Arc Length", ylab = "Perpendicular Distance")

   # Different strategies for initialising alpha
   # 	mest <- matrix(model$estimate, nrow = length(model$eval.points[[1]]))
   #    alpha  <- model$eval.points[[2]][apply(mest, 1, which.max)]
	alpha  <- matrix(rep(0, ngrid))
   # points(model$eval.points[[1]], alpha, col = "red", pch = 16)

   change <- 1
   i      <- 0
   while(change > 0.001 & i < 50) { 
      # cat(change, "")
      i      <- i + 1
      alpha0 <- alpha
      drv1   <- drv(alpha, model, 1, lambda)
      drv2   <- drv(alpha, model, 2, lambda)
      if (any(endpoints.fixed)) {
         ind    <- c(1, length(alpha))[endpoints.fixed]
         drv1   <- drv1[-ind]
         drv2   <- drv2[-ind, -ind]
         alpha  <- alpha[-ind]
         alpha0 <- alpha0[-ind]
      }
      alpha  <- alpha - solve(drv2) %*% drv1
      change <- max(abs(alpha0 - alpha) / abs(alpha0))
      if (endpoints.fixed[1]) alpha <- c(0, alpha)
      if (endpoints.fixed[2]) alpha <- c(alpha, 0)
      if (monitor) {
         par(mar = c(3, 3, 1, 1) + 0.1, mgp = c(1.5, 0.2, 0), tcl = -0.2)
         image(model$eval.points[[1]], model$eval.point[[2]],
               matrix(model$estimate, nrow = length(model$eval.points[[1]])),
               col = topo.colors(20), xlab = "Arc Length", ylab = "Perpendicular Distance")
         points(model$eval.points[[1]], alpha, col = "red", pch = 16)
      }
   }
   if (i == 50) warning("iterations have not converged.")
   
   invisible(list(x = model$eval.points[[1]], y = alpha, n.iter = i))
}
