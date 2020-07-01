pcurve.face3d <- function (shape, df = rep(5, ncol(shape$coords)),
                        weights = rep(1, nrow(shape$coords)),
                        fixed = rep(NA, 2), increasing = rep(FALSE, 2),
                        initial = NULL, start = 1, finish = 2,
                        thresh = 0.001, maxit = 10, stretch = 2, ...) {

    if (class(shape) == "face3d") x <- shape$coords
    else                          x <- shape
    if (all(!is.na(fixed))) {
       x <- rbind(x, fixed)
       weights <- c(weights, rep(0, nrow(fixed)))
    }
    n         <- nrow(x)
    p         <- ncol(x)
    dist.old  <- sum(diag(var(x)))

    if (is.null(initial))
       initial <- planepath.face3d(shape, fixed[start, ], fixed[finish, ])$path
    pcurve <- get.lam(x, initial, stretch = stretch)
    # pcurve <- project_to_curve(x, initial, stretch = stretch)

    sm.reg <- function(x, y, df, weights, fixed, increasing) {
       model <- smooth.face3d(x, y, df = df, fixed = fixed, increasing = increasing, weights = weights,
                   eval.points = x, display = "none")
       model$estimate
    }

    it           <- 0
    s            <- matrix(NA, nrow = n, ncol = p)
    hasConverged <- (abs((dist.old - pcurve$dist)/dist.old) <= thresh)
    while (!hasConverged && it < maxit) {
        it           <- it + 1
        for (j in 1:p)
             s[ , j] <- sm.reg(pcurve$lambda, x[ , j], df = df[j], increasing = FALSE,
                               weights = weights, fixed = matrix(NA, ncol = p))
        dist.old     <- pcurve$dist
        pcurve       <- get.lam(x, s = s, pcurve$tag, stretch = stretch)
        if (all(!is.na(fixed))) {
           if (is.vector(fixed)) fixed <- matrix(fixed, nrow = 1)
           lam.fxd   <- pcurve$lambda[(n + 1 - nrow(fixed)):n]
            for (j in 1:p)
             s[ , j] <- sm.reg(pcurve$lambda, x[ , j], df = df[j], increasing = increasing[j],
                               weights = weights, fixed = cbind(lam.fxd, fixed[ , j]))
           dist.old  <- pcurve$dist
           pcurve    <- get.lam(x, s = s, pcurve$tag, stretch = stretch)
        }
        hasConverged <- (abs((dist.old - pcurve$dist)/dist.old) <= thresh)
    }

    structure(list(s = pcurve$s, tag = pcurve$tag, lambda = pcurve$lambda,
        dist = pcurve$dist))
}
