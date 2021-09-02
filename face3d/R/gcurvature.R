gcurvature <- function(curve, df, direction, n = 200, monitor = FALSE) {

    if (missing(direction)) direction <- ""
    if (missing(df)) stop("df must be specified.")
    
    arc.length         <- cumsum(c(0, sqrt(rowSums(diff(curve)^2))))
    x.coords           <- approx(arc.length, curve[ , 1], n = n)$y
    y.coords           <- approx(arc.length, curve[ , 2], n = n)$y
    z.coords           <- approx(arc.length, curve[ , 3], n = n)$y
    resampled.curve    <- cbind(x.coords, y.coords, z.coords)
    arc.length         <- cumsum(c(0, sqrt(rowSums(diff(resampled.curve)^2))))
 
    # model              <- smooth.face3d(arc.length, x.coords, df = df, display = "none")
    model              <- smooth.spline(arc.length, x.coords, df = df)
    x0                 <- predict(model, arc.length, deriv = 0)$y
    x1                 <- predict(model, arc.length, deriv = 1)$y
    x2                 <- predict(model, arc.length, deriv = 2)$y
    # model              <- smooth.face3d(arc.length, y.coords, df = df, display = "none")
    model              <- smooth.spline(arc.length, y.coords, df = df)
    y0                 <- predict(model, arc.length, deriv = 0)$y
    y1                 <- predict(model, arc.length, deriv = 1)$y
    y2                 <- predict(model, arc.length, deriv = 2)$y
    # model              <- smooth.face3d(arc.length, y.coords, df = df, display = "none")
    model              <- smooth.spline(arc.length, z.coords, df = df)
    z0                 <- predict(model, arc.length, deriv = 0)$y
    z1                 <- predict(model, arc.length, deriv = 1)$y
    z2                 <- predict(model, arc.length, deriv = 2)$y
    gcurvature         <- sqrt((x1*y2 - x2*y1)^2 + (x1*z2 - x2*z1)^2 +
                               (y1*z2 - y2*z1)^2) / (x1^2 + y1^2 + z1^2)^1.5
    gcrv         <- gcurvature * switch(direction, 1, x = x2, y = y2, z = z2)
    ind.localmax <- which(c(0, diff(sign(diff(gcrv))), 0) < 0)
    pos.localmax <- resampled.curve[ind.localmax, ]
    ind.max      <- which.max(gcrv)
    pos.max      <- resampled.curve[ind.max, ]
    
    if (monitor) {
       plot(gcrv ~ arc.length)
       points(arc.length[ind.localmax], gcrv[ind.localmax], pch = 16, col = "blue")
       points(arc.length[ind.max], gcrv[ind.max], pch = 16, col = "red")
    }

    invisible(list(gcurvature = gcurvature, arclength = arc.length,
                   pos.localmax = pos.localmax, ind.localmax = ind.localmax,
                   pos.max = pos.max, ind.max = ind.max,
                   resampled.curve = resampled.curve, direction,
                   d2.x = x2, d2.y = y2, d2.z = z2))
}
