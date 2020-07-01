"gcurvature.face3d" <- function(curve, df, n = 200) {

    arc.length         <- cumsum(c(0, sqrt(rowSums(diff(curve)^2))))
    x.coords           <- approx(arc.length, curve[ , 1], n = n)$y
    y.coords           <- approx(arc.length, curve[ , 2], n = n)$y
    z.coords           <- approx(arc.length, curve[ , 3], n = n)$y
    resampled.curve    <- cbind(x.coords, y.coords, z.coords)
    arc.length         <- cumsum(c(0, sqrt(rowSums(diff(resampled.curve)^2))))
 
    # model              <- sm.psp(arc.length, x.coords, df = df, display = "none")
    model              <- smooth.spline(arc.length, x.coords, df = df)
    x0                 <- predict(model, arc.length, deriv = 0)$y
    x1                 <- predict(model, arc.length, deriv = 1)$y
    x2                 <- predict(model, arc.length, deriv = 2)$y
    # model              <- sm.psp(arc.length, y.coords, df = df, display = "none")
    model              <- smooth.spline(arc.length, y.coords, df = df)
    y0                 <- predict(model, arc.length, deriv = 0)$y
    y1                 <- predict(model, arc.length, deriv = 1)$y
    y2                 <- predict(model, arc.length, deriv = 2)$y
    # model              <- sm.psp(arc.length, y.coords, df = df, display = "none")
    model              <- smooth.spline(arc.length, z.coords, df = df)
    z0                 <- predict(model, arc.length, deriv = 0)$y
    z1                 <- predict(model, arc.length, deriv = 1)$y
    z2                 <- predict(model, arc.length, deriv = 2)$y
    gcurvature         <- sqrt((x1*y2 - x2*y1)^2 + (x1*z2 - x2*z1)^2 +
                               (y1*z2 - y2*z1)^2) / (x1^2 + y1^2 + z1^2)^1.5

    # xyz0               <- cbind(x0, y0, z0)
    # ind                <- 1 + which(diff(-sign(diff(gcurvature))) == 2)
    # ind.min            <- which.min(xyz0[ind, 1])
    # ind.max            <- which.max(xyz0[ind, 1])
    # ind                <- c(ind[ind.min],ind[ind.max])

    invisible(list(resampled.curve = resampled.curve, gcurvature = gcurvature, arc.length = arc.length))
}
