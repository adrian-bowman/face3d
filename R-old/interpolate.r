interpolate.face3d <- function(shape, y, id) {
    if (!require(MASS)) stop("the MASS package is required.")
    coords     <- shape$coords
    k1         <- length(id[id == TRUE])
    k2         <- length(id[id == FALSE])
    Sdx        <- outer(coords[id == TRUE, 1], coords[id == TRUE, 1], "-")
    Sdy        <- outer(coords[id == TRUE, 2], coords[id == TRUE, 2], "-")
    Sdz        <- outer(coords[id == TRUE, 3], coords[id == TRUE, 3], "-")
    Sphi       <- sqrt(Sdx * Sdx + Sdy * Sdy + Sdz * Sdz)
    ONE        <- matrix(1, k1, 1)
    ZERO       <- matrix(0, 4, 4)
    Q          <- rbind(Sphi, t(ONE), t(coords[id == TRUE, ]))
    O          <- cbind(ONE, coords[id == TRUE, ])
    O          <- as.matrix(O)
    U          <- rbind(O, ZERO)
    GAMMA      <- cbind(Q, U)
    Gi         <- ginv(GAMMA)
    tps.coeffs <- Gi %*% matrix(c(y, rep(0, 4)), k1 + 4, 1)
    # affine part
    affine     <- rep(1, k2) %o% tps.coeffs[k1 + 1, ] +
                  coords[id == FALSE, 1] %o% tps.coeffs[k1 + 2, ] +
                  coords[id == FALSE, 2] %o% tps.coeffs[k1 + 3, ] +
                  coords[id == FALSE, 3] %o% tps.coeffs[k1 + 4, ]
    # nonaffine part
    dx         <- outer(coords[id == TRUE, 1], coords[id == FALSE, 1], "-")
    dy         <- outer(coords[id == TRUE, 2], coords[id == FALSE, 2], "-")
    dz         <- outer(coords[id == TRUE, 3], coords[id == FALSE, 3], "-")
    phi        <- sqrt(dx * dx + dy * dy + dz * dz)
    nonaff     <- t(phi) %*% tps.coeffs[1:k1, ]
    TPS        <- c(affine + nonaff)
    invisible(TPS)
}
