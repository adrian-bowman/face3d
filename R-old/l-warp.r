"warp.face3d" <- function(template.pts, face.pts, subset, project = TRUE, face) {
	
	# template: lmks, curves or meshes from the template
	# face:     lmks, curves or meshes from the face
	# subset:   indices of the section required
	
    if (!require(MASS)) stop("the MASS package is required.")
    
    sbst       <- subset
    if (any(is.logical(sbst))) sbst <- which(sbst)
    
    shape      <- template.pts[subset, ]
    lmks1      <- template.pts[-subset, ]
    lmks2      <- face.pts[-subset, ]
    
    k1         <- dim(lmks1)[1]
    k2         <- dim(shape)[1]
    Sdx        <- outer(lmks1[,1],lmks1[,1], "-")
    Sdy        <- outer(lmks1[,2],lmks1[,2], "-")
    Sdz        <- outer(lmks1[,3],lmks1[,3], "-")
    Sphi       <- sqrt(Sdx * Sdx + Sdy * Sdy + Sdz * Sdz)
    ONE        <- matrix(1, k1, 1)
    ZERO       <- matrix(0, 4, 4)
    Q          <- rbind(Sphi, t(ONE), t(lmks1))
    O          <- cbind(ONE, lmks1)
    O          <- as.matrix(O)
    U          <- rbind(O, ZERO)
    GAMMA      <- cbind(Q, U)
    Gi         <- ginv(GAMMA)
    tps.coeffs <- Gi %*% rbind(lmks2, matrix(0, 4, 3))
    # affine part
    affine     <- rep(1, k2) %o% tps.coeffs[k1 + 1,] +
                  shape[,1]  %o% tps.coeffs[k1 + 2,] +
                  shape[,2]  %o% tps.coeffs[k1 + 3,] +
                  shape[,3]  %o% tps.coeffs[k1 + 4,]
    # nonaffine part
    dx         <- outer(lmks1[ , 1], shape[ , 1], "-")
    dy         <- outer(lmks1[ , 2], shape[ , 2], "-")
    dz         <- outer(lmks1[ , 3], shape[ , 3], "-")
    phi        <- sqrt(dx * dx + dy * dy + dz * dz)
    nonaff     <- t(phi) %*% tps.coeffs[1:k1, ]
    TPS        <- affine + nonaff
    
    if (project) TPS <- closest.face3d(TPS, face)$point
    
    invisible(TPS)
}
