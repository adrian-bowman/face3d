smoothpath.face3d <- function(shape, x1, x2, direction,   penalty = 0.002,
                              bothways = FALSE, ppath,
                              boundary = c(0.2, 0.5), si.target, distance = 5,
                              rotation.range = pi/2, monitor = 0) {

  if (missing(ppath))
    ppath <- planepath(shape, x1, x2, boundary = boundary, si.target = si.target,
                              distance = distance, monitor = monitor)
  sbst  <- ppath$shape
  ppath <- ppath$path
  
  perp.dist.bound <- 10
  
  if (monitor > 1) {
    plot(sbst)
    spheres3d(ppath)
  }
  
  values  <- if (si.target > 0) pmax(-sbst$kappa2, 0) else pmax(sbst$kappa1, 0)
  ccurve  <- closestcurve.face3d(sbst, ppath)
  lambda  <- cumsum(c(0, sqrt(apply((diff(ppath))^2, 1, sum))))
  xcoord  <- lambda[ccurve$closest.curvept] # distance along the curve
  ycoord  <- ccurve$closest.distance
  si      <- sbst$shape.index
  values  <- values * as.numeric(sign(si) == sign(si.target))
  ind     <- (abs(ycoord) < perp.dist.bound) & (xcoord > 0) & (xcoord < max(xcoord))
  ycoord  <- c(0, ycoord[ind], 0)
  xcoord  <- c(0, xcoord[ind], max(xcoord))
  values  <- c(0, values[ind], 0)
  ridge   <- ridge2d.face3d(xcoord, ycoord, values, lambda = penalty, monitor = monitor)

  plot(xcoord, ycoord, col = topo.colors(20)[cut(values, 20, labels = FALSE)])
  points(ridge$x, ridge$y, col = "red", pch = 16)
       
# Project back onto the surface full pts (interpbary)
  sbst     <- subset.face3d(sbst, ind, remove.singles = FALSE)
  X        <- cbind(xcoord, ycoord)[-c(1, length(xcoord)), ]
  ind      <- (ridge$x > min(X[ , 1])) & (ridge$x < max(X[ , 1]))
  Xnew     <- cbind(ridge$x, ridge$y)[ind, ]
  trngs    <- sbst$triangles
  interp.x <- interpbary(X, sbst$vertices[ , 1], trngs, Xnew)$fnew
  interp.y <- interpbary(X, sbst$vertices[ , 2], trngs, Xnew)$fnew
  interp.z <- interpbary(X, sbst$vertices[ , 3], trngs, Xnew)$fnew

  return(invisible(cbind(interp.x, interp.y, interp.z)))
}
