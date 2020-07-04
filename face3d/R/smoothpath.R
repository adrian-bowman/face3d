smoothpath <- function(shape, x1, x2, direction, bothways = FALSE, ppath,
                       boundary = c(0.2, 0.5), si.target, distance = 5,
                       rotation.range = pi/2, monitor = 0) {

  if (missing(ppath))
    ppath <- planepath.face3d(shape, x1, x2, boundary = boundary, si.target = si.target,
                              distance = distance, monitor = monitor)
  sbst  <- ppath$shape
  ppath <- ppath$path
  
  if (monitor > 1) {
    plot(sbst)
    spheres3d(ppath)
  }
    
  values  <- pmax(-sbst$kappa2, 0) #ridge
  ccurve  <- closestcurve.face3d(shape.smooth, curve)
  lambda  <- cumsum(c(0, sqrt(apply((diff(curve))^2, 1, sum))))
  xcoord  <- lambda[ccurve$closest.curvept] # distance along the curve
  ycoord  <- ccurve$closest.distance
  si      <- shape.smooth$shape.index
  values[sign(si)!= sign(si.target)] <- 0 
  ind      <- (abs(ycoord) < perp.dist.bound) & (xcoord > 0.005 * max(xcoord)) & (xcoord < 0.995 * max(xcoord))
  ycoord   <- c(0, ycoord[ind], 0)
  xcoord   <- c(0, xcoord[ind], max(xcoord))
  values   <- c(0, values[ind], 0)
  ridge           <- ridge2d.face3d(xcoord, ycoord, values, lambda = penalty, monitor = monitor)
  shape.smooth    <- subset.face3d(shape.smooth, ind, remove.singles=FALSE)
  arclength       <- ridge$x 
  new.pdist       <- ridge$y

# Project back onto the surface full pts (interp bary)
  full.ccurve   <- closestcurve.face3d(shape, curve)
  lambda        <- cumsum(c(0, sqrt(apply((diff(curve))^2, 1, sum))))
  full.xc       <- lambda[full.ccurve$closest.curvept] # distance along the curve
  full.xc   <- c(0, full.xc, max(lambda))
  full.yc       <- full.ccurve$closest.distance
  full.yc   <- c(0, full.yc, 0)
  full.va1      <- c(0, pmax(abs(shape$kappa1), abs(shape$kappa2)),0)

  X             <- cbind(full.xc, full.yc)
  # X             <- cbind(xc,yc)
  Xnew          <- cbind(arclength, new.pdist)
  interp.x      <- interp.barycentric(X, f = c(lmk1[1], shape$vertices[ , 1], lmk2[1]), Xnew)$fnew
  interp.y      <- interp.barycentric(X, f = c(lmk1[2], shape$vertices[ , 2], lmk2[2]), Xnew)$fnew
  interp.z      <- interp.barycentric(X, f = c(lmk1[3], shape$vertices[ , 3], lmk2[3]), Xnew)$fnew
  interp.values <- interp.barycentric(X, f = full.va1, Xnew)$fnew
  smooth.path   <- cbind(interp.x, interp.y, interp.z)
  aver.path     <- smooth.path
}
