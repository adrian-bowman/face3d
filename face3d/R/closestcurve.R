closestcurve.face3d <- function(shape, curve) {
   if (!("normals" %in% names(shape))) shape <- normals.face3d(shape)
   Distance         <- rdist(curve, shape$coords)
   closest.distance <- apply(Distance, 2, min)
   closest.curvept  <- apply(Distance, 2, which.min)
   closest.meshpt   <- apply(Distance, 1, which.min)
   normals          <- shape$normal[closest.meshpt,]
   ncurve           <- nrow(curve)
   curvevec         <- curve[pmin(closest.curvept + 1, ncurve), ] -  
                       curve[pmax(closest.curvept - 1, 1), ]
   meshvec          <- shape$coords - curve[closest.curvept, ]
   axis1            <- crossproduct(normals[closest.curvept, ], curvevec)
   ind              <- apply(meshvec * axis1, 1, sum)
   closest.distance <- closest.distance * sign(ind)
   invisible(list(closest.distance = closest.distance,
                  closest.curvept  = closest.curvept))
}
