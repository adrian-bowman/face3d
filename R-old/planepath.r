planepath.face3d <- function(shape, x1, x2, pts1, pts2, direction, bothways = TRUE,
                         distance = 5, ngrid = 50, boundary = c(0.2, 0.5),
                         rotation = "optimise", rotation.range = pi/2,
                         bridge.gaps = FALSE, si.target, directions = FALSE, graphics = FALSE) {

   if (missing(x1)) stop("x1 must be specified.")
   if (missing(pts1))      pts1      <- NA
   if (missing(pts2))      pts2      <- NA
   if (missing(x2))        x2        <- NA
   if (missing(direction)) direction <- NA
   if (any(is.na(x2)) & any(is.na(direction))) stop("one of x2 or direction must be specified.")
   si.target.present <- !missing(si.target)
   
   # This is an internal variable which, when switched on, shows the intermediate steps of the algorithm
   monitor <- TRUE
   
   if (!any(is.na(x2)) && sqrt(sum((x1 - x2)^2)) < .Machine$double.eps)
      return(invisible(list(path = c(x1, x2), arclength = c(0, 0))))

   #  Subset the object by a tube around x1 and the relevant direction
   
   if (all(!is.na(boundary))) {
      if (any(is.na(x2))) {
      	  rng      <- 0
      	  unit     <- direction / sqrt(sum(direction^2))
      }
      else {
      	  rng      <- sqrt(sum((x2 - x1)^2))
      	  boundary <- boundary * rng
     	  unit     <- (x2 - x1) / rng
      }
      prjn  <- c(sweep(shape$coords, 2, x1) %*% unit)
      near  <- function(x, pts, distance)
                      sqrt(.rowSums((sweep(pts, 2, x))^2, nrow(pts), 3)) < distance
      if (any(is.na(x2)))
      	  ind1 <- near(x1, shape$coords, boundary[1])
      else {
         ind1 <- (prjn > - boundary[1]) & (prjn < rng + boundary[1])
         prjn  <- outer(prjn, unit)
         ind1  <- ind1 & apply((sweep(shape$coords, 2, x1) - prjn)^2, 1,
                            function(x) sqrt(sum(x)) < boundary[2])
      }
      shp   <- subset.face3d(shape, ind1)
      shape <- shp
   }
   
   # Calculate the surface curvatures
   if (si.target.present | directions) {
   	  if (!all(c("shape.index", "kappa1", "kappa2") %in% names(shape)))
         shape <- index.face3d(shape, distance = distance, directions = directions)
   }
   if (si.target.present) {
      si       <- shape$shape.index
      k1       <- shape$kappa1
      k2       <- shape$kappa2
      ind      <- (sign(si) == sign(si.target))
      si[!ind] <- 0
      k1[!ind] <- 0
      k2[!ind] <- 0
      values   <- pmax(abs(k1), abs(k2))
      # values   <- -pmin(k1, 0)
   }
   if (directions) {
      drns             <- array(NA, dim = c(3, 3, dim(shape$directions)[3]))
      drns[ , 1:2, ]   <- shape$directions[ , 1:2, ]
      drns[ ,   3, ]   <- t(shape$normals)
      shape$directions <- drns
   }

   
   if (graphics) {
      if (si.target.present) 
     	 display.face3d(shape, colour = values + 2)
      else
     	 display.face3d(shape)
      spheres3d(rbind(x1, x2), radius = 0.5, col = "red")
  	}
  
   triangles <- matrix(shape$triples, ncol = 3, byrow = TRUE)
   
   if (any(is.na(pts1)))                         pts1 <- closest.face3d(x1, shape)$id
   if (any(is.na(direction)) & any(is.na(pts2))) pts2 <- closest.face3d(x2, shape)$id
   triangle1  <- as.numeric(triangles[ , 1] %in% pts1) +
                 as.numeric(triangles[ , 2] %in% pts1) +
                 as.numeric(triangles[ , 3] %in% pts1)
   triangles1 <- which(triangle1 == length(pts1))
   triangle1  <- triangles1[1]
   pts1       <- triangles[triangle1, ]
   triangle2  <- as.numeric(triangles[ , 1] %in% pts2) +
                 as.numeric(triangles[ , 2] %in% pts2) +
                 as.numeric(triangles[ , 3] %in% pts2)
   triangles2 <- which(triangle2 == length(pts2))
   triangle2  <- triangles2[1]
   pts2       <- triangles[triangle2, ]
   # triangle1 <- which((triangles[ , 1] %in% pts1) | (triangles[ , 2] %in% pts1) |
                      # (triangles[ , 3] %in% pts1))
   # triangle2 <- which((triangles[ , 1] %in% pts2) | (triangles[ , 2] %in% pts2) |
                      # (triangles[ , 3] %in% pts2))

   # spheres3d(shape$coords[triangles[triangle1, ], ], radius = 0.5, col = "green")
   # spheres3d(shape$coords[triangles[triangle2, ], ], radius = 0.5, col = "green")

   if (any(triangles1 %in% triangles2))
      return(invisible(list(path = c(x1, x2), arclength = c(0, sqrt(sum((x1 - x2)^2))))))

   if (is.vector(pts1)) pts1 <- matrix(pts1, ncol = 3)
   if (any(is.na(direction))) {
      if (is.vector(pts2)) pts2 <- matrix(pts2, ncol = 3)
      if ("normals" %in% names(shape))
         normals <- shape$normals[unique(c(pts1, pts2)), ]
      else
         normals <- rbind(crossproduct(shape$coords[pts1[, 2], ] - shape$coords[pts1[, 1], ],
                                       shape$coords[pts1[, 3], ] - shape$coords[pts1[, 1], ]),
                          crossproduct(shape$coords[pts2[, 2], ] - shape$coords[pts2[, 1], ],
                                       shape$coords[pts2[, 3], ] - shape$coords[pts2[, 1], ]))
      normal <- apply(normals, 2, mean)
      if (sum(normal^2) < sqrt(.Machine$double.eps))
         normal <- apply(matrix(normals[1:length(pts1), ], ncol = 3), 2, mean)
      direction <- x2 - x1
   }
   else {
      if ("normals" %in% names(shape))
         normals <- shape$normals[unique(c(pts1)), ]
      else
         normals <- matrix(crossproduct(shape$coords[pts1[, 2], ] - shape$coords[pts1[, 1], ],
                                       (shape$coords[pts1[, 3], ] - shape$coords[pts1[, 1], ])),
                            ncol = 3)
      normal <- apply(normals, 2, mean)
   }
   cross  <- c(crossproduct(direction, normal))

   if ((length(rotation) == 1) & (rotation %in% c("optimise", "optimize")))
      angle.values <- seq(-rotation.range, rotation.range, length = ngrid)[-1]
   else if ((length(rotation) == 1) & is.numeric(rotation))
      angle.values <- rotation
   else
      stop("an invalid value of rotation has been supplied.")
   criterion <- numeric(0)
   
   for (angle in angle.values) {

   rcross     <- rotate3d(cross, angle, direction[1], direction[2], direction[3])
   # normal     <- c(crossproduct(cross, direction))
   # rnormal    <- rotate3d(normal,    angle, direction[1], direction[2], direction[3])
   # rdirection <- rotate3d(direction, angle, direction[1], direction[2], direction[3])
   # segments3d(rbind(x1, x1 + 10*rnormal, x1, x1 + 10*rcross, x1, x1 + 10*rdirection), col = "blue")

   ind       <- t(apply(triangles, 1,
                      function(a) as.numeric(c(shape$coords[a, ] %*% rcross <= sum(x1 * rcross)))))
   indsum    <- apply(ind, 1, sum)
   # print(indsum)
   indtrng   <- which(indsum > 0 & indsum < 3)
   trngls    <- triangles[indtrng, ]
   triangle1 <- triangles1[which(triangles1 %in% indtrng)][1]
   triangle2 <- triangles2[which(triangles2 %in% indtrng)][1]
   ind       <- ind[indtrng, ]
   j         <- which(apply(ind, 1, sum) == 2)
   ind[j, ]  <- 1 - ind[j, ]
   ind1      <- c(ind %*% 1:3)
   nind1     <- length(ind1)
   indnot1   <- rep(1:3, nind1)[-((0:(nind1 - 1)) * 3 + ind1)]
   indnot1   <- matrix(indnot1, ncol = 2, byrow = TRUE)
   j         <- 1:nind1
   edges     <- rbind(cbind(trngls[cbind(j, ind1)], trngls[cbind(j, indnot1[ , 1])]),
                      cbind(trngls[cbind(j, ind1)], trngls[cbind(j, indnot1[ , 2])]))
   e1        <- c(shape$coords[edges[, 1], ] %*% rcross)
   e2        <- c(shape$coords[edges[, 2], ] %*% rcross)
   wts       <- (sum(x1 * rcross) - e1) / (e2 - e1)
   crossings <- (1 - wts) * shape$coords[edges[ , 1], ] +
                     wts  * shape$coords[edges[ , 2], ]
   if (si.target.present)
      cvals <- (1 - wts) * values[edges[ , 1]] + wts  * values[edges[ , 2]]
   if ("directions" %in% names(shape)) {
      drns  <- sweep(shape$directions[ , , edges[ , 1]], 3, 1 - wts, "*") + 
               sweep(shape$directions[ , , edges[ , 1]], 3,     wts, "*")
      kp1   <- (1 - wts) * shape$kappa1[edges[ , 1]] + wts  * shape$kappa1[edges[ , 2]]
      kp2   <- (1 - wts) * shape$kappa2[edges[ , 1]] + wts  * shape$kappa2[edges[ , 2]]
   }
   
   edges <- t(apply(edges, 1, sort))
   e     <- paste(edges[ , 1], edges[ , 2])
   tbl   <- table(e)
   ends  <- which(e %in% names(which(tbl == 1)))
   if (length(ends) == 0) ends <- 1

   # if (graphics) spheres3d(crossings, col = "green", radius = 0.3)

   ind  <- integer(0)
   idx1 <- ends[1]
   ends <- ends[-1]
   trng <- integer()
   
   # Find the section which contains x1 (and x2, if required)
   # flag counts how many of x1 and x2 have been found.
   flag <- 0
   while (flag < 2) {
      idx     <- if (idx1 <= nind1) idx1 + nind1 else idx1 - nind1
      trngnew <- if (idx  <= nind1) idx else idx - nind1
      trng    <- c(trng, trngnew)
      idx1    <- which(e == e[idx])
      idx1    <- idx1[idx1 != idx]
      # Deal with circular paths
      flag1   <- length(which(triangle1 == indtrng[trng]))
      flag2   <- if (all(!is.na(x2))) length(which(triangle2 == indtrng[trng])) else 0
      if (flag1 > 1 | flag2 > 1) idx1 <- NULL
      flag1   <- min(flag1, 1)
      flag2   <- min(flag2, 1)
      # flag1   <- as.numeric(any(triangle1 %in% indtrng[trng]))
      # flag2   <- if (all(!is.na(x2))) as.numeric(any(triangle2 %in% indtrng[trng])) else 0
      # # flag2   <- if (all(!is.na(x2))) as.numeric(any(pts2 %in% indtrng[trngnew])) else 0
      flag    <- flag1 + flag2
      if (flag == 1 | all(is.na(x2))) ind <- c(ind, idx)
      # If an edge is encountered and bridge.gaps is FALSE, jump to another edge 
      # if no points have been found.  Give up if one point has been found.
      if (length(idx1) == 0) {
         if (flag == 1 & bridge.gaps) {
            # segments3d(shape$coords[c(t(edges[ends, ])), ], col = "blue", lwd = 3)
            # spheres3d(crossings[ends, ], col = "green", radius = 1, alpha = 1)
            ends <- ends[ends != idx]
            # Find the end point which is closer to the missing one of x1 or x2
            # but of those is the closest to idx.
            xx    <- if (flag1) x2 else x1
            ed0   <- sqrt(sum((crossings[idx, ] - xx)^2))
            ed1   <- edist.face3d(crossings[ends, ], xx)
            ed2   <- edist.face3d(crossings[ends, ], crossings[idx, ])
            ed    <- ends[ed1 <= ed0]
            idx1  <- ed[which.min(ed2[ed1 <= ed0])]
         	# cat("ends", ends, "\n")
         	# cat("ed", ed, "\n")
         	# cat("idx1", idx1, "\n")
            ind   <- c(ind, idx1)
            ends <- ends[ends != idx1]
         }
         else if (flag == 0) {
            if (length(ends) > 0) ends <- ends[-which(ends == idx)]
            trng <- integer(0)
            if (length(ends) > 0) {
               idx1 <- ends[1]
               ends <- ends[-1]
               flag <- 0
               ind  <- integer(0)
               }
            else
               flag <- 3
         	 # else if (length(idx) < 2 * nind1) {
          	 # 	idx1 <- which(!(1:(2*nind1) %in% ind))[1]
         	 # 	ind  <- c(ind, idx1)
         	 # }
         	 # else stop("a path cannot be identified.")
         }
         else if (flag == 1)
            flag <- if (all(is.na(x2))) 2 else 3
      }
   }

   path  <- crossings[ind, ]
   if (!is.matrix(path)) path <- matrix(path, nrow = 1)
   if (si.target.present) cvals <- cvals[ind]
   if ("directions" %in% names(shape)) {
      drns <- drns[ , , ind]
      kp1  <- kp1[ind]
      kp2  <- kp2[ind]
   }
   
   # Remove repeated points
   if (flag < 3 & nrow(path) > 1) {
      tol  <- sqrt(.Machine$double.eps)
      eps  <- apply(diff(path), 1, function(x) sqrt(sum(x^2)))
      ind  <- 1 + which(eps < tol)
      if (length(ind) > 0) {
         path <- path[-ind, ]
         if (si.target.present) cvals <- cvals[-ind]
         if ("directions" %in% names(shape)) {
            drns <- drns[ , , -ind]
            kp1  <- kp1[-ind]
            kp2  <- kp2[-ind]
         }
      }
   }
  
   if (graphics) {
      display.face3d(shape)
      spheres3d(matrix(c(path), ncol = 3), col = "green", radius = 0.5)
   }
   
   # If x1 and x2 are both present, restrict the path to the section between x1 and x2.
   if (all(!is.na(x2)) & (flag == 2)) {
      dist1 <- sqrt(sum((path[1, ] - x1)^2))
      distn <- sqrt(sum((path[nrow(path), ] - x1)^2))
      if (distn < dist1) {
         path  <- path[nrow(path):1, ]
         if (si.target.present) cvals <- cvals[nrow(path):1]
         if ("directions" %in% names(shape)) {
            drns <- drns[ , , nrow(path):1]
            kp1  <- kp1[nrow(path):1]
            kp2  <- kp2[nrow(path):1]
         }
      }
      if (sqrt(sum((x1 - path[1, ])^2)) < 1e-8) {
         path  <- path[-1, ]
         if (si.target.present) cvals <- cvals[-1]
         if ("directions" %in% names(shape)) {
            drns <- drns[ , , -1]
            kp1  <- kp1[-1]
            kp2  <- kp2[-1]
         }
         if (!is.matrix(path)) path <- matrix(path, nrow = 1)
      }
      if (sqrt(sum((x2 - path[nrow(path), ])^2)) < 1e-8) {
         path <- path[-nrow(path), ]
         if (si.target.present) cvals <- cvals[-nrow(path)]
         if ("directions" %in% names(shape)) {
            drns <- drns[ , , -nrow(path)]
            kp1  <- kp1[-nrow(path)]
            kp2  <- kp2[-nrow(path)]
         }
      }
      path   <- rbind(x1, path, x2)
      lngths <- apply(diff(path)^2, 1, function(x) sqrt(sum(x)))
      if (si.target.present)
         criterion <- c(criterion, -sum(cvals * lngths[-length(lngths)]) / sum(lngths[-length(lngths)]))
      else
         criterion  <- c(criterion, sum(lngths))
   }
   # If only x1 is present, restrict the path to the section from x1 in the appropriate direction.   
   if (any(is.na(x2))) {
      dst  <- edist.face3d(path, x1)
      proj <- c(sweep(path, 2, x1) %*% direction)
      ind  <- which(proj > 0)
      ind1 <- ind[which.min(dst[ind])]
      ind2 <- if (proj[ind1 - 1] > 0) 1 else nrow(path)
      if (bothways) {
      	 ind          <- (nrow(path) + 1 - ind2):ind1
      	 ind          <- ind[-length(ind)]
      	 path1        <- rbind(path[ind, ], x1)
         x1.arclength <- sum(apply(diff(path1)^2, 1, function(x) sqrt(sum(x))))
         path         <- rbind(path1[-nrow(path1), ], path[ind1:ind2, ])
      }
      else {
         x1.arclength <- 0
         path         <- rbind(x1, path[ind1:ind2, ])
      }
   }

   if (all(!is.na(x2)) & (flag == 3))
      criterion <- c(criterion, Inf)
   if (all(is.na(x2))  & (flag == 2))
      criterion <- c(criterion, sum(apply(diff(path)^2, 1, function(x) sqrt(sum(x)))))
   
   # print(criterion[length(criterion)])
   ind <- which.min(criterion)
   if (ind == length(criterion)) {
      criterion.min <- criterion[ind]
      path.min      <- path
      angle.min     <- angle
      if (si.target.present) cvals.min <- cvals
      if ("directions" %in% names(shape)) {
         drns.min <- drns
         kp1.min  <- kp1
         kp2.min  <- kp2
      }
   }
   
   # print(c(angle, flag, criterion[length(criterion)]))
   # spheres3d(path, radius = 2, col = "green")
   # scan()
   # pop3d()
   }
   
   if (graphics) {
   	  if (si.target.present) display.face3d(shape, colour = 2 + values)
   	  else                   display.face3d(shape)
      spheres3d(path.min, col = "green", radius = 0.5)
      spheres3d(rbind(x1, x2), col = "red", radius = 1)
   }
   
   arclength <- apply(diff(path.min)^2, 1, function(x) sqrt(sum(x)))
   arclength <- cumsum(c(0, arclength))
  
   result <- list(path = path.min, arclength = arclength, criterion = criterion.min,
                  pts1 = pts1, shape = shape, angle = angle.min, normal = normal)
   if (!any(is.na(x2)))
      result$pts2         <- pts2
   else
      result$x1.arclength <- x1.arclength
   if (si.target.present)  result$values <- values
   if ("directions" %in% names(shape)) {
      result$directions <- drns.min
      result$kappa1    <- kp1.min
      result$kappa2    <- kp2.min
   }
   invisible(result)
   
}
