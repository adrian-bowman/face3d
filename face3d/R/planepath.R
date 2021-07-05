planepath.face3d <- function(shape, x1, x2, pts1, pts2, direction, normal,
                         bothways = FALSE,
                         distance = 5, ngrid = 50, boundary, 
                         rotation = "optimise", rotation.range = 0.8 * pi/2,
                         bridge.gaps = FALSE, si.target, directions = FALSE,
								 monitor = FALSE, monitor.prompt = monitor) {

   if (missing(x1)) stop("x1 must be specified.")
   if (is.null(x1) | any(is.na(x1))) stop("x1 cannot be null or contain missing values.")
   if (missing(pts1))      pts1      <- NA
   if (missing(pts2))      pts2      <- NA
   if (missing(x2))        x2        <- NA
   if (missing(direction)) direction <- NA
   if (missing(boundary)) {
      boundary   <- if (all(is.na(x2))) NA else c(1, 1)
   }
   if (any(is.na(x2)) & any(is.na(direction))) stop("one of x2 or direction must be specified.")
   si.target.present <- !missing(si.target)

   if (!any(is.na(x2)) && sqrt(sum((x1 - x2)^2)) < .Machine$double.eps)
      return(invisible(list(path = rbind(x1, x2), arclength = c(0, 0))))
   
   # x1 and x2 must lie on the surface
   x1 <- closest.face3d(x1, shape)$points
   if (!any(is.na(x2))) x2 <- closest.face3d(x2, shape)$points

   #  Subset the object by a tube around x1 and the relevant direction
   
   if (all(!is.na(boundary))) {
      if (any(is.na(x2)))
      	ind1 <- (apply(sweep(shape$vertices, 2, x1), 1, function(x) sqrt(sum(x^2))) < boundary[1])
      else {
         rng      <- sqrt(sum((x2 - x1)^2))
     	   unit     <- (x2 - x1) / rng
     	   ind1     <- TRUE
     	   boundary <- boundary * rng / 2
     	   j <- 1
     	   while(length(which(ind1)) < 12 & j < 6) {
      	   boundary <- boundary * 2
     	      prjn     <- c(sweep(shape$vertices, 2, x1) %*% unit)
            ind1     <- (apply(outer(prjn, unit) - sweep(shape$vertices, 2, x1), 1, function(x) sqrt(sum(x^2))) < boundary[2]) &
                        (prjn > - boundary[1]) & (prjn < rng + boundary[1])
            j        <- j + 1
     	   }
      }
      if (length(which(ind1)) < 12)
         warning("the subset defined by boundary is very small.  Boundary restrictions will not be applied.")
      else
         shape <- subset.face3d(shape, ind1, retain.indices = TRUE)
   }
   
   if (monitor > 0) plot(shape, display = c("mesh", "lines"))

      # Check whether the boundary setting has split the shape into parts with x1 and x2 in different parts
   if (!any(is.na(x2)) & !bridge.gaps) {
      parts <- connected.face3d(shape)
      nparts <- length(unique(parts))
      # spheres3d(shape$vertices[parts ==2, ], col = "red")
      if (nparts > 1) {
         d1 <- numeric(nparts)
         d2 <- numeric(nparts)
         for (i in 1:nparts) {
            shape.i <- subset(shape, parts == i)
            d1[i]   <- closest.face3d(x1, shape.i)$distances
            d2[i]   <- closest.face3d(x2, shape.i)$distances
         }
         if (which.min(d1) != which.min(d2)) {
            cat("The path cannot be computed with the current setting of the boundary parameter.\n")
            cat("Try setting boundary = NA.\n")
            return()
         }
      }
   }
   gspheres <- max(apply(shape$vertices, 2, function(x) diff(range(x)))) / 64

   # Calculate the surface curvatures
   
   if (si.target.present | directions) {
   	  if (!all(c("shape.index", "kappa1", "kappa2", "directions") %in% names(shape)))
         shape <- index.face3d(shape, distance = distance, directions = directions,
                               overwrite = TRUE)
   }
   if (si.target.present) {
      si       <- shape$shape.index
      k1       <- shape$kappa1
      k2       <- shape$kappa2
      ind      <- (sign(si) == sign(si.target))
      si[!ind] <- 0
      k1[!ind] <- 0
      k2[!ind] <- 0
      if(sign(si.target) == -1) 
         values <- pmax(shape$kappa1, 0) # valley
      else
      	values <- pmax(-shape$kappa2, 0) # ridge
      # values   <- -pmin(k1, 0)
   }
   if (directions) {
      drns             <- array(NA, dim = c(3, 3, dim(shape$directions)[3]))
      drns[ , 1:2, ]   <- shape$directions[ , 1:2, ]
      drns[ ,   3, ]   <- t(shape$normals)
      shape$directions <- drns
   }
   
   if (monitor) {
       if (si.target.present)
     	    plot(shape, display = c("mesh", "lines"), colour = values + 2)
     	 else
     	    plot(shape, display = c("mesh", "lines"))
     if (monitor.prompt) cat("Press the Return key to see successive images.\n")
     rgl::spheres3d(rbind(x1, x2), radius = gspheres, col = "red")
     # Create a dummy object for removal when iterations begin
     rgl::spheres3d(rbind(x1, x2), radius = gspheres, alpha = 0)
   }
   
   triangles <- shape$triangles
   
   if (any(is.na(pts1)))                         pts1 <- closest.face3d(x1, shape)$ids
   if (any(is.na(direction)) & any(is.na(pts2))) pts2 <- closest.face3d(x2, shape)$ids
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

   # If the triangles containing x1 and x2 are adjacent, return the path c(x1, x2).
   if (any(triangles1 %in% triangles2))
      return(invisible(list(path = rbind(x1, x2), arclength = c(0, sqrt(sum((x1 - x2)^2))))))
   
   if (is.vector(pts1)) pts1 <- matrix(pts1, ncol = 3)
   if (missing(normal)) {
      if (any(is.na(direction))) {
         if (is.vector(pts2)) pts2 <- matrix(pts2, ncol = 3)
         if ("normals" %in% names(shape))
            normals <- shape$normals[unique(c(pts1, pts2)), ]
         else
            normals <- rbind(crossproduct(shape$vertices[pts1[, 2], ] - shape$vertices[pts1[, 1], ],
                                          shape$vertices[pts1[, 3], ] - shape$vertices[pts1[, 1], ]),
                             crossproduct(shape$vertices[pts2[, 2], ] - shape$vertices[pts2[, 1], ],
                                          shape$vertices[pts2[, 3], ] - shape$vertices[pts2[, 1], ]))
         normal <- apply(normals, 2, mean)
         if (sum(normal^2) < sqrt(.Machine$double.eps))
            normal <- apply(matrix(normals[1:length(pts1), ], ncol = 3), 2, mean)
      }
      else {
         if ("normals" %in% names(shape))
            normals <- shape$normals[unique(c(pts1)), ]
         else
            normals <- matrix(crossproduct(shape$vertices[pts1[, 2], ] - shape$vertices[pts1[, 1], ],
                                          (shape$vertices[pts1[, 3], ] - shape$vertices[pts1[, 1], ])),
                               ncol = 3)
         normal <- apply(normals, 2, mean)
      }
   }
   if (any(is.na(direction))) direction <- x2 - x1
   cross  <- c(crossproduct(direction, normal))

   if ((length(rotation) == 1) &&
              (rotation %in% c("optimise", "optimize", "maximise", "maximize")))
      angle.values <- seq(-rotation.range, rotation.range, length = ngrid)[-1]
   else if (is.numeric(rotation))
      angle.values <- rotation
   else
      stop("an invalid value of rotation has been supplied.")
   criterion <- numeric(0)
   
   for (angle in angle.values) {

   rcross     <- c(rgl::rotate3d(cross, angle, direction[1], direction[2], direction[3]))
   # normal     <- c(crossproduct(cross, direction))
   # rnormal    <- rotate3d(normal,    angle, direction[1], direction[2], direction[3])
   # rdirection <- rotate3d(direction, angle, direction[1], direction[2], direction[3])
   # segments3d(rbind(x1, x1 + 10*rnormal, x1, x1 + 10*rcross, x1, x1 + 10*rdirection), col = "blue")

   ind       <- t(apply(triangles, 1,
                      function(a) as.numeric(c(shape$vertices[a, ] %*% rcross <= sum(x1 * rcross)))))
   indsum    <- apply(ind, 1, sum)
   indtrng   <- which(indsum > 0 & indsum < 3)
   # Proceed with this angle only if crossings points have been identified.
   if (length(indtrng) == 0) {
      if (monitor > 0) cat("No crossing points identified for angle", angle, "\n")
   }
   else {
      trngls    <- triangles[indtrng, , drop = FALSE]
      triangle1 <- triangles1[which(triangles1 %in% indtrng)][1]
      triangle2 <- triangles2[which(triangles2 %in% indtrng)][1]
      ind       <- ind[indtrng, , drop = FALSE]
      j         <- which(apply(ind, 1, sum) == 2)
      ind[j, ]  <- 1 - ind[j, ]
      ind1      <- c(ind %*% 1:3)
      nind1     <- length(ind1)
      indnot1   <- rep(1:3, nind1)[-((0:(nind1 - 1)) * 3 + ind1)]
      indnot1   <- matrix(indnot1, ncol = 2, byrow = TRUE)
      j         <- 1:nind1
      edges     <- rbind(cbind(trngls[cbind(j, ind1)], trngls[cbind(j, indnot1[ , 1])]),
                         cbind(trngls[cbind(j, ind1)], trngls[cbind(j, indnot1[ , 2])]))
      e1        <- c(shape$vertices[edges[, 1], ] %*% rcross)
      e2        <- c(shape$vertices[edges[, 2], ] %*% rcross)
      wts       <- (sum(x1 * rcross) - e1) / (e2 - e1)
      crossings <- (1 - wts) * shape$vertices[edges[ , 1], ] +
                        wts  * shape$vertices[edges[ , 2], ]
      if (si.target.present)
         cvals <- (1 - wts) * values[edges[ , 1]] + wts  * values[edges[ , 2]]
      if ("directions" %in% names(shape)) {
         drns  <- sweep(shape$directions[ , , edges[ , 1]], 3, 1 - wts, "*") + 
                  sweep(shape$directions[ , , edges[ , 2]], 3,     wts, "*")
         lngth <- apply(drns, 2:3, function(x) sqrt(sum(x^2)))
         drns  <- sweep(drns, 2:3, lngth, "/")
         kp1   <- (1 - wts) * shape$kappa1[edges[ , 1]] + wts  * shape$kappa1[edges[ , 2]]
         kp2   <- (1 - wts) * shape$kappa2[edges[ , 1]] + wts  * shape$kappa2[edges[ , 2]]
      }
      
      if (as.numeric(monitor) == 2) rgl::spheres3d(crossings, radius = gspheres/5)

      edges <- t(apply(edges, 1, sort))
         e     <- paste(edges[ , 1], edges[ , 2])
      tbl   <- table(e)
      ends  <- which(e %in% names(which(tbl == 1)))
      if (length(ends) == 0) ends <- 1

      # Find the edge point which is closest to x1
      ind  <- which.min(apply(sweep(matrix(c(crossings[ends, ]), ncol = 3), 2, x1)^2, 1, sum))
      idx1 <- ends[ind]
      ends <- ends[-ind]
      trng <- integer()
      ind  <- integer(0)
      # ind  <- idx1

      # Find the section which contains x1 (and x2, if required)
      # The variable flag counts how many of x1 and x2 have been found.
      flag    <- 0
      idx.vec <- numeric(0)
      while (flag < 2) {
         idx     <- if (idx1 <= nind1) idx1 + nind1 else idx1 - nind1
         trngnew <- if (idx  <= nind1) idx else idx - nind1
         trng    <- c(trng, trngnew)
         # spheres3d(shape$vertices[c(triangles[indtrng[trngnew], ]), ], radius = 0.5)
         if (!any(e == e[idx])) stop("e problem.")
         idx1    <- which(e == e[idx])
         idx1    <- idx1[idx1 != idx]
         # flag1   <- length(which(triangle1 == indtrng[trng]))
         # flag2   <- if (all(!is.na(x2))) length(which(triangle2 == indtrng[trng])) else 0
         flag1   <- length(which(triangle1 %in% indtrng[trng]))
         flag2   <- if (all(!is.na(x2))) length(which(triangle2 %in% indtrng[trng])) else 0
         idx1    <- if (flag1 > 1 | flag2 > 1) NULL else idx1
         flag1   <- min(flag1, 1)
         flag2   <- min(flag2, 1)
         # flag1   <- as.numeric(any(triangle1 %in% indtrng[trng]))
         # flag2   <- if (all(!is.na(x2))) as.numeric(any(triangle2 %in% indtrng[trng])) else 0
         # # flag2   <- if (all(!is.na(x2))) as.numeric(any(pts2 %in% indtrng[trngnew])) else 0
         flag    <- flag1 + flag2
         # Deal with circular paths
         if (any(idx == idx.vec)) flag <- 3 else idx.vec <- c(idx, idx.vec)
         if (flag == 1 | all(is.na(x2))) ind <- c(ind, idx)
         # print(ends)
         # print(flag)
         # print(idx1)
         # print(matrix(c(crossings[idx, ]), ncol = 3))
         # spheres3d(matrix(c(crossings[idx, ]), ncol = 3), radius = 1.2 * gspheres, col = "blue")
         # spheres3d(matrix(c(crossings[idx, ]), ncol = 3), radius = 0.21, col = "blue")
         # If an edge is encountered and bridge.gaps is FALSE, jump to another edge 
         # if no points have been found.  Give up if one point has been found.
         if (length(idx1) == 0) {
            if (flag == 1 & bridge.gaps) {
               # segments3d(shape$vertices[c(t(edges[ends, ])), ], col = "blue", lwd = 3)
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
   
      # Keep a record of the triangles in the path (only works in the case where x1 and x2 are both present)
      if (flag != 3 & all(!is.na(x2))) {
         # i1 <- which(triangle1 == indtrng[trng])[1]
         # i2 <- which(triangle2 == indtrng[trng])[1]
         # print(flag)
         # print(indtrng[trng])
         # print(triangle2)
         i1 <- which(indtrng[trng] %in% triangle1)[1]
         i2 <- which(indtrng[trng] %in% triangle2)[1]
         # print(c(i1, i2))
         trngls <- indtrng[trng][i1:i2]
      }
   
      path  <- crossings[ind, , drop = FALSE]
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
         if (si.target.present) cvals <- c(NA, cvals, NA)
         if ("directions" %in% names(shape)) {
            dm    <- dim(drns)
            dm[3] <- dm[3] + 2
            matd  <- array(NA, dim = dm)
            matd[ , , 2:(dm[3]-1)] <- drns
            drns  <- matd
            kp1   <- c(NA, kp1, NA)
            kp2   <- c(NA, kp2, NA)
         }
      }
   
      # If only x1 is present, restrict the path to the section from x1 in the appropriate direction.   
      if (any(is.na(x2)) & flag < 3) {
         dst  <- edist.face3d(path, x1)
         proj <- c(sweep(path, 2, x1) %*% direction)
         ind  <- which(proj > 0)
         ind1 <- ind[which.min(dst[ind])]
         ind2 <- if (ind1 > 1 && proj[ind1 - 1] > 0) 1 else nrow(path)
         if (bothways) {
            ind          <- ind2:ind1
         	path1        <- rbind(path[ind, ], x1)
         	x1.arclength <- sum(apply(diff(path1)^2, 1, function(x) sqrt(sum(x))))
         	ind3         <- if (ind2 ==1) ind1 + 1 else ind1 - 1
         	ind4         <- nrow(path) + 1 - ind2
            path         <- rbind(path1, path[ind3:ind4, ])
            ind5         <- c(ind, NA, ind3:ind4)
            if (si.target.present) cvals <- c(cvals[ind], NA, cvals[ind3:ind4])
            if ("directions" %in% names(shape)) {
               dm        <- dim(drns)
               drns1     <- array(dim = c(dm[1:2], length(ind) + 2 + abs(ind4 - ind3)))
               drns1[ , , 1:length(ind)] <- drns[ , , ind]
               drns1[ , , length(ind) + 2 + 0:(ind4 - ind3)] <- drns[ , , ind3:ind4]
               drns      <- drns1
               kp1       <- kp1[ind5]
               kp2       <- kp2[ind5]
            }
         }
         else {
            x1.arclength <- 0
            path         <- rbind(x1, path[ind1:ind2, ])
            if (si.target.present) cvals <- c(NA, cvals[ind1:ind2])
            if ("directions" %in% names(shape)) {
               dm        <- dim(drns)
               dm[3]     <- length(ind1:ind2) + 1
               matd      <- array(NA, dim = dm)
               matd[ , , 2:dm[3]] <- drns[ , , ind1:ind2]
               drns      <- matd
               kp1       <- c(NA, kp1[ind1:ind2])
               kp2       <- c(NA, kp2[ind1:ind2])
            }
         }
      }

      # Compute the criterion
      if (flag == 3)
         crit <- Inf
      else {
         lngths <- apply(diff(path)^2, 1, function(x) sqrt(sum(x)))
         if (si.target.present)
            crit <- -sum(cvals[-1] * lngths, na.rm = TRUE) / sum(lngths)
            # Old version
            # crit <- -sum(cvals[-1] * lngths[-length(lngths)], na.rm = TRUE) / 
            #                            sum(lngths[-length(lngths)]))
         else if (rotation[1] %in% c("maximise", "maxmize"))
            crit <- -sum(lngths)
         else
            crit <-  sum(lngths)
      }
      # if (all(is.na(x2)))
      #    criterion <- c(criterion, sum(apply(diff(path)^2, 1, function(x) sqrt(sum(x)))))
      criterion <- c(criterion, crit)
   
      # Identify whether this is the minimum so far
      ind <- which.min(criterion)
      if (ind == length(criterion)) {
         criterion.min <- criterion[ind]
         path.min      <- path
         angle.min     <- angle
         if (all(!is.na(x2))) trngs.min <- trngls
         if (si.target.present) cvals.min <- cvals
         if ("directions" %in% names(shape)) {
            drns.min <- drns
            kp1.min  <- kp1
            kp2.min  <- kp2
         }
      }
   
      if (monitor) {
      	  rgl::pop3d()
      	  rgl::spheres3d(matrix(c(path), ncol = 3), radius = gspheres)
      	  cat("angle:", angle, "criterion:", criterion[length(criterion)])
      	  if (monitor.prompt) scan(quiet = TRUE)
      }
      
   }  # End of if block after test for the presence of crossings
   
   }  # End of the angles for loop.
   
   if (nrow(path.min) > 0) {
      arclength <- apply(diff(path.min)^2, 1, function(x) sqrt(sum(x)))
      arclength <- cumsum(c(0, arclength))
   }
   else
      arclength <- NULL
   
   if (criterion.min == Inf) {
      path.min <- NULL
      cat("A path cannot be identified.\n")
      return(invisible(NULL))
   }
   
   if (monitor & !is.null(path.min)) {
      rgl::pop3d()
      if (si.target.present) plot(shape, colour = 2 + values)
      else                   plot(shape)
      rgl::spheres3d(path.min, col = "green", radius = gspheres)
      rgl::spheres3d(rbind(x1, x2), col = "red", radius = 1.2 * gspheres)
   }
   
   result <- list(path = path.min, arclength = arclength, criterion = criterion.min,
                  pts1 = c(pts1), shape = shape, angle = angle.min, normal = normal)
   if (!any(is.na(x2))) {
      result$pts2      <- c(pts2)
      if ("subset" %in% names(shape)) {
         sb.ind           <- shape$subset
         result$triangles <- matrix(sb.ind[c(triangles[trngs.min, ])], ncol = 3)
      }
   }
   else
      result$x1.arclength <- x1.arclength
   if (si.target.present)  result$values <- values
   if ("directions" %in% names(shape)) {
      result$directions  <- drns.min
      result$kappa1      <- kp1.min
      result$kappa2      <- kp2.min
      result$shape.index <- 2 / pi * (atan((kp2.min + kp1.min) / (kp2.min - kp1.min)))
   }
   invisible(result)
   
}
