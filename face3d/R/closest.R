closest.face3d    <- function(x, shape, nearest = 30, tol = 1e-8) {

   closestone <- function(x, shape, nearest, tol) {

      # Find the closest vertex
      dst <- apply(shape$coords, 1, function(y) sqrt(sum((x - y)^2)))
      ind <- which.min(dst)
      id  <- c(ind, NA, NA)
      pts <- shape$coords[ind, ]
      if (dst[ind] < tol)
         return(invisible(c(pts, id, dst[ind])))
      
      # Create a subset consisting of the nearest vertices
      ind    <- order(dst)[1:min(nearest, length(dst))]
      sbst   <- subset.face3d(shape, ind, remove.singles = FALSE, retain.indices = TRUE)
   
      # Find the smallest of the perpendicular projections onto edges
      edges0 <- matrix(sbst$triples, ncol = 3, byrow = TRUE)[ , c(1, 2, 2, 3, 3, 1)]
      edges0 <- matrix(c(t(edges0)), ncol = 2, byrow = TRUE)
      origin <- sbst$coords[edges0[ , 1], ]
      x1     <- -sweep(origin, 2, x)
      edges  <- sbst$coords[edges0[ , 2], ] - sbst$coords[edges0[ , 1], ]
      lngth  <- apply(edges, 1, function(x) sqrt(sum(x^2)))
      edges  <- edges / matrix(rep(lngth, 3), ncol = 3)
      prjn   <- apply(x1 * edges, 1, sum)
      dst    <- apply(x1 - sweep(edges, 1, prjn, "*"), 1, function(x) sqrt(sum(x^2)))
      ind    <- which((prjn >= 0) & (prjn <= lngth))
      if (length(ind) > 0) {
         ind1 <- ind[which.min(dst[ind])]
         x1   <- prjn[ind1] * edges[ind1, ] + origin[ind1, ]
   	   pts  <- rbind(pts, x1)
   	   id   <- rbind(id, c(edges0[ind1, ], NA))
   	   dimnames(id) <- list(NULL, NULL)
   	   # edg <- matrix(as.numeric(rownames(sbst$coords))[c(edges0[ind, ])], ncol = 2)
   	   # id  <- rbind(id, cbind(edg, rep(NA, nrow(edg))))
         if (dst[ind1] < tol) {
            # return(invisible(list(point = pts[2, ], id.sub = id[2, 1:2],
            #                       id = as.numeric(rownames(sbst$coords))[id[2, 1:2]],
            #                       distance = dst[ind1], subset = sbst)))
            return(invisible(c(pts[2, ], id = c(sbst$subset[id[2, 1:2]], NA), dst[ind1])))
         }
      }
   
      # Find the smallest of the perpendicular projections onto triangles
      trpls <- matrix(sbst$triples, ncol = 3, byrow = TRUE)
      PX    <- sbst$coords[trpls[ , 2], ] - sbst$coords[trpls[ , 1], ]
      if (!is.matrix(PX)) PX <- matrix(PX, ncol = 3)
      lngth <- apply(PX, 1, function(x) sqrt(sum(x^2)))
      PX    <- PX / matrix(rep(lngth, 3), ncol = 3)
      PZ    <- crossproduct(PX, sbst$coords[trpls[ , 3], ] - sbst$coords[trpls[ , 1], ])
      PY    <- crossproduct(PZ, PX)
      sc1   <- matrix(c(sbst$coords[trpls[ , 1], ]), ncol = 3)
      sc2   <- matrix(c(sbst$coords[trpls[ , 2], ]), ncol = 3)
      sc3   <- matrix(c(sbst$coords[trpls[ , 3], ]), ncol = 3)
      X     <- -sweep(sc1, 2, x)
      XX    <- cbind(apply(PX * X, 1, sum), apply(PY * X, 1, sum))
      P1    <- matrix(0, ncol = 2, nrow = nrow(X))
      P2    <- sc2 - sc1
      P2    <- cbind(apply(PX * P2, 1, sum), apply(PY * P2, 1, sum))
      P3    <- sc3 - sc1
      P3    <- cbind(apply(PX * P3, 1, sum), apply(PY * P3, 1, sum))
      D     <- P2 - P1
      N     <- cbind(-D[ , 2], D[ , 1])
      L1    <- (apply(XX * N, 1, sum) >= 0)
      D     <- P3 - P2
      N     <- cbind(-D[ , 2], D[ , 1])
      L2    <- (apply((XX - P2) * N, 1, sum) >= 0)
      D     <- P1 - P3
      N     <- cbind(-D[ , 2], D[ , 1])
      L3    <- (apply((XX - P3) * N, 1, sum) >= 0)
      ind   <- which(L1 & L2 & L3)
      if (length(ind) > 0) {
         x1  <- XX[ind, 1] * PX[ind, ] + XX[ind, 2] * PY[ind, ] + matrix(c(sbst$coords[trpls[ind, 1], ]), ncol = 3)
   	   pts <- rbind(pts, x1)
   	   trp <- matrix(sbst$subset[c(trpls[ind, ])], ncol = 3)
   	   # trp <- matrix(as.numeric(rownames(sbst$coords))[c(trpls[ind, ])], ncol = 3)
   	   id  <- rbind(id, trp)
      }

      # Locate the closest value of pts
      pts   <- matrix(c(pts), ncol = 3)   # Ensure that pts is a matrix even if there is only a single point
      dst   <- apply(pts, 1, function(y) sqrt(sum((x - y)^2)))
      ind   <- which.min(dst)
      # Make sure id is a matrix, even if there is only a single row
      id    <- matrix(c(id), ncol = 3)
      id    <- id[ind, ]
      # id    <- id[!is.na(id)]
   
      # trpls <- matrix(shape$triples, ncol = 3, byrow = TRUE)
      # indtr <- which((trpls[ , 1] %in% id[ind, ]) | (trpls[ , 2] %in% id[ind, ]) |
      #                (trpls[ , 3] %in% id[ind, ]))

      # dst <- apply(shape$coords, 1, function(y) sqrt(sum((x - y)^2)))
      # ind <- which.min(dst)
      
      invisible(c(pts[ind, ], id, dst[ind]))

   }
   
   if (is.vector(x)) x <- matrix(x, nrow = 1)
   if (!is.matrix(x)) stop("x should be a vector or matrix.")
   if (ncol(x) != 3)  stop("x should be a vector of length 3 or a matrix with 3 columns.")
   result <- t(apply(x, 1, closestone, shape = shape, nearest = nearest, tol = tol))

   result <- list(points = result[ , 1:3], ids = result[ , 4:6], distances = result[ , 7])
   if (is.vector(result$ids)) result$ids <- result$ids[!is.na(result$ids)]
   invisible(result)
   # invisible(list(point = pts[ind, ], id = id, distance = dst[ind], subset = sbst))
   # invisible(list(point = pts[ind, ], id = id, distance = dst[ind], subset = sbst, id.num = names(ind)))
}
