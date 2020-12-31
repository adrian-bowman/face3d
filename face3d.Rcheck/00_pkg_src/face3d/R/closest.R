closest.face3d    <- function(x, shape, nearest = 30){

   dst    <- apply(shape$coords, 1, function(y) sqrt(sum((x - y)^2)))
   id     <- matrix(c(which.min(dst), NA, NA), nrow = 1)
   pts    <- shape$coords[which.min(dst), ]
   
   ind    <- order(dst)[1:min(nearest, length(dst))]
   sbst   <- subset.face3d(shape, ind, remove.singles = FALSE)
   
   edges0 <- matrix(sbst$triples, ncol = 3, byrow = TRUE)[ , c(1, 2, 2, 3, 3, 1)]
   edges0 <- matrix(c(t(edges0)), ncol = 2, byrow = TRUE)
   origin <- sbst$coords[edges0[ , 1], ]
   x1     <- -sweep(origin, 2, x)
   edges  <- sbst$coords[edges0[ , 2], ] - sbst$coords[edges0[ , 1], ]
   lngth  <- apply(edges, 1, function(x) sqrt(sum(x^2)))
   edges  <- edges / matrix(rep(lngth, 3), ncol = 3)
   prjn   <- apply(x1 * edges, 1, sum)
   ind    <- which((prjn >= 0) & (prjn <= lngth))
   if (length(ind) > 0) {
   	  x1  <- prjn[ind] * edges[ind, ] + origin[ind, ]
   	  pts <- rbind(pts, x1)
   	  edg <- matrix(as.numeric(rownames(sbst$coords))[c(edges0[ind, ])], ncol = 2)
   	  id  <- rbind(id, cbind(edg, rep(NA, nrow(edg))))
   }

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
   	trp <- matrix(as.numeric(rownames(sbst$coords))[c(trpls[ind, ])], ncol = 3)
   	id  <- rbind(id, trp)
   }
   
   # Locate the closest value of pts
   pts   <- matrix(c(pts), ncol = 3)   # Ensure that pts is a matrix even if there is only a single point
   dst   <- apply(pts, 1, function(y) sqrt(sum((x - y)^2)))
   ind   <- which.min(dst)
   id    <- id[ind, ]
   id    <- id[!is.na(id)]

   # trpls <- matrix(shape$triples, ncol = 3, byrow = TRUE)
   # indtr <- which((trpls[ , 1] %in% id[ind, ]) | (trpls[ , 2] %in% id[ind, ]) |
   #                (trpls[ , 3] %in% id[ind, ]))

   # dst <- apply(shape$coords, 1, function(y) sqrt(sum((x - y)^2)))
   # ind <- which.min(dst)
   # invisible(list(point = shape$coords[ind, ], id = ind, distance = dst[ind]))
   
   invisible(list(point = pts[ind, ], id = id, distance = dst[ind], subset = sbst, id.num = names(ind)))
}
