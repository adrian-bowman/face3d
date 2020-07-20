# Find acL/R from pn and se

findac <- function(face) {
   
   # Subset to area around the nose
   mn2.image <- apply(face$landmarks[c("pn", "se"), ], 2, mean)
   dst       <- c(rdist(t(mn2.image), face$vertices))
   sbst      <- subset(face, dst < 50)
   
   plot(sbst)
   spheres3d(face$landmarks, col = "blue",   radius = 2)
   spheres3d(sbst$lmks[landmark.names, ], col = "yellow", radius = 2)
   
   # Move the image to match the mid-point of pn and se and the direction to se.
   mn2.popn       <- apply(mn.popn[c("pn", "se"), ], 2, mean)
   sbst$landmarks <- sweep(sbst$landmarks, 2, mn2.image - mn2.popn)
   sbst$vertices  <- sweep(sbst$vertices,  2, mn2.image - mn2.popn)
   sbst$lmks      <- sweep(sbst$lmks,      2, mn2.image - mn2.popn)
   a1             <- mn.popn["se", ] - mn2.popn
   a2             <- face$landmarks["se", ] - mn2.image
   raxis          <- c(crossproduct(a1, a2))
   angle          <- acos(sum(a1 * a2) / sqrt(sum(a1^2) * sum(a2^2)))
   sbst$vertices  <- sweep(sbst$vertices,  2, mn2.popn)
   sbst$landmarks <- sweep(sbst$landmarks, 2, mn2.popn)
   sbst$lmks      <- sweep(sbst$lmks, 2, mn2.popn)
   sbst$vertices  <- rotate3d(sbst$vertices,  angle, raxis[1] , raxis[2], raxis[3])
   sbst$landmarks <- rotate3d(sbst$landmarks, angle, raxis[1] , raxis[2], raxis[3])
   sbst$lmks      <- rotate3d(sbst$lmks, angle, raxis[1] , raxis[2], raxis[3])
   sbst$vertices  <- sweep(sbst$vertices,  2, mn2.popn, "+")
   sbst$landmarks <- sweep(sbst$landmarks, 2, mn2.popn, "+")
   sbst$lmks      <- sweep(sbst$lmks, 2, mn2.popn, "+")
   
   plot(sbst)
   spheres3d(sbst$landmarks[c("pn", "se"), ], col = "blue", radius = 2)
   # spheres3d(mn.popn, col = "yellow", radius = 2)
   for (j in 1:n.lmks) {
      ind    <- c(j, j + n.lmks, j + 2 * n.lmks)
      plot3d(ellipse3d(covt[ind, ind], centre = mn.popn[j, ]),
             col = "lightblue", alpha = 0.5, add = TRUE)
   }
   spheres3d(sbst$lmks[landmark.names, ], col = "yellow", radius = 2)
   
   # Find starting points for acL and acR.
   # Move the distribution to locate the highest density at the two (currently)
   # fixed landmarks and the highest values at the other two locations
   
   ind   <- c(3, 3 + n.lmks, 3 + 2 * n.lmks)
   imax3 <- which.max(dmvnorm(sbst$vertices, c(mn.popn)[ind], covt[ind, ind], log = TRUE))
   ind   <- c(4, 4 + n.lmks, 4 + 2 * n.lmks)
   imax4 <- which.max(dmvnorm(sbst$vertices, c(mn.popn)[ind], covt[ind, ind], log = TRUE))
   lmks  <- sbst$vertices[c(imax3, imax4), ]
   lmks           <- rbind(sbst$landmarks[c("pn", "se"), ], lmks)
   sbst$vertices  <- opa.face3d(lmks, mn.popn, sbst$vertices,  scale = FALSE)
   sbst$lmks      <- opa.face3d(lmks, mn.popn, sbst$lmks,      scale = FALSE)
   lmks           <- opa.face3d(lmks, mn.popn, scale = FALSE)
   
   plot(sbst)
   spheres3d(lmks, col = "blue", radius = 2)
   for (j in 1:n.lmks) {
      ind    <- c(j, j + n.lmks, j + 2 * n.lmks)
      plot3d(ellipse3d(covt[ind, ind], centre = mn.popn[j, ]),
             col = "lightblue", alpha = 0.5, add = TRUE)
   }
   spheres3d(sbst$lmks[landmark.names, ], col = "yellow", radius = 2)
   
   inv <- ginv(covt)
   x   <- c(lmks) - c(mn.popn)
   mhd <- c(t(x) %*% inv %*% x)
   print(mhd)
   
   # Try neighbouring positions for acL/R
   threshold <- 3
   dst       <- c(rdist(t(sbst$lmks["acL", ]), sbst$vertices))
   nhd.L     <- which(dst < 3)
   dst       <- c(rdist(t(sbst$lmks["acR", ]), sbst$vertices))
   nhd.R     <- which(dst < 3)
   
   
}