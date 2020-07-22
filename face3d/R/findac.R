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
   sbst           <- translate.face3d(sbst, mn2.popn - mn2.image)
   a1             <- mn.popn["se", ] - mn2.popn
   a2             <- face$landmarks["se", ] - mn2.image
   raxis          <- c(crossproduct(a1, a2))
   angle          <- acos(sum(a1 * a2) / sqrt(sum(a1^2) * sum(a2^2)))
   sbst           <- rotate.face3d(sbst, angle, raxis, mn2.popn)

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
   imax3 <- which.max(mvtnorm::dmvnorm(sbst$vertices, c(mn.popn)[ind], covt[ind, ind], log = TRUE))
   ind   <- c(4, 4 + n.lmks, 4 + 2 * n.lmks)
   imax4 <- which.max(mvtnorm::dmvnorm(sbst$vertices, c(mn.popn)[ind], covt[ind, ind], log = TRUE))
   lmks  <- sbst$vertices[c(imax3, imax4), ]
   rownames(lmks) <- c("acR", "acL")
   lmks1 <- rbind(sbst$landmarks[c("pn", "se"), ], lmks)
   sbst  <- opa.face3d(lmks1, mn.popn, sbst,  scale = FALSE)
   lmks  <- opa.face3d(lmks1, mn.popn,        scale = FALSE)
   
   plot(sbst)
   spheres3d(lmks, col = "blue", radius = 2)
   for (j in 1:n.lmks) {
      ind    <- c(j, j + n.lmks, j + 2 * n.lmks)
      plot3d(ellipse3d(covt[ind, ind], centre = mn.popn[j, ]),
             col = "lightblue", alpha = 0.5, add = TRUE)
   }

   # inv <- ginv(covt)
   # x   <- c(lmks) - c(mn.popn)
   # mhd <- c(t(x) %*% inv %*% x)
   # print(mhd)
   
   sbst <- index.face3d(sbst, distance = 10, overwrite = TRUE, directions = TRUE)
   
   # plot(sbst, col = sbst$kappa1)
   # plot(sbst, col = pmax(sbst$kappa1, 0))
   # sbst1 <- subset(sbst, sbst$kappa1 > 0)
   # plot(sbst1, col = sbst1$kappa1)
   # plot(sbst, col = sbst$kappa2)
   # plot(sbst, col = pmin(sbst$kappa2, 0))
   # plot(sbst, col = sbst$kappa1 * sbst$kappa2)
   # spheres3d(face$landmarks, col = "yellow", radius = 2)
   # spheres3d(face$lmks[landmark.names, ], col = "red", radius = 2)
   
   sbst <- index.face3d(sbst, overwrite = TRUE, distance = 4,
                        subset = sbst$kappa2 > 0, extension = TRUE, directions = TRUE)
   # sbst <- index.face3d(sbst, overwrite = TRUE, distance = 10, subset = sbst$kappa1 < 0)
   plot(sbst, col = "kappa1")
   plot(sbst, col = "kappa2")
   spheres3d(lmks, col = "yellow", radius = 2)

   # dst  <- c(rdist(t(sbst$lmks["acR", ]), sbst$vertices))
   dst  <- c(rdist(t(lmks["acR", ]), sbst$vertices))
   boxplot(dst)
   ind  <- which.min(dst)
   acR  <- sbst$vertices[ind, ]
   # spheres3d(acR, col = "red", radius = 2)
   drn  <- sbst$directions[ , 1, ind]
   # dst1 <- rdist(t(sbst$lmks["pn", ]), rbind(acR, acR + drn))
   dst1 <- rdist(t(sbst$landmarks["pn", ]), rbind(acR, acR + drn))
   if (dst1[2] > dst1[1]) drn <- -drn
   # lines3d(rbind(acR, acR + 5 * drn))
   # plot(sbst, col = sbst$kappa2)
   # spheres3d(acR, col = "red", radius = 2)
   sbst1 <- subset(sbst, dst < 10)
   # plot(sbst1, col = sbst1$kappa2)
   # spheres3d(acR, col = "red", radius = 2)
   
   plot(sbst, col = "kappa1")
   plot(sbst, col = "kappa2")
   plot(sbst1, col = "kappa1")
   plot(sbst1, col = "kappa2")
   
   jlst <- which(sbst1$kappa1 > quantile(sbst1$kappa1, 0.5))
   crv  <- numeric(0)
   for (j in jlst) {
      acR   <- sbst1$vertices[j, ]
      dst   <- c(rdist(t(acR), sbst$vertices))
      sbst2 <- subset(sbst, dst < 10)
      drn   <- sbst1$directions[ , 1, j]
      # dst1  <- rdist(t(sbst1$lmks["pn", ]), rbind(acR, acR + drn))
      dst1  <- rdist(t(sbst1$landmarks["pn", ]), rbind(acR, acR + drn))
      if (dst1[2] > dst1[1]) drn <- -drn
      pnac  <- (sbst1$landmarks["pn", ] - acR)
      pnse  <- (sbst1$landmarks["se", ] - sbst1$landmarks["pn", ])
      nrm   <- sbst1$normals[j, ]
      print(pnac)
      print(pnse)
      print(nrm)
      print(acos(sum(drn * pnac) / sqrt(sum(pnac^2))))
      print(acos(sum(nrm * pnac) / sqrt(sum(pnac^2))))
      print(sum(nrm * pnac))
      if ((acos(sum(drn * pnac) / sqrt(sum(pnac^2))) < pi / 4) &
          (abs(acos(sum(nrm * pnac) / sqrt(sum(pnac^2))) - pi / 2) < pi / 4)) { 
         path  <- planepath.face3d(sbst2, acR, direction = drn, si.target = 1, rotation = 0)$path
         cdst   <- closestcurve.face3d(sbst2, path)
         area   <- area.face3d(sbst2)$points
         ind    <- (cdst$closest.curvept != 1) & (abs(cdst$closest.distance) < 4) 
         rdg    <- -sum(area[ind] * sbst2$kappa2[ind])
      }
      else
         rdg <- NA
      # ind1   <- ind & (cdst$closest.distance > 0)
      # ind2   <- ind & (cdst$closest.distance < 0)
      # rdg    <- sum(area[ind1] * sbst2$kappa2[ind1]) - sum(area[ind2] * sbst2$kappa2[ind2])
      
      # plot(sbst1, col = sbst1$kappa1)
      # plot(sbst1, col = sbst1$kappa2)
      # spheres3d(acR,  col = "red", radius = 2)
      # spheres3d(path, col = "yellow")
      # spheres3d(sbst1$vertices[ind, ])
      
      crv  <- c(crv, sbst1$kappa1[j] * rdg)
   }
   
   plot(sbst1, col = "kappa2")
   plot(sbst1, col = "kappa1")
   ind   <- which(!is.na(crv))
   sbstj <- subset(sbst1, jlst[ind], remove.singles = FALSE)
   plot(sbstj, col = crv[ind], display = "spheres", add = TRUE)
   acR <- sbstj$vertices[which.max(crv[ind]), ]
   spheres3d(acR, radius = 3)

   # Try neighbouring positions for acL/R
   threshold <- 3
   dst       <- c(rdist(t(sbst$lmks["acL", ]), sbst$vertices))
   nhd.L     <- which(dst < 3)
   dst       <- c(rdist(t(sbst$lmks["acR", ]), sbst$vertices))
   nhd.R     <- which(dst < 3)
   
   sbst$landmarks <- rbind(sbst$landmarks, lmks[3:4, ])
   rownames(sbst$landmarks)[5:6] <- c("acL", "acR")
   return(invisible(sbst))
}
