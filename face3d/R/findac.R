# Find acL/R from pn and se

findac <- function(face) {
   
   # Subset to area around the nose
   mn2.image <- apply(face$landmarks[c("pn", "se"), ], 2, mean)
   dst       <- c(rdist(t(mn2.image), face$vertices))
   sbst      <- subset(face, dst < 50)
   
   # Move the image to match the mid-point of pn and se and the direction to se.
   mn2.popn       <- apply(mn.popn[c("pn", "se"), ], 2, mean)
   sbst           <- translate.face3d(sbst, mn2.popn - mn2.image)
   a1             <- mn.popn["se", ] - mn2.popn
   a2             <- face$landmarks["se", ] - mn2.image
   raxis          <- c(crossproduct(a1, a2))
   angle          <- acos(sum(a1 * a2) / sqrt(sum(a1^2) * sum(a2^2)))
   sbst           <- rotate.face3d(sbst, angle, raxis, mn2.popn)
   
   # Rotate the image around the pn-se line so that enL/R are as close to the image surface as possible
   raxis     <- sbst$landmarks["se", ] - mn2.popn
   ang.vec   <- seq(0, 2 * pi, length = 360)[-1]
   dens.vec  <- numeric(0)
   ind.mat   <- matrix(nrow = 0, ncol = 2)
   for (angle in ang.vec) {
      sbst1    <- rotate.face3d(sbst, angle, raxis, mn2.popn)
      ind      <- c(3, 3 + n.lmks, 3 + 2 * n.lmks)
      dens3    <- mvtnorm::dmvnorm(sbst$vertices, c(mn.popn)[ind], covt[ind, ind], log = TRUE)
      ind      <- c(4, 4 + n.lmks, 4 + 2 * n.lmks)
      dens4    <- mvtnorm::dmvnorm(sbst$vertices, c(mn.popn)[ind], covt[ind, ind], log = TRUE)
      ind3     <- which.max(dens3)
      ind4     <- which.max(dens4)
      dens.vec <- c(dens.vec, dens3[ind3] + dens4[ind4])
      ind.mat  <- rbind(ind.mat, c(ind3, ind4))
   }
   ind  <- which.max(dens.vec)
   lmks <- sbst$vertices[ind.mat[ind, ], ]
   rownames(lmks) <- c("acL", "acR")
   sbst <- rotate.face3d(sbst, ang.vec[ind], raxis, mn2.popn)

   # Find starting points for acL and acR.
   # Move the distribution to locate the highest density at the two (currently)
   # fixed landmarks and the highest values at the other two locations
   ind   <- c(3, 3 + n.lmks, 3 + 2 * n.lmks)
   imax3 <- which.max(mvtnorm::dmvnorm(sbst$vertices, c(mn.popn)[ind], covt[ind, ind], log = TRUE))
   ind   <- c(4, 4 + n.lmks, 4 + 2 * n.lmks)
   imax4 <- which.max(mvtnorm::dmvnorm(sbst$vertices, c(mn.popn)[ind], covt[ind, ind], log = TRUE))
   lmks  <- sbst$vertices[c(imax3, imax4), ]
   rownames(lmks) <- c("acL", "acR")
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
   sbst <- index.face3d(sbst, overwrite = TRUE, distance = 4,
                        subset = sbst$kappa2 > 0, extension = TRUE, directions = TRUE)
   
   plot(sbst, col = "kappa2")
   spheres3d(lmks, col = "yellow", radius = 2)

   for (ac.nm in c("acL", "acR")) {
      
      dst  <- c(rdist(t(lmks[ac.nm, ]), sbst$vertices))
      ind  <- which.min(dst)
      ac   <- sbst$vertices[ind, ]
      # spheres3d(acR, col = "red", radius = 2)
      drn  <- sbst$directions[ , 1, ind]
      dst1 <- rdist(t(sbst$landmarks["pn", ]), rbind(ac, ac + drn))
      if (dst1[2] > dst1[1]) drn <- -drn
      sbst1 <- subset(sbst, dst < 10)
      
      sbst0 <- subset(sbst, dst < 20 & sbst$kappa2 < 0)
      wts   <- sbst0$kappa2 / sum(sbst0$kappa2)
      vec   <- apply(sbst0$directions[ , 1, ], 1, function(x) sum(x * wts))
      vec   <- vec / sqrt(sum(vec^2))
      # spheres3d(ac, radius = 5)
      # lines3d(rbind(ac + 10 * vec, ac - 10 * vec))

      # lines3d(rbind(acR, acR + 5 * drn))
      # plot(sbst, col = sbst$kappa2)
      # spheres3d(acR, col = "red", radius = 2)
      # plot(sbst1, col = sbst1$kappa2)
      # spheres3d(acR, col = "red", radius = 2)
   
      # plot(sbst, col = "kappa1")
      # plot(sbst, col = "kappa2")
      # plot(sbst1, col = "kappa1")
      # plot(sbst1, col = "kappa2")
   
      jlst <- which(sbst1$kappa1 > quantile(sbst1$kappa1, 0.5))
      crv  <- numeric(0)
      for (j in jlst) {
         ac    <- sbst1$vertices[j, ]
         dst   <- c(rdist(t(ac), sbst$vertices))
         sbst2 <- subset(sbst, dst < 10)
         # drn   <- sbst1$directions[ , 1, j]
         drn   <- vec
         # dst1  <- rdist(t(sbst1$lmks["pn", ]), rbind(acR, acR + drn))
         dst1  <- rdist(t(sbst1$landmarks["pn", ]), rbind(ac, ac + drn))
         if (dst1[2] > dst1[1]) drn <- -drn
         pnac  <- (sbst1$landmarks["pn", ] - ac)
         pnse  <- (sbst1$landmarks["se", ] - sbst1$landmarks["pn", ])
         nrm   <- sbst1$normals[j, ]
         # if ((acos(sum(drn * pnac) / sqrt(sum(pnac^2))) < pi / 4) &
         #     (abs(acos(sum(nrm * pnac) / sqrt(sum(pnac^2))) - pi / 2) < pi / 4)) {
            path <- planepath.face3d(sbst2, ac, direction = drn, si.target = 1, rotation = 0)$path
            if (is.null(path)) 
               rdg <- NA
            else {
               cdst <- closestcurve.face3d(sbst2, path)
               area <- area.face3d(sbst2)$points
               ind  <- (cdst$closest.curvept != 1) & (abs(cdst$closest.distance) < 4) 
               rdg  <- -sum(area[ind] * sbst2$kappa2[ind])
            }
         # }
         # else
         #    rdg <- NA
            
         # ind1   <- ind & (cdst$closest.distance > 0)
         # ind2   <- ind & (cdst$closest.distance < 0)
         # rdg    <- sum(area[ind1] * sbst2$kappa2[ind1]) - sum(area[ind2] * sbst2$kappa2[ind2])
         crv  <- c(crv, sbst1$kappa1[j] * rdg)
      }
   
      ind   <- which(!is.na(crv))
      sbstj <- subset(sbst1, jlst[ind], remove.singles = FALSE)
      ac   <- sbstj$vertices[which.max(crv[ind]), ]
      rnms  <- rownames(sbst$landmarks)
      sbst$landmarks <- rbind(sbst$landmarks, ac)
      rownames(sbst$landmarks) <- c(rnms, ac.nm)
   
      plot(sbstj, col = crv[ind], display = "spheres", add = TRUE)
      spheres3d(ac, radius = 3)
   }

   return(invisible(sbst))
}
