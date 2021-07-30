# Find acL/R from pn and se

findac <- function(face, monitor = 0) {
   
   # Subset to area around the nose
   mn2.image <- apply(face$landmarks[c("pn", "se"), ], 2, mean)
   dst       <- c(rdist(t(mn2.image), face$vertices))
   sbst      <- face
   
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
   
   dstL <- c(rdist(t(lmks["acL", ]), sbst$vertices))
   dstR <- c(rdist(t(lmks["acR", ]), sbst$vertices))
   ind  <- (dstL < 20) | (dstR < 20)
   sbst <- curvatures(sbst, distance = 10, overwrite = TRUE, directions = TRUE,
                        subset = ind)
   sbst <- curvatures(sbst, distance =  4, overwrite = TRUE, directions = TRUE,
                        subset = ind & sbst$kappa2 > 0)
   
   plot(sbst, col = "kappa2")
   plot(sbst, display = "principal 1", add = TRUE)
   spheres3d(lmks, col = "yellow", radius = 2)
   
   results <- matrix(nrow = 0, ncol = 6)

   for (ac.nm in c("acL", "acR")) {
      
      dst   <- c(rdist(t(lmks[ac.nm, ]), sbst$vertices))
      ind   <- which.min(dst)
      ac    <- sbst$vertices[ind, ]
      sbst1 <- subset(sbst, dst < 10, retain.indices = TRUE)
      # drn  <- sbst$directions[ , 1, ind]
      # dst1 <- rdist(t(sbst$landmarks["pn", ]), rbind(ac, ac + drn))
      # if (dst1[2] > dst1[1]) drn <- -drn

      # sbst0 <- subset(sbst, dst < 20 & sbst$kappa2 < 0)
      # wts   <- sbst0$kappa2 / sum(sbst0$kappa2)
      # vec   <- apply(sbst0$directions[ , 1, ], 1, function(x) sum(x * wts))
      # vec   <- vec / sqrt(sum(vec^2))
      # lines3d(rbind(ac + 10 * vec, ac - 10 * vec))

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
   
      jlst  <- which(sbst1$kappa1 > quantile(sbst1$kappa1, 0.5))
      crv   <- numeric(0)
      first <- TRUE
      for (j in jlst) {
         ac    <- sbst1$vertices[j, ]
         dst   <- c(rdist(t(ac), sbst$vertices))
         sbst2 <- subset(sbst, dst < 10)
         # Set the ppath direction, making sure we are heading towards the nose tip.
         drn   <- sbst1$directions[ , 1, j]
         # drn   <- vec
         # dst1  <- rdist(t(sbst1$lmks["pn", ]), rbind(acR, acR + drn))
         dst1  <- rdist(t(sbst1$landmarks["pn", ]), rbind(ac, ac + drn))
         if (dst1[2] > dst1[1]) drn <- -drn
         # pnac  <- (sbst1$landmarks["pn", ] - ac)
         # pnse  <- (sbst1$landmarks["se", ] - sbst1$landmarks["pn", ])
         # nrm   <- sbst1$normals[j, ]
         # if ((acos(sum(drn * pnac) / sqrt(sum(pnac^2))) < pi / 4) &
         #     (abs(acos(sum(nrm * pnac) / sqrt(sum(pnac^2))) - pi / 2) < pi / 4)) {
            path <- planepath(sbst2, ac, direction = drn, directions = TRUE, rotation = 0)
            if (is.null(path)) 
               rdg <- NA
            else {
               path$path       <- path$path[path$arclength < 10, ]
               path$directions <- path$directions[ , , path$arclength < 10]
               cdst <- closestcurve.face3d(sbst2, path$path)
               ind  <- !(cdst$closest.curvept %in% c(1, nrow(path$path))) &
                        (abs(cdst$closest.distance) < 4)
               area <- area.face3d(sbst2)$points
               rdg  <- -sum(area[ind] * sbst2$kappa2[ind])
               # ang  <- abs(apply(t(sbst2$directions[ , 2, ind]) * sbst1$directions[ , 2, j], 1, sum))
               ang1 <- mean(abs(c(t(path$directions[ , 2, -1]) %*% sbst1$directions[ , 2, j])))
               ang  <- exp(-(acos(ang1)^2 / 0.3))

               if (monitor > 2 & ac.nm == "acR" & sbst1$kappa1[j] * rdg * ang > 0.4) {
                  # print(range(apply(t(path$directions[ , 2, -1]), 1, function(x) sum(x^2))))
                  # print(sum(sbst1$directions[ , 2, j]^2))
                  cat(ang1, acos(ang1), ang, "\n")
                  cat(j, sbst1$kappa1[j], rdg , ang, sbst1$kappa1[j] * rdg * ang, "\n")
                  if (!first) for (jj in 1:5) pop3d()
                  spheres3d(path$path)
                  spheres3d(ac, col = "yellow", radius = 2)
                  spheres3d(sbst2$vertices[ind, ], col = "yellow", radius = 0.5)
                  vec1 <- apply(t(path$directions[ , 2, -1]), 2, mean)
                  vec2 <- sbst1$directions[ , 2, j]
                  lines3d(rbind(ac + 10 * vec1, ac - 10 * vec1), col = "yellow")
                  lines3d(rbind(ac + 10 * vec2, ac - 10 * vec2), col = "red")
                  first <- FALSE
                  invisible(readline(prompt = " Press [enter] to continue"))
               }
            }
         # }
         # else
         #    rdg <- NA
            
         # ind1   <- ind & (cdst$closest.distance > 0)
         # ind2   <- ind & (cdst$closest.distance < 0)
         # rdg    <- sum(area[ind1] * sbst2$kappa2[ind1]) - sum(area[ind2] * sbst2$kappa2[ind2])

         # if (j %in% c(29, 82)) {
         # if (sbst1$kappa1[j] * rdg * ang > 1) {
         #    print(j)
         #    cat(j, sbst1$kappa1[j], rdg , ang, sbst1$kappa1[j] * rdg * ang, "\n")
         #    print(dim(path$path))
         #    angs <- c(t(path$directions[ , 2, -1]) %*% sbst1$directions[ , 2, j])
         #    boxplot(abs(angs))
         #    print(length(angs))
         #    if (j == 82) stop()
         # }
         results <- rbind(results, c(ac, sbst1$kappa1[j], rdg , ang1))
         crv  <- c(crv, sbst1$kappa1[j] * rdg * ang)
      }
   
      ind   <- which(!is.na(crv))
      sbstj <- subset(sbst1, jlst[ind], remove.singles = FALSE)
      ind1  <- ind[which.max(crv[ind])]
      # cat(jlst[ind1], crv[ind1], "\n")
      ac    <- sbstj$vertices[ind1, ]
      rnms  <- rownames(sbst$landmarks)
      sbst$landmarks <- rbind(sbst$landmarks, ac)
      rownames(sbst$landmarks) <- c(rnms, ac.nm)
      
      # Store the curvature of the local region
      nm <- paste(ac.nm, ".ind", sep = "")
      sbst[[nm]] <- sbst1$subset[jlst[ind]]
      nm <- paste(ac.nm, ".crv", sep = "")
      sbst[[nm]] <- crv[ind]
      
      plot(sbstj, col = crv[ind], display = "spheres", add = TRUE)
      spheres3d(ac, radius = 3)
   }

   # return(invisible(list(sbst = sbst, results = results)))
   return(invisible(sbst))
}
