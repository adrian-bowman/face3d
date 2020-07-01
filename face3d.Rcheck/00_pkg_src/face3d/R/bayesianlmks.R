bayesian.lmks <- function(lmks.initial, lmks.prior, graphics = FALSE, monitor = FALSE) {
   
   lmks    <- lmks.initial
   n.lmks  <- nrow(lmks)
   lmk.nms <- rownames(lmks)
   
   # Create prior distributions

   gpa    <- procGPA(lmks.prior[lmk.nms, , ])
   rownames(gpa$mshape) <- lmk.nms
   mnt    <- apply(gpa$tan, 1, mean)
   covt   <- cov(t(gpa$tan))

   if (graphics) {
      clear3d()
      mfrow3d(1, 2)
      view3d()
      plot(face, new = FALSE, col = "shape index")
      lmks.id <- spheres3d(lmks, col = "blue")
      scene.f <- currentSubscene3d()
      next3d()
      view3d()
      # for (i in 1:gpa$n) spheres3d(gpa$rotated[ , , i], radius = 0.2)
      for (i in 1:n.lmks) {
         ind <- c(i, i + n.lmks, i + 2 * n.lmks)
         plot3d(ellipse3d(covt[ind, ind], centre = gpa$mshape[i, ]), col = "yellow",
                alpha = 0.7, add = TRUE)
      }
      opa      <- procOPA(gpa$mshape, lmks)
      lmksp.id <- spheres3d(opa$Bhat, col = "blue")
      scene.p  <- currentSubscene3d()
   }

   del.lmks <- 1
   del.Bhat <- 1
   opa      <- procOPA(gpa$mshape, lmks)
   while (del.lmks > 0.1 | del.Bhat > 0.1) {
      lmks.old <- lmks
      opaB.old <- opa$Bhat
      for (lmk.nm in rownames(lmks)) {
         ind   <- which(edist.face3d(face$coords, lmks[lmk.nm, ]) < 2)
         post  <- numeric()
         for (i in ind) {
            lmks[lmk.nm, ] <- face$coords[i, ]
            opa    <- procOPA(gpa$mshape, lmks)
            args   <- c(c(opa$Bhat - gpa$mshape) %*% gpa$pcar)
            lprior <- sum(dnorm(args, 0, gpa$pcasd, log = TRUE))
            ll     <- 0
            for (j in 1:nrow(lmks)) {
               dst   <- edist.face3d(face$coords, lmks[j, ])
               nhd   <- which(dst < 10)
               x1    <- dst[nhd]
               y     <- (face$kappa1 * face$kappa2)[nhd]
   	         fit   <- lsfit(x1^2, y)
   	         beta  <- fit$coef["Intercept"]
   	         ll    <- ll - sum((fit$residuals)^2) / (2 * 1)
            }
   	      post  <- c(post, ll + lprior)
         }
         lmks[lmk.nm, ] <- face$coords[ind[which.max(post)], ]
         opa            <- procOPA(gpa$mshape, lmks)
         if (graphics) {
            useSubscene3d(scene.f)
            delFromSubscene3d(lmks.id)
            lmks.id  <- spheres3d(lmks, col = "blue")
            useSubscene3d(scene.p)
            delFromSubscene3d(lmksp.id)
            lmksp.id <- spheres3d(opa$Bhat, col = "blue")
         }
      }
      del.lmks <- sqrt(sum((lmks - lmks.old)^2))
      del.Bhat <- sqrt(sum((opa$Bhat - opaB.old)^2))
      if (monitor)
         cat("Change in lmks and matched lmks:", c(del.lmks, del.Bhat), "\n")
   }
   
   lmks
}
