findlandmarks.face3d <- function(face, lmks, orient = TRUE, monitor = FALSE, graphics = monitor,
                                 overwrite = FALSE, niter = 10) {

   if (missing(lmks)) lmks <- c("pn",  "acL",  "acR",  "sn",  "se",  "n", "exL", "exR", "enL", "enR",
                                "ls", "cphL", "cphR", "chL", "chR", "st",  "li",  "sl", "gn")

   if ("lmks" %in% names(face)) {
      ind <- which(is.na(match(lmks, rownames(face$lmks))))
      if (length(ind) > 0)
         face$lmks <- rbind(face$lmks, matrix(nrow = length(ind), ncol = 3, dimnames = list(lmks[ind])))
   }
   else
      face$lmks <- matrix(nrow = length(lmks), ncol = 3, dimnames = list(lmks))

   # lmks4 <- matrix(nrow = 4, ncol = 3, dimnames = list(c("pn", "enL", "enR", "se"), NULL))

   if (orient) {
      if (monitor) cat("orienting ... ")
      face <- orient.face3d(face)
      orient <- FALSE
   }
   else if (!("nearest" %in% names(face)))
      stop("orient is set to FALSE but face$nearest is not present.")

   quadlocate <- function(face, constraint, si, si.within, pair = FALSE, graphics = FALSE) {
      sbst    <- subset.face3d(face, constraint)
      indices <- as.numeric(rownames(sbst$coords))
      sbst    <- index.face3d(sbst)
      if (monitor.extra) display.face3d(sbst, colour = "shape index", new = FALSE)
      if (monitor.extra) scan()
       if (monitor.extra) display.face3d(sbst, colour = 2 + sbst$kappa1 * sbst$kappa2, new = FALSE)
      if (monitor.extra) scan()
      sbst    <- subset.face3d(sbst, abs(sbst$shape.index - si) < si.within)
      indices <- indices[as.numeric(rownames(sbst$coords))]
      # crv   <- sbst$kappa1 * sbst$kappa2
      # sbst  <- subset.face3d(sbst, crv > quantile(crv, 0.9))
      crv   <- sbst$kappa1 * sbst$kappa2
      if (monitor.extra) display.face3d(sbst, colour = 2 + crv, new = FALSE)
      if (monitor.extra) scan()
      parts <- connected.face3d(sbst)
      # ord   <- order(tapply(crv, parts, mean), decreasing = TRUE)
      # parts <- ord[parts]
      indices.all <- indices
      for (i in 1:(1 + as.numeric(pair))) {
         sbst1     <- subset.face3d(sbst, parts == i) 
         indices.i <- indices.all[as.numeric(rownames(sbst1$coords))]
         # sbst1 <- index.face3d(sbst1, overwrite = TRUE)
         # sbst1  <- subset.face3d(sbst1, abs(sbst1$shape.index - si) < si.within)
         dst   <- as.matrix(dist(sbst1$coords))
         crv   <- sbst1$kappa1
         crv   <- sbst1$kappa1 + sbst1$kappa2
         crv   <- sbst1$kappa1 * sbst1$kappa2
         rss.i <- apply(dst, 1, function(x) sum(lsfit(x^2, crv)$residuals^2))
         beta  <- apply(dst, 1, function(x) lsfit(x^2, crv)$coeff[2])
         # ind   <- (abs(si) > si.within) & (sign(beta) == sign(si))
         # rss.i[ind] <- max(rss.i)
         neg   <- (abs(si) > si.within)
         ind   <- if (neg) beta > 0 else beta < 0
         rss.i[ind] <- max(rss.i)
         # rss.i[ind] <- NA
         lmk.i <- sbst1$coords[which.min(rss.i), ]
         lmk   <- if (i == 1) lmk.i else rbind(lmk, lmk.i)
         rss   <- if (i == 1) rss.i else c(rss, rss.i)
         indices <- if (i == 1) indices.i else c(indices, indices.i)
         if (graphics) {
            # display.face3d(sbst1, new = FALSE, colour = "shape index", add = (i == 2))
            # scan()
            display.face3d(sbst1, colour = 2 + crv, new = FALSE, add = (i == 2))
            spheres3d(lmk, radius = 1.5, col = "green")
            if (monitor.extra) scan()
            # ind <- which.min(rss)
            # plot(dst[ind, ], crv)
         }
      }
      if (pair && (lmk[1, 1] < lmk[2, 1])) lmk <- lmk[2:1, ]
      invisible(list(lmk = lmk, indices = indices, rss = rss))
   }
   

   #---------------------------------------------------------------------------------
   #                       Initial location of pn, enL, enR, se
   #---------------------------------------------------------------------------------

      monitor.extra <- graphics
   
      if (any(is.na(face$lmks["pn", ]))) {
      if (monitor) cat("Locating initial pn ... ")
      nrst <- face$nearest
      nrst <- nrst[!is.na(nrst)] 
      nrst <- face$coords[nrst, ]
      ord  <- order(nrst[ , 3], decreasing = TRUE)
      nrst <- nrst[ord, ]
      flag <- FALSE
      i    <- 0
      while (!flag & (i < nrow(nrst))) {
         i    <- i + 1
         ind  <- apply(sweep(face$coords, 2, nrst[i, ]), 1, function(x) sqrt(sum(x^2))) < 20
         face <- index.face3d(face, subset = ind)
         crv  <- pmin(-face$kappa1[ind], -face$kappa2[ind])
         crv[face$shape.index[ind] < 0.8] <- 0
         if (length(which(crv > 0.1)) > 10) flag <- TRUE
         # display.face3d(face, colour = "shape index", new = FALSE)
      }
      if (!flag) stop("pn could not be identified.")                        
      face$lmks["pn", ] <- quadlocate(face, edist.face3d(face$coords, nrst[i, ]) < 30,
                                      1, 0.25, graphics = graphics)$lmk
      }
      
      if (any(is.na(c(face$lmks[c("enL", "enR"), ])))) {
      if (monitor) cat("enL/R ... ")
      face$lmks[c("enL", "enR"), ] <- quadlocate(face,
              (edist.face3d(face$coords, face$lmks["pn", ]) < 70) & (face$coords[ , 2] > face$lmks["pn", 2] + 10),
              -1, 0.3, pair = TRUE, graphics = graphics)$lmk
      # face$lmks["enL", ] <- quadlocate(face, edist.face3d(face$coords, face$lmks["enL", ]) < 20,
              # -1, 0.25, graphics = graphics)$lmk
      # face$lmks["enR", ] <- quadlocate(face, edist.face3d(face$coords, face$lmks["enR", ]) < 20,
              # -1, 0.25, graphics = graphics)$lmk
      }
      
      if (any(is.na(face$lmks["se", ]))) {
      if (monitor) cat("se ... ")
      curve <- planepath.face3d(face, face$lmks["enL", ], face$lmks["enR", ], boundary = c(0.2, 2))$path
      gcrv  <- gcurvature.face3d(curve, 4)
      ind   <- which.max(gcrv$gcurvature)
      face$lmks["se", ] <- gcrv$resampled.curve[ind, ]
      if (monitor.extra) {
         display.face3d(face, new = FALSE)
         spheres3d(gcrv$resampled.curve, col = "green")
         spheres3d(face$lmks["se", ], radius = 1.5, col = "blue")
         scan()
      }
      # se is a saddlepoint but the negative curvature will be stronger so go from -0.25 to ridge (0.5).
      face$lmks["se", ] <- quadlocate(face, edist.face3d(face$coords, face$lmks["se", ]) < 10,
                                                 0.25, 0.35, graphics = graphics)$lmk
      # face$lmks["se", ] <- quadlocate(face, edist.face3d(face$coords, face$lmks["se", ]) < 10,
                                                 # 0.125, 0.375, graphics = graphics)$lmk
      }
      
   # load("/Volumes/adrian/research/shape/Face3D/bayesian-faces/lmks.dmp")
   load("~/Desktop/bayesian-faces/lmks.dmp")
   ind        <- match(c("oiL", "oiR", "tL", "tR"), dimnames(lmks)[[1]])
   lmks       <- lmks[-ind, , ]
   # ind        <- match(c("pn", "enL", "enR", "se"), dimnames(lmks)[[1]])
   # lmks       <- lmks[ind, , ]
   if (!require(shapes)) stop("the shapes package is required.")
   gpa        <- procGPA(lmks, scale = FALSE)
   # gpa        <- procrustes.face3d(lmks, scale = FALSE)
   # gmean      <- gpa$mean
   gmean      <- gpa$mshape
   rownames(gmean) <- dimnames(lmks)[[1]]
   glmks      <- gpa$rotated
   k          <- dim(lmks)[1]
   n          <- dim(lmks)[3]
   X          <- matrix(c(aperm(sweep(glmks, 1:2, gmean), c(2, 1, 3))), nrow = n, byrow = TRUE)
   covmat     <- cov(X)


      
   if (monitor) cat("completed.\n")
      

   #---------------------------------------------------------------------------------
   #                       Find a suitable registration
   #---------------------------------------------------------------------------------
   
   lmks4       <- face$lmks[ c("pn", "enL", "enR", "se"), ]
   g4          <- gmean[c("pn", "enL", "enR", "se"), ]
   g4mean      <- apply(g4,    2, mean)
   lmean       <- apply(lmks4, 2, mean)
   opa         <- procOPA(g4, lmks4, scale = FALSE)
   face$coords <- sweep(sweep(face$coords, 2, lmean) %*% opa$R, 2, g4mean, "+")
   face$lmks   <- sweep(sweep(face$lmks,   2, lmean) %*% opa$R, 2, g4mean, "+")
   
   qd1 <- quadlocate(face, edist.face3d(face$coords, face$lmks["pn",  ]) < 10,    1, 0.25, graphics = FALSE)
   qd2 <- quadlocate(face, edist.face3d(face$coords, face$lmks["enL", ]) < 10,   -1, 0.30, graphics = FALSE)
   qd3 <- quadlocate(face, edist.face3d(face$coords, face$lmks["enR", ]) < 10,   -1, 0.30, graphics = FALSE)
   qd4 <- quadlocate(face, edist.face3d(face$coords, face$lmks["se",  ]) < 10, 0.25, 0.35, graphics = FALSE)
   face$lmks["pn",  ] <- qd1$lmk
   face$lmks["enL", ] <- qd2$lmk
   face$lmks["enR", ] <- qd3$lmk
   face$lmks["se",  ] <- qd4$lmk
   face$rss <- rep(NA, nrow(face$coords))
   face$rss[qd1$indices] <- qd1$rss   
   face$rss[qd2$indices] <- qd2$rss   
   face$rss[qd3$indices] <- qd3$rss   
   face$rss[qd4$indices] <- qd4$rss

   lmks4       <- face$lmks[c("pn", "enL", "enR", "se"), ]
   lmean       <- apply(lmks4, 2, mean)
   opa         <- procOPA(g4, lmks4, scale = FALSE)
   face$coords <- sweep(sweep(face$coords, 2, lmean) %*% opa$R, 2, g4mean, "+")
   face$lmks   <- sweep(sweep(face$lmks,   2, lmean) %*% opa$R, 2, g4mean, "+")
   
   ind1 <- (edist.face3d(face$coords, face$lmks["pn", ])  <= 40)
   ind2 <- (edist.face3d(face$coords, face$lmks["enL", ]) <= 20)
   ind3 <- (edist.face3d(face$coords, face$lmks["enR", ]) <= 20)
   ind4 <- (edist.face3d(face$coords, face$lmks["se", ])  <= 20)
   face <- index.face3d(face, subset = ind1 | ind2 | ind3 | ind4)
   
   for (nm in c("acL", "acR")) {
   	  i        <- match(nm, rownames(gmean))
      ind      <- (i - 1) * 3 + 1:3
      covmat.i <- covmat[ind, ind]
      crds.i   <- sweep(face$coords, 2, gmean[i, ])
      prior    <- exp(-rowSums((crds.i %*% solve(covmat.i)) * crds.i))
      face$lmks[nm, ] <- face$coords[which.max(prior), ]
   }
   
   qd5 <- quadlocate(face, edist.face3d(face$coords, face$lmks["acL", ]) < 10,
                                    -0.75, 0.25, graphics = TRUE)
   qd6 <- quadlocate(face, edist.face3d(face$coords, face$lmks["acR", ]) < 10,
                                    -0.75, 0.25, graphics = TRUE)

   face$lmks["acL", ]    <- qd5$lmk
   face$lmks["acR", ]    <- qd6$lmk
   face$rss[qd5$indices] <- qd5$rss   
   face$rss[qd6$indices] <- qd6$rss   

   # display.face3d(face, colour = "shape index", new = FALSE)
   # display.face3d(face, colour = 1 + face$kappa1 * face$kappa2, new = FALSE)
   # spheres3d(face$lmks, col = "blue", radius = 1.5)
   # print(face$lmks[c("acL", "acR"), ])

   # display.face3d(face, colour = "shape index", new = FALSE)
   ind  <- !is.na(face$rss)
   sbst <- subset.face3d(face, ind)
   sbst$rss <- face$rss[ind]
   # display.face3d(sbst, colour = 2 + sbst$rss, new = FALSE)
   # spheres3d(face$lmks, col = "yellow")
   # # scan()
   # display.face3d(face, colour = 1 + face$rss, new = FALSE)
   # spheres3d(face$lmks, col = "yellow")
   # for (nm in c("pn", "enL", "enR", "se")) {
   	  # i        <- match(nm, rownames(gmean))
      # ind      <- (i - 1) * 3 + 1:3
      # covmat.i <- covmat[ind, ind]
      # ell      <- ellipse3d(covmat.i, centre = gmean[i, ], level = 0.95)
      # plot3d(ell, add = TRUE, col = i, alpha = 1)
   # }
   
   rotate <- function(par, face, lmks.ref, gmean, covmat, graphics = FALSE, optim = TRUE) {
   	  mn.ref   <- apply(gmean[lmks.ref, ], 2, mean)
   	  crds     <- face$coords
   	  crds     <- sweep(crds, 2, mn.ref)
      crds     <- rotate3d(crds, par[1], 1, 0, 0)
      crds     <- rotate3d(crds, par[2], 0, 1, 0)
      crds     <- rotate3d(crds, par[3], 0, 0, 1)
   	  crds     <- sweep(crds, 2, mn.ref, "+")
   	  crds     <- sweep(crds, 2, par[4:6], "+")
      face$coords <- crds
      if (graphics) display.face3d(face, colour = 1 + face$rss, new = FALSE)
   	  logpost <- 0
   	  
      for (nm in lmks.ref) {
   	     i        <- match(nm, rownames(gmean))
         ind      <- (i - 1) * 3 + 1:3
         covmat.i <- covmat[ind, ind]
         crds.i   <- sweep(crds, 2, gmean[nm, ])
      	 # logprior <- rowSums((crds.i %*% solve(covmat.i)) * crds.i)
      	 # rss      <- face$rss
      	 # rss[is.na(rss)] <- 1e8
      	 # logpost  <- logpost + sum(logprior + rss, na.rm = TRUE)
      	 prior    <- exp(-rowSums((crds.i %*% solve(covmat.i)) * crds.i))
      	 like     <- exp(-face$rss)
      	 like[is.na(like)] <- 0
      	 logpost  <- logpost + log(sum(prior * like))
         if (graphics) {
            ell      <- ellipse3d(covmat.i, centre = gmean[nm, ], level = 0.95)
            plot3d(ell, add = TRUE, col = i, alpha = 1)
         }
      }
      # print(c(par, logpost))
      result <- if (optim) logpost else face
      invisible(result)
   }
   
   lmks.ref <- c("pn", "enL", "enR", "se", "acL", "acR")
   rotate(rep(0, 6), face, lmks.ref, gmean, covmat, graphics = TRUE)
   opt <- optim(rep(0, 6), rotate, face = face, lmks.ref = lmks.ref, gmean = gmean, covmat = covmat,
                method = "L-BFGS-B", lower = c(rep(-0.2, 3), rep(-5, 3)), upper = c(rep(0.2, 3), rep(5, 3)),
                control = list(trace = 1, fnscale = -1, maxit = 50))
   # print(opt$par)
   face <- rotate(opt$par, face, lmks.ref, gmean, covmat, graphics = TRUE, optim = FALSE)
   # for (i in 1:k) {
      # ind      <- (i - 1) * 3 + 1:3
      # covmat.i <- covmat[ind, ind]
      # ell      <- ellipse3d(covmat.i, centre = gmean[i, ], level = 0.95)
      # plot3d(ell, add = TRUE, col = i, alpha = 1)
   # }
   
   #---------------------------------------------------------------------------------
   #                       Find the four initial landmarks
   #---------------------------------------------------------------------------------
   
   for (nm in lmks.ref) {
   	  i        <- match(nm, rownames(gmean))
      ind      <- (i - 1) * 3 + 1:3
      covmat.i <- covmat[ind, ind]
      crds.i   <- sweep(face$coords, 2, gmean[nm, ])
   	  prior    <- exp(-rowSums((crds.i %*% solve(covmat.i)) * crds.i))
      like     <- exp(-face$rss)
      like[is.na(like)] <- 0
   	  logpost  <- log(prior * like)
   	  ind      <- which.max(logpost)
      face$lmks[nm, ] <- face$coords[ind, ]
   }

   

   
   #---------------------------------------------------------------------------------
   #                       Find the other landmarks
   #---------------------------------------------------------------------------------
   
   
      
   return(invisible(face))
   
   #---------------------------------------------------------------------------------
   #                       locate the others - older code
   #---------------------------------------------------------------------------------
   
   # load("/Volumes/adrian/research/shape/Face3D/bayesian-faces/lmks.dmp")
   ind        <- match(c("oiL", "oiR", "tL", "tR"), dimnames(lmks)[[1]])
   lmks       <- lmks[-ind, , ]
   gpa        <- procrustes.face3d(lmks, scale = FALSE)
   glmks      <- gpa$rotated
   gmean      <- gpa$mean
   gmean.mean <- apply(gmean, 2, mean)
   k          <- dim(lmks)[1]
   n          <- dim(lmks)[3]
   X          <- matrix(c(aperm(sweep(glmks, 1:2, gmean), c(2, 1, 3))), nrow = n, byrow = TRUE)
   covmat     <- cov(X)

   if (!require(shapes)) stop("the shapes package is required.")
   g4          <- gpa$mean[c("pn", "enL", "enR", "se"), ]
   g4mean      <- apply(g4, 2, mean)
   lmean       <- apply(lmks4, 2, mean)
   opa         <- procOPA(g4, lmks4, scale = FALSE)
   face$coords <- sweep(sweep(face$coords, 2, lmean) %*% opa$R, 2, g4mean, "+")

   if (graphics) display.face3d(face, new = FALSE)
   crds        <- face$coords
   lmks        <- gmean
   # nrmls       <- face$normals
   # spheres3d(lmks, radius = 1, col = "red")
   for (iter in 1:niter) {
      R          <- procOPA(gmean, lmks, scale = FALSE)$R
      lmks.mean  <- apply(lmks, 2, mean)
      lmks       <- sweep(sweep(lmks, 2, lmks.mean) %*% R, 2, gmean.mean, "+")
      crds       <- sweep(sweep(crds, 2, lmks.mean) %*% R, 2, gmean.mean, "+")
      # nrmls      <- nrmls %*% R
      logpost    <- 0
      for (i in 1:k) {
         ind       <- (i - 1) * 3 + 1:3
         covmat.i  <- covmat[ind, ind]
         crds      <- sweep(crds, 2, gmean[i, ])
         dst       <- rowSums((crds %*% solve(covmat.i)) * crds)
         crds      <- sweep(crds, 2, gmean[i, ], "+")
         ind       <- (dst <= 3)
         cat("iteration", iter, rownames(lmks)[i])
         face      <- index.face3d(face, subset = ind)
         cat(" finished\n")
         logprior  <- -dst
         loglik    <- log(pmax(abs(face$kappa1), abs(face$kappa2)))
         logpost.i <- loglik + logprior
         lmks[i, ] <- crds[which.max(logpost.i), ]
         logpost   <- logpost + logpost.i
      }
      if (graphics) {
         face$coords <- crds
         display.face3d(face, new = FALSE)
         spheres3d(lmks, radius = 1, col = "red")
      }
   }

   face$coords <- crds
   face$lmks   <- lmks
   invisible(face)
}


   #---------------------------------------------------------------------------------
   #                       locate acL, acR
   #---------------------------------------------------------------------------------
   
         # ind   <- apply(face$coords, 1, function(x) sqrt(sum((x - face$lmks["pn", ])^2))) < 70
         # face  <- index.face3d(face, distance = 10, subset = ind)
         # ind   <- (face$shape.index < -0.25) & abs(face$coords[ , 2] - face$lmks["pn", 2]) < 20 &
                       # apply(face$coords, 1, function(x) sqrt(sum((x - face$lmks["pn", ])^2))) < 50
         # sbst  <- subset.face3d(face, ind)
         # parts <- connected.face3d(sbst)
         # ind   <- as.numeric(rownames(sbst$coords)[parts %in% 1:2])
         # sbst  <- subset.face3d(sbst, parts %in% 1:2)
         # crv   <- pmin(sbst$kappa1, sbst$kappa2)
         # brks  <- seq(0.0, 0.2, length = 21)
         # brks  <- c(min(crv) - 1, quantile(crv, seq(0.05, 0.95, 0.05)), max(crv) + 1)
         # parts <- parts[parts %in% 1:2]
         # sbst1 <- subset.face3d(sbst, parts == 1)
         # crv1  <- crv[parts == 1]
         # alL   <- sbst1$coords[which.max(crv1), ]
         # sbst2 <- subset.face3d(sbst, parts == 2)
         # crv2  <- crv[parts == 2]
         # alR   <- sbst2$coords[which.max(crv2), ]
         # al    <- if (alL[1] > alR[1]) cbind(alL, alR) else cbind(alR, alL)
         
         # ind   <- apply(face$coords, 1, function(x) sqrt(sum((x - face$lmks["pn", ])^2))) < 40
         # face  <- index.face3d(face, distance = 10, subset = ind)
         # ind   <- ind & (face$shape.index > 0.25)
         # sbst  <- subset.face3d(face, ind)
         # sbst  <- subset(sbst, connected.face3d(sbst) == 1)

         # if ("acL" %in% lmks) face$lmks["acL", ] <- al[ , 1]
         # if ("acR" %in% lmks) face$lmks["acR", ] <- al[ , 2]


