initiallandmarks.face3d <- function(face, lmk.names, orient = TRUE, monitor, graphics = FALSE,
                                    overwrite = FALSE, niter = 10) {

   if (missing(lmk.names)) lmk.names <- c("pn", "enL", "enR", "se")
   if (all(lmk.names %in% rownames(face$lmks)) & !overwrite) {
      cat("All requested landmarks are already present - use the overwrite argument if required.")
      return(invisible(face))
   }
   lmks <- matrix(nrow = 0, ncol = 3)
   if ("lmks" %in% names(face)) {
      ind <- rownames(face$lmks) %in% lmk.names
      if (any(!ind))             lmks <- rbind(lmks, face$lmks[!ind, ])
      if (any(ind) & !overwrite) lmks <- rbind(lmks, face$lmks[ ind, ])
   }
   if (missing(monitor)) monitor <- 0

   # if ("lmks" %in% names(face)) {
   #    ind <- which(is.na(match(lmks, rownames(face$lmks))))
   #    if (length(ind) > 0)
   #       face$lmks <- rbind(face$lmks, matrix(nrow = length(ind), ncol = 3, dimnames = list(lmks[ind])))
   # }
   # else
      # face$lmks <- matrix(nrow = length(lmks), ncol = 3, dimnames = list(lmks))

   # lmks4 <- matrix(nrow = 4, ncol = 3, dimnames = list(c("pn", "enL", "enR", "se"), NULL))

   if (graphics != FALSE & graphics != "add") plot(face)
   if (graphics == "add") graphics <- TRUE

   if (monitor > 0) cat("Locating approximate landmarks: ...\n")
   if (graphics) plot(face, new = FALSE)

   #---------------------------------------------------------------------------------
   #                       locate pn
   #---------------------------------------------------------------------------------
   
   if (("pn" %in% lmk.names)) {
      if (!(("lmks" %in% names(face)) && ("pn" %in% rownames(face$lmks)) && !overwrite)) {

      if (monitor > 0) cat("  pn ...\n")

      if (monitor > 1) cat("    Finding convex hull points ...\n")
      chull <- convhulln(face$coords)
      cpts  <- unique(c(chull))
      if (graphics) spheres3d(face$coords[cpts, ], col = "green")
      
      # if (monitor > 1) cat("    Removing points near the edges ...\n")
      # epts  <- unique(c(edges.face3d(face)))
      # dst   <- rdist(face$coords[cpts, ], face$coords[epts, ])
      # ind   <- apply(dst, 1, function(x) all(x > 10))
      # if (graphics) {
      #    p3d <- par3d(skipRedraw = TRUE)
      #    pop3d()
      #    spheres3d(face$coords[cpts, ], col = "green")
      #    par3d(p3d)
      # }

      if (monitor > 1) cat("    Removing isolated points ...\n")
      dst  <- rdist(face$coords[cpts, ])
      ind  <- apply(dst, 1, function(x) length(which(x < 10)) > 10)
      cpts <- cpts[ind]
      if (graphics) {
         p3d <- par3d(skipRedraw = TRUE)
         pop3d()
         spheres3d(face$coords[cpts, ], col = "green")
         par3d(p3d)
      }

      if (monitor > 1) cat("    Locating point of maximum density and Gaussian curvature ...\n")
      face <- index.face3d(face, distance = 5, subset = cpts)
      gc   <- face$kappa1[cpts] * face$kappa2[cpts]
      clr  <- topo.colors(20)[cut(gc, 20, labels = FALSE)]
      if (graphics) {
         p3d <- par3d(skipRedraw = TRUE)
         pop3d()
         spheres3d(face$coords[cpts, ], col = clr)
         par3d(p3d)
      }

      dst   <- rdist(face$coords[cpts, ], face$coords)
      val   <- apply(dst, 1, function(x) length(which(x < 5)))
      dst   <- rdist(face$coords[cpts, ], face$coords[cpts, ])
      val   <- apply(dst, 1, function(x) length(which(x < 5))) / val
      val   <- apply(dst, 1, function(x) mean(gc[x < 5])) * val
      pn    <- face$coords[cpts, ][which.max(val), ]
      lmks  <- rbind(lmks, pn = pn)
      if (graphics) spheres3d(pn, col = "red", radius = 2)
   }
   else cat("\npn is already present: use the overwrite argument if required.\n")
   }

   
   #---------------------------------------------------------------------------------
   #               Trim the image to focus on the face
   #---------------------------------------------------------------------------------
   
   # tface <- subset(face, edist.face3d(face$coords, lmks["pn", ]) < 100)
   
   #---------------------------------------------------------------------------------
   #                       locate enL, enR
   #---------------------------------------------------------------------------------
   
   if (all(c("enR", "enL") %in% lmk.names)) {
      if (("lmks" %in% names(face)) && all(c("enR", "enL") %in% rownames(face$lmks)) && !overwrite)
         cat("\nenL/R are already present: use the overwrite argument if required.\n")
      else {
         if ("pn" %in% rownames(face$lmks)) lmks <- rbind(lmks, pn = face$lmks["pn", ])
   	   if (!("pn" %in% rownames(lmks)))
   	      cat("enR/L cannot be identified because pn is needed.\n")
   	   else {
   	      if (monitor > 0) cat("  enL/R ...\n")
   	      ed    <- edist.face3d(face$coords, lmks["pn", ])
   	      # ind   <- (ed < 70) & (ed > 40) & (face$coords[ , 2] > lmks["pn", 2])
   	      ind   <- (ed < 70)
   	      face  <- index.face3d(face, subset = ind, distance = 5, overwrite = TRUE)
            sbst  <- subset.face3d(face, ind & face$shape.index < -0.7)
            parts <- connected.face3d(sbst)
            sbst1 <- subset.face3d(sbst, parts == 1)
            sbst2 <- subset.face3d(sbst, parts == 2)
            gc1   <- sbst1$kappa1 * sbst1$kappa2
            gc2   <- sbst2$kappa1 * sbst2$kappa2
            # dst    <- as.matrix(dist(sbst1$coords))
            # rss    <- apply(dst, 1, function(x) sum(lsfit(x^2, gc1)$residuals^2))
            # enL    <- sbst1$coords[which.min(rss), ]
            # dst    <- as.matrix(dist(sbst2$coords))
            # rss    <- apply(dst, 1, function(x) sum(lsfit(x^2, gc2)$residuals^2))
            # enR    <- sbst2$coords[which.min(rss), ]
            # en     <- if (enL[1] < enR[1]) rbind(enR, enL) else rbind(enL, enR)
            h     <- 2
            dst   <- rdist(sbst1$coords)
            val   <- apply(dst, 1, function(x) sum(gc1 * exp(-0.5 *x^2 / h^2)) / sum(exp(-0.5 *x^2 / h^2)))
            enL   <- sbst1$coords[which.max(val), ]
            dst   <- rdist(sbst2$coords)
            val   <- apply(dst, 1, function(x) sum(gc2 * exp(-0.5 *x^2 / h^2)) / sum(exp(-0.5 *x^2 / h^2)))
            enR   <- sbst2$coords[which.max(val), ]
            en    <- rbind(enL, enR)
            if (graphics) {
               plot(sbst1, display = "spheres", col = gc1, add = TRUE)
               plot(sbst2, display = "spheres", col =gc2, add = TRUE)
               spheres3d(en, radius = 2, col = "red")
            }
            lmks <- rbind(lmks, enL = en[1, ], enR = en[2, ])
   	   }
      }
   }

   #---------------------------------------------------------------------------------
   #                       locate se
   #---------------------------------------------------------------------------------
   
   if ("se" %in% lmk.names) {
      if (!(("lmks" %in% names(face)) && ("se" %in% rownames(face$lmks)) && !overwrite)) {

      if (monitor > 0) cat("  se ... ")
   	if (!all(c("enL", "enR") %in% rownames(lmks)))
   	   cat("se cannot be identified because enL/R are both needed.\n")
      else {
         curve <- planepath.face3d(face, lmks["enL", ], lmks["enR", ], boundary = c(0.2, 2))$path
         gc    <- gcurvature.face3d(curve, 3)
         lmks  <- rbind(lmks, se = gc$pos.max)

         if (graphics) {
            plot(gc$arclength, gc$gcurvature)
            sbst <- subset(face, edist.face3d(face$coords, lmks["se", ]) < 15)
            sbst <- index.face3d(sbst, distance = 10, overwrite = TRUE)
            plot(sbst, colour = "shape index", add = TRUE)
            spheres3d(curve, radius = 0.7)
            spheres3d(lmks["se", ], col = "red", radius = 2)
         }
      }
      
      # clr   <- -sbst$kappa1 * sbst$kappa2
      # # clr    <- sbst$kappa1 - sbst$kappa2
      # # clr    <- sbst$kappa1
      # # clr    <- -sbst$kappa2
      # # clr    <- sbst$shape.index
      # clr[sbst$shape.index < -0.25] <- 0
      # clr[edist.face3d(sbst$coords, face$lmks[lmk.name, ]) > 20] <- 0
      # # clr    <- pmax(clr, quantile(clr, 0.98))
      # # display.face3d(sbst, colour = 100 * clr, new = FALSE)
      # # spheres3d(sbst$lmks[lmk.name, ], col = "red")
      # ind   <- (clr > 0)
      # sbst1 <- subset(sbst, ind, remove.singles = FALSE)
      # clr   <- clr[ind]
      # ind   <- (clr > quantile(clr, 0.75))
      # sbst1 <- subset(sbst1, ind, remove.singles = FALSE)
      # clr   <- clr[ind]
      # se    <- mean.face3d(sbst1, clr)
      # lmks <- rbind(lmks, se)
      # rownames(lmks)[nrow(lmks)] <- "se"
      # if (graphics) {
      #    display.face3d(face, new = FALSE)
      #    spheres3d(sbst$lmks[lmk.name, ], radius = 3, col = "red")
      #    spheres3d(lmk, radius = 3)
      # }
      
   }
   else cat("\nse already present: use the overwrite argument if required.\n")
   }

   
   #---------------------------------------------------------------------------------
   #                       locate sn
   #---------------------------------------------------------------------------------
   
   if ("sn" %in% lmk.names) {
      if (!(("lmks" %in% names(face)) && ("pn" %in% rownames(face$lmks)) && !overwrite)) {
      
      if (monitor > 0) cat("sn ... ")
         
      # Warp the template curves onto the face using a few landmarks
      template   <- template_male
      nms        <- c("pn", "se", "enL", "enR")
      ind        <- match(nms, rownames(template$lmks))
      mat1       <- rbind(template$lmks[ind, ], template$curves, template$curves + template$cnormals)
      wp         <- warp.face3d(mat1, lmks[nms, ], subset = (length(ind) + 1):nrow(mat1), general = TRUE, project = FALSE)
      wp.curves  <- wp[1:nrow(template$curves), ]
      wp.normals <- wp[nrow(template$curves) + (1:nrow(template$curves)), ] - wp.curves
      wp.normals <- wp.normals / apply(wp.normals, 1, function(x) sqrt(sum(x^2)))
   
      # Project the template onto the image by moving along the curve normals
      cname <- "mid-line columella"
      ind   <- which(substr(rownames(template$curves), 1, nchar(cname)) == cname)
      wp.c  <- wp.curves[ind, ]
      wp.n  <- wp.normals[ind, ]
      proj  <- project.face3d(wp.c, wp.n, face)
      lmks  <- rbind(lmks, sn = proj[nrow(proj), ])

      if (graphics) {
         spheres3d(wp.curves, radius = 1.5, col = "red")
         lns <- matrix(c(t(cbind(wp.curves, wp.curves + 5 * wp.normals))), ncol = 3, byrow = TRUE)
         segments3d(lns, col = "red")
         spheres3d(proj, col = "blue")
         spheres3d(lmks["sn", ], col = "yellow", radius = 1.5)
      }
   
   }
   else cat("\nse already present: use the overwrite argument if required.\n")
   }
   
   # Old code
   # if (("sn" %in% lmk.names) && (!("sn" %in% rownames(face$lmks)) | overwrite)) {
   # 
   #    ind    <- edist.face3d(face$coords, lmks["pn", ]) < 50
   #    face   <- index.face3d(face, distance = 10, subset = ind)
   #    ind    <- (face$shape.index < -0.65) &
   #                apply(face$coords, 1, function(x) sqrt(sum((x - lmks["pn", ])^2))) > 30
   #    sbst   <- subset.face3d(face, ind)
   #    parts  <- connected.face3d(sbst)
   #    sbst   <- subset.face3d(sbst, parts == 1)
   #    if (graphics) {
   #       display.face3d(sbst, colour = "shape index", add = TRUE)
   #       spheres3d(lmks["se", ], radius = 1.5)
   #    }
   #    
   # }
   
   #---------------------------------------------------------------------------------
   #                       locate acL, acR
   #---------------------------------------------------------------------------------
   
   if (any(c("acR", "acL") %in% lmk.names) && (!any(c("acR", "acL") %in% rownames(face$lmks)) | overwrite)) {
      lmk.pn <- if ("pn" %in% rownames(face$lmks)) face$lmks["pn", ] else
                if ("pn" %in% rownames(lmks)) lmks["pn", ] else NULL
      lmk.se <- if ("se" %in% rownames(face$lmks)) face$lmks["se", ] else
                if ("se" %in% rownames(lmks)) lmks["se", ] else NULL
      if (is.null(lmk.pn) | is.null(lmk.se))
   	   cat("acR/L cannot be identified because pn and se are needed.\n")
   	else {
   	   k1     <- face$kappa1
   	   mp     <- (lmk.se + lmk.pn) / 2
   	   proj   <- sweep(face$coords, 2, mp) %*% (lmk.se - mp)
   	   k1[proj > 0] <- NA
   	   p      <- 0.8
   	   parts  <- 1:2
   	   while (length(table(parts)) < 3) {
   	      ind    <- k1 > quantile(k1, p, na.rm = TRUE)
   	      sbst   <- subset(face, ind)
   	      parts  <- connected.face3d(sbst)
   	      p      <- p - 0.05
   	   }
   	   sbst1  <- subset(sbst, parts == 1)
   	   sbst2  <- subset(sbst, parts == 2)
   	   sbst3  <- subset(sbst, parts == 3)
   	   dist1  <- min(rdist(sbst1$coords, t(lmk.pn)))
   	   dist2  <- min(rdist(sbst2$coords, t(lmk.pn)))
   	   dist3  <- min(rdist(sbst3$coords, t(lmk.pn)))
   	   ord    <- order(c(dist1, dist2, dist3))
   	   nose   <- switch(match(1, ord), sbst1, sbst2, sbst3)
   	   mouth  <- switch(match(2, ord), sbst1, sbst2, sbst3)
   	   sulcus <- switch(match(3, ord), sbst1, sbst2, sbst3)
   	   nose   <- index.face3d(nose,   distance = 3, directions = TRUE, overwrite = TRUE)
   	   mouth  <- index.face3d(mouth,  distance = 3, directions = TRUE, overwrite = TRUE)
   	   sulcus <- index.face3d(sulcus, distance = 3, directions = TRUE, overwrite = TRUE)
   	   noseR  <- subset(nose, edist.face3d(nose$coords, face$lmks["acR", ]) < 10)
   	   noseL  <- subset(nose, edist.face3d(nose$coords, face$lmks["acL", ]) < 10)
   	   plot(nose,   display = "spheres", add = TRUE, col = nose$kappa1)
   	   plot(mouth,  display = "spheres", add = TRUE, col = mouth$kappa1)
   	   # plot(sulcus, display = "spheres", add = TRUE, col = sulcus$kappa1)
   	   # plot(noseL,   display = "spheres", add = TRUE, col = noseL$kappa1)
   	   # plot(noseR,   display = "spheres", add = TRUE, col = noseR$kappa1)
   	   return()
   	   
   	   for (side in c("R", "L")) {
   	      
   	   noseRL <- switch(side, R = noseR, L = noseL)
   	      
   	   # Focus on areas of high kappa1 which are not on edges.  Order these points.
   	   # plot(noseRL, col = "shape index", display = "mesh", new = FALSE)
   	   edge.pts <- unique(c(edges.face3d(noseRL)))
   	   ind      <- (1:nrow(noseRL$coords))[-c(edge.pts)]
   	   ind      <- ind[noseRL$kappa1[ind] > quantile(noseRL$kappa1, 0.7, na.rm = TRUE)]
   	   ord      <- order(noseRL$kappa1[ind], decreasing = TRUE)
   	   if (graphics) {
   	      plot(noseRL, col = noseRL$kappa1, display = "mesh", new = FALSE)
   	      spheres3d(noseRL$coords[ind, ], radius = 0.15)
   	   }
   	   # Find a regularly spaced set of points with high kappa1
   	   # # Use rdist to decide on a threshold separation distance (currently 2)
   	   # # rd       <- rdist()
   	   # tops     <- ind[ord[1]]
   	   # for (j in ind[ord[-1]])
   	   #    if (all(rdist(t(noseRL$coords[j, ]), noseRL$coords[tops, ]) > 2)) tops <- c(tops, j)
   	   # plot(noseRL, col = noseRL$kappa1, display = "mesh", new = FALSE)
   	   # spheres3d(noseRL$coords[tops, ], col = "red", radius =0.15)
   	   
   	   # Search over these points using planepath to maximise the integral of kappa1
   	   # The planepaths are too variable
   	   # plot(noseRL, col = noseRL$kappa1, display = "mesh", new = FALSE)
   	   # for (j in tops) {
   	   #    drn <- c(crossproduct(noseRL$normals[j, ], noseRL$directions[ , 1, j]))
   	   #    spheres3d(noseRL$coords[j, ], radius = 0.2, col = "blue")
   	   #    lines3d(rbind(noseRL$coords[j, ], noseRL$coords[j, ] + drn), col = "blue")
   	   #    pp  <- planepath.face3d(noseRL, noseRL$coords[j, ], direction = drn, si.target = -1,
   	   #                            rotation.range = pi/4, boundary = c(100, 100), graphics = FALSE)
   	   #    spheres3d(pp$path, radius = 0.15)
   	   #    scan()
   	   #    for (k in 1:2) pop3d()
   	   # }
   	   
   	   # Move each point to the maximum of kappa1 along the principal direction
   	   # save(noseRL, ind, ord, file = "~/Desktop/noseRL.Rdata")
   	   if (graphics) plot(noseRL, col = noseRL$kappa1, display = "mesh", new = FALSE)
   	   roll <- function(starts, dirns) {
   	      valley <- matrix(nrow = 0, ncol = 3)
   	      for (j in 1:nrow(starts)) {
   	         # cat(j, "")
   	         pp      <- planepath.face3d(noseRL, starts[j, ], rotation = 0,
                                           direction = dirns[ , j], directions = TRUE,
                                           boundary = c(100, 100), graphics = FALSE)
   	         crv     <- gcurvature.face3d(pp$path, 6)
   	         al      <- cumsum(c(0, sqrt(rowSums(diff(pp$path)^2))))
   	         si      <- approx(al, pp$shape.index, n = length(crv$gcurvature))$y
   	         ind.crv <- which.max(crv$gcurvature * (-sign(si)))
   	         valley  <- rbind(valley, crv$resampled.curve[ind.crv, ])
   	         
   	         # apply(rdist(pp$path), 1, function(x) sum(x == 0))
   	         # col <- topo.colors(20)[cut(crv$gcurvature, 20, labels = FALSE)]
   	         # open3d()
   	         # spheres3d(crv$resampled.curve, radius = 0.01, col = gcol)
   	         # for (jj in 1:200) spheres3d(crv$resampled.curve[jj, ], radius = 0.01)
   	         if (graphics) {
   	            plot(crv$arc.length, crv$gcurvature * (-sign(si)))
   	            points(crv$arc.length[ind.crv], crv$gcurvature[ind.crv], pch = 16, col = "red")
   	            # spheres3d(crv$resampled.curve[ind.crv, ], col = "red", radius = 0.15)
   	            spheres3d(valley[nrow(valley), ], col = "green", radius = 0.15)
   	            spheres3d(pp$path, radius = 0.1)
   	            spheres3d(noseRL$coords[j, ], radius = 0.12, col = "blue")
   	            lines3d(rbind(noseRL$coords[j, ], noseRL$coords[j, ] + noseRL$directions[ , 1, j]), col = "blue")
   	            # scan()
   	            for (k in 1:3) pop3d()
   	         }
   	      }
   	      invisible(valley)
   	   }
   	   valley <- roll(noseRL$coords[ind[ord], ], noseRL$directions[ , 1, ind[ord]])

   	   # Find the order of the points in arclength simply by the next closest
   	   ind    <- which.max(valley[ , 2])
   	   pts    <- matrix(valley[ind, ], nrow = 1)
   	   valley <- valley[-ind, ]
   	   dst    <- c(rdist(t(pts[nrow(pts), ]), valley))
   	   # Find a starting point which is not isolated
   	   while(min(dst) > 2) {
   	      ind    <- which.min(dst)
   	      pts    <- matrix(valley[ind, ], nrow = 1)
   	      valley <- valley[-ind, ]
   	      dst    <- c(rdist(t(pts[nrow(pts), ]), valley))
   	   }
   	   # spheres3d(pts[nrow(pts), ], col = "blue", radius = 0.2)
   	   for (j in 2:nrow(valley)) {
   	      if (!is.matrix(valley)) valley <- t(valley)
   	      dst      <- c(rdist(t(pts[nrow(pts), ]), valley))
   	      ind      <- which.min(dst)
   	      if (dst[ind] < 2) {
   	         pts <- rbind(pts, valley[ind, ])
   	         # spheres3d(pts[nrow(pts), ], col = "blue", radius = 0.2)
   	      }
   	      valley <- valley[-ind, ]
   	   }
   	   gc  <- gcurvature.face3d(pts, df = 4)
   	   ind <- which(c(0, diff(sign(diff(gc$gcurvature))), 0) < 0)[1]
   	   ac  <- gc$resampled.curve[ind, ]
   	   if (graphics) {
   	      plot(noseRL, col = noseRL$kappa1, display = "mesh", new = FALSE)
   	      spheres3d(valley, col = "green", radius = 0.15)
   	      spheres3d(pts, col = "blue", radius = 0.2)
   	      plot(gc$arc.length, gc$gcurvature, col = "blue")
   	      points(gc$arc.length[ind], gc$gcurvature[ind], pch = 16, col = "red")
   	      spheres3d(ac, col = "red", radius = 0.3)
   	   }
   	   
         # ind   <- apply(face$coords, 1, function(x) sqrt(sum((x - face$lmks["pn", ])^2))) < 30
         # face  <- index.face3d(face, distance = 10, subset = ind)
         # ind   <- (face$shape.index < -0.25) & abs(face$coords[ , 2] - face$lmks["pn", 2]) < 10 &
         #             apply(face$coords, 1, function(x) sqrt(sum((x - face$lmks["pn", ])^2))) < 30
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
         # acL   <- sbst1$coords[which.max(crv1), ]
         # sbst2 <- subset.face3d(sbst, parts == 2)
         # crv2  <- crv[parts == 2]
         # ac   <- sbst2$coords[which.max(crv2), ]

         # ind   <- apply(face$coords, 1, function(x) sqrt(sum((x - face$lmks["pn", ])^2))) < 40
         # face  <- index.face3d(face, distance = 10, subset = ind)
         # ind   <- ind & (face$shape.index > 0.25)
         # sbst  <- subset.face3d(face, ind)
         # sbst  <- subset(sbst, connected.face3d(sbst) == 1)

         if (graphics) {
         	 # clr <- topo.colors(20)[cut(crv1, brks, labels = FALSE)]
             # spheres3d(sbst1$coords, col = clr)
         	 # clr <- topo.colors(20)[cut(crv2, brks, labels = FALSE)]
             # spheres3d(sbst2$coords, col = clr)
            # plot(noseRL, new = FALSE, colour = noseRL$kappa1)
            spheres3d(ac, radius = 0.3, col ="red")
            # spheres3d(rbind(acL, ac), radius = 1.5)
         }
         # lmks <- rbind(lmks, acL = acL)
   	   lmks <- if (side == "R") rbind(lmks, acR = ac) else rbind(lmks, acL = ac)
   	   if (side == "R") face$curve.acR <- pts else face$curve.acL <- pts
   	   }
      }
   }
   
   # Put the located landmarks into the face object
   ind <- rownames(lmks) %in% rownames(face$lmks)
   nms <- rownames(lmks)[ind]
   if (length(which( ind)) > 0) face$lmks[nms, ] <- lmks[ind, ]
   if (length(which(!ind)) > 0)
      face$lmks <- rbind(face$lmks, lmks[!ind, , drop = FALSE])

   invisible(face)
}

#---------------------------------------------------------------------------------
#                                    Old code
#---------------------------------------------------------------------------------

# pn using the orientation idea
#------------------------------

# if (orient) {
#    if (monitor) cat("Orienting: ... ")
#    face <- orient.face3d(face)
#    if (monitor) cat("complete\n")
# }
# else if (!("nearest" %in% names(face)))
#    stop("orient is set to FALSE but face$nearest is not present.")

# nrst <- face$nearest
# nrst <- nrst[!is.na(nrst)] 
# nrst <- face$coords[nrst, ]
# ord  <- order(nrst[ , 3], decreasing = TRUE)
# nrst <- nrst[ord, ]
# flag <- FALSE
# i <- 0
# while (!flag & (i < nrow(nrst))) {
#    i    <- i + 1
#    ind  <- edist.face3d(face$coords, nrst[i, ]) < 20
#    face <- index.face3d(face, subset = ind)
#    crv  <- pmin(-face$kappa1[ind], -face$kappa2[ind])
#    # crv[face$shape.index[ind] < 0.8] <- 0
#    # if (graphics) {
#    #    if (i > 1) for (j in 1:2) pop3d()
#    #    plot(subset(face, ind), display = "spheres", col = crv, add = TRUE)
#    #    spheres3d(nrst[i, ], radius = 1.5, col = "red")
#    # }
#    if (length(which(crv > 0.1)) > 10) flag <- TRUE
# }
# 
# if (flag) {
#    sbst  <- subset.face3d(face, edist.face3d(face$coords, nrst[i, ]) < 15)
#    sbst  <- index.face3d(sbst)
#    # sbst  <- subset.face3d(sbst, sbst$shape.index > 0.8)
#    # parts <- connected.face3d(sbst)
#    # sbst  <- subset.face3d(sbst, parts == 1)
#    dst   <- as.matrix(dist(sbst$coords))
#    crv   <- sbst$kappa1 * sbst$kappa2
#    # crv   <- pmax(sbst$kappa1 * sbst$kappa2, 0)
#    rss   <- apply(dst, 1, function(x) sum(lsfit(x^2, crv)$residuals^2))
#    pn    <- sbst$coords[which.min(rss), ]
#    lmks  <- rbind(lmks, pn = pn)
#    if (graphics) {
#       plot(sbst, display = "spheres", radius = 1.1, colour = crv, add = TRUE)
#       spheres3d(lmks["pn", ], radius = 1.5, col = "yellow")
#       ind <- which.min(rss)
#       plot(dst[ind, ], crv)
#    }
# }
# else {
#    cat("pn could not be identified.")
# }



# This function doesn't work
# Checking whether points are in the margins of the object

# margins.face3d <- function(shape, x = shape$coords, threshold = 10) {
#    
#    if (missing(x)) x <- shape$coords
#    if (is.vector(x) && length(x) ==3) x <- matrix(x, ncol = 3)
#    
#    fn <- function(x) {
#       dst  <- c(rdist(t(x), shape$coords))
#       nbrs <- shape$coords[dst <= threshold, ]
#       nbrs <- sweep(nbrs, 2, x)
#       vec  <- apply(nbrs, 2, mean)
#       sgn  <- c(sign(nbrs %*% vec))
#       sgn  <- length(which(sgn > 0)) / length(sgn)
#       (sgn > 0.8) | (sgn < 0.2)
#    }
#    invisible(apply(x, 1, fn))
# }
# 
