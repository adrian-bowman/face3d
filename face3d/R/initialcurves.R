#     Initial location of curves
   
initialcurves.face3d <- function(face, curve.names,
                                 monitor = FALSE,new.window = FALSE) {

   if (missing(curve.names))
      curve.names <- c("mid-line nasal profile", "mid-line columella",
                       "nasal bridge left",      "nasal bridge right",
                       "nasal base left",        "nasal base right")
   curves <- list()
   face1  <- subset(face, edist.face3d(face$vertices, face$landmarks["pn", ]) < 100)

   if (monitor) cat("Computing curvatures ... ")
   face1  <- curvatures(face1)
   if (monitor) cat("complete.\n")
   
   # Identify the nasal region by thresholding kappa2

   thresh <- quantile(face1$kappa2[face1$kappa2 < 0], 0.35, na.rm = TRUE)
   sbst1  <- subset(face1,  face1$kappa2 < thresh)
   sbst   <- subset(sbst1, connected.face3d(sbst1) == 1)
   ind1   <- as.numeric(rownames(sbst1$vertices))
   ind2   <- as.numeric(rownames(sbst$vertices))
   ind    <- ind1[ind2]
   # Possible adjustment of distance for curvature computation
   # sbst  <- curvatures(sbst, distance = 5, overwrite = TRUE)
   
   if (monitor) {
      # Possibly plot face and add the nose
      # face2  <- subset(face1, -ind)
      # plot(face2, new = FALSE)
      plot(sbst, col = sbst$kappa2, new = new.window)
      spheres3d(face$landmarks["pn", ], radius = 1.2, col = "yellow")
   }
   
   #--------------------------------------------------------------------------
   #     mid-line nasal profile (landmark: se)
   #--------------------------------------------------------------------------

   if ("mid-line nasal profile" %in% curve.names) {
      se  <- face$landmarks["se", ]
      pp  <- planepath.face3d(sbst, face$landmarks["se", ], face$landmarks["pn", ],
                              rotation.range = pi/8, bridge.gaps = TRUE)
      drn <- face$landmarks["se", ] - face$landmarks["pn", ]
      pp  <- planepath.face3d(sbst, face$landmarks["se", ], direction = drn, rotation = 0)
      gc  <- gcurvature.face3d(pp$path, 4)
      face$landmarks["se", ] <- gc$resampled.curve[which.max(gc$gcurvature), ]
      pp  <- planepath.face3d(sbst, face$landmarks["se", ], face$landmarks["pn", ],
                              rotation.range = pi/8, bridge.gaps = TRUE)
      curves$"mid-line nasal profile" <- pp$path
      if (monitor) {
         spheres3d(pp$path, col = "red", radius = 0.5)
         spheres3d(face$landmarks["se", ], radius = 1.2, col = "yellow")
      }
   }
   
   #--------------------------------------------------------------------------
   #     nasal bridge L/R (landmarks: acL/R)
   #--------------------------------------------------------------------------

   curve.nms <- c("nasal bridge left", "nasal bridge right")
   flag      <- curve.nms %in% curve.names
   if (any(flag)) {
      for (curve.nm in curve.nms[flag]) {
         lmk.nm   <- if (curve.nm == "nasal bridge left") "acL" else "acR"
         # First approximation through a planepath
         # ids  <- closest.face3d(face$landmarks["se", ], sbst)$ids
         ids  <- closest.face3d(face$landmarks["pn", ], sbst)$ids
         # Allow the possibility that only one id is returned.
         nrms <- matrix(c(sbst$normals[ids, ]), ncol = 3)
         nrm  <- apply(nrms, 2, mean)
         drn  <- face$landmarks["se", ] - face$landmarks["pn", ]
         drn  <- if (curve.nm == "nasal bridge right") -drn else drn
         drn  <- c(crossproduct(drn, nrm))
         pph  <- planepath.face3d(sbst, face$landmarks["pn", ], si.target = 1,
                                   direction = drn, normal = nrm, rotation.range = pi/12)
         crv  <- pph$path
         # Now adjust to a smooth path
         # proj  <- projectline.face3d(sbst, face$landmarks["pn", ], ac)
         # sbst3 <- subset(sbst, (proj$x >= 0) & (proj$y < 10))
         # hist(proj$x)
         # hist(proj$y)
         # plot(sbst3, new = FALSE, col = sbst3$kappa2)
         cc    <- closestcurve.face3d(sbst, crv)
         arcl  <- arclength.face3d(crv)[cc$closest.curvept]
         dst   <- cc$closest.distance
         ind   <- (abs(dst) < 10) & (arcl > 0.4 * max(arcl)) & (arcl < max(arcl))
         # ind   <- (abs(dst) < 5) & (arcl > 0) & (arcl < max(arcl))
         arcl  <- arcl[ind]
         dst   <- dst[ind]
         sbst2 <- subset(sbst, ind, remove.singles = FALSE)
         ridge <- ridge2d.face3d(arcl, dst, -sbst2$kappa2, lambda = 0.02,
                                 monitor = FALSE,
                                 endpoints.fixed = c(TRUE, FALSE), ngrid = 20)
         crv.x <- interp.barycentric(cbind(arcl, dst), sbst2$vertices[ , 1],
                                      cbind(ridge$x, ridge$y))$fnew
         crv.y <- interp.barycentric(cbind(arcl, dst), sbst2$vertices[ , 2],
                                      cbind(ridge$x, ridge$y))$fnew
         crv.z <- interp.barycentric(cbind(arcl, dst), sbst2$vertices[ , 3],
                                      cbind(ridge$x, ridge$y))$fnew
         crv   <- rbind(crv[pph$arclength <= 0.4 * max(pph$arclength), ],
                        cbind(crv.x, crv.y, crv.z))
         # Extend the new curve and find the point of maximum curvature
         ac    <- c(tail(crv, 1))
         sbst3 <- subset(face1, edist.face3d(face1$vertices, ac) < 10)
         sbst3 <- curvatures(sbst3, distance = 3, directions = TRUE, overwrite = TRUE)
         drn   <- c(diff(tail(crv, 2)))
         pp    <- planepath.face3d(sbst3, ac, direction = drn, rotation = 0)
         gc    <- gcurvature.face3d(pp$path, 4)
         ac    <- gc$resampled.curve[gc$ind.max, ]
         al    <- arclength.face3d(pp$path)
         crv   <- rbind(crv, pp$path[al < gc$arclength[gc$ind.max], ], ac)
         if (monitor) {
            # spheres3d(pph$path, radius = 0.5)
            spheres3d(crv, radius = 0.5, col = "red")
            # spheres3d(pp$path, col = "green", radius = 0.2)
            spheres3d(ac, col = "yellow")
	         model  <- smooth.face3d(cbind(arcl, dst), -sbst2$kappa2, df = 12, ngrid = 100)
	         image(model$eval.points[[1]], model$eval.point[[2]],
                  matrix(model$estimate, nrow = length(model$eval.points[[1]])),
                  col = topo.colors(20),
	               xlab = "Arc Length", ylab = "Perpendicular Distance")
            points(ridge, pch = 16, col = "red")
         }
         if (lmk.nm %in% names(face$landmarks))
            face$landmarks[lmk.nm, ] <- ac
         else {
            face$landmarks <- rbind(face$landmarks, ac)
            rownames(face$landmarks)[nrow(face$landmarks)] <- lmk.nm
         }
         curves[[curve.nm]] <- crv
      }
   }

   #--------------------------------------------------------------------------
   #     mid-line columella (landmark: sn)
   #--------------------------------------------------------------------------

   if ("mid-line columella" %in% curve.names) {
      pp0 <- planepath.face3d(face1, face$landmarks["acL", ], face$landmarks["acR", ], boundary = NA)
      gc  <- gcurvature.face3d(pp0$path, 3, monitor = monitor)
      sn  <- gc$pos.max
      pp1 <- planepath.face3d(face1, face$landmarks["pn", ], sn, rotation = 0)
      pp2 <- planepath.face3d(face1, sn,
                              direction = c(diff(tail(pp1$path, 2))),
                              normal = pp1$normal, rotation = 0)
      pp2 <- rbind(pp1$path, pp2$path[pp2$arclength < 10, ])
      gc  <- gcurvature.face3d(pp2, 3, monitor = monitor)
      crv <- gc$gcurvature * gc$d2.z
      ind.max <- which.max(crv)
      pos.max <- gc$resampled.curve[ind.max, ]
      if ("sn" %in% rownames(face$landmarks))
         face$landmarks["sn", ] <- pos.max
      else {
         face$landmarks <- rbind(face$landmarks, pos.max)
         rownames(face$landmarks)[nrow(face$landmarks)] <- "sn"
      }
      pp2 <- planepath.face3d(face1, face$landmarks["pn", ], face$landmarks["sn", ],
                              normal = pp1$normal, rotation = 0)
      curves$"mid-line columella" <- pp2$path
      if (monitor) {
         # spheres3d(pp0$path, col = "yellow", radius = 0.5)
         spheres3d(gc$pos.max, col = "yellow")
         spheres3d(pp0$path, radius = 0.5, col = "black")
         spheres3d(pp2$path, radius = 0.5, col = "red")
      }
   }
   
   #--------------------------------------------------------------------------
   #     nasal base L/R
   #--------------------------------------------------------------------------

   if ("nasal base left" %in% curve.names) {
      curves$"nasal base left" <-
         planepath.face3d(face1, face$landmarks["acL", ], face$landmarks["sn",])$path
      if (monitor) spheres3d(curves$"nasal base left", col = "red", radius = 0.5)
   }
   if ("nasal base right" %in% curve.names) {
      curves$"nasal base right" <-
         planepath.face3d(face1, face$landmarks["acR", ], face$landmarks["sn",])$path
      if (monitor) spheres3d(curves$"nasal base right", col = "red", radius = 0.5)
   }

   if ("mid-line philtral" %in% curve.names) {
      sbst <- subset(face, edist.face3d(face$vertices, face$landmarks["sn", ]) < 40)
      plot(sbst, new = FALSE, col = sbst$kappa2)
      spheres3d(face$landmarks[c("pn", "sn"), ], col = "yellow")
   }

   if ("curves" %in% names(face))
      for (i in curve.names) face$curves[[i]] <- curves[[i]]
   else
      face$curves <- curves
   
   invisible(face)
}
