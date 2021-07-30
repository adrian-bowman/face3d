# Likelihood function for a landmark at the end of a ridge

lmklik <- function(face, lmk.name, reflmk.name, monitor = 0) {
   
   # Restrict attention to the area around approximate se, using sbst0 as a larger buffer
   
   dst   <- c(rdist(t(face$landmarks[lmk.name, ]), face$vertices))
   if (!("directions" %in% names(face)))
      face  <- curvatures(face, subset = dst < 20, directions = TRUE, overwrite = TRUE)
   sbst0 <- subset(face, dst < 20)
   sbst  <- subset(face, dst < 10, retain.indices = TRUE)
   drn1  <- sbst0$directions[ , 1, which.min(dst[dst < 20])]
   
   if (monitor > 0) {
      plot(sbst0, col = "kappa2")
      plot(sbst0, display = "principal 1", add = TRUE)
      lmk.initial <- face$landmarks[lmk.name, ]
      spheres3d(lmk.initial, radius = 2)
   }
   
   # Evaluate the curvature at each vertex
   nvert      <- nrow(sbst$vertices)
   jlst       <- 1:nvert
   crv        <- numeric(0)
   sbst0$area <- area.face3d(sbst0)$points
   for (j in jlst) {
      if (monitor > 1 & j != jlst[1]) for (i in 1:3) pop3d()
      lmk    <- sbst$vertices[j, ]
      drn1   <- sbst$directions[ , 1, j]
      dst    <- rdist(t(sbst$landmarks[reflmk.name, ]), rbind(lmk, lmk + drn1))
      drn1   <- if (dst[2] > dst[1]) -drn1 else drn1
      dst    <- c(rdist(t(lmk), sbst0$vertices))
      drn2   <- c(crossproduct(drn1, sbst$normals[j, ]))
      sbst2  <- sbst0
      sbst2$directions <- NULL
      path1   <- planepath(sbst2, lmk, direction = drn1, bothways = TRUE, rotation = 0)
      path2   <- planepath(sbst2, lmk, direction = drn2, bothways = TRUE, rotation = 0)
      if (!is.null(path1) & !is.null(path2)) {
         lmkarc  <- path1$arclength[which.min(c(rdist(t(lmk), path1$path)))]
         gcrv    <- gcurvature.face3d(path1$path, 3)
         gcrvlmk <- approx(gcrv$arclength, gcrv$gcurvature, lmkarc, ties = median)$y
         crv1    <- gcrvlmk / gcrv$gcurvature[gcrv$ind.max]
         lmkarc  <- path2$arclength[which.min(c(rdist(t(lmk), path2$path)))]
         gcrv    <- gcurvature.face3d(path2$path, 3)
         gcrvlmk <- approx(gcrv$arclength, gcrv$gcurvature, lmkarc, ties = median)$y
         crv2    <- gcrvlmk / gcrv$gcurvature[gcrv$ind.max]
         # if (!is.null(path)) {
         #    cdst <- closestcurve.face3d(sbst1, path)
         #    area <- sbst1$area
         #    ind  <- (cdst$closest.curvept != 1) & (abs(cdst$closest.distance) < 4) 
         #    rdg  <- -sum(area[ind] * sbst1$kappa2[ind])
         # }
         # else
         #    rdg <- NA
         # crv  <- c(crv, sbst$kappa1[j] * rdg)
         # crv   <- c(crv, -sbst$kappa1[j] * sbst$kappa2[j])
         crv   <- c(crv, crv1 * crv2)
      }
      else
         crv <- c(crv, NA)
      
      if (monitor > 1) {
         spheres3d(path1$path)
         spheres3d(path2$path)
         spheres3d(lmk, col = "yellow", radius = 2)
         # spheres3d(sbst1$vertices[ind, ], col = "yellow", radius = 0.5)
         if (monitor > 2) invisible(readline(prompt = " Press [enter] to continue"))
      }
   }
   lmk <- sbst$vertices[which.max(crv), ]

   if (monitor > 0) {
      plot(sbst, col = crv)
      spheres3d(lmk.initial)
      spheres3d(lmk, col = "yellow")
   }
   
   # Add the landmark and local curvature to the face object
   if (lmk.name %in% rownames(face$landmarks))
      face$landmarks[lmk.name, ] <- lmk
   else {
      face$landmarks <- rbind(face$landmarks, lmk)
      rownames(face$landmarks)[nrow(face$landmarks)] <- lmk.name
   }
   nm <- paste(lmk.name, ".ind", sep = "")
   face[[nm]] <- sbst$subset
   nm <- paste(lmk.name, ".crv", sep = "")
   face[[nm]] <- crv
   
   return(invisible(face))
}
