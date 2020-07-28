# Likelihood function for a landmark at the end of a ridge

lmklik <- function(face, lmk.name, reflmk.name, monitor = 0) {
   
   # Restrict attention to the area around approximate se, using sbst0 as a larger buffer
   
   dst   <- c(rdist(t(face$landmarks[lmk.name, ]), face$vertices))
   if (!("directions" %in% names(face)))
      face  <- index.face3d(face, subset = dst < 20, directions = TRUE, overwrite = TRUE)
   sbst0 <- subset(face, dst < 20)
   sbst  <- subset(face, dst < 10)
   
   if (monitor > 0) {
      plot(sbst, col = "kappa2")
      plot(sbst, display = "direction 1", add = TRUE)
      lmk.initial <- face$landmarks[lmk.name, ]
      spheres3d(lmk.initial, radius = 2)
   }
   
   # Evaluate the curvature at each vertex
   nvert <- nrow(sbst$vertices)
   jlst  <- 1:nvert
   crv   <- numeric(0)
   for (j in jlst) {
      if (monitor > 1 & j != jlst[1]) for (i in 1:2) pop3d()
      lmk   <- sbst$vertices[j, ]
      drn   <- sbst$directions[ , 1, j]
      dst   <- rdist(t(sbst$landmarks[reflmk.name, ]), rbind(lmk, lmk + drn))
      drn   <- if (dst[2] > dst[1]) -drn else drn
      dst   <- c(rdist(t(lmk), sbst0$vertices))
      sbst1 <- subset(sbst0, dst < 10)
      path  <- planepath.face3d(sbst1, lmk, direction = drn, si.target = 1, rotation = 0)$path
      if (!is.null(path)) {
         cdst <- closestcurve.face3d(sbst, path)
         area <- area.face3d(sbst)$points
         ind  <- (cdst$closest.curvept != 1) & (abs(cdst$closest.distance) < 4) 
         rdg  <- -sum(area[ind] * sbst$kappa2[ind])
      }
      else
         rdg <- NA
      crv  <- c(crv, sbst$kappa1[j] * rdg)
      if (monitor > 1) {
         spheres3d(path)
         spheres3d(lmk, col = "yellow", radius = 2)
      }
   }
   
   ind   <- which(!is.na(crv))
   sbstj <- subset(sbst, jlst[ind], remove.singles = FALSE)
   lmk   <- sbstj$vertices[which.max(crv[ind]), ]

   if (monitor > 0) {
      plot(sbstj, col = crv[ind])
      spheres3d(lmk, col = "yellow", radius = 2)
   }
   
   # Add the landmark to the face object
   if (lmk.name %in% rownames(face$landmarks))
      face$landmarks[lmk.name, ] <- lmk
   else {
      face$landmarks <- rbind(face$landmarks, lmk)
      rownames(face$landmarks)[nrow(face$landmarks)] <- lmk.name
   }
   
   return(face)
}