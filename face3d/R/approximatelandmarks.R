approximatelandmarks.face3d <- function(face, landmark.names = c("pn", "enL", "enR", "se"),
                                   sample.spacing, trim = 2 * sample.spacing,
                                   distance = 10, sample.distance = 5 * distance,
                                   monitor = 1, overwrite = FALSE) {

   # Check whether any requested landmarks are already present
   if ("landmarks" %in% names(face)) {
      ind <- match(rownames(face$landmarks), landmark.names)
      ind <- ind[!is.na(ind)]
      if (length(ind) > 0 & !overwrite) {
         multiple <- 1 + as.numeric(length(ind) > 1)
         cat(c("Landmark", "Landmarks")[multiple],
             landmark.names[ind], c("is", "are")[multiple],
             "already present. Use the overwrite argument if required.\n")
         landmark.names <- landmark.names[-ind]      }
   }
   if (length(landmark.names) == 0) return(invisible(face))

   lmks <- matrix(nrow = length(landmark.names), ncol = 3, dimnames = list(landmark.names))

   #---------------------------------------------------------------------------------
   #        Estimate global curvature
   #---------------------------------------------------------------------------------
   
   if (any(c("pn", "enL", "enR") %in% landmark.names) |all(landmark.names == "none")) {
      
      if (monitor > 0) cat("Sampling ... ")
   
      ncoord   <- nrow(face$coords)
      selected <- 1
      mindist  <- c(rdist(t(face$coords[selected, ]), face$coords))
      while(max(mindist) > sample.spacing & length(selected) < ncoord) {
         iselected <- which.max(mindist)
         selected  <- c(selected, iselected)
         idist     <- rdist(t(face$coords[iselected, ]), face$coords)
         mindist   <- pmin(idist, mindist)
      }
      if (monitor > 0) cat(length(selected), "points ... ")
      if (monitor > 1) {
         plot(face)
         spheres3d(face$coords[selected, ], radius = sample.spacing / 5)
      }

      if (monitor > 0) cat("interpolating curvature ... ")
   
      face     <- index.face3d(face, subset = selected, distance = sample.distance, overwrite = TRUE)
      ind      <- !is.na(face$shape.index[selected])
      selected <- selected[ind]
      to       <- cbind(face$shape.index, face$kappa1, face$kappa2)[selected, ]
      wp       <- warp.face3d(face$coords[selected, ], to, face$coords)
      face$shape.index <- pmax(pmin(wp[ , 1], 1), -1)
      face$kappa1      <- wp[ , 2]
      face$kappa2      <- wp[ , 3]
      face$gc          <- face$kappa1 * face$kappa2

      if (monitor > 0) cat("completed.\n")
      if (monitor > 2) {
         plot(face, col = "shape index")
         invisible(readline(prompt = "      Press [enter] to continue"))
      }
   
      if (all(landmark.names == "none")) return(invisible(face))
      
      sbst.pos  <- subset(face, face$shape.index >  0, retain.indices = TRUE)
      sbst.neg  <- subset(face, face$shape.index <= 0, retain.indices = TRUE)
      ind.neg   <- face$shape.index <= 0
   }

   if (monitor > 0) cat("Approximating location of")
   
   #---------------------------------------------------------------------------------
   #        pn
   #---------------------------------------------------------------------------------
   
   if ("pn" %in% landmark.names) {
      
      if (monitor == 1) cat(" ... pn ...")
      if (monitor > 1)  cat(":\n  pn: main area of positive shape index ... ")
      
      parts.pos <- connected.face3d(sbst.pos)
      parts.neg <- connected.face3d(sbst.neg)
      sbst.pos  <- subset(sbst.pos, parts.pos == 1)
      edges.pos <- edges.face3d(sbst.pos)
      edge.ind  <- which.max(sapply(edges.pos, function(x) max(rdist(sbst.pos$coords[x, ]))))
      
      if (monitor > 1) {
         plot(sbst.pos, col = "shape index")
         if (monitor > 2) {
            invisible(readline(prompt = "      Press [enter] to continue"))
            cat("      ")
         }
         cat("strongest Gaussian curvature ...")
         plot(sbst.pos, col = sbst.pos$gc)
      }

      # Find the patches with very high Gaussian curvature
      sbst1    <- subset(sbst.pos, sbst.pos$gc > quantile(sbst.pos$gc, 0.90), retain.indices = TRUE)

      # Remove those which are close to the largest edge
      p1       <- connected.face3d(sbst1)
      edge.pts <- sbst.pos$coords[edges.pos[[edge.ind]], ]
      ind      <- which(tapply(1:length(p1), p1, function(x) 
                               min(rdist(sbst1$coords[x, ], edge.pts)) > trim))
      if (length(ind) > 0) sbst1 <- subset(sbst1, p1 %in% ind)
      p1       <- p1[p1 %in% ind]
      
      # Ensure there is a patch of negative curvature (the eyes) very close
      ind      <- logical(length = length(unique(p1)))
      for (j in 1:length(unique(p1)))
         ind[j] <- (min(rdist(subset(sbst1, p1 == unique(p1)[j])$coords, sbst.neg$coords)) < 20)
      sbst1 <- subset(sbst1, p1 %in% unique(p1)[ind])
      p1 <- p1[p1 %in% unique(p1)[ind]]
      
      # Integrate the scores over each patch and find the largest
      scores <- numeric(0)
      for (j in unique(p1))
         scores <- c(scores, sum(area.face3d(subset(sbst1, p1 == j))$points * sbst1$gc[p1 == j]))
      ind  <- order(scores, decreasing = TRUE)
      sbs  <- subset(sbst1, p1 == unique(p1)[ind[1]])
      mode <- mode.face3d(sbs, sbs$gc, 10)
      if (length(ind) > 1) {
         sbs   <- subset(sbst1, p1 == unique(p1)[ind[2]])
         mode2 <- mode.face3d(sbs, sbs$gc, 10)
         if (mode2$value > mode$value) mode <- mode2
      }
      lmks["pn", ] <- mode$mode

      if (monitor > 1) {
         plot(subset(sbst.pos, -sbst1$subset))
         plot(sbst1, col = sbst1$gc, add = TRUE)
         spheres3d(lmks["pn", ], radius = 3, col = "red")
         if (monitor > 2) {
            invisible(readline(prompt = "      Press [enter] to continue"))
            cat("      ")
         }
         cat("recompute curvatures with a smaller distance ...")
      }
      
      # Refine the estimate by computing shape index with a smaller value of distance
      face$shape.index <- NULL
      dst  <- c(rdist(t(lmks["pn", ]), face$coords))
      face <- index.face3d(face, subset = dst < 30, distance = distance)
      sbst <- subset(face, dst < 30)
      
      # Focus only on curvature where shape index is positive
      crv          <- sbst$kappa1 * sbst$kappa2
      ind          <- (sbst$kappa1 > 0) & (sbst$kappa2 > 0)
      crv[ind]     <- -crv[ind]
      lmks["pn", ] <- mode.face3d(sbst, crv, 10)$mode

      if (monitor > 1) {
         dst  <- c(rdist(t(lmks["pn", ]), sbst.pos$coords))
         plot(subset(sbst.pos, dst > 30))
         plot(sbst, col = crv, add = TRUE)
         spheres3d(lmks["pn", ], radius = 3, col = "red")
         cat("completed.")
         if (monitor > 2) invisible(readline(prompt = "      Press [enter] to continue"))
      }

   }
   
   #---------------------------------------------------------------------------------
   #        enL/R
   #---------------------------------------------------------------------------------
   
   if (all(c("enR", "enL") %in% landmark.names)) {

      pn <- if ("pn" %in% rownames(lmks)) lmks["pn", ] else
            if ("pn" %in% rownames(face$landmarks)) face$landmarks["pn", ] else
               stop("pn is needed to find enL/R.")
      
      if (monitor == 1) cat("en ... ")
      if (monitor == 2) cat("\n")
      if (monitor >  1) cat("  en: main area of negative shape index ...")
      
      dst         <- c(rdist(t(pn), face$coords))
      
      # Calculate the ind.neg for large distance if "pn" has not been estimated.
      # This is not the same as the earlier calculation because that one interpolated.
      # Remove this?
      if (!("pn" %in% landmark.names)) {
         face     <- index.face3d(face, subset = dst < 70, distance = sample.distance)
         sbst.neg <- subset(face, dst < 70)
         ind.neg  <- (sbst.neg$shape.index <= 0)
         sbst.pos <- subset(sbst.neg, !ind.neg)
         face$shape.index <- NULL
      }
      
      if (monitor > 1) {
         plot(sbst.pos, col = "shape index")
         plot(sbst.neg, col = "shape index", add = TRUE)
      }
      
      # Compute shape index with a smaller distance
      face <- index.face3d(face, subset = dst < 70, distance = distance)
      if (("pn" %in% landmark.names))
         sbst.neg   <- subset(face, ind.neg & (dst < 70))
      else
         sbst.neg    <- subset(sbst.neg, ind.neg)
      dst            <- c(rdist(t(pn), sbst.neg$coords))
      
      # Use the new shape index to define the areas with positive Gaussian curvature
      sbst.neg$gc    <- sbst.neg$kappa1 * sbst.neg$kappa2
      sbst.neg       <- subset(sbst.neg, sbst.neg$gc > 0 & dst > 40)
      # sbst.neg       <- subset(sbst.neg, sbst.neg$shape.index > 0 & dst > 40)
      sbst.neg$parts <- connected.face3d(sbst.neg)
      sbst.neg       <- subset(sbst.neg, sbst.neg$parts %in% 1:2)
      
      mode           <- mode.face3d(sbst.neg, sbst.neg$gc, 5)
      lmks["enL", ]  <- mode$mode
      dst            <- rdist(t(lmks["enL", ]), sbst.neg$coords)
      ind            <- which.min(dst)
      sbst.neg1      <- subset(sbst.neg, sbst.neg$parts != sbst.neg$parts[ind])
      mode           <- mode.face3d(sbst.neg1, sbst.neg1$gc, 5)
      lmks["enR", ]  <- mode$mode

      if (monitor > 1) {
         cat("completed.")
         plot(sbst.neg, col = sbst.neg$gc)
         spheres3d(lmks[c("enL", "enR"), ], radius = 5, col = "red")
         if (monitor > 2) invisible(readline(prompt = "      Press [enter] to continue"))
      }

   }
   
   #---------------------------------------------------------------------------------
   #        se (requires enL and enR)
   #---------------------------------------------------------------------------------
   
   if ("se" %in% landmark.names) {
      
      pn  <- if ("pn" %in% rownames(lmks)) lmks["pn", ] else
             if ("pn" %in% rownames(face$landmarks)) face$landmarks["pn", ] else
                stop("pn is needed to find se.")
      enL <- if ("enL" %in% rownames(lmks)) lmks["enL", ] else
             if ("enL" %in% rownames(face$landmarks)) face$landmarks["enL", ] else
             stop("enL is needed to find se.")
      enR <- if ("enR" %in% rownames(lmks)) lmks["enR", ] else
             if ("enR" %in% rownames(face$landmarks)) face$landmarks["enR", ] else
             stop("enR is needed to find se.")
      
      if (monitor == 1) cat("se ... ")
      if (monitor == 2) cat("\n")
      if (monitor >  1) cat("  se: shortest path between enL and enR ...")

      # First approximation through the curve from enL to enR
      curve  <- planepath.face3d(face, enL, enR, boundary = c(0.2, 2))$path
      gcrv   <- gcurvature.face3d(curve, 3)
      se     <- gcrv$pos.max
      
      dst   <- rdist(t(se), face$coords)
      sbst  <- subset(face, dst < 10)
      ppath <- planepath.face3d(sbst, se, direction = se - pn, rotation = 0,
                                bothways = TRUE, distance = distance)$path
      gcrv  <- gcurvature.face3d(ppath, 3)
      se    <- gcrv$pos.max
      
      # Searching for the mode of -gc does not seem to work very well
      # dst          <- c(rdist(t(se), face$coords))
      # face         <- index.face3d(face, subset = dst < 30)
      # sbst         <- subset(face, dst < 30)
      # crv          <- -sbst$kappa1 * sbst$kappa2
      # lmks["se", ] <- mode.face3d(sbst, crv, 5)$mode
      
      if (monitor > 0) cat("completed.\n")
      if (monitor > 1) {
         plot(sbst)
         spheres3d(ppath)
         spheres3d(se,  radius = 2, col = "red")
         if (monitor > 2) invisible(readline(prompt = "      Press [enter] to continue"))
      }
      
      lmks["se", ] <- se
   }
   
   #---------------------------------------------------------------------------------
   #        Assign L/R to enL and enR (requires se)
   #---------------------------------------------------------------------------------

   if (all(c("enL", "enR") %in% rownames(lmks))) {
      se  <- if ("se" %in% rownames(lmks)) lmks["se", ] else
             if ("se" %in% rownames(face$landmarks)) face$landmarks["se", ] else
                NULL
      if (is.null(se))
         cat("se is required to find L and R for en.")
      else {
         dst <- c(rdist(face$coords, t(se)))
         nrm <- face$normals[which.min(dst), ]
         axs <- c(crossproduct(lmks["pn", ] - se, nrm))
         en  <- c("enL", "enR")
         en1 <- if (sum((lmks["enL", ] - se) * axs) < 0) en else rev(en)
         lmks[en, ] <- lmks[en1, ]
      }
   }

   #---------------------------------------------------------------------------------
   #        acL and acR (requires pn)
   #---------------------------------------------------------------------------------
   
   if (all(c("acR", "acL") %in% landmark.names)) {
      
      pn <- if ("pn" %in% rownames(lmks)) lmks["pn", ] else
         if ("pn" %in% rownames(face$landmarks)) face$landmarks["pn", ] else
            stop("pn is needed to find acL/R.")
      
      dst       <- c(rdist(face$coords, t(pn)))
      ind       <- which(dst < 50)
      overwrite <- (face$si.distance != distance)
      face      <- index.face3d(face, subset = ind, distance = distance, overwrite = overwrite)
      sbst      <- subset(face, ind)

      if (monitor == 1) cat("acL/acR ... ")
      if (monitor > 1) {
         cat("\n  acL/acR: ...")
         plot(sbst, col = sbst$kappa2)
         spheres3d(pn, col = "red", radius = 2)
         if (monitor > 2) invisible(readline(prompt = "      Press [enter] to continue"))
      }
      
      nrm           <- face$normals[which.min(dst), ]
      drn           <- c(crossproduct(face$landmarks["se", ] - pn, nrm))
      pathL         <- planepath.face3d(sbst, pn, direction =  drn,
                                        rotation.range = pi / 4, si.target = 1)$path
      pathR         <- planepath.face3d(sbst, pn, direction = -drn,
                                        rotation.range = pi / 4, si.target = 1)$path
      lmks["acL", ] <- gcurvature.face3d(pathL, 10)$pos.max
      lmks["acR", ] <- gcurvature.face3d(pathR, 10)$pos.max
   
      # drn   <- pn -face$landmarks["se", ]
      # pathM <- planepath.face3d(sbst, pn, direction =  drn, rotation.range = pi / 4, si.target = 1)
      # sn    <- gcurvature.face3d(pathM$path, 10)$pos.max
   
      if (monitor > 0) cat("completed.\n")
      if (monitor > 1) {
         spheres3d(pathL)
         spheres3d(pathR)
         spheres3d(lmks["acL", ], radius = 3, col = "orange")
         spheres3d(lmks["acR", ], radius = 3, col = "green")
         if (monitor > 2) invisible(readline(prompt = "      Press [enter] to continue"))
      }
   }
   
   #---------------------------------------------------------------------------------
   #        Return the landmarks
   #---------------------------------------------------------------------------------
   
   # if (monitor > 1) {
   #    plot(face)
   #    clr                            <- rep("blue", nrow(lmks))
   #    clr[grep("R", rownames(lmks))] <- "green"
   #    clr[grep("L", rownames(lmks))] <- "orange"
   #    spheres3d(lmks, radius = 3, col = clr)
   # }
   # if (monitor == 1) cat(" ... completed.\n")
   
   
   # Place lmks in face$landmarks
   for (nm in rownames(lmks)) {
      if (nm %in% rownames(face$landmarks))
         face$landmarks[nm, ] <- lmks[nm, ]
      else {
         rnms                     <- c(rownames(face$landmarks), nm)
         face$landmarks           <- rbind(face$landmarks, lmks[nm, ])
         rownames(face$landmarks) <- rnms
      }
   }
   
   return(invisible(face))
}
