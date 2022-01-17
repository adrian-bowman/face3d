facialmodel <- function(face, pca, npc, pn.id, se.id,
                        distance = 10, sample.spacing, sample.distance = 30,
                        reltol = 0.01, gc.quantile = 0.9, gcint.minpropmax = 0.1,
                        area.min = 500, trim = 0,
                        monitor = 1, overwrite = FALSE) {

   if (monitor > 0) cat("Fitting facial model ...\n")
   
   if (!("shape.index" %in% names(face))) {
      if (!missing(sample.spacing)) {
         if (monitor > 0) cat("  sampling ... ")
         sampled <- sample.face3d(face, spacing = sample.spacing)
         if (monitor > 0) cat(length(sampled), "points ... ")
         if (monitor > 1) {
            plot(face)
            spheres3d(face$vertices[sampled, ], radius = sample.spacing / 5)
         }
         if (monitor > 0) cat("interpolating curvatures ...\n")
         face <- curvatures(face, distance = sample.distance, subset = sampled, interpolate = TRUE,
                              monitor = monitor)
      }
      else {
         if (monitor > 0) cat("estimating curvatures ...\n")
         face <- curvatures(face, distance = sample.distance)
      }
   }
   face$gc <- face$kappa1 * face$kappa2

   if (monitor > 1) {
      plot(face, col = "shape index")
      if (monitor > 2) monitor3()
   }

   if (monitor > 0) cat("  locating main convex area ...\n")
   sbst.pos  <- subset(face, face$shape.index > 0, retain.indices = TRUE)
   parts.pos <- connected.face3d(sbst.pos)
   sbst.pos  <- subset(sbst.pos, parts.pos == 1)

   if (monitor > 1) {
      plot(sbst.pos, col = "shape index")
      if (monitor > 2) monitor3()
   }

   if (monitor > 0) cat("  locating areas of high curvature ...\n")
   sbst.high <- subset(sbst.pos, sbst.pos$gc > quantile(sbst.pos$gc, gc.quantile), retain.indices = TRUE)
   p1        <- connected.face3d(sbst.high)

   if (monitor > 0) cat("  excluding areas which are ... small ... ")
   areas     <- sapply(unique(p1), function(x) areas(subset(sbst.high, p1 == x))$area)
   ind       <- unique(p1)[areas > area.min]
   sbst.high <- subset(sbst.high, p1 %in% ind)
   p1        <- p1[p1 %in% ind]
   
   if (monitor > 0) cat("low in curvature ... ")
   scores    <- sapply(unique(p1), function(x)
                       sum(areas(subset(sbst.high, p1 == x))$points * sbst.high$gc[p1 == x]))
   ind       <- unique(p1)[scores > gcint.minpropmax * max(scores)]
   sbst.high <- subset(sbst.high, p1 %in% ind)
   p1        <- p1[p1 %in% ind]
   
   if (trim > 0) {
      if (monitor > 0) cat("close to the edge ...\n")
      edge      <- edges(face)[[1]]
      edge      <- face$vertices[edge, ]
      edge.dst  <- sapply(unique(p1), function(x)
                          edist(subset(sbst.high, p1 == x), edge, minsum = TRUE) / length(p1[p1 == x]))
      ind       <- unique(p1)[edge.dst > trim]
      if (length(ind) > 0) sbst.high <- subset(sbst.high, p1 %in% ind)
      p1        <- p1[p1 %in% ind]
   }
   
   if (monitor > 1) {
      plot(sbst.high, col = sbst.high$gc, display = "spheres", add = TRUE)
      if (monitor > 2) monitor3()
   }
   
   fn.rt <- function(pars, angle1.start, angle2.start, axis1.start, axis2.start, pn) {
      # Transform the curves by translation rotation and PCs
      crvs <- pca$mean
      for (i in 1:(length(pars) - 6))
         crvs <- crvs + pars[i + 6] * pca$sd[i] * matrix(pca$evecs[ , i], ncol = 3)
      crvs <- rotate3d(crvs, angle1.start, axis1.start[1], axis1.start[2], axis1.start[3])
      crvs <- rotate3d(crvs, angle2.start, axis2.start[1], axis2.start[2], axis2.start[3])
      crvs <- rotate3d(crvs, pars[1], 1, 0, 0)
      crvs <- rotate3d(crvs, pars[2], 0, 1, 0)
      crvs <- rotate3d(crvs, pars[3], 0, 0, 1)
      crvs <- sweep(crvs, 2, pn - crvs[pn.id, ] + pars[4:6], "+")
      invisible(crvs)
   }
   
   distance.fn <- function(pars, angle1.start, angle2.start, axis1.start, axis2.start, pn, face.pn,
                  distance.type, monitor = 0) {
      # if (any(abs(pars[-(1:6)]) > 3)) return(Inf)
      crvs1 <- fn.rt(pars, angle1.start, angle2.start, axis1.start, axis2.start, pn)
      
      if (distance.type == "sampled")
         dst   <- edist(crvs1, face$vertices[sampled, ], minsum = TRUE) / nrow(crvs1)
      else if (distance.type == "all")
         dst   <- edist(crvs1, face.pn$vertices, minsum = TRUE) / nrow(crvs1)
      else if (distance.type == "projection") {
         # Distance based on projection using the tangent to the curve.
         # This needs diffmat from create_priors_curves.
         tangents <- crvs1[diffmat[ , 1], ] - crvs1[diffmat[ , 2], ]
         deltas   <- edist(tangents)
         tangents <- tangents / deltas
         deltas   <- deltas / 2
         dist.fn <- function(i) {
            ind      <- which(edist(crvs1[i, ], face.pn) < distance)
            if (length(ind) < 5) return(nrow(crvs1) * 100 * distance)
            vertices <- sweep(face.pn$vertices[ind, ], 2, crvs1[i, ])
            proj     <- c(vertices %*% tangents[i, ])
            ind1     <- which(abs(proj) < deltas[i])
            if (length(ind1) < 5) return(nrow(crvs1) * 100 * distance)
            min(edist(crvs1[i, ], face.pn$vertices[ind[ind1], ]))
         }
         dst <- sapply(1:nrow(crvs1), dist.fn)
         dst <- sum(dst) / nrow(crvs1)
      }
      
      if (monitor > 1) {
         pop3d()
         spheres3d(crvs1, col = "yellow", radius = 2)
      }
      invisible(dst)
   }
   
   curvefit.fn <- function(x, distance.type) {
      
      if (monitor > 0) cat("   ", match(x, unique(p1)), "... ")
      
      nose    <- subset(sbst.high, p1 == x)
      amax    <- areamax(nose, nose$gc, distance)
      pn      <- amax$point
      face.pn <- subset(face, edist(face, pn) < 120)
      
      # Try fitting a nose shape from pca$mean
      # ind <- grepl("mid-line nasal profile", rownames(pca$mean)) |
      #    grepl("mid-line nasal root",    rownames(pca$mean)) |
      #    grepl("mid-line columella",    rownames(pca$mean)) |
      #    grepl("nasal base",    rownames(pca$mean)) |
      #    grepl("nasal bridge",    rownames(pca$mean))
      # spheres3d(pca$mean[ind, ])
      # scan()
      
      if (monitor > 1) {
         plot(nose, display = "spheres", col = nose$gc, add = TRUE)
         if (monitor > 2) monitor3()
      }
      
      if (monitor > 0) cat("enL/R ... ")
      sbst <- subset(face, edist(face, pn) < 8 * distance & face$shape.index < 0)
      eye1 <- areamax(sbst, -sbst$shape.index, distance)$point
      sbst <- subset(sbst, edist(sbst, eye1) > distance)
      eye2 <- areamax(sbst, -sbst$shape.index, distance)$point
      eyes.angle <- acos(sum((eye1 - pn) * (eye2 - pn)) / (edist(eye1 - pn) * edist(eye2 - pn)))
      if (monitor > 1) spheres3d(rbind(pn, eye1, eye2), col = "green", radius = 3)


      # The angle between the eyes from pn should be < pi/2 and could be used as a filter
      # if (monitor > 0) cat("angle between the eyes from pn:", round(eyes.angle, 2), " ...")

      if (monitor > 0) cat("se ... ")
      crv  <- planepath(face, eye1, eye2, boundary = c(0.2, 2), monitor = 0)$path
      if (is.null(crv)) {
         if (monitor > 0) cat("model mismatch: null\n")
         return(invisible(c(Inf, rep(0, 20))))
      }
      gcrv <- gcurvature(crv, 3)
      se   <- gcrv$pos.max
       
      if (monitor > 1) {
         spheres3d(se, col = "yellow", radius = 3)
         if (monitor > 2) monitor3()
      }
      
      # Find a good starting point by rotating around the pn-se line
      ft         <- numeric(0)
      ftmin      <- Inf
      a          <- se - pn
      b          <- pca$mean[se.id, ] - pca$mean[pn.id, ]
      angle1     <- acos(sum(a * b) / (edist(a) * edist(b)))
      axis1      <- crossproduct(a, b)
      axis2      <- a
      angle.grid <- seq(-pi, pi, length = 32)[-1]
      fta        <- numeric(0)
      for (angle in angle.grid)
         fta <- c(fta, distance.fn(rep(0, 6 + npc), angle1, angle, axis1, axis2, pn,
                                   distance.type = distance.type, monitor = 0))
      ind        <- which.min(fta)
      ftmin      <- fta[ind]
      angle2     <- angle.grid[ind]
      if (is.infinite(ftmin)) return(invisible(c(Inf, rep(0, 18 + npc))))
      opt        <- list(value = ftmin, par = rep(0, 6 + npc))

      if (monitor > 1) {
         distance.fn(rep(0, 6 + npc), angle1, angle2, axis1, axis2, pn,
                     distance.type = distance.type, monitor = monitor)
         if (monitor > 2) monitor3() 
      }

      # Fast search using selected vertices
      # opt  <- optim(rep(0, 9), distance.fn, distance.type = distance.type, monitor = monitor,
      #               angle1.start = angle1, angle2.start = angle2,
      #               axis1.start = axis1, axis2.start = axis2, pn = pn,
      #               control = list(reltol = reltol))

     # Use projection with all vertices
      # opt  <- optim(rep(0, 9), fn, monitor = monitor,
      #               angle1.start = angle1, angle2.start = angle2,
      #               axis1.start = axis1, axis2.start = axis2, pn = pn,
      #               face.pn = face.pn,
      #               control = list(reltol = reltol))
      
      # Refine by using all vertices
      # crvs1 <- rotate3d(pca$mean, opt$par[1], 1, 0, 0)
      # crvs1 <- rotate3d(crvs1,    opt$par[2], 0, 1, 0)
      # crvs1 <- rotate3d(crvs1,    opt$par[3], 0, 0, 1)
      # pars   <- c(opt$par, pn - crvs1[pn.id, ], rep(0, 3))
      # opt    <- optim(pars, fn, monitor = monitor, face.pn = face.pn, control = list(reltol = reltol))
      
      if (monitor > 0) {
         cat("model mismatch:", opt$value, "\n")
         if (monitor > 1) for (j in 1:3) pop3d()
      }

      invisible(c(opt$value, opt$par, angle1, angle2, axis1, axis2, pn))
   }
   
   sampled <- if ("sampled" %in% names(face)) face$sampled
              else sample.face3d(face, spacing = sample.spacing)

   if (monitor > 0) {
      cand <- if (length(unique(p1)) > 1) "candidates" else "candidate"
      cat("  examining", length(unique(p1)), cand, "for the nose ...\n")
      if (monitor > 1) {
         plot(face)
         spheres3d(face$vertices[sampled, ])
      }
   }
   
   values <- sapply(unique(p1), curvefit.fn, distance.type = "sampled")
   if (!all(is.null(values[1, ]))) {
      ind   <- which.min(values[1, ])
      first <- values[ , ind]  
   }
   else {
      if (monitor > 0) cat("A facial model cannot be located.\n")
      return(invisible(NULL))
   }
   
   # curves <- fn.rt(values[-1, ind])
   # sbs    <- subset(sbst.high, p1 == unique(p1)[ind])
   # pn     <- areamax(sbs, sbs$gc, distance)$mode
   # curves <- fn.rt(values[2:10, ind], values[11, ind], values[12, ind],
   #                 values[13:15, ind], values[16:18, ind], values[19:21, ind])
   
   if (monitor > 0) cat("  refining by using all vertices ....\n")
   # Refine by using all vertices
   pn      <- first[16:18 + npc]
   angle1  <- first[8 + npc]
   angle2  <- first[9 + npc]
   axis1   <- first[10:12 + npc]
   axis2   <- first[13:15 + npc]
   face.pn <- subset(face, edist(face, pn) < 120)
   opt     <- optim(first[2:(7 + npc)], distance.fn, distance.type = "all",
                    angle1 = angle1, angle2 = angle2, axis1 = axis1, axis2 = axis2,
                    pn = pn, face.pn = face.pn,
                    monitor = monitor, control = list(reltol = reltol))
   curves  <- fn.rt(opt$par, angle1.start = angle1, angle2.start = angle2,
                     axis1.start = axis1, axis2.start = axis2, pn = pn)
   if (monitor > 0) cat("  completed.\n")
      
   if (monitor > 0) {
      cat("  location parameters:", round(opt$par[4:6], 2), "\n")
      cat("  rotation parameters:", round(opt$par[1:3] * 360 / (2 * pi), 2), "\n")
      cat("        pc parameters:", round(opt$par[-(1:6)], 2), "\n")
   }
   if (monitor > 1) {
      plot(face)
      spheres3d(curves, col = "yellow", radius = 2)
   }

   invisible(curves)
}   
