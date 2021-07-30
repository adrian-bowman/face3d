facialmodel.face3d <- function(face, pca, npc, pn.id,
                               distance = 10, sample.spacing,
                               reltol = 0.01, gc.quantile = 0.9, gcint.minpropmax = 0.1, area.min = 500,
                               monitor = 1, overwrite = FALSE) {

   if (monitor > 0) cat("Fitting shape model ...\n")
   
   if (!("shape.index" %in% names(face))) {
      if (!missing(sample.spacing)) {
         if (monitor > 0) cat("sampling ... ")
         sampled <- sample.face3d(face, spacing = sample.spacing)
         if (monitor > 0) cat(length(sampled), "points ... ")
         if (monitor > 1) {
            plot(face)
            spheres3d(face$vertices[sampled, ], radius = sample.spacing / 5)
         }
         if (monitor > 0) cat("interpolating curvatures ... ")
         face <- curvatures(face, distance = distance, subset = sampled, interpolate = TRUE,
                              monitor = monitor)
      }
      else {
         if (monitor > 0) cat("estimating curvatures ... ")
         face <- curvatures(face, distance = distance)
      }
   }
   else {
      if (monitor > 1) cat("  plotting shape index ...\n")
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

   # Exclude small areas
   areas     <- sapply(unique(p1), function(x) areas(subset(sbst.high, p1 == x))$area)
   ind       <- unique(p1)[areas > area.min]
   sbst.high <- subset(sbst.high, p1 %in% ind)
   p1        <- p1[p1 %in% ind]
   
   # Exclude areas with low integrated curvature
   scores    <- sapply(unique(p1), function(x)
      sum(areas(subset(sbst.high, p1 == x))$points * sbst.high$gc[p1 == x]))
   ind       <- unique(p1)[scores > gcint.minpropmax * max(scores)]
   sbst.high <- subset(sbst.high, p1 %in% ind)
   p1        <- p1[p1 %in% ind]
   
   if (monitor > 1) {
      plot(sbst.high, col = sbst.high$gc, display = "spheres", add = TRUE)
      if (monitor > 2) monitor3()
   }
   
   fn.rt <- function(pars, angle1.start, angle2.start, axis1.start, axis2.start, pn) {
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
   
   fn <- function(pars, angle1.start, angle2.start, axis1.start, axis2.start, pn, face.pn, monitor = 0) {
      if (any(abs(pars[7:9]) > 3)) return(Inf)
      crvs1 <- fn.rt(pars, angle1.start, angle2.start, axis1.start, axis2.start, pn)
      
      # Distance computed from closest point to sampled vertices
      # dst   <- edist(crvs1, face$vertices[sampled, ], minsum = TRUE) / nrow(crvs1)
      # Distance computed from closest point to all vertices
      # dst   <- edist(crvs1, face.pn$vertices, minsum = TRUE) / nrow(crvs1)
      
      # Distance based on projection using the tangent to the curve
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
      
      if (monitor > 1) {
         pop3d()
         spheres3d(crvs1, col = "yellow", radius = 2)
      }
      invisible(dst)
   }
      
   fn.selected <- function(pars, angle1.start, angle2.start, axis1.start, axis2.start, pn, monitor = 0) {
      if (any(abs(pars[7:9]) > 3)) return(Inf)
      crvs1 <- fn.rt(pars, angle1.start, angle2.start, axis1.start, axis2.start, pn)
      if (monitor > 1) {
         pop3d()
         spheres3d(crvs1, col = "yellow", radius = 2)
      }
      edist(crvs1, face$vertices[sampled, ], minsum = TRUE) / nrow(crvs1)
   }
   
   
   curvefit.fn <- function(nose) {
      amax    <- areamax(nose, nose$gc, distance)
      pn      <- amax$point
      nrm     <- nose$normals[amax$id, ]
      pdir1   <- nose$axes[[amax$id]][ , 2]
      face.pn <- subset(face, edist(pn, face) < 12 * distance)
      nose    <- subset(face, edist(pn, face) <  5 * distance)
      nose    <- curvatures(nose, distance = distance, overwrite = TRUE)
      nose    <- subset(nose, nose$shape.index > 0.25)

      if (monitor > 1) {
         plot(nose, display = "spheres", col = "shape index", add = TRUE)
         # Plotting pn gives something for pop3d to remove in code below
         spheres3d(pn, radius = 1.5, col = "blue")
         if (monitor > 2) monitor3()
      }
      
      # Find the nose ridge by locating the longest planepath from pn
      angle.grid <- seq(0, 2 * pi, length = 32)[-1]
      lngth      <- numeric(0)
      for (angle in angle.grid) {
         dirn  <- c(rotate3d(pdir1, angle = angle, nrm[1], nrm[2], nrm[3]))
         path  <- planepath(nose, pn, direction = dirn, rotation = 0)$path
         lng   <- if (!is.null(path)) max(arclength(path)) else NA
         lngth <- c(lngth, lng)
         if (monitor > 1) {
            pop3d()
            spheres3d(path, radius = 2)
         }
      }
      angle1.start <- angle.grid[which.max(lngth)]
      if (length(angle1.start) == 0) return(invisible(c(Inf, rep(0, 20))))
      
      dirn  <- c(rotate3d(pdir1, angle1.start, nrm[1], nrm[2], nrm[3]))
      path  <- planepath(nose, pn, direction = dirn, rotation = 0)$path

      if (monitor > 1) {
         pop3d()
         spheres3d(path, radius = 2)
         if (monitor > 2) monitor3() 
      }

      # Locate a starting point for the curve template
      # The information on the pca$mean curve needs to be passed as an argument of curvefit.fn
      path.start   <- pca$mean[pn.id - 1:14, ]
      ind.start    <- which.min(abs(arclength(path.start) - 20))
      ind.path     <- which.min(abs(arclength(path) - 20))
      a            <- path[ind.path, ] - path[1, ]
      b            <- path.start[ind.start, ] - pca$mean[pn.id, ]
      angle1.start <- acos(sum(a * b) / (edist(a) * edist(b)))
      axis1.start  <- crossproduct(a, b)
      axis2.start  <- a
      
      ft <- numeric(0)
      for (angle in angle.grid) {
         ft <- c(ft, fn.selected(rep(0, 9), angle1.start, angle, axis1.start, axis2.start, pn, monitor = 0))
      }
      angle2.start <- angle.grid[which.min(ft)]
      
      if (monitor > 1) {
         plot(face)
         spheres3d(pn + a, radius = 2, col = "red")
         fn.selected(rep(0, 9), angle1.start, angle2.start, axis1.start, axis2.start, pn, monitor = monitor)
         if (monitor > 2) monitor3()
      }

      # Fast search using selected vertices
      opt  <- optim(rep(0, 9), fn.selected, monitor = monitor,
                    angle1.start = angle1.start, angle2.start = angle2.start,
                    axis1.start = axis1.start, axis2.start = axis2.start, pn = pn,
                    control = list(reltol = reltol))
      
      # Refine by using all vertices
      # crvs1 <- rotate3d(pca$mean, opt$par[1], 1, 0, 0)
      # crvs1 <- rotate3d(crvs1,    opt$par[2], 0, 1, 0)
      # crvs1 <- rotate3d(crvs1,    opt$par[3], 0, 0, 1)
      # pars   <- c(opt$par, pn - crvs1[pn.id, ], rep(0, 3))
      # opt    <- optim(pars, fn, monitor = monitor, face.pn = face.pn, control = list(reltol = reltol))
      
      if (monitor > 1) {
         pop3d()
         pop3d()
      }

      invisible(c(opt$value, opt$par, angle1.start, angle2.start, axis1.start, axis2.start, pn))
   }
   
   sampled <- if ("sampled" %in% names(face)) face$sampled
              else sample.face3d(face, spacing = sample.spacing)

   if (monitor > 0) {
      cand <- if (length(unique(p1)) > 1) "candidates" else "candidate"
      cat("  examining", length(unique(p1)), cand, "for pn ... \n")
   }
   if (monitor > 1) {
      plot(face)
      spheres3d(face$vertices[sampled, ])
   }
   
   values <- sapply(unique(p1), function(x) curvefit.fn(subset(sbst.high, p1 == x)))
   ind    <- which.min(values[1, ])
   # curves <- fn.rt(values[-1, ind])
   sbs    <- subset(sbst.high, p1 == unique(p1)[ind])
   # pn     <- areamax(sbs, sbs$gc, distance)$mode
   curves <- fn.rt(values[2:10, ind], values[11, ind], values[12, ind],
                   values[13:15, ind], values[16:18, ind], values[19:21, ind])

   if (monitor > 0) cat("npc parameters:", values[8:10, ind], "\ncompleted.")
   if (monitor > 1) {
      plot(face)
      spheres3d(curves, col = "yellow", radius = 2)
   }

   invisible(curves)
}   
