facialmodel <- function(face, pca, npc, pn.id,
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
                          edist(subset(sbst.high, p1 == x), edge, "minsum") / length(p1[p1 == x]))
      ind       <- unique(p1)[edge.dst > trim]
      if (length(ind) > 0) sbst.high <- subset(sbst.high, p1 %in% ind)
      p1        <- p1[p1 %in% ind]
   }
   
   if (monitor > 1) {
      plot(sbst.high, col = sbst.high$gc, display = "spheres", add = TRUE)
      if (monitor > 2) monitor3()
   }
   
   fn.rt <- function(pars, reference) {
      # Transform the curves by translation rotation and PCs
      crvs <- pca$mean
      if (length(pars) > 6) {
         for (i in 1:(length(pars) - 6))
            crvs <- crvs + pars[i + 6] * pca$sd[i] * matrix(pca$evecs[ , i], ncol = 3)
      }
      crvs <- rotate3d(crvs, reference$angle1, reference$axis1[1], reference$axis1[2], reference$axis1[3])
      crvs <- rotate3d(crvs, reference$angle2, reference$axis2[1], reference$axis2[2], reference$axis2[3])
      crvs <- rotate3d(crvs, pars[1], 1, 0, 0)
      crvs <- rotate3d(crvs, pars[2], 0, 1, 0)
      crvs <- rotate3d(crvs, pars[3], 0, 0, 1)
      crvs <- sweep(crvs, 2, reference$centre + pars[4:6], "+")
      invisible(crvs)
   }
   
   distance.fn <- function(pars, reference, distance.type, monitor = 0) {
      # if (any(abs(pars[-(1:6)]) > 3)) return(Inf)
      crvs1 <- fn.rt(pars, reference)
      
      if (distance.type == "sampled")
         dst   <- edist(crvs1, face$vertices[sampled, ], "minsum") / nrow(crvs1)
      else if (distance.type == "all")
         dst   <- edist(crvs1, reference$face.pn$vertices, "minsum") / nrow(crvs1)
      else if (distance.type == "projection") {
         # Distance based on projection using the tangent to the curve.
         # This needs diffmat from create_priors_curves.
         tangents <- crvs1[diffmat[ , 1], ] - crvs1[diffmat[ , 2], ]
         deltas   <- edist(tangents)
         tangents <- tangents / deltas
         deltas   <- deltas / 2
         dist.fn <- function(i) {
            ind      <- which(edist(crvs1[i, ], reference$face.pn) < distance)
            if (length(ind) < 5) return(nrow(crvs1) * 100 * distance)
            vertices <- sweep(reference$face.pn$vertices[ind, ], 2, crvs1[i, ])
            proj     <- c(vertices %*% tangents[i, ])
            ind1     <- which(abs(proj) < deltas[i])
            if (length(ind1) < 5) return(nrow(crvs1) * 100 * distance)
            min(edist(crvs1[i, ], reference$face.pn$vertices[ind[ind1], ]))
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
   
   initialfit.fn <- function(x, distance.type) {
      
      if (monitor > 0) cat("   ", match(x, unique(p1)), "... ")
      
      # Find approximate pn
      nose    <- subset(sbst.high, p1 == x)
      amax    <- areamax(nose, nose$gc, distance)
      pn      <- amax$point
      face.pn <- subset(face, edist(face, pn) < 120)
      
      if (monitor > 1) {
         plot(nose, display = "spheres", col = nose$gc, add = TRUE)
         if (monitor > 2) monitor3()
      }
      
      # Rotate the template (pca$mean) around pn to find the best fit
      reference <- list(pn = pn, centre = pn - pca$mean[pn.id, ],
                        angle1 = 0, angle2 = 0, axis1 = c(1, 0, 0), axis2 = c(0, 1, 0))
      opt  <- optim(rep(0, 6), distance.fn, distance.type = distance.type, monitor = monitor,
                    reference = reference, control = list(reltol = reltol))

      if (monitor > 1) {
         distance.fn(c(opt$par, rep(0, npc)), reference = reference,
                     distance.type = distance.type, monitor = monitor)
         if (monitor > 2) monitor3() 
      }

      if (monitor > 0) {
         cat("model mismatch:", opt$value, "\n")
         if (monitor > 1) pop3d()
      }

      invisible(c(opt$value, opt$par,
                  reference$angle1, reference$angle2, reference$axis1, reference$axis2, reference$pn))
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
   
   values <- sapply(unique(p1), initialfit.fn, distance.type = "sampled")
   if (!all(is.null(values[1, ]))) {
      ind       <- which.min(values[1, ])
      first     <- values[ , ind]  
      reference <- list(angle1 = first[8], angle2 = first[9], axis1 = first[10:12], axis2 = first[13:15],
                        pn = first[16:18], centre = first[16:18] - pca$mean[pn.id, ])
      curves    <- fn.rt(first[2:7], reference)
   }
   else {
      if (monitor > 0) cat("A facial model cannot be located.\n")
      return(invisible(NULL))
   }
   
   if (monitor > 0) cat("  refining by using all vertices ....\n")
   # reference$centre  <- apply(curves - pca$mean, 2, mean)
   reference$centre  <- reference$pn - pca$mean[pn.id, ]
   # reference$face.pn <- subset(face, edist(face, reference$pn) < 120)
   reference$face.pn <- subset(face, edist(face, curves, "min") < 10)
   opt    <- optim(c(first[2:7], rep(0, npc)), distance.fn, distance.type = "all",
                    reference = reference, monitor = monitor, control = list(reltol = reltol))
   curves <- fn.rt(opt$par, reference)
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
