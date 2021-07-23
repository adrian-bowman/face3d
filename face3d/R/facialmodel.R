facialmodel.face3d <- function(face, pca, npc, pn.id,
                                   reltol = 0.01, sample.spacing, trim = 2 * sample.spacing,
                                   distance = 10,
                                   monitor = 1, overwrite = FALSE) {

   if (!("shape.index" %in% names(face))) {
      if (!missing(sample.spacing)) {
         if (monitor > 0) cat("Sampling ... ")
         sampled <- sample.face3d(face, spacing = sample.spacing)
         if (monitor > 0) cat(length(sampled), "points ... ")
         if (monitor > 1) {
            plot(face)
            spheres3d(face$vertices[sampled, ], radius = sample.spacing / 5)
         }
         if (monitor > 0) cat("interpolating curvature ... ")
         face <- index.face3d(face, distance = distance, subset = sampled, interpolate = TRUE,
                              monitor = monitor)
      }
      else
         face <- index.face3d(face, distance = distance)
   }
   face$gc <- face$kappa1 * face$kappa2

   if (monitor > 1) {
      plot(face, col = "shape index")
      # plot(face, col = face$gc, key = TRUE)
      if (monitor > 2) invisible(readline(prompt = "      Press [enter] to continue"))
   }

   if (monitor > 0) cat("Finding main convex area ... ")
   sbst.pos  <- subset(face, face$shape.index > 0, retain.indices = TRUE)
   # sbst.neg  <- subset(face, face$shape.index <= 0, retain.indices = TRUE)
   # ind.neg   <- face$shape.index <= 0
   parts.pos <- connected.face3d(sbst.pos)
   sbst.pos  <- subset(sbst.pos, parts.pos == 1)

   if (monitor > 0) cat("completed.\n")
   if (monitor > 1) {
      # plot(sbst.pos, col = "shape index")
      if (monitor > 2) {
         invisible(readline(prompt = "      Press [enter] to continue"))
         cat("      ")
      }
      plot(sbst.pos, col = sbst.pos$gc)
   }

   if (monitor > 0) cat("Restricting to areas with high curvature ... ")
   sbst.high <- subset(sbst.pos, sbst.pos$gc > quantile(sbst.pos$gc, 0.90), retain.indices = TRUE)
   p1        <- connected.face3d(sbst.high)

   # Exclude small areas and low integrated curvature
   areas     <- sapply(unique(p1), function(x) area.face3d(subset(sbst.high, p1 == x))$area)
   ind       <- unique(p1)[areas > 1000]
   sbst.high <- subset(sbst.high, p1 %in% ind)
   p1        <- p1[p1 %in% ind]
   
   # Exclude areas with low integrated curvature
   scores    <- sapply(unique(p1), function(x)
      sum(area.face3d(subset(sbst.high, p1 == x))$points * sbst.high$gc[p1 == x]))
   ind       <- unique(p1)[scores > 0.1 * max(scores)]
   sbst.high <- subset(sbst.high, p1 %in% ind)
   p1        <- p1[p1 %in% ind]
   
   if (monitor > 1) {
      plot(sbst.high, col = sbst.high$gc, display = "spheres", add = TRUE)
      if (monitor > 2) {
         invisible(readline(prompt = "      Press [enter] to continue"))
         cat("      ")
      }
   }

   # Remove those which are close to the largest edge
   # edges.pos <- edges.face3d(sbst.pos)
   # # lapply(edges.pos, function(x) rgl::lines3d(sbst.pos$vertices[x, ], lwd = 2, col = "yellow"))
   # edge.ind  <- which.max(sapply(edges.pos, function(x) max(rdist(sbst.pos$vertices[x, ]))))
   # edge.pts  <- sbst.pos$vertices[edges.pos[[edge.ind]], ]
   # ind       <- which(tapply(1:length(p1), p1, function(x) 
   #                           min(rdist(sbst.high$vertices[x, ], edge.pts)) > trim))
   # if (length(ind) > 0) sbst.high <- subset(sbst.high, p1 %in% ind)
   # p1        <- p1[p1 %in% ind]
   
   # Attempt to reconfigure the edges issue
   # sampled <- sample.face3d(sbst.pos, spacing = 10)
   # mn      <- apply(sbst.pos$vertices[sampled, ], 2, mean)
   # spheres3d(mn, col = "red", radius = 3)
   # vec.fn <- function(x) {
   #    sbs    <- subset(sbst.high, p1 == x)
   #    pn     <- mode.face3d(sbs, sbs$gc, 10)$mode
   #    # mean(edist(pn, sbst.pos$vertices[sampled, ]))
   #    # mn.vec <- apply(sweep(sbst.pos$vertices[sampled, ], 2, pn), 2, mean)
   #    # edist(mn.vec)
   #    edist(mn, pn)
   # }
   # sapply(unique(p1), vec.fn)
      
   # Ensure there is a patch of negative curvature (the eyes) very close
   # ind      <- logical(length = length(unique(p1)))
   # for (j in 1:length(unique(p1)))
   #    ind[j] <- (min(rdist(subset(sbst.high, p1 == unique(p1)[j])$vertices, sbst.neg$vertices)) < 30)
   # # plot(sbst.neg, col = "shape.index", add = TRUE)
   # ind       <- unique(p1)[ind]
   # sbst.high <- subset(sbst.high, p1 %in% ind)
   # p1        <- p1[p1 %in% ind]
      
   fn.rt <- function(pars) {
      crvs <- pca$mean
      for (i in 1:(length(pars) - 6))
         crvs <- crvs + pars[i + 6] * pca$sd[i] * matrix(pca$evecs[ , i], ncol = 3)
      crvs <- rotate3d(crvs, pars[1], 1, 0, 0)
      crvs <- rotate3d(crvs, pars[2], 0, 1, 0)
      crvs <- rotate3d(crvs, pars[3], 0, 0, 1)
      crvs <- sweep(crvs, 2, pars[4:6], "+")
      invisible(crvs)
   }
   
   fn <- function(pars, monitor = 0, face.pn) {
      if (any(abs(pars[7:9]) > 3)) return(Inf)
      crvs1 <- fn.rt(pars)
      
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
      
   fn.pn <- function(pars, monitor = 0) {
      crvs <- pca$mean
      crvs <- sweep(crvs, 2, crvs[pn.id, ])
      crvs <- rotate3d(crvs, pars[1], 1, 0, 0)
      crvs <- rotate3d(crvs, pars[2], 0, 1, 0)
      crvs <- rotate3d(crvs, pars[3], 0, 0, 1)
      crvs <- sweep(crvs, 2, pn, "+")
      dst  <- edist(crvs, face$vertices[sampled, ], minsum = TRUE) / nrow(crvs)
      if (monitor > 1) {
         pop3d()
         spheres3d(crvs, col = "yellow", radius = 2)
      }
      invisible(dst)
   }
   
   fn.selected <- function(pars, monitor = 0) {
      if (any(abs(pars[7:9]) > 3)) return(Inf)
      crvs1 <- fn.rt(pars)
      if (monitor > 1) {
         pop3d()
         spheres3d(crvs1, col = "yellow", radius = 2)
      }
      edist(crvs1, face$vertices[sampled, ], minsum = TRUE) / nrow(crvs1)
   }
   
   
   curvefit.fn <- function(sbs) {
      pn      <- mode.face3d(sbs, sbs$gc, distance)$mode
      face.pn <- subset(face, edist(pn, face) < 12 * distance)
      
      if (monitor > 1) {
         plot(face.pn)
         spheres3d(face$vertices[sampled, ])
         # spheres3d(pn, col = "red", radius = 1.5)
      }
      
      # if (monitor > 0) cat("Initial grid search ... ")
      # ngrid  <- 8
      # pgrid  <- seq(-pi, pi, length = ngrid + 1)[-1]
      # pgrid  <- as.matrix(expand.grid(pgrid, pgrid, pgrid))
      # shift  <- pn - pca$mean[pn.id, ]
      # pgrid  <- cbind(pgrid, shift[1], shift[2], shift[3], 0, 0, 0)
      # values <- apply(pgrid, 1, fn.selected, monitor = 0)
      # rots   <- pgrid[which.min(values), 1:3]
      # if(monitor > 0) cat("completed.\n")

      # Optimise while fixing pn
      # pars   <- c(rep(0, 3))
      # opt    <- optim(pars, fn.pn, monitor = monitor, pn = pn,
      #                 control = list(reltol = reltol))
      # return(invisible(c(opt$value, opt$par)))
      
      # pgrid  <- cbind(pgrid, parmn[1], parmn[2], parmn[3], 0, 0, 0)
      # values <- apply(pgrid, 1, fn, graphics = FALSE)
      # pars   <- pgrid[which.min(values), ]
      # fn(pars, graphics = TRUE)
      
      # Choose starting parameters carefully
      nose    <- subset(face, edist(pn, face) < 50)
      nose    <- index.face3d(nose, distance = distance, overwrite = TRUE)
      nose$gc <- nose$kappa1 * nose$kappa2
      nose    <- subset(nose, nose$shape.index > 0.25)
      parts   <- connected.face3d(nose)
      nose    <- subset(nose, parts == 1)
      print(area.face3d(nose)$area)
      pn      <- mode.face3d(nose, nose$gc, distance)$mode
      
      plot(nose, col = nose$gc)
      spheres3d(pn, radius = 2.5, col = "black")
      
      angle.grid <- seq(0, 2 * pi, length = 32)
      lngth <- numeric(0)
      for (angle in angle.grid) {
         # dirn need to be chosen more carefully
         dirn <- c(rotate3d(c(0,1,0), angle = angle, 0, 0, 1))
         path <- planepath.face3d(nose, pn, direction = dirn, rotation = 0)$path
         # pop3d()
         # spheres3d(path)
         lngth <- c(lngth, max(arclength.face3d(path)))
      }
      angle <- angle.grid[which.max(lngth)]
      dirn  <- c(rotate3d(c(0,1,0), angle = angle, 0, 0, 1))
      path  <- planepath.face3d(nose, pn, direction = dirn, rotation = 0)$path
      
      spheres3d(path, col = "green")
      scan()
      
      # Fast search using selected vertices
      rots <- rep(0, 3)
      pars <- c(rots, pn - pca$mean[pn.id, ], rep(0, 3))
      opt  <- optim(pars, fn.selected, monitor = monitor, control = list(reltol = reltol))
      
      # Refine by using all vertices
      # crvs1 <- rotate3d(pca$mean, opt$par[1], 1, 0, 0)
      # crvs1 <- rotate3d(crvs1,    opt$par[2], 0, 1, 0)
      # crvs1 <- rotate3d(crvs1,    opt$par[3], 0, 0, 1)
      # pars   <- c(opt$par, pn - crvs1[pn.id, ], rep(0, 3))
      # opt    <- optim(pars, fn, monitor = monitor, face.pn = face.pn, control = list(reltol = reltol))

      invisible(c(opt$value, opt$par))
   }
   
   if (monitor > 0) cat("examining", length(unique(p1)), "candidates for pn ... ")
   sampled <- if ("sampled" %in% names(face)) face$sampled
              else sample.face3d(face, spacing = sample.spacing)
   values <- sapply(unique(p1), function(x) curvefit.fn(subset(sbst.high, p1 == x)))
   ind    <- which.min(values[1, ])
   # curves <- fn.rt(values[-1, ind])
   sbs    <- subset(sbst.high, p1 == unique(p1)[ind])
   pn     <- mode.face3d(sbs, sbs$gc, distance)$mode
   curves <- fn.rt(values[-1, ind])
   if (monitor > 0) cat("completed.\npc parameters:", values[-(1:7), ind], "\n")

   invisible(curves)
}   
