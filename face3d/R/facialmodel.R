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
   
   sbst.pos  <- subset(face, face$shape.index >  0, retain.indices = TRUE)
   sbst.neg  <- subset(face, face$shape.index <= 0, retain.indices = TRUE)
   ind.neg   <- face$shape.index <= 0
   parts.pos <- connected.face3d(sbst.pos)
   parts.neg <- connected.face3d(sbst.neg)
   sbst.pos  <- subset(sbst.pos, parts.pos == 1)

   if (monitor > 0) cat("completed.\n")
   if (monitor > 1) {
      plot(sbst.pos, col = "shape index")
      if (monitor > 2) {
         invisible(readline(prompt = "      Press [enter] to continue"))
         cat("      ")
      }
      plot(sbst.pos, col = sbst.pos$gc)
   }
   if (monitor > 0) cat("Examining pn candidates ... ")

   # Find the patches with very high Gaussian curvature
   sbst.high  <- subset(sbst.pos, sbst.pos$gc > quantile(sbst.pos$gc, 0.90), retain.indices = TRUE)

   if (monitor > 1) {
      plot(sbst.high, col = sbst.high$gc, display = "spheres", add = TRUE)
      if (monitor > 2) {
         invisible(readline(prompt = "      Press [enter] to continue"))
         cat("      ")
      }
   }

   # sbst2    <- subset(sbst.pos, sbst.pos$gc > quantile(sbst.pos$gc, 0.80), retain.indices = TRUE)
   # plot(sbst2, col = sbst2$gc)
   # return()
      
   # Remove those which are close to the largest edge
   p1        <- connected.face3d(sbst.high)
   edges.pos <- edges.face3d(sbst.pos)
   # lapply(edges.pos, function(x) rgl::lines3d(sbst.pos$vertices[x, ], lwd = 2, col = "yellow"))
   edge.ind  <- which.max(sapply(edges.pos, function(x) max(rdist(sbst.pos$vertices[x, ]))))
   edge.pts  <- sbst.pos$vertices[edges.pos[[edge.ind]], ]
   ind       <- which(tapply(1:length(p1), p1, function(x) 
                             min(rdist(sbst.high$vertices[x, ], edge.pts)) > trim))
   if (length(ind) > 0) sbst.high <- subset(sbst.high, p1 %in% ind)
   p1        <- p1[p1 %in% ind]
      
   # Ensure there is a patch of negative curvature (the eyes) very close
   ind      <- logical(length = length(unique(p1)))
   for (j in 1:length(unique(p1)))
      ind[j] <- (min(rdist(subset(sbst.high, p1 == unique(p1)[j])$vertices, sbst.neg$vertices)) < 20)
   ind       <- unique(p1)[ind]
   sbst.high <- subset(sbst.high, p1 %in% ind)
   p1        <- p1[p1 %in% ind]
      
   # Select those patches with high integrated curvature
   scores    <- sapply(unique(p1), function(x)
                       sum(area.face3d(subset(sbst.high, p1 == x))$points * sbst.high$gc[p1 == x]))
   ind       <- unique(p1)[scores > 0.1 * max(scores)]
   sbst.high <- subset(sbst.high, p1 %in% ind)
   p1        <- p1[p1 %in% ind]
   
   if (monitor > 1) {
      pop3d()
      plot(sbst.high, col = sbst.high$gc, display = "spheres", add = TRUE)
      if (monitor > 2) {
         invisible(readline(prompt = "      Press [enter] to continue"))
         cat("      ")
      }
   }
   
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
   
   fn <- function(pars, graphics = FALSE) {
      if (any(abs(pars[7:9]) > 3)) return(Inf)
      crvs1 <- fn.rt(pars)
      dst   <- distance.face3d(crvs1, face$vertices[sampled, ], minsum = TRUE) / nrow(crvs1)
      if (graphics) {
         pop3d()
         spheres3d(crvs1, col = "yellow", radius = 2)
      }
      invisible(dst)
   }
      
   fn.pn <- function(pars, pn, graphics = FALSE) {
      crvs <- pca$mean
      crvs <- sweep(crvs, 2, crvs[pn.id, ])
      crvs <- rotate3d(crvs, pars[1], 1, 0, 0)
      crvs <- rotate3d(crvs, pars[2], 0, 1, 0)
      crvs <- rotate3d(crvs, pars[3], 0, 0, 1)
      crvs <- sweep(crvs, 2, pn, "+")
      dst  <- distance.face3d(crvs, face$vertices[sampled, ], minsum = TRUE) / nrow(crvs)
      if (graphics) {
         pop3d()
         spheres3d(crvs, col = "yellow", radius = 2)
      }
      invisible(dst)
   }
   
   curvefit.fn <- function(sbs) {
      pn <- mode.face3d(sbs, sbs$gc, 10)$mode

      if (monitor > 1) {
         face.pn <- subset(face, distance.face3d(pn, face) < 120)
         plot(face.pn)
         spheres3d(face$vertices[sampled, ])
         spheres3d(pn, col = "red", radius = 1.5)
      }
      
      # ngrid  <- 8
      # pgrid  <- seq(-pi, pi, length = ngrid + 1)[-1]
      # pgrid  <- as.matrix(expand.grid(pgrid, pgrid, pgrid))
      # values <- apply(pgrid, 1, fn.pn, graphics = TRUE)
      # pars   <- pgrid[which.min(values), ]
      # fn.pn(pars, graphics = TRUE)
      
      pars   <- c(rep(0, 3))
      opt    <- optim(pars, fn.pn, graphics = monitor > 1, pn = pn,
                      control = list(reltol = reltol))
      
      # pgrid  <- cbind(pgrid, parmn[1], parmn[2], parmn[3], 0, 0, 0)
      # values <- apply(pgrid, 1, fn, graphics = FALSE)
      # pars   <- pgrid[which.min(values), ]
      # fn(pars, graphics = TRUE)
      
      crvs1 <- rotate3d(pca$mean, opt$par[1], 1, 0, 0)
      crvs1 <- rotate3d(crvs1,    opt$par[2], 0, 1, 0)
      crvs1 <- rotate3d(crvs1,    opt$par[3], 0, 0, 1)
      pars   <- c(opt$par, pn - crvs1[pn.id, ], rep(0, 3))
      opt    <- optim(pars, fn, graphics = monitor > 1, control = list(reltol = reltol))

      invisible(c(opt$value, opt$par))
   }
   
   # sbs <- subset(sbst.high, p1 == unique(p1)[1])
   
   # Use a sample of vertices to make the calculation of distances efficient
   sampled <- sample.face3d(face, spacing = sample.spacing)
      
   # Fit the curves to each candidate location and choose the best
   values <- sapply(unique(p1), function(x) curvefit.fn(subset(sbst.high, p1 == x)))
   ind    <- which.min(values[1, ])
   curves <- fn.rt(values[-1, ind])
   
   if (monitor > 0) cat("completed.\npc parameters:", values[-(1:7), ind], "\n")
   if (monitor > 1) fn(values[-1, ind])

   invisible(curves)
}   
