facialmodel.face3d <- function(face, pca, npc, pn.id,
                                   reltol = 0.01, sample.spacing, trim = 2 * sample.spacing,
                                   distance = 10, sample.distance = 5 * distance,
                                   monitor = 1, overwrite = FALSE) {

   if (monitor > 0) cat("Sampling ... ")

   ncoord   <- nrow(face$vertices)
   selected <- as.integer(1)
   mindist  <- c(rdist(t(face$vertices[selected, ]), face$vertices))
   while(max(mindist) > sample.spacing & length(selected) < ncoord) {
      iselected <- which.max(mindist)
      selected  <- c(selected, iselected)
      idist     <- rdist(t(face$vertices[iselected, ]), face$vertices)
      mindist   <- pmin(idist, mindist)
   }
      
   # if (monitor > 0) cat(length(selected), "points ... ")
   if (monitor > 1) {
      plot(face)
      spheres3d(face$vertices[selected, ], radius = sample.spacing / 5)
   }
   if (monitor > 0) cat("interpolating curvature ... ")
   
   face     <- index.face3d(face, subset = selected, distance = sample.distance, overwrite = TRUE)
   ind      <- !is.na(face$shape.index[selected])
   selected <- selected[ind]
   to       <- cbind(face$shape.index, face$kappa1, face$kappa2)[selected, ]
   wp       <- warp.face3d(face$vertices[selected, ], to, face$vertices)
   face$shape.index <- pmax(pmin(wp[ , 1], 1), -1)      # Warping may extrapolate beyond [-1,1]
   face$kappa1      <- wp[ , 2]
   face$kappa2      <- wp[ , 3]
   face$gc          <- face$kappa1 * face$kappa2

   if (monitor > 1) {
      plot(face, col = "shape index")
      # plot(face, col = face$gc, key = TRUE)
      if (monitor > 2) invisible(readline(prompt = "      Press [enter] to continue"))
   }
   if (monitor > 0) cat("finding main convex area  ... ")
   
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
   
   # pop3d()
   # plot(sbst.high, col = sbst1$gc, display = "spheres", add = TRUE)

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
      
   fn <- function(pars) {
      if (any(abs(pars[7:9]) > 3)) return(Inf)
      crvs1 <- fn.rt(pars)
      dst   <- distance.face3d(crvs1, face$vertices[selected, ], minsum = TRUE) / nrow(crvs1)
      if (monitor > 1) {
         pop3d()
         spheres3d(crvs1, col = "yellow", radius = 2)
      }
      invisible(dst)
   }
      
   curvefit.fn <- function(sbs) {
      pn <- mode.face3d(sbs, sbs$gc, 10)$mode

      if (monitor > 1) {
         face.pn <- subset(face, distance.face3d(pn, face) < 120)
         plot(face.pn)
         spheres3d(face$vertices[selected, ])
         spheres3d(pn, col = "red", radius = 1.5)
      }
      
      # Simple grid search to locate the curves approximately
      # ngrid <- 4
      # angle.grid <- seq(-pi, pi, length = ngrid + 1)[-1]
      # angle.grid <- as.matrix(expand.grid(angle.grid, angle.grid, angle.grid))
      # values     <- apply(angle.grid, 1, fn, pn = pn)
      # pars       <- angle.grid[which.min(values), ]
      # fn(pars, graphics = TRUE, pn = pn)
      
      crvs  <- pca$mean
      crvs1 <- sweep(crvs, 2, crvs[pn.id, ] - pn)
      pars  <- c(rep(0, 3), apply(crvs1, 2, mean), rep(0, npc))
      opt   <- optim(pars, fn,
                     # lower = c(-pi, -pi, -pi, -Inf, -Inf, -Inf, -3, -3, -3),
                     # upper = c( pi,  pi,  pi,  Inf,  Inf,  Inf,  3,  3,  3),
                     control = list(reltol = reltol))

      invisible(c(opt$value, opt$par))
   }
      
   # Fit the curves to each candidate location and choose the best
   crvs   <- pca$mean
   values <- sapply(unique(p1), function(x) curvefit.fn(subset(sbst.high, p1 == x)))
   ind    <- which.min(values[1, ])
   curves <- fn.rt(values[-1, ind])
   
   # curvefit.fn(subset(sbst.high, p1 == 1))
   # 
   # plot(subset(sbst.high, p1 == 1))
   # summary(subset(sbst.high, p1 == 1))
   # summary(subset(sbst.high, p1 == 2))
   # area.face3d(subset(sbst.high, p1 == 1))$area
   # area.face3d(subset(sbst.high, p1 == 2))$area
   # 
   # plot(sbst.pos, col = sbst.pos$gc)
   # plot(sbst.high, col = sbst.high$gc, display = "spheres", add = TRUE)
   # fn(values[-1, ind])
   
   if (monitor > 0) cat("completed.\n")
   if (monitor > 1) fn(values[-1, ind])

   invisible(curves)
}   
