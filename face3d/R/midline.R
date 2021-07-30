#     Estimate a mid-line based on local asymmetry

midline.face3d <- function(shape, x1, x2, reference, d.asym = 20, d.search = 20, nv = 19, nh = 19,
                           extension = c(0, 0), df = 24, lambda = 1, quantile = 0.1, monitor = FALSE) {

   reference.missing <- missing(reference)
   if (reference.missing)
      reference  <- planepath(shape, x1, x2, si.target = 1, directions = TRUE)
   kappa1.ref <- reference$kappa1
   kappa2.ref <- reference$kappa2
   si.ref     <- reference$shape.index
   drns.ref   <- reference$directions
   reference  <- reference$path
   if (reference.missing) {
      if (extension[1] > 0) {
         pp.x1      <- planepath(shape, x1, direction = x1 - x2, rotation = 0, bothways = FALSE, directions = TRUE)
         ind        <- which(pp.x1$arclength < extension[1])
         reference  <- rbind(pp.x1$path[rev(ind), ], reference)
         kappa1.ref <- c(pp.x1$kappa1, kappa1.ref)
         kappa2.ref <- c(pp.x1$kappa1, kappa2.ref)
      }
      if (extension[2] > 0) {
         pp.x1     <- planepath(shape, x2, direction = x2 - x1, rotation = 0, bothways = FALSE, directions = TRUE)
         ind       <- which(pp.x1$arclength < extension[2])
         reference <- rbind(reference, pp.x1$path[ind, ])
         kappa1.ref <- c(kappa1.ref, pp.x1$kappa1)
         kappa2.ref <- c(kappa2.ref, pp.x1$kappa2)
      }
   } 
   arclength  <- function(path) c(0, cumsum(apply(diff(path)^2, 1, function(x) sqrt(sum(x)))))
   al         <- arclength(reference)
   kappa1.ref <- approx(al, kappa1.ref, n = nv)$y
   kappa2.ref <- approx(al, kappa2.ref, n = nv)$y
   si.ref     <- approx(al, si.ref,     n = nv)$y
   reference  <- resample.face3d(reference, nv)
   rd         <- rdist(shape$vertices, reference)
   rd         <- apply(rd, 1, min)
   ind.sbst   <- (rd < (d.search + d.asym) * 1.1 * 2)
   sbst       <- subset.face3d(shape, ind.sbst)

   if (monitor) {
      plot(sbst)
      if (reference.missing) spheres3d(rbind(x1, x2), radius = 0.5, col = "red")
      spheres3d(reference, radius = 0.4, col = "black")
      # if (monitor.prompt) invisible(readline(prompt="Press [enter] to continue:"))
   }
   
   # Find the average kappa2 direction for points with a ridge signal
   ind     <- !is.na(si.ref) & (si.ref > 1/8) & (si.ref < 5/8)
   drns    <- drns.ref[ , 2, ind]
   dotp    <- apply(drns, 2, function(x) sum(x * drns[ , 2]))
   ind1    <- as.numeric(dotp > 0)
   drns    <- sweep(drns, 2, (2 * as.numeric(dotp > 0) - 1), "*")
   direct  <- apply(drns, 1, mean, na.rm = TRUE)
   direct  <- direct / sqrt(sum(direct^2))

   # Create a grid of locations and compute asymmetry in kappa2
   ndel  <- floor((nh - 1) / 2)
   nasy  <- floor((ndel - 1) / 2)
   grd   <- matrix(nrow = 0, ncol = 3)
   asym  <- numeric(0)

   y.name <- "kappa2"
   # y.name <- "shape.index"
   y.ref  <- switch(y.name, kappa1 = kappa2.ref, kappa2 = kappa2.ref, shape.index = si.ref)
   
   for (i in 1:nv) {
      path.r <- planepath(sbst, reference[i, ], direction =  direct, rotation = 0, boundary = d.search * 1.1,
                                 bothways = FALSE, directions = TRUE)
      path.l <- planepath(sbst, reference[i, ], direction = -direct, rotation = 0, boundary = d.search * 1.1,
                                 bothways = FALSE, directions = TRUE)
      y.r    <- path.r[[y.name]]
      y.l    <- path.l[[y.name]]
      al.r   <- arclength(path.r$path)
      al.l   <- arclength(path.l$path)
      drns.r <- path.r$directions
      drns.l <- path.l$directions

      y.r    <- y.r[al.r < d.search]
      y.l    <- y.l[al.l < d.search]
      path.r <- path.r$path[al.r < d.search, ]
      path.l <- path.l$path[al.l < d.search, ]
      drns.r <- drns.r[ , , al.r < d.search]
      drns.l <- drns.l[ , , al.l < d.search]
      # path   <- rbind(resample.face3d(path.l, ndel)[ndel:1, ], resample.face3d(path.r, ndel)[-1, ])
      path   <- rbind(path.l[nrow(path.l):1, ], path.r[-1, ])
      al     <- arclength(path)
      mn     <- min(max(al.r[al.r < d.search]), max(al.l[al.l < d.search]))
      del    <- mn / ndel
      path   <- resample.face3d(path, 2 * ndel + 1)
      y      <- c(rev(y.l[-1]), y.ref[i], y.r[-1])
      y      <- approx(al, y, n = nh)$y
      drns   <- abind(drns.l[ , , rev(1:dim(drns.l)[3])], drns.r[ , , -1])
      drns   <- aperm(apply(drns, 1:2, function(x) approx(al, x, n = nh)$y), c(2, 3, 1))
      drns   <- drns[ , , (nasy + 1):(nh - nasy)]
      path   <- path[(nasy + 1):(nh - nasy), ]
      grd    <- rbind(grd, path)
      
      # asym.i <- numeric(0)
      # for (j in (nasy + 1):(nh - nasy))
      #    asym.i <- c(asym.i, sum((y[(j - nasy):(j + nasy)] - y[(j + nasy):(j - nasy)])^2))
      # if (si.ref[i] >= 5/8) asym.i <- asym.i * NA
      
      if ((abs(si.ref[i]) < -1) | (abs(si.ref[i]) > 1))
      # if ((si.ref[i] < 0) | (si.ref[i] > 5/8))
            asym.i <- rep(NA, nrow(path))
      else {
         asym.i <- numeric(0)
         for (j in 1:nrow(path)) {
            path.r <- planepath(sbst, path[j, ], direction =  drns[ , 2, j], rotation = 0, boundary = d.asym * 1.1,
                                       bothways = FALSE, directions = TRUE)
            path.l <- planepath(sbst, path[j, ], direction = -drns[ , 2, j], rotation = 0, boundary = d.asym * 1.1,
                                       bothways = FALSE, directions = TRUE)
            al.r   <- arclength(path.r$path)
            al.l   <- arclength(path.l$path)
            y.r    <- path.r[[y.name]]
            y.l    <- path.l[[y.name]]
            y.r    <- approx(al.r, y.r, (1:ndel) * d.asym / ndel)$y
            y.l    <- approx(al.l, y.l, (1:ndel) * d.asym / ndel)$y
            nmiss  <- length(which(is.na(y.r - y.l)))
            if (nmiss > 0) cat("Number of missing entries:", nmiss, "\n")
            # asym.i <- c(asym.i, sum((y.r - y.l)^2))
            asym.ij <- sum((y.r - y.l)^2, na.rm = TRUE)
            asym.i  <- c(asym.i, asym.ij)
            # cat("asym:", asym.i, "\n")
            if (any(is.na(asym.i))) {
               spheres3d(path.r$path, col = "red")
               spheres3d(path.l$path, col = "yellow")
               print(y.r)
               print(al.r)
               stop()
            }
         }
      }
      asym <- c(asym, asym.i)
      
      if (monitor) {
         if (!all(is.na(asym.i))) {
            spheres3d(path, col = topo.colors(20)[cut(asym.i, 20, labels = FALSE)], radius = 0.6)
            lines3d(path, col = topo.colors(20)[cut(asym.i, 20, labels = FALSE)], lwd = 5)
         }
         spheres3d(path[which.min(asym.i), ], col = "red", radius = 0.5)
      }
      # print(asym.i)
      # if (i == 1) stop()
   }
   # if (monitor.prompt) invisible(readline(prompt="Press [enter] to continue:"))
   
   xcoord  <- rep(0:(nv - 1), each = nh - 2 * nasy) * max(al) / (nv - 1)
   ycoord  <- rep((nasy + 1):(nh - nasy) - ndel, nv) * del
   values  <- -asym
   ind     <- which(is.na(values))
   xcoord0 <- xcoord
   ycoord0 <- ycoord
   grd0    <- grd
   values0 <- values
   if (length(ind) > 0) {
      xcoord  <- xcoord[-ind]
      ycoord  <- ycoord[-ind]
      values  <- values[-ind]
      grd     <- grd[-ind, ]
   }
   # if (monitor) {
   #    clr <- topo.colors(20)[cut(values, 20, labels = FALSE)]
   #    plot(xcoord, ycoord, col = clr, xlab = "arc length", ylab = "adjustment")
   # }
   values <- pmax(values, quantile(values, quantile))
   ridge <- ridge2d.face3d(xcoord, ycoord, values, lambda = lambda, df = df,
                           endpoints.fixed = FALSE, monitor = monitor)
   if (monitor) points(ridge$x, ridge$y, col = "red", pch = 16)
   
   # Project back onto the surface
   X              <- cbind(xcoord0, ycoord0)
   Xnew           <- cbind(ridge$x, ridge$y)
   interp.x       <- interp.barycentric(X, f = grd0[ , 1], Xnew)$fnew
   interp.y       <- interp.barycentric(X, f = grd0[ , 2], Xnew)$fnew
   interp.z       <- interp.barycentric(X, f = grd0[ , 3], Xnew)$fnew
   interp.values  <- interp.barycentric(X, f = values0, Xnew)$fnew
   path           <- cbind(interp.x, interp.y, interp.z)

   if (monitor) {
      # clr <- sbst[[y.name]]
      # if (y.name != "shape.index") clr <- clr + 1
      # plot(sbst, new = FALSE, col = clr)
      pop3d()
      spheres3d(path, col = "red", radius = 0.7)
      # spheres3d(reference, radius = 0.5)
   }

   invisible(list(path = path, x = xcoord, y = ycoord, z = values, subset = ind.sbst))
}
