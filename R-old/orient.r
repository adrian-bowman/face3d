orient.face3d <- function(face, ngrid = c(16, 4), angle.max = c(pi/4, pi/8),
                          nlevels = 1, graphics = FALSE) {

   if (!require(rgl)) stop("orient.face3d requires the rgl package.")

   if (!("normals" %in% names(face))) face <- normals.face3d(face)
   centroid    <- apply(face$coords, 2, mean)
   face$coords <- sweep(face$coords, 2, centroid)
   if ("lmks" %in% names(face)) {
      if (is.vector(face$lmks) && length(face$lmks) == 3) face$lmks <- as.matrix(face$lmks, ncol = 3)
      face$lmks    <- sweep(face$lmks, 2, centroid)
   }
   # parts       <- connected.face3d(face)
   # if (length(table(parts)) > 1) {
      # face1 <- face
      # face  <- subset(face, parts == 1)
   # }
   
   fn.nearest <- function(angles, dim = 1) {
   	  coords <- face$coords
   	  if (abs(angles[1]) > 0) coords  <- rotate3d(coords,  angles[1], 1, 0, 0)
   	  if (abs(angles[2]) > 0) coords  <- rotate3d(coords,  angles[2], 0, 1, 0)
   	  if (abs(angles[3]) > 0) coords  <- rotate3d(coords,  angles[3], 0, 0, 1)
      rngy    <- range(coords[ , 2])
      ind     <- (coords[ , 2] > 0.75 * rngy[1] + 0.25 * rngy[2]) &
                 (coords[ , 2] < 0.25 * rngy[1] + 0.75 * rngy[2])
      coords1 <- coords[ind, ]
      # face1   <- face
      # face1$coords <- coords
      # display.face3d(subset.face3d(face1, ind))
      nstrip  <- 8
      strip   <- cut(coords1[ , 2], nstrip, labels = FALSE)
      # brks    <- seq(0.80 * rngy[1] + 0.20 * rngy[2], 
                     # 0.20 * rngy[1] + 0.80 * rngy[2], length = nstrip + 1)
      # strip   <- cut(coords[ , 2], brks, labels = FALSE)
      nearest <- tapply(1:nrow(coords1), strip, function(x) x[which.max(coords1[x, 3])])
      if (graphics) {
         face1   <- face
         face1$coords <- coords
         display.face3d(face1, new = FALSE)
         spheres3d(coords1[nearest, ], col = "green", radius = 2)
         print(c(angles, sd(coords1[nearest, dim])))
      }
      if (length(nearest) < nstrip) nearest <- c(nearest, rep(NA, nstrip - length(nearest)))
      stdev <- sd(coords1[nearest, 1]) # - sd(coords1[nearest, 3])
      nearest <- (1:nrow(coords))[which(ind)[nearest]]
      c(angles, stdev, nearest)
   }
   
   fn.z <- function(angles) {
   	  normals <- face$normals
      if (abs(angles[1]) > 0) normals <- rotate3d(normals, angles[1], 1, 0, 0)
      if (abs(angles[2]) > 0) normals <- rotate3d(normals, angles[2], 0, 1, 0)
      if (abs(angles[3]) > 0) normals <- rotate3d(normals, angles[3], 0, 0, 1)
      # print(c(angles, sum(normals[ , 3])))
      c(angles, sum(normals[ , 3]))
   }
   
   fn.asymmetry <- function(angles) {
   	  normals <- face$normals
   	  coords  <- face$coords
      if (abs(angles[1]) > 0) normals <- rotate3d(normals, angles[1], 1, 0, 0)
      if (abs(angles[2]) > 0) normals <- rotate3d(normals, angles[2], 0, 1, 0)
      if (abs(angles[3]) > 0) normals <- rotate3d(normals, angles[3], 0, 0, 1)
      if (abs(angles[1]) > 0) coords  <- rotate3d(coords,  angles[1], 1, 0, 0)
      if (abs(angles[2]) > 0) coords  <- rotate3d(coords,  angles[2], 0, 1, 0)
      if (abs(angles[3]) > 0) coords  <- rotate3d(coords,  angles[3], 0, 0, 1)
      rngy               <- range(coords[ , 2])
      ind                <- (coords[ , 3] > mean(coords[ , 3])) &
                            (coords[ , 2] > 0.80 * rngy[1] + 0.20 * rngy[2]) &
                            (coords[ , 2] < 0.20 * rngy[1] + 0.80 * rngy[2])
      coords             <- coords[ind, ]
      normals            <- normals[ind, ]
      brks               <- seq(min(coords[ , 1]), max(coords[ , 1]), by = 2)
      brks[1]            <- brks[1] - 1
      brks[length(brks)] <- brks[length(brks)] + 1
      nstrip             <- length(brks) - 1
      vstrip             <- cut(coords[ , 1], brks)
      brks               <- seq(min(coords[ , 2]), max(coords[ , 2]), by = 6)
      brks[1]            <- brks[1] - 1
      brks[length(brks)] <- brks[length(brks)] + 1
      hstrip             <- cut(coords[ , 2], brks)
      grid1              <- tapply(abs(normals[ , 1]), list(vstrip, hstrip), mean)
      nw                 <- 15
      fn <- function(i, grid) {
      	       gdiff  <- abs(grid[i + (-nw:nw), ] - grid[i - (-nw:nw), ])
      	       gdiffm <-  mean(colMeans(gdiff), na.rm = TRUE)
      	       gdiffm
               }
      fval1 <- apply(as.matrix((1 + nw):(nstrip - nw)), 1, fn, grid1)
      if (graphics) {
         face1         <- face
         face1$coords  <- coords
         face1$normals <- normals
      	 print(summary.face3d(face1))
         display.face3d(face1, type = "points", new = FALSE, colour = "normal-x")
         plot((1 + nw):(nstrip - nw), fval1, ylim = c(0, 0.6))
         title(paste(angles, collapse = " "))
         # scan()
      }
      # print(fval1)
      # print(c(angles, min(fval1, na.rm = TRUE)))
      c(angles, min(fval1, na.rm = TRUE))
   }
   
   # Search over x and y rotations
   angle1 <- 0
   angle2 <- 0
   angle3 <- 0
   for (i in 1:2) {
   	  scl     <- if (i == 1) 1 else ngrid[1]
      angles  <- seq(-angle.max[1], angle.max[1], length = ngrid[i]) / scl
      angles  <- as.matrix(expand.grid(angle1 + angles, angle2 + angles, 0))
      results <- t(apply(angles, 1, fn.z))
      ind     <- which.max(results[ , 4])
      angle1  <- results[ind, 1]
      angle2  <- results[ind, 2]
   }
   face$coords    <- rotate3d(face$coords,  angle1, 1, 0, 0)
   face$coords    <- rotate3d(face$coords,  angle2, 0, 1, 0)
   face$normals   <- rotate3d(face$normals, angle1, 1, 0, 0)
   face$normals   <- rotate3d(face$normals, angle2, 0, 1, 0)
   if ("lmks" %in% names(face)) {
   	  face$lmks <- rotate3d(face$lmks, angle1, 1, 0, 0)
   	  face$lmks <- rotate3d(face$lmks, angle2, 0, 1, 0)
   }
   
   if (graphics) {
      display.face3d(face, new = FALSE)
      }

   # Search over y and z rotations
   if (nlevels > 1) {
      angle2 <- 0
      angle3 <- 0
      for (i in 1:2) {
      	  scl     <- if (i == 1) 1 else ngrid[1]
         anglesy <- seq(-angle.max[2], angle.max[2], length = ngrid[i]) / scl
         anglesz <- seq(-angle.max[1], angle.max[1], length = ngrid[i]) / scl
         angles  <- as.matrix(expand.grid(0, angle2 + anglesy, angle3 + anglesz))
         results <- t(apply(angles, 1, fn.nearest))
         # print(results)
         ind     <- which.min(results[ , 4])
         angle2  <- results[ind, 2]
         angle3  <- results[ind, 3]
      }
      face$coords    <- rotate3d(face$coords,  angle2, 0, 1, 0)
      face$coords    <- rotate3d(face$coords,  angle3, 0, 0, 1)
      face$normals   <- rotate3d(face$normals, angle2, 0, 1, 0)
      face$normals   <- rotate3d(face$normals, angle3, 0, 0, 1)
      if ("lmks" %in% names(face)) {
   	     face$lmks <- rotate3d(face$lmks, angle2, 0, 1, 0)
   	     face$lmks <- rotate3d(face$lmks, angle3, 0, 0, 1)
      }
   }

   face$nearest   <- round(fn.nearest(rep(0, 3))[-(1:4)])
   face$coords    <- sweep(face$coords, 2, centroid, "+")
   if ("lmks" %in% names(face)) {
   	  if (!is.matrix(face$lmks)) face$lmks <- t(face$lmks)
      face$lmks <- sweep(face$lmks, 2, centroid, "+")
   }
   face$rotations <- c("x-axis" = angle1, "y-axis" = angle2, "z-axis" = angle3)
   
   invisible(face)
}

# orientnorm.face3d <- function(face, ngrid = 21, grid = 5, threshold = 0.4) {

   # results <- matrix(nrow = 0, ncol = 4)
   # for (angle1 in seq(-pi/4, pi/4, length = ngrid)) {
   # for (angle2 in seq(-pi/4, pi/4, length = ngrid)) {
   	  # coords  <- rotate3d(face$coords,  angle1, 0, 1, 0)
   	  # coords  <- rotate3d(     coords,  angle2, 0, 0, 1)      
      # normals <- rotate3d(face$normals, angle1, 0, 1, 0)
      # normals <- rotate3d(     normals, angle2, 0, 0, 1)
      # dir.x   <- abs(normals[ , 1])
      # dir.z   <-     normals[ , 3]
      # # clr     <- topo.colors(10)[cut(dir.x, 10, labels = FALSE)]
      # # clr     <- rep("grey", length(dir.x))
      # # clr[dir < threshold] <- "blue"
      # # eqscplot(coords[ , 1:2], col = clr, pch = ".")
      # # title(angle)
      # brks    <- seq(min(coords[ , 1]), max(coords[ , 1]), by = grid)
      # brks    <- c(brks, max(coords[ , 1]) + grid)
      # strip   <- cut(coords[ , 1], brks, labels = FALSE, include.lowest = TRUE)
      # # hist(abs(blue[ , 1] - median(blue[ , 1])), breaks = 20)
      # blue    <- (dir.x < threshold & dir.z > 0)
      # tbl     <- table(strip, blue)
      # ind     <- as.numeric(rownames(tbl)[which.max(tbl[ , 2])])
      # cat(angle1, angle2, tbl[ind, 2], ind, "\n")
      # results <- rbind(results, c(angle1, angle2, tbl[ind, 2], ind))
      # # abline(v = (brks[ind] + brks[ind + 1]) / 2, col = "red")
      # # abline(v = median(blue[ , 1]), col = "green")
      # # print(c(angle, mean(abs(blue[ , 1] - median(blue[ , 1])))))
      # # plot(abs(blue[ , 1] - median(blue[ , 1])), blue[ , 2])
      # # plot(blue[ , 1], blue[ , 2], pch = ".")
      # # scan()
# }
# }
   # image(seq(-pi/4, pi/4, length = ngrid),
         # seq(-pi/4, pi/4, length = ngrid),
         # t(matrix(results[ , 3], ncol = ngrid)))

   # ind            <- which.max(results[ , 3])
   # angles         <- results[ind, 1:2]
   # face$coords    <- rotate3d(face$coords,  angles[1], 0, 1, 0)
   # face$coords    <- rotate3d(face$coords,  angles[2], 0, 0, 1)
   # face$normals   <- rotate3d(face$normals, angles[1], 0, 1, 0)
   # face$normals   <- rotate3d(face$normals, angles[2], 0, 0, 1)
   # face$rotations <- c("y-axis" = angles[1], "z-axis" = angles[2], "x-axis" = 0)
   # brks           <- seq(min(face$coords[ , 1]), max(face$coords[ , 1]), by = grid)
   # brks           <- c(brks, max(face$coords[ , 1]) + grid)
   # ind1           <- results[ind, 4]
   # face$midline   <- (brks[ind1] + brks[ind1 + 1]) / 2

   # invisible(face)
# }

