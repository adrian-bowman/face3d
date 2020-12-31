orient.face3d <- function(face, ngrid = c(16, 4), angle.max = c(pi/4, pi/8),
                          nlevels = 1, graphics = FALSE) {

   if (!requireNamespace("rgl", quietly = TRUE)) stop("orient.face3d requires the rgl package.")

   if (!("normals" %in% names(face))) face <- normals.face3d(face)

   centroid <- apply(face$vertices, 2, mean)
   vertices   <- sweep(face$vertices, 2, centroid)
   normals  <- face$normals

   fn.nearest <- function(angles, dim = 1) {
   	if (abs(angles[1]) > 0) vertices  <- rgl::rotate3d(vertices,  angles[1], 1, 0, 0)
   	if (abs(angles[2]) > 0) vertices  <- rgl::rotate3d(vertices,  angles[2], 0, 1, 0)
   	if (abs(angles[3]) > 0) vertices  <- rgl::rotate3d(vertices,  angles[3], 0, 0, 1)
      rngy    <- range(vertices[ , 2])
      ind     <- (vertices[ , 2] > 0.75 * rngy[1] + 0.25 * rngy[2]) &
                 (vertices[ , 2] < 0.25 * rngy[1] + 0.75 * rngy[2])
      vertices1 <- vertices[ind, ]
      nstrip  <- 8
      strip   <- cut(vertices1[ , 2], nstrip, labels = FALSE)
      nearest <- tapply(1:nrow(vertices1), strip, function(x) x[which.max(vertices1[x, 3])])
      if (graphics) {
         face1   <- face
         face1$vertices <- vertices
         plot(face1, new = FALSE)
         rgl::spheres3d(vertices1[nearest, ], col = "green", radius = 2)
         print(c(angles, sd(vertices1[nearest, dim])))
      }
      if (length(nearest) < nstrip) nearest <- c(nearest, rep(NA, nstrip - length(nearest)))
      stdev   <- sd(vertices1[nearest, 1]) # - sd(vertices1[nearest, 3])
      nearest <- (1:nrow(vertices))[which(ind)[nearest]]
      c(angles, stdev, nearest)
   }
   
   fn.z <- function(angles) {
      normals <- normals + vertices
      if (abs(angles[1]) > 0) vertices  <- rgl::rotate3d(vertices,  angles[1], 1, 0, 0)
   	if (abs(angles[2]) > 0) vertices  <- rgl::rotate3d(vertices,  angles[2], 0, 1, 0)
   	if (abs(angles[3]) > 0) vertices  <- rgl::rotate3d(vertices,  angles[3], 0, 0, 1)
   	if (abs(angles[1]) > 0) normals <- rgl::rotate3d(normals, angles[1], 1, 0, 0)
   	if (abs(angles[2]) > 0) normals <- rgl::rotate3d(normals, angles[2], 0, 1, 0)
   	if (abs(angles[3]) > 0) normals <- rgl::rotate3d(normals, angles[3], 0, 0, 1)
   	normals <- normals - vertices
      # print(c(angles, sum(normals[ , 3])))
      c(angles, sum(normals[ , 3]))
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
   normals <- normals + vertices
   vertices  <- rgl::rotate3d(vertices,  angle1, 1, 0, 0)
   vertices  <- rgl::rotate3d(vertices,  angle2, 0, 1, 0)
   normals <- rgl::rotate3d(normals, angle1, 1, 0, 0)
   normals <- rgl::rotate3d(normals, angle2, 0, 1, 0)
   normals <- normals - vertices
   
   if (graphics) {
      face1 <- list(vertices = vertices, triangles = face$triangles)
      plot(face, new = FALSE)
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
         ind     <- which.min(results[ , 4])
         angle2  <- results[ind, 2]
         angle3  <- results[ind, 3]
      }
      normals <- normals + vertices
      vertices  <- rgl::rotate3d(vertices,  angle2, 0, 1, 0)
      vertices  <- rgl::rotate3d(vertices,  angle3, 0, 0, 1)
      normals <- rgl::rotate3d(normals, angle2, 0, 1, 0)
      normals <- rgl::rotate3d(normals, angle3, 0, 0, 1)
      normals <- normals - vertices
   }

   face$nearest         <- round(fn.nearest(rep(0, 3))[-(1:4)])
   face$centroid        <- centroid
   face$rotation.angles <- c("x-axis" = angle1, "y-axis" = angle2, "z-axis" = angle3)
   face$orient <- function(x, centroid, rotation.angles) {
      x <- sweep(x, 2, centroid)
      x <- rgl::rotate3d(x,  rotation.angles[1], 1, 0, 0)
      x <- rgl::rotate3d(x,  rotation.angles[2], 0, 1, 0)
      x <- rgl::rotate3d(x,  rotation.angles[3], 0, 0, 1)
      sweep(x, 2, centroid, "+")
   }
   
   invisible(face)
}

# orientnorm.face3d <- function(face, ngrid = 21, grid = 5, threshold = 0.4) {

   # results <- matrix(nrow = 0, ncol = 4)
   # for (angle1 in seq(-pi/4, pi/4, length = ngrid)) {
   # for (angle2 in seq(-pi/4, pi/4, length = ngrid)) {
   	  # vertices  <- rotate3d(face$vertices,  angle1, 0, 1, 0)
   	  # vertices  <- rotate3d(     vertices,  angle2, 0, 0, 1)      
      # normals <- rotate3d(face$normals, angle1, 0, 1, 0)
      # normals <- rotate3d(     normals, angle2, 0, 0, 1)
      # dir.x   <- abs(normals[ , 1])
      # dir.z   <-     normals[ , 3]
      # # clr     <- topo.colors(10)[cut(dir.x, 10, labels = FALSE)]
      # # clr     <- rep("grey", length(dir.x))
      # # clr[dir < threshold] <- "blue"
      # # eqscplot(vertices[ , 1:2], col = clr, pch = ".")
      # # title(angle)
      # brks    <- seq(min(vertices[ , 1]), max(vertices[ , 1]), by = grid)
      # brks    <- c(brks, max(vertices[ , 1]) + grid)
      # strip   <- cut(vertices[ , 1], brks, labels = FALSE, include.lowest = TRUE)
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
   # face$vertices    <- rotate3d(face$vertices,  angles[1], 0, 1, 0)
   # face$vertices    <- rotate3d(face$vertices,  angles[2], 0, 0, 1)
   # face$normals   <- rotate3d(face$normals, angles[1], 0, 1, 0)
   # face$normals   <- rotate3d(face$normals, angles[2], 0, 0, 1)
   # face$rotations <- c("y-axis" = angles[1], "z-axis" = angles[2], "x-axis" = 0)
   # brks           <- seq(min(face$vertices[ , 1]), max(face$vertices[ , 1]), by = grid)
   # brks           <- c(brks, max(face$vertices[ , 1]) + grid)
   # ind1           <- results[ind, 4]
   # face$midline   <- (brks[ind1] + brks[ind1 + 1]) / 2

   # invisible(face)
# }

