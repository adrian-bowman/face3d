"distance.face3d" <- function(x1, x2, minsum = FALSE) {
  
   if (is.vector(x1)) x1 <- t(as.matrix(x1))
   if (is.vector(x2)) x2 <- t(as.matrix(x2))
   
   if (is.matrix(x1) | is.matrix(x2)) {
      if (is.face3d(x1)) x1 <- x1$vertices
      if (is.face3d(x2)) x2 <- x2$vertices
      if (ncol(x1) != 3) stop("the dimension of x1 do not correspond to 3d co-ordinates.")
      if (ncol(x2) != 3) stop("the dimension of x2 do not correspond to 3d co-ordinates.")
      if (minsum)
         dst <- sum(c(apply(x1, 1, function(x) min(fields::rdist(t(x), x2)))))
      else {
         dst <- fields::rdist(x1, x2)
         if (any(dim(dst) == 1)) dst <- c(dst)
      }
      return(dst)
   }
   
   if (!(is.face3d(x1) & is.face3d(x2)))
     stop("x1 and x2 must both be face3d objects.")
   
  if (!(nrow(x1$vertices) == nrow(x2$vertices)))
     stop("x1 and x2 have different numbers of vertices.")

   if (!all("normals" %in% names(x1)))
     x1 <- normals.face3d(x1)

  ## area of a triangle
  # log.area.ratio   <- log(area.face3d(x2)$triangles / area.face3d(x1)$triangles)
  
  ## Distance in normal direction
  dist.mat         <- x2$vertices - x1$vertices
  dist.euclid      <- sqrt(apply(dist.mat^2, 1, sum))
  dist.normal      <- apply(x1$normals * dist.mat, 1, sum)

  # # Old version
  # ## 3d distance and sign distances with respect to x,y,z direction
  # ## 3d sign distance with respect to x direction
  # sign.lr.x        <- sign(dist.mat[,1]) # left-right on the screen
  # k2               <- dim(x1$vertices)[1]
  # index.dir        <- rep(1, k2)
  # index            <- which((x1$vertices[,1] + x2$vertices[,1]) / 2 <= 0)
  # index.dir[index] <- -1
  # sign.lr.x        <- sign.lr.x * index.dir
  # dist.x           <- sign.lr.x * abs(dist.mat[,1])
  # ## 3d sign distance with respect to y direction
  # sign.lr.y        <- sign(dist.mat[,2]) # up-down on the screen
  # dist.y           <- sign.lr.y * abs(dist.mat[,2])
  # ## 3d sign distance with respect to z direction
  # sign.lr.z        <- sign(dist.mat[,3]) # in-out of the screen
  # dist.z           <- sign.lr.z * abs(dist.mat[,3])

  # New version
  #sign.lr.x <- sign(dist.mat[, 1])
  k2 <- dim(x1$vertices)[1]
  #index.dir <- rep(1, k2)
  #index <- which((x1$vertices[, 1] + x2$vertices[, 1])/2 <= 0)
  #index.dir[index] <- -1
  #sign.lr.x <- sign.lr.x * index.dir
  # dist.x <- sign.lr.x * abs(dist.mat[, 1])
  sign.lr.x <- sign(dist.mat[,1])
  dist.x <- sign.lr.x * abs(dist.mat[, 1])
  # dist.x <- dist.mat[, 1]
  sign.lr.y <- sign(dist.mat[, 2])
  dist.y <- sign.lr.y * abs(dist.mat[, 2])
  sign.lr.z <- sign(dist.mat[, 3])
  dist.z <- sign.lr.z * abs(dist.mat[, 3])
  
  result           <- list(x = dist.x, y = dist.y, z = dist.z, xyz = dist.euclid,
                           # log.area.ratio = log.area.ratio,
                           normal = dist.normal)
  return(result)
}
