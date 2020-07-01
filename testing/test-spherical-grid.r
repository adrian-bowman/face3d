# Identify features from a low resolution grid

setwd("/Volumes/adrian/research/shape/Face3D/Face3D-package/Face3D_0.1-1")

library(rgl)

source("spherical-grid.r")
source("Face3D/R/crossproduct.r")
grid <- spherical.grid()
points3d(grid$points)

ind  <- 1
pt   <- grid$points[ind, ]
spheres3d(matrix(pt, ncol = 3), radius = 0.03, col = "green")
ind1 <- grid$nbr[1, ]
ind1 <- ind1[!is.na(ind1)]
pts1 <- grid$points[ind1, ]
apply(pts1, 1, function(x) sqrt(sum(x^2)))
for (i in 1:nrow(pts1)) {
   spheres3d(matrix(pts1[i, ], ncol = 3), radius = 0.03, col = "red")
   scan()
}
pts1 <- sweep(pts1, 2, pt)
pt1  <- pts1[ 1, ]
pts1 <- pts1[-1, ]
dp   <- apply(pts1, 1, function(x) sum(x * pt1))
cp   <- crossproduct(pts1, matrix(pt1, ncol = 3, nrow = nrow(pts1),
              byrow = TRUE))
dir  <- apply(cp, 1, function(x) sum(x * pt))


