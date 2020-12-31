#     Finding the nearest neighbours of each point on a shape

setwd("/Volumes/adrian/research/shape/Face3D/Face3D-package/Face3D_0.1-1")

library(Face3D)

source("Face3D/R/crossproduct.r")

crossproduct(1:3, 4:6, scale = FALSE)
a <- matrix(1:12, ncol = 3)
b <- matrix(12:1, ncol = 3)
crossproduct(a, b)
