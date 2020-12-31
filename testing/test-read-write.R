setwd("~/ownCloud/Face3D_0.1-1/")

library(Face3D)

#       ply files
face <- read.face3d("test-data/mesh.ply")
summary(face)
plot(face)

face <- read.face3d("test-data/mesh-binary.ply", monitor = TRUE)
summary(face)
plot(face)

