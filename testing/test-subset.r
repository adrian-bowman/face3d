#     Test code for subset.face3d

setwd("/Volumes/adrian/research/shape/Face3D/Face3D-package/Face3D_0.1-1/Face3D/R")
# setwd("~/Desktop/Face3D-package/Face3D_0.1-1/Face3D/R")

library(Face3D)

source("subset.r")

data(face)
display.face3d(face)

display.face3d(face, subset = face$coords[ , 1] < 0, new = FALSE)
display.face3d(face, subset = -(1:500), new = FALSE)


face1 <- subset(face, face$coords[ , 1] > 0)
display.face3d(face1)

display.face3d(face)
fn    <- select3d()
ind   <- which(fn(face$coords))
face2 <- subset.face3d(face, ind)
display.face3d(face2)

detach(package:Face3D, unload = TRUE)
