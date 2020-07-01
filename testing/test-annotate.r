#     Finding the nearest neighbours of each point on a shape

setwd("/Volumes/adrian/research/shape/Face3D/Face3D-package/Face3D_0.1-1")

library(Face3D)
data(face)

nose <- subset.face3d(face, face$coords[ , 3] > 85 & 
           face$coords[ , 2] > -60 & face$coords[ , 2] < -15)
nose <- index.face3d(nose, extent = 2)

display.face3d(nose, type = "mesh", new = FALSE)
source("Face3D/R/annotate.r")
annotate.face3d(nose)
