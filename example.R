library(face3d)
library(rgl)
library(fields)
library(geometry)

load("example.dmp")
plot(shape)
# plot(shape, colour = "shape index")
spheres3d(shape$lmks[14,],col=2)
spheres3d(shape$lmks[16,],col=1)
lip_curve <- planepath(shape,shape$lmks[16,],shape$lmks[14,])
spheres3d(lip_curve$path)
