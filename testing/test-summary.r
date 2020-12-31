#     Test of the summary.face3d function

install.packages("~/ownCloud/Face3D_0.1-1/Face3D", repos = NULL, type = "source")
library(Face3D)

summary.face3d(template)

# Check for isolated points

load("test-data/isolated.RData")
plot(shape,colour = "grey")
spheres3d(shape$coords[outer.idf, ], col = "red")
newshape <- subset(shape,-outer.idf, remove.singles = FALSE)
smy <- summary(newshape, checks = TRUE)
smy$isolated
nrow(shape$coords)
nrow(newshape$coords) + length(outer.idf)
plot(newshape)
spheres3d(newshape$coords, col = "blue")
newpath  <- connected.face3d(newshape)
subshape <- subset.face3d(newshape, newpath==2)
points3d(subshape$coords)
spheres3d(shape$coords[lostid,],col=2,radius = 0.2)

