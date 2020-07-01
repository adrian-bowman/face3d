# Locate the face by approximating the shape index from a random sample of locations
# and using warping to spread this out across the surface.  Using a large distance
# for shape index allows the ball shape of the face to be identified at large scale.

# Try this with a wider selection of faces.

# Put a check on -1 <= si <= 1 in plot function.

library(Face3D)
library(rgl)

plot(face)
face <- normals.face3d(face)

ind    <- sample(1:nrow(face$coords), 500)
face1  <- face
face1  <- index.face3d(face1, subset = ind, distance = 50, overwrite = TRUE)
ind    <- ind[which(!is.na(face1$shape.index[ind]))]
to     <- cbind(face1$kappa1[ind], face1$kappa2[ind], face1$shape.index[ind])
result <- warp.face3d(face1$coords[ind, ], to, face1$coords)
face1$kappa1 <- result[ , 1]
face1$kappa2 <- result[ , 2]
face1$shape.index <- result[ , 3]
face1$shape.index <- pmax(face1$shape.index, -1)
face1$shape.index <- pmin(face1$shape.index,  1)
plot(face1, new = FALSE, col = face1$kappa1)
plot(face1, new = FALSE, col = face1$kappa2)
plot(face1, new = FALSE, col = "shape index")
plot(face1, new = FALSE, col = face1$kappa1 * face1$kappa2)
spheres3d(face1$coords[ind, ])


sbst <- subset(face1, face1$shape.index > 0.5)
plot(face1, new = FALSE, col = "shape index")
plot(sbst, add= TRUE, display = "spheres", col = "red")

parts <- connected.face3d(sbst)
sbst1 <- subset(sbst, parts == 1)

ind   <- sample(1:nrow(sbst1$coords), 500)
sbst1 <- index.face3d(sbst1, subset = ind, distance = 20, overwrite = TRUE)
ind   <- ind[which(!is.na(sbst1$shape.index[ind]))]
to    <- sbst1$shape.index[ind]
to    <- cbind(to, to, to)
sbst1$shape.index <- warp.face3d(sbst1$coords[ind, ], to, sbst1$coords)[[1]][ , 1]
sbst1$shape.index <- pmax(sbst1$shape.index, -1)
sbst1$shape.index <- pmin(sbst1$shape.index,  1)
plot(sbst1, new = FALSE, col = "shape index")
