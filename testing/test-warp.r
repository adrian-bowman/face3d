#     Test of the warp.face3d function

library(Face3D)
library(rgl)

# Warp a template to an image

plot(template.female)
mtof <- warp.face3d(template.male$curves, template.female$curves, template.male)
plot(mtof)
# Why doesn't the step below work?
dst <- distance.face3d(template.female, mtof)
plot(template.female, col = dst$normal)

# Warp a scalar value

plot(face)
face <- normals.face3d(face)

ind   <- sample(1:nrow(face$coords), 500)
face1 <- face
face1 <- index.face3d(face1, subset = ind, distance = 50, overwrite = TRUE)
ind   <- ind[which(!is.na(face1$shape.index[ind]))]
to    <- face1$shape.index[ind]
face1$shape.index <- warp.face3d(face1$coords[ind, ], to, face1$coords)
face1$shape.index <- pmax(face1$shape.index, -1)
face1$shape.index <- pmin(face1$shape.index,  1)
plot(face1, new = FALSE, col = "shape index")
spheres3d(face1$coords[ind, ])

# Warp two scalar values

to  <- cbind(face1$kappa1[ind], face1$kappa2[ind])
k12 <- warp.face3d(face1$coords[ind, ], to, face1$coords)
plot(face1, new = FALSE, col = k12[ , 1])
plot(face1, new = FALSE, col = k12[ , 2])

