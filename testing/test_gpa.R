#     Tests for the gpa.face3d function

library(face3d)
if (reinstall) devtools::install("face3d")

library(rgl)
library(shapes)
library(abind)

test_label("gpa", test.prompt)
X      <- abind(template_male$mesh, template_female$mesh, along = 3)
trngls <- as.face3d(template_male$mesh, model.mesh = TRUE)$triangles
ind    <- which(apply(trngls, 1, function(x) any(is.na(x))))
trngls[ind, ]
gp     <- gpa.face3d(X, trngls)
spheres3d(template_male$mesh)
spheres3d(template_male$curves, radius = 1.2, col = "blue")

nvert <- nrow(template_male$vertices)
X <- array(template_male$vertices, dim = c(nvert, 3, 3))
X <- X + array(rnorm(nvert * 3 * 3), dim = c(nvert, 3, 3))
gp <- gpa.face3d(X, template_male$triangles)
points3d(gp$aligned[ , , 1])
points3d(gp$aligned[ , , 2], col = 2)
points3d(gp$aligned[ , , 3], col = 3)
apply(sweep(X, c(1, 3), gp$weights), 2:3, mean)
