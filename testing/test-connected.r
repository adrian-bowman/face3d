#     Tests for connected.face3d

install.packages("~/OneDrive - University of Glasgow/research/face3d_0.1-1/face3d",
                 repos = NULL, type = "source")
library(face3d)
library(rgl)


sbst  <- subset(template_male, template_male$shape.index > 0)
parts <- connected.face3d(sbst)
sbst  <- subset(sbst, parts == 1)
plot(sbst, col = "kappa2")



load("testing/test-data/connected.RData")
plot(face, display = "points", new = FALSE)
plot(face, display = "mesh", add = TRUE)

parts <- connected.face3d(face)
table(parts)

spheres3d(matrix(face$coords[38027, ], ncol = 3),
      radius = 5, col = "green")

face1 <- subset.face3d(face, parts == 1)
plot(face1, colour = "blue")
face2 <- subset.face3d(face, parts == 2)
plot(face2, colour = "red", add = TRUE)
spheres3d(face$coords[parts == 0, ], col = "yellow", radius = 5)
