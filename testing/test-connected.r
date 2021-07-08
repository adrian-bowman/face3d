#     Tests for connected.face3d

install.packages("~/research/face3d/face3d", repos = NULL, type = "source")
library(face3d)
library(rgl)

sbst  <- subset(template_male,
                template_male$vertices[ , 1] > 10 | template_male$vertices[ , 1] < 5)
plot(sbst)
parts <- connected.face3d(sbst)
summary(sbst)
table(parts)
plot(subset(sbst, parts == 1), col = "grey")
plot(subset(sbst, parts == 2), col = "green", add = TRUE)




load("~/research/face3d/testing/test-data/connected.RData")
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
