# Checks on the plot function

library(Face3D)

liberty.data <- "~/Desktop/Data-Liberty/"
# liberty.data <- "/Volumes/Face3D/Glasgow-controls/Liberty/Liberty-Controls-dmp/Liberty-Controls-dmp-final/Data-Liberty/"
fls <- list.files(liberty.data, full.names = TRUE)

load(fls[1])
summary(face)

plot(face)
plot(face, new = FALSE, type = "points")
plot(face, new = FALSE, type = "mesh")
plot(face, new = FALSE, windowRect = c(500, 500, 1000, 1000))
plot(face, new = FALSE, col = "green")
plot(face, new = FALSE, col = "normal-x")
plot(face, new = FALSE, col = "normal-y")
plot(face, new = FALSE, col = "normal-z")

plot(face, new = FALSE, col = face$kappa1)
plot(face, new = FALSE, col = "kappa1")
plot(face, new = FALSE, col = face$kappa1 + 1)

# lines
load ("analysis/nose.RData")
nose$vertices  <- nose$coords
nose$triangles <- matrix(nose$triples, ncol = 3, byrow = TRUE)
nose$coords    <- NULL
nose$triples   <- NULL
nose$lines     <- matrix(c(55, rep(130:135, each = 2), 136,
                           109, rep(118:123, each = 2), rep(136, 2),
                           rep(129:124, each = 2), 1),
                         ncol = 2, byrow = TRUE)
plot(nose)
plot(nose, display = c("surface", "lines", "spheres"))

