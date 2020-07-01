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
