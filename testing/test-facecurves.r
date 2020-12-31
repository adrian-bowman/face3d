library(Face3D)
library(fields)
# creating geodesics for liberty controls
setwd("/Volumes/Face3D/Glasgow-controls/Liberty/liberty-dmp")


#load("controls-liberty-065.dmp")
#load("controls-liberty-067.dmp")
#load("controls-liberty-069.dmp")

#ind <- which(rownames(face$lmks) == "pm")
#rownames(face$lmks)[ind] <- "pn"

setwd("/Volumes/adrian/research/shape/Face3D/Face3D-package/Face3D_0.1-1/Face3D/R")
for (file in list.files()) source(file)


display.face3d(face)
id <- rgl.cur()

rgl.set(id)
display.face3d(face, new = FALSE)
spheres3d(face$lmks, col = "red", radius = 1.5)

face <- facecurves.face3d(face, "mid-line lip", graphics = FALSE)
face <- facecurves.face3d(face, "upper lip right", graphics = TRUE)
face <- facecurves.face3d(face, "upper lip left", graphics = TRUE)
face <- facecurves.face3d(face, "mandible right", graphics = TRUE)
face <- facecurves.face3d(face, "brow ridge right", graphics = TRUE)

spheres3d(face$lmks[c("cphL", "chL"), ], col = "red", radius = 1.5)

rgl.set(id)
open3d()
spheres3d(face$curves, col = "blue", radius = 0.5)
spheres3d(face$curves, col = "green", radius = 0.5)


pop3d()
x1 <- closest.face3d(c(0, -125,  90), face)$point
x1 <- closest.face3d(x1 + c(15, 5, 0), face)$point
x2 <- closest.face3d(c(-30, -60, -5), face)$point
spheres3d(rbind(x1, x2), col = "red", radius = 1.5)

m <- matrix(c(rbind(x1, x2)), ncol = 3)
rownames(m) <- c("gn", "oiR")
face$lmks <- m
face <- facecurves.face3d(face, "mandible right", graphics = TRUE)
spheres3d(face$curves, col = "blue", radius = 0.5)
spheres3d(face$curves, col = "green", radius = 1)

face1 <- facecurves.face3d(face, monitor = TRUE)
spheres3d(face1$curves, col = "green", radius = 1)

n.vec <- c(10,10,10,10,15,15,5,20,10,10,5,5,5,10,10,10,20,20,20,20,15,15,10,10)
geodesics.liberty <- resamplecurves.face3d(face, n.vec = n.vec)
