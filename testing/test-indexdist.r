## Face Data ##

# Load libraries and data #

library(LPCM)
library(rgl)
library(Face3D)
library(princurve)
library(sm)

#  Create a template from kathryn001

load("/Volumes/Face3D/Glasgow-controls/Kathryn/glasgow-controls-dmp-final/kathryn001.dmp")
# display.face3d(face)
nms <- rownames(face$curves)
ind <- c(grep("upper lip", nms), grep("lower lip", nms),
         grep("mid lip", nms))
ind <- ind[-grep("mid-line", nms[ind])]
# points3d(face$curves[ind, ])
template.curves.mouth <- face$curves[ind, ]

select.lips <- (face$coords[,2] < -60) & (face$coords[,2] > -95) & (face$coords[,1] > -10) & (face$coords[,1] < 50)
lips <- subset.face3d(face, select.lips)
# display.face3d(lips)

source("/Volumes/adrian/research/shape/Face3D/Face3D-package/Face3D_0.1-1/Face3D/R/index-dist.r")
lips.si <- indexdist.face3d(lips, distance = 4)$shape.index
display.face3d(lips, colour = lips.si)
template.curves.mouth <-
   template.curves.mouth[-c(1,11,21,31,41,51, 10, 20,30,40,50,60), ]
spheres3d(template.curves.mouth, radius = 0.5)
lmks.mouth <- c("chR", "chL", "ls", "cphL", "cphR", "li", "st")
template.lmks.mouth <- face$lmks[lmks.mouth, ]
spheres3d(template.lmks.mouth, radius = 0.7, col = "red")


#   Warp the template to another case

load("/Volumes/Face3D/Glasgow-controls/Kathryn/glasgow-controls-dmp-final/kathryn002.dmp")
display.face3d(face)

select.lips <- (face$coords[,2] < -60) & (face$coords[,2] > -95) & (face$coords[,1] > -50) & (face$coords[,1] < 25)
lips <- subset.face3d(face, select.lips)
face.lmks.mouth <- face$lmks[lmks.mouth, ]
lips.si <- indexdist.face3d(lips, distance = 4)$shape.index
display.face3d(lips, colour = lips.si, type = "mesh")
spheres3d(face.lmks.mouth, radius = 0.7, col = "red")

face.curves.mouth <- warp.face3d(template.curves.mouth,
                      template.lmks.mouth, face.lmks.mouth)
spheres3d(face.curves.mouth, radius = 0.5)
for (i in 1:nrow(face.curves.mouth))
   face.curves.mouth[i, ] <- closest.face3d(face.curves.mouth[i, ],
                                     lips)$point
spheres3d(face.curves.mouth, radius = 0.5)


#   Planepath algorithm

load("/Volumes/Face3D/Glasgow-controls/Kathryn/glasgow-controls-dmp-final/kathryn002.dmp")
display.face3d(face)

select.lips <- (face$coords[,2] < -60) & (face$coords[,2] > -95) & (face$coords[,1] > -50) & (face$coords[,1] < 25)
lips <- subset.face3d(face, select.lips)
source("/Volumes/adrian/research/shape/Face3D/Face3D-package/Face3D_0.1-1/Face3D/R/index-dist.r")
lips.si <- indexdist.face3d(lips, distance = 4)$shape.index
display.face3d(lips, colour = lips.si, type = "mesh", new = FALSE)
lmks.mouth <- c("chR", "chL", "ls", "cphL", "cphR", "li", "st")
face.lmks.mouth <- face$lmks[lmks.mouth, ]
spheres3d(face.lmks.mouth, radius = 0.7, col = "red")

source("/Volumes/adrian/research/shape/Face3D/Face3D-package/Face3D_0.1-1/Face3D/R/planepath.r")
source("/Volumes/adrian/research/shape/Face3D/Face3D-package/Face3D_0.1-1/Face3D/R/subset.r")
cl1  <- closest.face3d(face.lmks.mouth["chR", ], lips)
cl2  <- closest.face3d(face.lmks.mouth["st", ],  lips)
x1   <- cl1$point
x2   <- cl2$point
nrml <- (lips$normal[cl1$id, ] + lips$normal[cl2$id, ]) / 2
path <- planepath.face3d(lips, x1, x2, rotation = 0)
# spheres3d(path$path, radius = 0.3)
midline <- matrix(nrow = 0, ncol = 3)
for (i in 2:(nrow(path$path) - 1)) {
  drn   <- c(crossproduct(path$path[i+1, ] - path$path[i-1, ], nrml))
  opath <- planepath.face3d(lips, path$path[i, ], values = lips.si,
    direction = drn, boundary = c(10, 2), rotation = 0)
  spheres3d(opath$path, radius = 0.7)
  print(opath$values)
  plot(opath$values ~ opath$arclength, type = "l")
  ind     <- which.min(opath$values)
  midline <- c(midline, opath$path[ind, ])
}
spheres3d(midline, radius = 0.3)


load("~/Desktop/temp.dmp")
shp <- subset.face3d(lips, temp)
dim(shp$coords)

###############
# Select lips #

select.lips <- (face$coords[,2] < -60) & (face$coords[,2] > -95) & (face$coords[,1] > -10) & (face$coords[,1] < 50)
lips <- subset.face3d(face, select.lips)
display.face3d(lips, new=FALSE)

# Calculate shape index for the lips #

source("/Volumes/adrian/research/shape/Face3D/Face3D-package/Face3D_0.1-1/Face3D/R/index-dist.r")
lips.si <- indexdist.face3d(lips, distance = 4)$shape.index
display.face3d(lips, colour = lips.si)


ind <- (lips$shape.index < -3/8)
lips.trough <- subset.face3d(lips, ind)
display.face3d(lips.trough, colour = "shape index", extent = 10)


