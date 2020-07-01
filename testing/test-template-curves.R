#              Template curve fitting

setwd("~/ownCloud/Face3D_0.1-1")

library(Face3D)
library(fields)
library(rgl)
library(shapes)
library(geometry)
library(colorspace)

source("Face3D/R/a-midline.R")
source("Face3D/R/a-ridge2d.R")
source("Face3D/R/a-display.R")
source("Face3D/R/a-initiallandmarks.R")
source("Face3D/R/a-planepath.R")

liberty.data <- "~/Desktop/Data-Liberty/"
# liberty.data <- "/Volumes/Face3D/Glasgow-controls/Liberty/Liberty-Controls-dmp/Liberty-Controls-dmp-final/Data-Liberty/"
fls          <- list.files(liberty.data, full.names = TRUE)

# Inspect curves on template
template <- template.male
ind <- grep("nasal", rownames(template$curves), value = TRUE)
ind <- c(ind, grep("columella", rownames(template$curves), value = TRUE))
ind <- grep("boundary", ind, value = TRUE, invert = TRUE)
tcurves <- template$curves[ind, ]
plot(template, new = FALSE)
spheres3d(tcurves, col = "red")
curves.ind <- ind

warp.curves <- function(lmks, template, curves.ind) {
   ind   <- match(rownames(lmks), rownames(template$lmks))
   mat1  <- rbind(template$lmks[ind, ], template$lmks, template$curves[curves.ind, ])
   wp    <- warp.face3d(mat1, lmks, subset = (nrow(lmks) + 1):nrow(mat1),
                        general = TRUE, project = FALSE)
   wp[nrow(template$lmks) + 1:nrow(template$curves[curves.ind, ]), ]
}

pcmodel <- function(scores, lmk.names, gpa, face, template) {
   lmks <- gpa$mshape + sweep(matrix(gpa$pcar[ , 1], ncol = 3), 1, scores * gpa$pcasd[1], "*") 
   rownames(lmks) <- lmk.names
   warp.template(lmks[lmk.names, ], template)
}

i <- 85
load(fls[i])
face$lmks <- NULL
face <- initiallandmarks.face3d(face, orient = TRUE, monitor = FALSE, graphics = FALSE, overwrite = TRUE)
face <- subset(face, edist.face3d(face$coords, face$lmks["pn", ]) < 100)
plot(face, type = "mesh", col = "shape index", new = FALSE)
spheres3d(face$lmks, col = "yellow")

wp <- warp.curves(face$lmks[c("pn", "se", "enL", "enR"), ], template, curves.ind)
spheres3d(wp, col = "red")

model <- warp.template(face$lmks[c("pn", "se", "enL", "enR"), ], template)

i <- 85
load(fls[i])
centroid    <- apply(face$coords, 2, mean)
face$curves <- sweep(face$curves, 2, centroid)
face$curves <- rgl::rotate3d(face$curves, face$rotations[1], 1, 0, 0)
face$curves <- rgl::rotate3d(face$curves, face$rotations[2], 0, 1, 0)
face$curves <- rgl::rotate3d(face$curves, face$rotations[3], 0, 0, 1)
face$curves <- sweep(face$curves, 2, centroid, "+")
plot(face, new = FALSE)
spheres3d(face$curves[curves.ind, ], col = "red")

gpa   <- procGPA(lmks.liberty[lmk.names, , ])

model <- subset(model, edist.face3d(model$coords, model$lmks["pn", ]) < 100)
plot(model, type = "mesh", colour = clr[2], add= TRUE)
pcm <- pcmodel(c(-3, rep(0, 12)), lmk.names, gpa, face, model)
pcm <- warp.template(face$lmks[c("pn", "se", "enL", "enR"), ], pcm)
pop3d()
plot(pcm, type = "mesh", colour = clr[2], add= TRUE)


open3d()
spheres3d(gpa$rotated)
spheres3d(gpa$mshape)
shp <- gpa$mshape - 3 * gpa$pcasd[1] * matrix(gpa$pcar[ , 1], ncol = 3)

pop3d()
spheres3d(shp, col = "red", radius = 1.1)

