#           Test asymmetry-path.face3d function

library(Face3D)
library(fields)
library(rgl)
library(geometry)

source("~/ownCloud/Face3D_0.1-1/Face3D/R/a-midline.R")
source("~/ownCloud/Face3D_0.1-1/Face3D/R/a-ridge2d.R")
source("~/ownCloud/Face3D_0.1-1/Face3D/R/a-display.R")
source("~/ownCloud/Face3D_0.1-1/Face3D/R/a-initiallandmarks.R")

x1 <- c(12, -7, 85)
x1 <- closest.face3d(x1, face)$point
x2 <- c(17.47, -75, 96.70)
x2 <- closest.face3d(x2, face)$point

face <- orient.face3d(face)

face$lmks <- NULL
rgl.close()
face <- initiallandmarks.face3d(face, orient = FALSE, monitor = TRUE, graphics = TRUE)

plot(template)
spheres3d(template$lmks, radius = 3)

plot(face)
spheres3d(face$lmks, col = "red")

ind.sbst  <- (abs(face$coords[ , 1] - x1[1]) < 40) & (face$coords[ , 2] < x1[2] + 40) & (face$coords[ , 2] > x2[2] - 40) 
sbst      <- subset.face3d(face, ind.sbst)
sbst      <- index.face3d(sbst, directions = TRUE, overwrite = TRUE)
plot(sbst, new = FALSE)
spheres3d(rbind(x1, x2), col = "red")

reference  <- planepath.face3d(sbst, x1, x2, si.target = 1, directions = TRUE)

# plot(sbst, col = "shape index", new = FALSE)
# plot(sbst, col = sbst$kappa2 + 1, new = FALSE)
# spheres3d(rbind(x1, x2), col = "red")
# breaks     <- unique(quantile(sbst$kappa2, seq(0, 1, by = 0.05), na.rm = TRUE))
# breaks[length(breaks)] <- breaks[length(breaks)] + 1
# clr        <- cut(reference$kappa2, breaks, labels = FALSE, include.lowest = TRUE)
# clr        <- topo.colors(length(breaks) - 1)[clr]
# spheres3d(reference$path, col = clr, radius = 0.5)
# d1 <- matrix(c(t(cbind(reference$path, reference$path + 0.7 * t(reference$directions[ , 1, ])))), ncol = 3, byrow = TRUE)
# d2 <- matrix(c(t(cbind(reference$path, reference$path + 0.7 * t(reference$directions[ , 2, ])))), ncol = 3, byrow = TRUE)
# segments3d(d1)
# segments3d(d2)
# plot(sbst, col = "shape index", new = FALSE)
# ind <- (reference$shape.index > 1/8) & (reference$shape.index < 5/8)
# ind[is.na(ind)] <- FALSE
# spheres3d(reference$path[ind, ], col = "yellow", radius = 0.5)
# d2 <- matrix(c(t(cbind(reference$path[ind, ], reference$path[ind, ] + 0.7 * t(reference$directions[ , 2, ind])))), ncol = 3, byrow = TRUE)
# segments3d(d2)
# drns <- reference$directions[ , 2, ind]
# dotp <- apply(drns, 2, function(x) sum(x * drns[ , 2]))
# ind1 <- (dotp > 0)
# drns[ , !ind1] <- -drns[ , !ind1]
# drns <- t(drns)
# d2 <- matrix(c(t(cbind(reference$path[ind, ], reference$path[ind, ] + 0.7 * drns))), ncol = 3, byrow = TRUE)
# segments3d(d2)
# direct <- apply(drns, 2, mean)
# direct <- direct / sqrt(sum(direct^2))
# dd <- matrix(c(t(cbind(reference$path[ind, ], reference$path[ind, ] + 0.7 * direct))), ncol = 3, byrow = TRUE)
# segments3d(d2, col = "red")
# for (i in 1:20) {
#    # print(i)
#    path.temp <- planepath.face3d(sbst, reference$path[i, ], direction =  sbst$directions[ , 2, i], rotation = 0, boundary = c(10, 10), bothways = FALSE)
# }

plot(sbst, new = FALSE, col = "shape index")
plot(sbst, new = FALSE, col = sbst$kappa2 + 1)
plot(sbst, new = FALSE, col = sbst$kappa1 + 1)
plot(sbst, new = FALSE, col = sbst$kappa1 * sbst$kappa2 + 1)
plot(sbst, new = FALSE, col = sbst$kappa1 + sbst$kappa2 + 1)
plot(sbst, add = TRUE, type = "direction 1")
plot(sbst, add = TRUE, type = "direction 2")

rgl.close()
mline  <- midline.face3d(sbst, reference = reference, lambda = 0.005, df = 24, d.asym = 25, nv = 31,
                         quantile = 0.25, extension = c(0, 0), monitor = TRUE)

# clr <- topo.colors(20)[cut(mline$z, 20, labels = FALSE)]
# plot(mline$x, mline$y, col = clr, xlab = "arc length", ylab = "adjustment")
zz <- mline$z
zz <- pmax(zz, quantile(zz, 0.25))
ridge <- ridge2d.face3d(mline$x, mline$y, zz, lambda = 0.1, df = 24, endpoints.fixed = FALSE, monitor = TRUE)
points(ridge$x, ridge$y, col = "red", pch = 16)

sbst1 <- subset(sbst, mline$subset)
plot(sbst1, col = "shape index")

plot(sbst1, new = FALSE, col = sbst1$kappa2 + 1)
plot(sbst1, add = TRUE, type = "direction 2")
plot(sbst1, new = FALSE, col = sbst1$kappa1 + 1)
plot(sbst1, add = TRUE, type = "direction 1")
plot(sbst1, new = FALSE, type = c("surface", "direction 2"), col = sbst1$kappa2 + 1)
spheres3d(mline$path, radius = 0.5)

plot(sbst, new = FALSE, type = c("surface", "directions 1"))

rgl.close()
mline  <- midline.face3d(sbst, reference = mline$path, lambda = 0.01, dh = 20, extension = c(0, 0),
                         monitor = TRUE, monitor.prompt = FALSE)
spheres3d(reference$path, radius = 0.5)
# mline  <- midline.face3d(face, reference = mline$path, lambda = 1, extension = c(30, 10), monitor = TRUE)
sbst   <- mline$subset
mline  <- mline$path

subst <- subset(face, sbst)
plot(subst, new = FALSE)
subst <- index.face3d(subst, overwrite = TRUE, directions = TRUE)
plot(subst, colour = "shape index", new = FALSE)
plot(subst, colour = 1 + subst$shape.index, new = FALSE)
plot(subst, colour = 1 + subst$kappa1, new = FALSE)
plot(subst, colour = 1 + subst$kappa2, new = FALSE)
plot(subst, type = "points", new = FALSE)
d1 <- matrix(c(t(cbind(subst$coords, subst$coords + 0.7 * t(subst$directions[ , 1, ])))), ncol = 3, byrow = TRUE)
d2 <- matrix(c(t(cbind(subst$coords, subst$coords + 0.7 * t(subst$directions[ , 2, ])))), ncol = 3, byrow = TRUE)
segments3d(d1)
segments3d(d2)

gc     <- gcurvature.face3d(mline, 15)
gc.sgn <- gc$gcurvature * sign(gc$d2.z)
prn    <- gc$resampled.curve[which.min(gc.sgn), ]
sn     <- gc$resampled.curve[which.max(gc.sgn), ]
gcrv   <- gc.sgn[gc$arc.length < gc$arc.length[which.min(gc.sgn)]]
se     <- gc$resampled.curve[which.max(gcrv), ]
plot(gc$arc.length, gc.sgn)

# face <- index.face3d(face, subset = sbst, directions = TRUE)
# plot(face, colour = "shape index", new = FALSE)
plot(face, new = FALSE)
spheres3d(rbind(x1, x2), col = "red")
spheres3d(mline, radius = 0.5, col = "green")
spheres3d(prn, radius = 2, col = "yellow")
spheres3d(sn, radius = 2, col = "yellow")
spheres3d(se, radius = 2, col = "yellow")

plot(face)
spheres3d(face$TPS, col = "green")
spheres3d(face$lmk, col = "red", radius = 1.5)
spheres3d(mline[1, ], col = "red", radius = 2.5)


plot(face, new = FALSE)
