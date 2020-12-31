setwd("~/ownCloud/Face3D_0.1-1/")

library(Face3D)
library(rgl)
library(rpanel)
library(fields)

load("liberty.dmp")
x1  <- face$lmks["chR", ]
x2  <- face$lmks["cphR", ]
drn <- x2 - x1

path  <- planepath.face3d(face, x1, x2)
dst   <- closestcurve.face3d(face, path$path)$closest.distance
ind   <- abs(dst) < 10
dst   <- dst[ind]
shape <- subset(face, ind)
shape$lmks <- face$lmks[c("chR", "cphR"), ]
shape <- index.face3d(shape, directions = TRUE)

display.face3d(shape)
segs <- matrix(c(t(cbind(shape$coords, shape$coords + shape$normals))), ncol = 3, byrow = TRUE)
segments3d(segs)

editlandmarks.face3d(shape, lmk.name = "chR")
editlandmarks.face3d(shape, panel = FALSE, lmk.name = "chR", directions = TRUE)


display.face3d(shape)
spheres3d(rbind(x1, x2), col = "blue", radius = 0.5)
i <- 420
x <- shape$coords[i, ]
spheres3d(t(x), col = "green", radius = 0.5)

shape <- index.face3d(shape, subset = (dst > -4), distance = 10, overwrite = TRUE,
                  directions = TRUE)

segments3d(rbind(x, x + 5 * shape$axes[, 1, i]))
segments3d(rbind(x, x + 5 * shape$axes[, 2, i]), col = "green")
segments3d(rbind(x, x + 5 * shape$axes[, 3, i]), col = "blue")

segments3d(rbind(x, x + 5 * shape$directions[ , 1, i]), col = "blue")
segments3d(rbind(x, x + 5 * shape$directions[ , 2, i]), col = "red")

xgrid <- seq(-shape$si.distance, shape$si.distance, length = 20)
dgrid <- xgrid %o% shape$directions[ , 1, i]
zgrid <- sweep(dgrid + shape$kappa1[i] * xgrid^2 %o% shape$axes[, 1, i], 2, x, "+")
lines3d(zgrid, col = "blue")
dgrid <- xgrid %o% shape$directions[ , 2, i]
zgrid <- sweep(dgrid + shape$kappa2[i] * xgrid^2 %o% shape$axes[, 1, i], 2, x, "+")
lines3d(zgrid, col = "red")

shape <- index.face3d(shape, subset = (dst > -4), distance = 10, overwrite = TRUE)
shape$shape.index[dst <= -4] <- 0
display.face3d(shape, colour = "shape index", new = FALSE)
spheres3d(rbind(x1, x2), col = "blue", radius = 0.5)

path     <- planepath.face3d(shape, x1, x2, rotation = 0)
normal   <- path$normal
ind      <- (sign(shape$shape.index) == 1)
k1       <- shape$kappa1
k2       <- shape$kappa2
k1[!ind] <- 0
k2[!ind] <- 0
clr      <- pmax(abs(k1), abs(k2))
clr      <- pmax(-k2, 0)
clr      <- pmax(-k1, 0)
clr[is.na(clr)] <- 0
display.face3d(shape, colour = clr + 2, new = FALSE)
spheres3d(rbind(x1, x2), col = "red")

face <- index.face3d(face, subset = rdist(face$coords, t(face$lmks["st", ])) < 25, distance = 10, overwrite = TRUE)
sbst <- subset.face3d(face, rdist(face$coords, t(face$lmks["st", ])) < 25)
display.face3d(sbst, colour = -sbst$kappa2 + 2, new = FALSE)


spheres3d(t(x1 + 25 * path$normal), alpha = 0)
rot.grid <- seq(-0.65*pi/2, 0.65*pi/2, length = 100)
for (i in 1:length(rot.grid)) {
  rotation <- rot.grid[i]
  path     <- planepath.face3d(shape, x1, x2, rotation = rotation)
  par3d(skipRedraw = TRUE)
  if (i > 1) for (j in 1:2) pop3d()
  lines3d(path$path, lwd = 4)
  normal <- rotate3d(path$normal, rotation, drn[1], drn[2], drn[3])
  quads3d(rbind(x1, x1 + 10 * normal, x2 + 10 * normal, x2), alpha = 0.5, col = "red")
  par3d(skipRedraw = FALSE)
}

path <- planepath.face3d(shape, x1, x2, si.target = 0.5)
display.face3d(shape, new = FALSE, colour = clr + 2)
spheres3d(rbind(x1, x2), col = "red")
lines3d(path$path, lwd = 4)
snapshot3d("lip4.png")

display.face3d(shape, new = FALSE)
spheres3d(rbind(x1, x2), col = "red")
lines3d(path$path, lwd = 4)
snapshot3d("lip5.png")

g <- seq(-1, 1, length = 50)
grid <- as.matrix(expand.grid(g, g))
kk1  <- grid[ , 1]
kk2  <- grid[ , 2]
si <- matrix(atan((kk2 + kk1) / (kk2 - kk1)) * 2 / pi, ncol = 50)
filled.contour(g, g, si)



x1        <- face$coords[which.max(face$coords[ , 3]), ]
x2        <- closest.face3d(x1 + c(0, 20, 0), face)$point
nose      <- subset(face, c(rdist(face$coords, t(x1))) < 20)
nose$lmks <- matrix(c(x1, x2), nrow = 2, byrow = TRUE,
                    dimnames = list(c("lmk1", "lmk2"), c("x", "y", "z")))

nose <- index.face3d(nose, subset = abs(nose$coords[ , 1] - x1[1]) <  5 &
                                    nose$coords[ , 2] - x1[2] > -5)




setwd("Face3D/R")

load("face.Rda")
library(Face3D)

library(rpanel)
library(rgl)
library(fields)
source("display.r")
source("subset.r")
source("crossproduct.r")
source("closest.r")
source("planepath.r")

x1        <- face$coords[which.max(face$coords[ , 3]), ]
x2        <- closest.face3d(x1 - c(2, 0, 0), face)$point
nose      <- subset(face, c(rdist(face$coords, t(x1))) < 20)
nose$lmks <- matrix(c(x1, x2), nrow = 2, byrow = TRUE,
                    dimnames = list(c("lmk1", "lmk2"), c("x", "y", "z")))

source("edit.r")
nose.new <- editlandmarks.face3d(nose, lmk.names = c(rownames(nose$lmks), "lmk3"))

view3d(0, 0)
origin <- apply(summary.face3d(nose)$ranges, 2, mean)
edist  <- function(x, x0) sqrt(apply(sweep(x, 2, x0)^2, 1, sum))
ind    <- which(edist(nose$coords[ , 1:2], origin[1:2]) < 1)
ind1   <- which.max(nose$coords[ind, 3])
pt     <- nose$coords[ind[ind1], ]
# pt     <- rgl.user2window(0.5, 0.5, 0.5)
spheres3d(t(pt), col = "blue", radius = 0.5)

display.face3d(nose)
view3d(0, 0)
umat <- par3d("userMatrix")
crds <- nose$coords %*% solve(umat[1:3, 1:3])
open3d()
view3d(0, 0)
points3d(crds)


display.face3d(nose)
spheres3d(nose$lmks, col = "green", radius = 0.5)
path    <- planepath.face3d(nose, x1,
                            direction = c(0, 1, 0),
                            rotation = 0, graphics = FALSE, boundary = c(Inf, Inf))
lines3d(path$path, col = "green", lwd = 2)

source("edit.r")
editlandmarks.face3d(face, lmk.names = c(rownames(shape.lmks), "lmk3"))
