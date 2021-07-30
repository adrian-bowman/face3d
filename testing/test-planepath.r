#     Planepath algorithm

install.packages("~/research/face3d/face3d", repos = NULL, type = "source")
library(face3d)
library(rgl)

# Interpolate the principal directions and normals
dst   <- rdist(t(template_male$landmarks["pn", ]), template_male$vertices)
nose  <- subset(template_male, dst < 30)
ppath <- planepath(nose, nose$landmarks["pn", ], nose$landmarks["acL", ],
                          si.target = 1, directions = TRUE)
plot(nose)
spheres3d(ppath$path, radius = 0.1)
dm <- 3
segs <- rbind(t(ppath$path) + 5 * ppath$directions[ , dm, ],
              t(ppath$path) - 5 * ppath$directions[ , dm, ])
segments3d(matrix(c(segs), ncol = 3, byrow = TRUE), col = "blue")

# Dealing with points on the edges

load("nose.Rda")

i <- 1494
j <- 1491

plot(nose)
spheres3d(nose$coords[c(i, j), ], radius = 0.5)

AUX <- planepath(nose,nose$coords[i,], nose$coords[j,])
if (is.null(AUX$path)) {
   edges <- edges.face3d(nose)
   ind   <- which(unlist(lapply(edges, function(x) all(c(i, j) %in% x))))
   if (length(ind) > 0) {
      edge <- edges[[ind]]
      path <- nose$coords[edge, ]
      arc  <- c(apply(diff(path)^2, 1, function(x) sqrt(sum(x))), 0)
      ind  <- which(edge %in% c(i, j))
      ind  <- sort(ind)
      d1   <- sum(arc[ind[1]:(ind[2] - 1)])
      d2   <- sum(arc[1:(ind[1] - 1)]) + sum(arc[ind[2]:arc[length(arc)]])
      dst  <- min(d1, d2)
   }
}

# -----------------------------

setwd("/Volumes/Face3D/Glasgow-controls/Kathryn/Kat-Controls-dmp/Kathryn-Controls-dmp-final")
load("kathryn004.dmp")
class(face) <- "face3d"

# Identify the triangles which lie along the planepath between x1 and x2
ind  <- (edist(face$coords, face$lmks["pn", ]) < 100)
face <- subset(face, ind)
plot(face, new = FALSE)
pp <- planepath(face, face$lmks["pn", ], face$lmks["se", ], rotation = 0)
triangles3d(face$coords[c(t(pp$triangles)), ], front= "line", back = "line")
spheres3d(face$lmks, col = "red", radius = 2)

lmks <- face$lndms
rownames(lmks) <- c("pn","acL","acR","sn","se","n","exL","exR","enL","enR","tL","tR","ls","cphL","cphR","chL","chR","st","li","sl","gn","oiL","oiR")
face$lmks <- lmks
face <- curves.face3d(face)
n.vec <- c(10,10,10,10,15,15,5,20,10,10,5,5,5,10,10,10,20,20,20,20,15,15,10,10)
rcrvs <- resamplecurves.face3d(face, n.vec = n.vec)
points3d(rcrvs)


x1 <- closest.face3d(face$lndms[11,], face)$point
x2 <- closest.face3d(face$lndms[ 2,], face)$point
x1 <- face$lndms[ 9,]
x2 <- face$lndms[ 5,]
plot(face, type = "mesh", new = FALSE)
spheres3d(rbind(x1, x2), radius = 1, col = "blue")
p <- planepath(face, x1, x2)
lines3d(p$path, col = "green", lwd = 3)

triangles <- matrix(face$triples, ncol = 3, byrow = TRUE)
triangles3d(face$coords[triangles[1516, ], ], col = "green")


planepath(face, x1 = face$lndms[9,], x2 = face$lndms[5,])$path

setwd("/Volumes/Face3D/reliability-project/reliability-pts/glasgow-pts/")
lmks <- read.face3d("SKPLY01.pts")$lmks
setwd("/Volumes/Face3D/reliability-project/reliability-images/images-glasgow")
face <- read.face3d("001.obj", jpgfile.addition = "")
   
x1 <- closest.face3d(lmks[ 2,], face)$point
x2 <- closest.face3d(lmks[ 5,], face)$point
display.face3d(face, type = "mesh")
spheres3d(rbind(x1, x2), radius = 1, col = "blue")
source("Face3D/R/planepath.r")
p <- planepath(face, x1, x2)
lines3d(p$path, col = "green")


files <- list.files("/Volumes/Face3D/Glasgow-controls/Kathryn/kathryn-dmp", pattern = "kathryn*", full.names = TRUE)
load(files[1])
display.face3d(face)
i <- 12
for (i in 1:length(files)) {
   print(i)
   load(files[i])
   x1 <- closest.face3d(face$lndms[ 2,], face)$point
   x2 <- closest.face3d(face$lndms[16,], face)$point
   display.face3d(face, type = "mesh", new = FALSE)
   spheres3d(rbind(x1, x2), radius = 0.5, col = "blue")
   source("Face3D/R/planepath.r")
   p <- planepath.face3d(face, x1, x2)
   lines3d(p$path, col = "green")
}

data(face)
nose <- subset(face, face$coords[ , 3] > 85 & 
           face$coords[ , 2] > -60 & face$coords[ , 2] < -15)
plot(nose, type = "mesh")

plot(nose, type = "mesh", new = FALSE)
x1 <- nose$coords[2737, ]
x1 <- (nose$coords[2736, ] + nose$coords[2737, ]) / 2
x1 <- (nose$coords[2736, ] + nose$coords[2737, ] +
       nose$coords[2704, ]) / 3
x2 <- nose$coords[146, ]
spheres3d(t(matrix(x1)), col = "red", radius = 0.2)
spheres3d(t(matrix(x2)), col = "red", radius = 0.2)

source("Face3D/R/planepath.r")
p <- planepath(nose, x1, x2, rotation = "optimise")
lines3d(p$path, col = "green")
p$length

nose0 <- subset(nose, 
              !((abs(nose$coords[ , 1] - 20) < 3) &
                (abs(nose$coords[ , 2] + 35) < 2)))
plot(nose0, type = "mesh", new = FALSE)
spheres3d(t(matrix(x1)), col = "red", radius = 0.2)
spheres3d(t(matrix(x2)), col = "red", radius = 0.2)

source("Face3D/R/planepath.r")
p <- planepath(nose0, x1, x2, bridge.gaps = TRUE)
lines3d(p$path, col = "green")



for (i in 1:nrow(p$path)) {
   spheres3d(matrix(p$path[i,], nrow = 1), radius = 0.3, col = "red")
   scan()
}

plot(face, type = "mesh", new = FALSE)
spheres3d(t(matrix(x1)), col = "red", radius = 0.2)
spheres3d(t(matrix(x2)), col = "red", radius = 0.2)

pop3d()
source("Face3D/R/planepath.r")
p <- planepath(face, x1, x2, rotation = "optimise")
#             boundary = c(1, 1))
#             boundary = c(100, 100))
lines3d(p$path, col = "green")
p$length

up <- c(crossproduct(rcross, x2 - x1))
quads3d(rbind(x1, x1 + 5*up, x2 + 5*up, x2), col = "red", alpha = 0.5)

triangles <- matrix(nose$triples, ncol = 3, byrow = TRUE)
ind <- which(apply(triangles, 1, function(x) all(c(2507, 2508) %in% x)))
triangles3d(nose$coords[c(t(triangles[4720, ])), ], col = "green")
rcross <- c(0.5579822, -0.5587915, -0.6135210)
a <- triangles[4720, ]
c(nose$coords[a, ] %*% rcross)
sum(x1 * rcross)


pop3d()
source("Face3D/R/planepath.r")
p <- planepath(nose, x1, direction = c(1, 1, 0))
lines3d(p$path, col = "green")

surface <- subset.face3d(face, 
   face$coords[ , 3] > 35 & face$coords[ , 2] < 50)
display.face3d(surface, type = "mesh")

source("Face3D/R/planepath.r")
crossings <- planepath(surface, x1, x2)
lines3d(crossings, col = "green")

# Try this instead of indexing to find the neighbours and shape index
# Could also be used to reduce the search area for find.triangle.
ed <- function(x, coords) sqrt(apply(sweep(coords, 2, x)^2, 1, sum))
sort(ed(x1, face$coords))
which.min(ed(x1, face$coords))

# x1 only

setwd("/Volumes/Face3D/Glasgow-controls/Liberty/Liberty-Controls-dmp/Liberty-Controls-dmp-final/")

library(Face3D)
library(rgl)

load("controls-liberty-041.dmp")
summary(face)
x1 <- face$lmks["pn", ]
shape <- subset(face, edist.face3d(face$coords, x1) < 20)
shape <- curvatures(shape, directions = TRUE, overwrite = TRUE)
plot(shape, colour = "shape index")
spheres3d(x1)

pop3d()

pp      <- planepath(shape, x1, direction = c(0, 1, 0), boundary = 1, rotation = 0, directions = TRUE)
si      <- 2 / pi * (atan((pp$kappa2 + pp$kappa1) / (pp$kappa2 - pp$kappa1)))
clr.r   <- c(rep(0, 3), 0.5, rep(1, 5))
clr.g   <- c(rep(1, 7), 0.5, 0)
clr.b   <- c(0, 0.5, rep(1, 3), 0.5, rep(0, 3))
clr.pts <- seq(-1, 1, by = 0.25)
clr     <- rep("white", length(shape$shape.index))
clr     <- rgb(approx(clr.pts, clr.r, xout = si)$y,
               approx(clr.pts, clr.g, xout = si)$y,
               approx(clr.pts, clr.b, xout = si)$y)

spheres3d(pp$path, radius = 0.5, col = clr)
