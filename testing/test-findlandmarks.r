library(Face3D)

setwd("/Volumes/adrian/research/shape/Face3D/Face3D-package/Face3D_0.1-1/Face3D/R")
source("normals.r")
source("index.r")
source("orient.r")
source("landmarks.r")
source("display.r")

display.face3d(face)

rm(face)
face <- orient.face3d(face)
display.face3d(face, new = FALSE)
spheres3d(face$coords[face$nearest, ], radius = 2, col = "blue")
face <- findlandmarks.face3d(face, lmks = c("pn"),         orient = FALSE, monitor = TRUE)
face <- findlandmarks.face3d(face, lmks = c("enL", "enR"), orient = FALSE, monitor = TRUE)
face <- findlandmarks.face3d(face, lmks =   "se",          orient = FALSE, monitor = TRUE)
display.face3d(face, new = FALSE)
spheres3d(face$lmks, radius = 2, col = "blue")

# setwd("/Volumes/Face3D/Glasgow-controls/Liberty/Liberty-controls-dmp/Liberty-controls-dmp-final")
# files <- list.files()
# for (fname in files) {
   # cat(match(fname, files), "of", length(files), ":", fname, "\n")
   # load(fname)
   # face <- orient.face3d(face)
   # ind  <- apply(face$coords, 1, function(x) sqrt(sum((x - face$lmks["pn", ])^2))) < 100
   # centre <- subset.face3d(face, ind)
   # centre <- orient.face3d(centre)
   # display.face3d(centre, new = FALSE)
   # ind  <- apply(centre$coords, 1, function(x) sqrt(sum((x - centre$lmks["pn", ])^2))) < 90
   # centre <- index.face3d(centre, subset = ind)
   # display.face3d(centre, colour = "shape index", new = FALSE)
   # save(centre, file = paste("/Users/abowman/Desktop/faces/", fname, sep = ""))
# }

setwd("/Users/abowman/Desktop/faces/")
files <- list.files()
load(files[1])
display.face3d(centre)

display.face3d(centre, colour = "shape index", new = FALSE)

lmks <- c("pn", "enL", "enR")
lmks <- c("pn", "enL", "enR", "alL", "alR")
lmks <- c("pn", "acL", "acR")
lmks <- c("pn", "enL", "enR", "exL", "exR")

for (fname in files) {
   print(fname)
   load(fname)
   mtch        <- match(lmks, rownames(centre$lmks))
   mtch        <- mtch[!is.na(mtch)]
   centre$lmks <- centre$lmks[-mtch, ]
   # display.face3d(centre, new = FALSE)
   display.face3d(centre, colour = "shape index", new = FALSE)
   centre      <- landmarks.face3d(centre, lmks, orient = FALSE, graphics = TRUE)
   spheres3d(centre$lmks[lmks, ], radius = 3, col = "red")
   title3d(fname)
   # snapshot3d(paste("/Users/abowman/Desktop/images/",
                    # substr(fname, 1, nchar(fname) - 4), ".png", sep = ""))
   scan()
}

pop3d()


ind    <- (centre$shape.index < -0.8)
endo   <- subset.face3d(centre, ind)
display.face3d(endo, colour = "shape index", new = FALSE)
parts  <- connected.face3d(endo)
endo   <- subset.face3d(endo, parts %in% 1:2)
display.face3d(endo, colour = "shape index", new = FALSE)
# endo   <- index.face3d(endo, distance = 30)
# display.face3d(endo, colour = "shape index", new = FALSE)
crv    <- pmin(endo$kappa1, endo$kappa2)
clr    <- topo.colors(20)[cut(crv, 20, labels = FALSE)]
display.face3d(endo, colour = clr, new = FALSE)
spheres3d(endo$coords, col = clr, radius = 0.5)


#     Nose ridges

ind  <- apply(face$coords, 1, function(x) sqrt(sum((x - face$lmks["pn", ])^2))) < 50
face <- index.face3d(face, distance = 10, subset = ind)
display.face3d(face, new = FALSE)
display.face3d(face, colour = "shape index", new = FALSE)

nose    <- subset.face3d(face, ind, remove.singles = FALSE)
display.face3d(nose, colour = "shape index", new = FALSE)
ridges  <- subset.face3d(nose, nose$shape.index > 0)
valleys <- subset.face3d(nose, nose$shape.index < 0)
display.face3d(ridges,  colour = "shape index", new = FALSE)
display.face3d(valleys, colour = "shape index", add = TRUE)
crv  <- abs(pmin(ridges$kappa1, ridges$kappa2))
clr   <- topo.colors(20)[cut(crv, 20, labels = FALSE)]
spheres3d(ridges$coords, col = clr)
crv  <- abs(pmin(valleys$kappa1, valleys$kappa2))
clr   <- topo.colors(20)[cut(crv, 20, labels = FALSE)]
spheres3d(valleys$coords, col = clr)

parts <- connected.face3d(ridges)
ridges <- subset.face3d(ridges, parts == 1)
display.face3d(ridges, colour = "shape index", new = FALSE)



spheres3d(t(face$lmks["pn", ]), radius = 2, col = "red")


valleys <- subset.face3d(nose, nose$shape.index < 0, remove.singles = FALSE)
crv2    <- crv[crv > 0.2]
clr2    <- clr[crv > 0.2]
display.face3d(valleys)
spheres3d(valleys$coords, col = clr2)
parts <- connected.face3d(valleys)
valley <- subset.face3d(valleys, parts == 1)
crv1    <- crv2[parts == 1]
clr1    <- clr2[parts == 1]
display.face3d(valley, new = FALSE)
spheres3d(valley$coords, col = clr1)
spheres3d(t(face$lmks["acL", ]), radius = 2, col = "red")


# Try radial strips

ind   <- apply(face$coords, 1, function(x) sqrt(sum((x - face$lmks["pn", ])^2))) < 50
nose  <- subset.face3d(face, ind)
display.face3d(nose)
angle  <- pi/2
coords <- sweep(nose$coords, 2, face$lmks["pn", ])
coords <- rotate3d(coords, angle, 0, 0, 1)
coords <- sweep(coords, 2, face$lmks["pn", ], "+")

ind    <- (abs(coords[ , 1] - face$lmks["pn", 1]) < 50) &
          (abs(coords[ , 2] - face$lmks["pn", 2]) < 1)
spheres3d(nose$coords[ind, ], col = "green")


setwd("/Volumes/Face3D/Glasgow-controls/Liberty/Liberty-controls-dmp/Liberty-controls-dmp-final")
files <- list.files()
# ind <- (as.numeric(substr(files, 18, 20)) %in%
           # c(51, 123, 192, 199, 200, 202, 207, 210, 217, 218, 222, 228))
# files <- files[ind]
for (fname in files) {
   print(fname)
   load(fname)
   face$lmks <- face$lmks[-which(rownames(face$lmks) == "pn"), ]
   face <- landmarks.face3d(face)
   display.face3d(face, new = FALSE)
   # spheres3d(face$coords[face$nearest, ], col = "green", radius = 2)
   if (!any(is.na(face$lmks["pn", ])))
      spheres3d(t(face$lmks["pn", ]), col = "blue", radius = 2)
   title3d(fname)
   snapshot3d(paste("/Users/abowman/Desktop/images/",
                    substr(fname, 1, nchar(fname) - 4), ".png", sep = ""))
}




fname <- files[9]
load(fname)
display.face3d(face, new = FALSE)
face <- orient.face3d(face, graphics = TRUE)
display.face3d(face, new = FALSE)
spheres3d(face$coords[face$nearest, ], col = "green", radius = 2)

face <- landmarks.face3d(face)
spheres3d(t(face$lmks["pn", ]), col = "blue", radius = 2)

fname <- files[1]
load(fname)
face$coords  <- rotate3d(face$coords,  pi/8, 0, 1, 0)
face$normals <- rotate3d(face$normals, pi/8, 0, 1, 0)
display.face3d(face, new = FALSE)
display.face3d(face, new = FALSE, colour = "normal-x")
display.face3d(face, new = FALSE, colour = "normal-y")
display.face3d(face, new = FALSE, colour = "normal-z")

ind <- (face$coords[ , 3] > mean(face$coords[ , 3]))
face1 <- subset(face, ind)
nstrip <- 50
vstrip <- cut(face1$coords[ , 1], nstrip)
hstrip <- cut(face1$coords[ , 2], nstrip)
grid1  <- tapply(abs(face1$normals[ , 1]), list(vstrip, hstrip), mean, na.rm = TRUE)
grid2  <- tapply(abs(face1$normals[ , 2]), list(vstrip, hstrip), mean, na.rm = TRUE)
grid3  <- tapply(abs(face1$normals[ , 3]), list(vstrip, hstrip), mean, na.rm = TRUE)
image(grid1, col = topo.colors(20))
nw <- 15
fn <- function(i, grid) mean(colMeans(abs(grid[i + (-nw:nw), ] - grid1[i - (-nw:nw), ])), na.rm = TRUE)
fval1 <- apply(as.matrix((1 + nw):(nstrip - nw)), 1, fn, grid1)
fval2 <- apply(as.matrix((1 + nw):(nstrip - nw)), 1, fn, grid2)
fval3 <- apply(as.matrix((1 + nw):(nstrip - nw)), 1, fn, grid3)
plot((1 + nw):(nstrip - nw), fval1)
plot((1 + nw):(nstrip - nw), fval2)
plot((1 + nw):(nstrip - nw), fval3)
plot((1 + nw):(nstrip - nw), fval1 + fval2 + fval3)

rm(face)
fname <- files[23]
load(fname)
display.face3d(face, new = FALSE)
face <- normals.face3d(face)

face1 <- orient.face3d(face, graphics = TRUE)
display.face3d(face1, new = FALSE)
spheres3d(face1$coords[face1$nearest, ], col = "green", radius = 2)
spheres3d(face$coords[face$nearest, ], col = "green", radius = 2)


# Old ideas

library(geometry)

x <- face$coords[ , 1] - 20
y <- face$coords[ , 2] + 45
ind <- abs(x) < 50 & abs(y) < 75
spheres3d(face$coords[ind, ], col = "blue")

centre <- subset.face3d(face, ind)
centre <- normals.face3d(centre)

display.face3d(centre, new = FALSE)
centre <- index.face3d(centre, subset = 2000)
ind    <- !is.na(centre$shape.index)
clr    <- topo.colors(20)[cut(centre$shape.index[ind], 20, labels = FALSE)]
spheres3d(centre$coords[ind, ], col = clr)

crv <- pmin(-centre$kappa1, -centre$kappa2)
crv[centre$shape.index < 0.8] <- 0
display.face3d(centre, colour = centre$shape.index, new = FALSE)
display.face3d(centre, colour = 100 * crv, new = FALSE)

crv  <- pmin(-centre$kappa1, -centre$kappa2)
crv1 <- interpolate.face3d(centre, crv[!is.na(crv)], ind)
crv[!ind] <- crv1
si   <- centre$shape.index
si1  <- interpolate.face3d(centre, si[!is.na(si)], ind)
si[!ind] <- si1
display.face3d(centre, colour = si, new = FALSE)
crv[si < 0.8] <- 0
display.face3d(centre, colour = 100 * crv, new = FALSE)

# nose tip

crv <- pmin(-centre$kappa1, -centre$kappa2)
crv[centre$shape.index < 0.8] <- 0
nose.tip <- centre$coords[which.max(crv), ]
display.face3d(centre, colour = 100 * crv, new = FALSE)
spheres3d(nose.tip, col = "blue")
display.face3d(centre, new = FALSE)
spheres3d(nose.tip, col = "blue")


# eye pit

crv <- pmin(centre$kappa1, centre$kappa2)
crv[centre$shape.index > -0.6] <- 0
eyes <- subset.face3d(centre, centre$shape.index < -0.6)
display.face3d(eyes)
parts <- connected.face3d(eyes)
eyes <- subset.face3d(eyes, parts <= 2)
crv <- pmin(eyes$kappa1, eyes$kappa2)
display.face3d(eyes, colour = 100 * crv)
parts <- parts[parts <= 2]
eye1 <- subset.face3d(eyes, parts == 1)
eye2 <- subset.face3d(eyes, parts == 2)
display.face3d(centre, new = FALSE)
spheres3d(eye1$coords[which.max(crv[parts == 1]), ], col = "blue")
spheres3d(eye2$coords[which.max(crv[parts == 2]), ], col = "blue")

# alare

display.face3d(centre, colour = "shape index", new = FALSE)
ind <- (centre$shape.index > -0.8) & (centre$shape.index < -0.2)
valleys <- subset.face3d(centre, ind)
display.face3d(valleys, colour = "shape index")




ch <- convhulln(centre$coords)
triangles3d(centre$coords[c(t(ch)), ], col = "green", alpha = 0.7)
pop3d()
spheres3d(centre$coords[unique(c(ch)), ], col = "blue")

x   <- centre$coords[ , 1]
y   <- centre$coords[ , 2]
ind <- abs(x) < 50 & abs(y) < 50
spheres3d(centre$coords[ind, ], col = "blue")

i <- match("shape.index", names(centre))
centre <- centre[-i]
eye <- subset.face3d(centre, ind)
display.face3d(eye)
eye$shape.index <- index.face3d(eye, distance = 20)$shape.index
display.face3d(eye, colour = eye$shape.index, new = FALSE)

