#     Computing the curvatures and shape index

setwd("~/research/face3d")

library(face3d)

load("~/research/face3d/testing/test-data/nose.RData")
shape <- subset(template_male, c(rdist(t(shape$landmarks["pn", ]), shape$vertices)) < 30)
plot(shape)

shape <- normals.face3d(shape)
shape <- index.face3d(shape, distance = 5, overwrite = TRUE)
shape <- index.face3d(shape, distance = 5)
cbind(shape$kappa1, shape$kappa2)

plot(shape$kappa1, shape$kappa2)

plot(shape, colour = "shape index")
plot(shape, colour = shape$kappa1)
plot(shape, colour = shape$kappa2)

ind <- which.max(shape$kappa1)
spheres3d(shape$vertices[ind, ])

normals.face3d(shape)$axes[[1]]

plot(shape, "normal")







# Old material

setwd("/Volumes/adrian/research/shape/Face3D/Face3D-package/Face3D_0.1-1")

library(Face3D)

data(face)

display.face3d(face, col = "grey")
view3d()

open3d()
points3d(face$coords[sample(1:nrow(face$coords), 5000), ])

face <- index.face3d(face, extent = 1)
display.face3d(face, colour = "shape index", extent = 1)

nose <- subset.face3d(face, face$coords[ , 3] > 85 & 
           face$coords[ , 2] > -60 & face$coords[ , 2] < -15)
nose <- index.face3d(nose)
           
source("Face3D/R/summary.r")
summary(nose)
source("Face3D/R/index.r")
nose  <- index.face3d(nose, extent = 2)
for (i in 1:10) nose <- index.face3d(nose, extent = i)

source("Face3D/R/display.r")
display.face3d(nose, colour = "shape index", extent = 1)
for (i in 2:10) {
   display.face3d(nose, colour = "shape index", extent = i, new = FALSE)
   scan()
}
display.face3d(nose, colour = "shape index", extent = 2, new = FALSE)
display.face3d(nose, colour = "shape index", extent = 3, new = FALSE)
display.face3d(nose, colour = "shape index", extent = 4, new = FALSE)
display.face3d(nose, colour = "shape index", extent = 5, new = FALSE)
display.face3d(nose, colour = "shape index", extent = 6, new = FALSE)

surface <- index.face3d(surface, extent = 5)

shape <- surface
indx  <- unlist(surface$shape.index[,1])
indy  <- indx
indy  <- sm.shape

triples  <- matrix(shape$triples, ncol = 3, byrow = TRUE)
nbrs     <- lapply(shape$trngs1, function(x) unique(c(triples[x, ])))
av.fn    <- function(x) mean(indy[x], na.rm = TRUE)
sm.shape <- sapply(nbrs, av.fn)
cbind(indx, sm.shape)

SI.levels <- c(-1,-7/8,-5/8,-3/8,-1/8,1/8,3/8,5/8,7/8,1)
clr    <- c("green", "cyan", "blue", "turquoise", "white", "khaki1",
            "yellow", "orange", "red")
si.clr <- as.character(cut(unlist(surface$shape.index[,"extent-5"]),
           breaks = SI.levels, labels = clr, include.lowest = TRUE))
display.face3d(surface, colour = si.clr)
view3d()


nose <- index.face3d(nose, extent = 3)

display.face3d(nose, type = "mesh")
view3d()
spheres3d(nose$coords[1454, ], radius = 0.5, col = "green")

source("Face3D/R/read.r")
template <- read.face3d("/Volumes/Face3D/RCSItemplate.obj")
display.face3d(template)
view3d()

source("Face3D/R/index.r")
template <- index.face3d(template, extent = 2)
surface  <- index.face3d(surface, extent = 2)
range(index)
SI.levels <- c(-1,-7/8,-5/8,-3/8,-1/8,1/8,3/8,5/8,7/8,1)
clr    <- c("green", "cyan", "blue", "turquoise", "white", "khaki1",
            "yellow", "orange", "red")
si.clr <- as.character(cut(index, breaks = SI.levels,
                 labels = clr, include.lowest = TRUE))
display.face3d(template, colour = si.clr)

display.face3d(template)
view3d()

source("tneighbours.r")
tnbrs2 <- tneighbours(template, extent = 2)
coords <- template$coords
rotations <- array(unlist(tnbrs2$rotations),
              dim = c(3, 3, nrow(coords)))
normals <- t(rotations[ , 1, ])
segs <- matrix(c(t(cbind(coords, coords + 5 * normals))),
               ncol = 3, byrow = TRUE)
segments3d(segs)

ind    <- (1:nrow(template$coords) - 1) * 9
nrmls  <- unlist(tnbrs2$rotations)
nrmls  <- cbind(nrmls[ind], nrmls[ind + 1], nrmls[ind + 2])
segments3d(matrix(c(t(cbind(coords, coords + 5 * nrmls))),
                  ncol = 3, byrow = TRUE))


source("Face3D/R/subset.r")
surface <- subset.face3d(face, 
   face$coords[ , 3] > 35 & face$coords[ , 2] < 50)
nose <- subset.face3d(face,
           face$coords[ , 3] > 85 & face$coords[ , 2] > -60 &
           face$coords[ , 2] < -15)
display.face3d(nose)
view3d()
for (i in 5:12) surface <- index.face3d(surface, extent = i)
for (i in 1:12) {
   display.face3d(surface, col = "shape index", extent = i, new = FALSE)
   scan()
}

# Can we identify where the major features are just from the shape index?


surface <- subset.face3d(face, 
   face$coords[ , 3] > 35 & face$coords[ , 2] < 50)
surface <- index.face3d(surface, extent = 10)
summary.face3d(surface)

source("Face3D/R/summary.r")
summary(features)
source("Face3D/R/display.r")
apply(surface$shape.index, 2, range)
for (i in 1:12) {
   nm       <- paste("extent", i, sep = "=")
   ind      <- abs(surface$shape.index[ , nm]) > 0.7
   ind      <- surface$shape.index[ , nm] < -0.2
   features <- subset.face3d(surface, ind)
   display.face3d(features, colour = "shape index", extent = i,
          new = FALSE)
   scan()
}
mn <- apply(features$coords, 2, mean)
spheres3d(matrix(mn, nrow = 1), col = "red", radius = 2)
axes <- princomp(features$coords)$loadings
segments3d(rbind(mn, mn + 10 * axes[ , 1], mn, mn + 10 * axes[ , 2],
                 mn, mn + 10 * axes[ , 3]))

source("Face3D/R/subset.r")
source("Face3D/R/index.r")
source("Face3D/R/display.r")

surface <- subset.face3d(face, 
              face$coords[ , 3] > 35 & face$coords[ , 2] < 50)
surface <- index.face3d(surface, extent = 2)
surface <- index.face3d(surface, extent = 10)

display.face3d(surface, new = FALSE)
display.face3d(surface, col = "shape index", extent = 2, new = FALSE)
display.face3d(surface, col = "shape index", extent = 10, new = FALSE)

display.face3d(surface, col = "grey", alpha = 0.2, new = FALSE)
si       <- surface$shape.index[ , "extent=10"]
ind      <- (si < -3/8)
ind      <- (si < -5/8)
ind      <- (si < -3/8 & si > -5/8)
ind      <- (si < -5/8 & si > -7/8)
ind      <- (si < -7/8)
ind      <- (si < quantile(si, prob = 0.1))
pop3d()
features <- subset.face3d(surface, ind)
display.face3d(features, col = "shape index", extent = 10, add = TRUE)

triples <- matrix(features$triples, ncol = 3, byrow = TRUE)
npts    <- nrow(features$coords)
dst     <- matrix(0, ncol = npts, nrow = npts)
dst[rbind(triples[ , 1:2], triples[ , 2:3], triples[ , c(1, 3)])] <- 1
dst     <- dst + t(dst)
dst[dst > 0] <- 1
dst     <- as.dist(1 - dst)
clusters <- cutree(hclust(dst, "single"), h = 0.5)
ind <- as.numeric(names(sort(table(clusters), decreasing = TRUE))[1:4])
ind <- which(clusters %in% ind)
clusters <- clusters[ind]
features <- subset.face3d(features, ind)
display.face3d(features, col = "shape index", extent = 10, new = FALSE, add = TRUE)

si   <- features$shape.index[ , "extent=10"]
ind  <- tapply(si, clusters, function(x) which.min(x))
lbls <- unique(clusters)
lmks <- matrix(nrow = 0, ncol = 3)
for (i in 1:5) lmks <- rbind(lmks,
                 features$coords[clusters == lbls[i], ][ind[i], ])
spheres3d(lmks, col = "red", radius = 1)

display.face3d(surface)
spheres3d(lmks, col = "red", radius = 1)


library(fastcluster)
hclust(dist(features$coords)

source("tneighbours.r")
tnbrs1 <- tneighbours(surface)
source("tneighbours.r")
tnbrs2 <- tneighbours(surface, extent = 2)
source("tneighbours.r")
tnbrs3 <- tneighbours(surface, extent = 3)
tnbrs5 <- tneighbours(surface, extent = 5)
tnbrs6 <- tneighbours(surface, extent = 6)
index  <- tnbrs2$index
index  <- tnbrs3$index
index  <- tnbrs5$index
index  <- tnbrs6$index
index  <- tnbrs7$index
index  <- tnbrs8$index
index  <- tnbrs9$index
index  <- tnbrs10$index
range(index)
hist(unlist(lapply(tnbrs6$pts, length)))

i <- 1000
spheres3d(surface$coords[tnbrs10$pts[[i]], ], radius = 1, col = "green")
pop3d()


SI.levels <- c(-1,-7/8,-5/8,-3/8,-1/8,1/8,3/8,5/8,7/8,1)
clr    <- c("green", "cyan", "blue", "turquoise", "white", "khaki1",
            "yellow", "orange", "red")
si.clr <- as.character(cut(index, breaks = SI.levels,
                 labels = clr, include.lowest = TRUE))
display.face3d(surface, colour = si.clr)
view3d()

display.face3d(surface, colour = "grey")


nrmls  <- tnbrs2$normals
coords <- surface$coords[as.numeric(rownames(nrmls)), ]
segments3d(matrix(c(t(cbind(coords, coords + 5 * nrmls))),
           ncol = 3, byrow = TRUE))

index <- lapply(tnbrs2$pts, function(x) 
           lsfit(surface$coords[x, 1:2], surface$coords[x, 3]))


tnbrs3  <- tneighbours(surface, extent = 3)

source("tneighbours.r")
nbrs1  <- tneighbours(face)$nbr.pts
nbrs2  <- tneighbours(face, extent = 2)$nbr.pts
nbrs3  <- tneighbours(face, extent = 3)$nbr.pts


i <- 100
spheres3d(matrix(surface$coords[i, ], nrow = 1), radius = 0.5, col = "red")
spheres3d(surface$coords[nbrs1[[i]], ], radius = 0.5, col = "green")
pop3d()
spheres3d(surface$coords[nbrs2[[i]], ], radius = 0.5, col = "blue")
pop3d()
spheres3d(surface$coords[nbrs3[[i]], ], radius = 0.5, col = "yellow")


nbrs1 <- neighbours(face)
coef  <- lapply(nbrs1, function(x)
             coef(lsfit(face$coords[x, 1:2], face$coords[x, 3])))

