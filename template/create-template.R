#     Create a template face

setwd("~/ownCloud/Face3D_0.1-1/template")

library(Face3D)
library(rgl)
library(MASS)

# Choose male or female
sex <- "male"

# Load Artec version of Liberty's face
load("pretemplate.Rda")
face      <- trimimage.face3d(face)
face$lmks <- face$lmks[-c(11, 12, 22, 23), ]
template  <- face
plot(template)

# Load control data and compute means (lmks + curves + mesh)
load("all.Rda")
all.m   <- gpa$rotated[ , , 1:61]
all.f   <- gpa$rotated[ , , 62:130]
mean.m  <- apply(all.m, 1:2, mean)
mean.f  <- apply(all.f, 1:2, mean)
ntotal  <- dim(mean.f)[1]

mn  <- if (sex == "male") mean.m else mean.f

# Warp template to mean
from  <- rbind(template$lmks, template$curves, template$mesh)
to    <- mn
carry <- template
Temp  <- warp.face3d(from, to, carry)

plot(Temp, new = FALSE)
spheres3d(Temp$lmks, col = "blue", radius = 2)
spheres3d(Temp$mesh)
spheres3d(Temp$curves, col = "red")

# Rotate to the standard x-y-z axes
rotation.lmks <- c("sn", "n", "exR", "exL")                               
angles        <- rotate.face3d(gpa$mshape, id.lndms = c(4, 6, 8, 7), rotation = "coronal")$angles
Temp$coords   <- rotate3d(Temp$coords, angles[1], 0, 0, 1)
Temp$coords   <- rotate3d(Temp$coords, angles[2], 1, 0, 0)
Temp$coords   <- rotate3d(Temp$coords, angles[3], 0, 1, 0)
Temp$mesh     <- rotate3d(Temp$mesh,   angles[1], 0, 0, 1)
Temp$mesh     <- rotate3d(Temp$mesh,   angles[2], 1, 0, 0)
Temp$mesh     <- rotate3d(Temp$mesh,   angles[3], 0, 1, 0)
Temp$lmks     <- rotate3d(Temp$lmks,   angles[1], 0, 0, 1)
Temp$lmks     <- rotate3d(Temp$lmks,   angles[2], 1, 0, 0)
Temp$lmks     <- rotate3d(Temp$lmks,   angles[3], 0, 1, 0)
Temp$curves   <- rotate3d(Temp$curves, angles[1], 0, 0, 1)
Temp$curves   <- rotate3d(Temp$curves, angles[2], 1, 0, 0)
Temp$curves   <- rotate3d(Temp$curves, angles[3], 0, 1, 0)

plot(Temp, new = FALSE)

# Cut the face in half to create the right hand side
ind                <- which(Temp$coords[,1] < 0) 
right.face         <- list(coords = matrix(Temp$coords[ind, ], ncol = 3,
                                           dimnames = list(as.character(ind))))
trpls              <- matrix(Temp$triples, ncol = 3, byrow = TRUE)
indt               <- (trpls[ , 1] %in% ind) & (trpls[ , 2] %in% ind) & (trpls[ , 3] %in% ind)
trpls              <- c(t(trpls[indt, ]))
vec                <- 1:length(ind)
nms                <- format(ind, scientific = FALSE)
nms                <- sub("^ +", "", nms)                   # remove leading zeroes 
names(vec)         <- nms
nms1               <- format(trpls, scientific = FALSE)
nms1               <- sub("^ +", "", nms1)

right.face[["triples"]]   <- vec[nms1]
right.face[["colour"]]    <- Temp$colour[ind]
class(right.face)         <- "face3d"

# Crop the right edge neatly
plot(right.face, new = FALSE)
ind <- grep("brow ridge right", rownames(Temp$curves))
spheres3d(Temp$curves[ind, ])
min.z <- min(Temp$curves[ind, 3])
right.face <- subset(right.face, right.face$coords[ , 3] >= min.z)
plot(right.face, new = FALSE)

#  Create the left side by reflection
left.face              <- right.face
left.face$coords[ , 1] <- -left.face$coords[ , 1]
left.face$triples      <- rev(left.face$triples)
left.face$colour       <- right.face$colour

plot(right.face, display = "mesh",    new = FALSE)
plot(right.face, display = "normals", add = TRUE)
plot( left.face, display = "mesh",    add = TRUE)
plot( left.face, display = "normals", add = TRUE)

#  Finding the ordered edge points
edge.pts <- edges.face3d(right.face)[[1]]
plot(right.face, display = "mesh", new = FALSE)
plot(left.face,  display = "mesh", add = TRUE)
# for (i in 1:length(edge.r[[1]])) {
#    spheres3d(right.face$coords[edge.pts[i], ], col = "green")
#    spheres3d( left.face$coords[edge.pts[i], ], col = "green")
#    scan()
# }
spheres3d(right.face$coords[edge.pts, ], col = "green", radius = 0.3)
spheres3d( left.face$coords[edge.pts, ], col = "blue",  radius = 0.3)

# Restrict to points where left and right are close
dst <- apply(right.face$coords[edge.pts, ] - left.face$coords[edge.pts, ], 1, function(x) sqrt(sum(x^2)))
ind <- (dst < 10)
edge.pts <- edge.pts[ind]
plot(right.face, display = "mesh", new = FALSE)
plot( left.face, display = "mesh", add = TRUE)
spheres3d(right.face$coords[edge.pts, ], col = "green", radius = 0.3)
spheres3d( left.face$coords[edge.pts, ], col = "blue", radius = 0.3)

# Manual adjustment to omit the first and last few, by inspection
edge.pts <- edge.pts[5:(length(edge.pts) - 4)]
plot(right.face, display = "mesh", new = FALSE)
plot( left.face, display = "mesh", add = TRUE)
spheres3d(right.face$coords[edge.pts, ], col = "green", radius = 0.3)
spheres3d( left.face$coords[edge.pts, ], col = "blue", radius = 0.3)

# Find the mid-point of each pair of edge points and triangulate the surrounding area
trp.new <- integer(0)
# mid.old  <- (right.face$coords[edge.pts[1], ] + left.face$coords[edge.pts[1], ]) / 2
# mid.pts  <- matrix(mid.old, nrow = 1)
mid.pts <- (right.face$coords[edge.pts, ] + left.face$coords[edge.pts, ]) / 2
mid.col <- right.face$colour[edge.pts]
mid.ind <- TRUE
ncoords  <- nrow(right.face$coords)
for (i in 2:length(edge.pts)) {
   # mid.new <- (right.face$coords[edge.pts[i], ] + left.face$coords[edge.pts[i], ]) / 2
   projn <- sum((right.face$coords[edge.pts[i], ] - mid.pts[i - 1, ]) *
                                    (mid.pts[i, ] - mid.pts[i - 1, ]))
   if (projn > 0) {
      imid    <- sum(as.numeric(mid.ind)) + 2 * ncoords
      trp.new <- c(trp.new, c(imid, imid + 1, edge.pts[i]),
                            c(imid, edge.pts[i] + ncoords, imid + 1),
                            c(imid, edge.pts[i], edge.pts[i - 1]),
                            c(imid, edge.pts[i - 1] + ncoords, edge.pts[i] + ncoords))
      mid.ind <- c(mid.ind, TRUE)
   }
   else {
      trp.new <- c(edge.pts[i - 1], edge.pts[i + 1], edge.pts[i],
                   edge.pts[i - 1] + nccords, edge.pts[i] + ncoords, edge.pts[i + 1] + ncoords)
      mid.ind <- c(mid.ind, FALSE)
   }
}

crds <- rbind(right.face$coords, left.face$coords, mid.pts[mid.ind, ])
rownames(crds) <- c(paste("right", 1:nrow(right.face$coords)),
                    paste("left",  1:nrow(left.face$coords)),
                    paste("mid",   1:nrow(mid.pts)))
template <- as.face3d(list(coords  = crds,
                           triples = c(right.face$triples, left.face$triples + ncoords, trp.new),
                           lmks = Temp$lmks, curves = Temp$curves, mesh = Temp$mesh,
                           colour = c(right.face$colour, left.face$colour,  mid.col[mid.ind])))

# Check that the mid-line triangles are correct
# plot(template, display = "mesh", new = FALSE)
# trpls <- template$triples[2 * length(right.face$triples) + 1:length(trp.new)]
# trpls <- matrix(trpls, ncol = 3, byrow = TRUE)
# for (i in 1:nrow(trpls)) {
#    spheres3d(t(template$coords[trpls[i, 1], ]), col = "blue", radius = 0.2)
#    scan()
#    spheres3d(t(template$coords[trpls[i, 2], ]), col = "blue", radius = 0.2)
#    scan()
#    spheres3d(t(template$coords[trpls[i, 3], ]), col = "blue", radius = 0.2)
#    scan()
#    pop3d()
#    pop3d()
#    pop3d()
# }

# plot(template, new = FALSE, display = "mesh")
# plot(template, new = FALSE)
# plot(template, new = FALSE, col = "grey")

if (sex == "male") template.male <- template else template.female <- template

# save(template.male,   file = "../Face3D/data/templateMale.Rda")
# save(template.female, file = "../Face3D/data/templateFemale.Rda")

