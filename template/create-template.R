#     Create a template face

setwd("~/OneDrive - University of Glasgow/research/face3d")
install.packages("face3d", repos = NULL, type = "source")
library(face3d)
library(face3d)
library(rgl)
library(MASS)

setwd("~/OneDrive - University of Glasgow/research/face3d/template")

# Choose male or female
sex <- "female"

# Load Artec version of Liberty's face
load("pretemplate.Rda")
face$vertices <- face$coords
face$triangles <- matrix(face$triples, ncol = 3, byrow = TRUE)
face$landmarks <- face$lmks
face      <- trimimage.face3d(face)
face$landmarks <- face$landmarks[-c(11, 12, 22, 23), ]
template  <- face
plot(template)

# Load control data and compute means (landmarks + curves + mesh)
load("all.Rda")
all.m   <- gpa$rotated[ , , 1:61]
all.f   <- gpa$rotated[ , , 62:130]
mean.m  <- apply(all.m, 1:2, mean)
mean.f  <- apply(all.f, 1:2, mean)
ntotal  <- dim(mean.f)[1]

mn  <- if (sex == "male") mean.m else mean.f

# Warp template to mean
from  <- rbind(template$landmarks, template$curves, template$mesh)
to    <- mn
carry <- template
Temp  <- warp.face3d(from, to, carry)

plot(Temp)
spheres3d(Temp$landmarks, col = "blue", radius = 2)
spheres3d(Temp$mesh)
spheres3d(Temp$curves, col = "red")

# Rotate to the standard x-y-z axes
rotation.landmarks <- c("sn", "n", "exR", "exL")                               
angles         <- rotate.face3d(gpa$mshape, id.lndms = c(4, 6, 8, 7), rotation = "coronal")$angles
Temp$vertices  <- rotate3d(Temp$vertices,  angles[1], 0, 0, 1)
Temp$vertices  <- rotate3d(Temp$vertices,  angles[2], 1, 0, 0)
Temp$vertices  <- rotate3d(Temp$vertices,  angles[3], 0, 1, 0)
Temp$mesh      <- rotate3d(Temp$mesh,      angles[1], 0, 0, 1)
Temp$mesh      <- rotate3d(Temp$mesh,      angles[2], 1, 0, 0)
Temp$mesh      <- rotate3d(Temp$mesh,      angles[3], 0, 1, 0)
Temp$landmarks <- rotate3d(Temp$landmarks, angles[1], 0, 0, 1)
Temp$landmarks <- rotate3d(Temp$landmarks, angles[2], 1, 0, 0)
Temp$landmarks <- rotate3d(Temp$landmarks, angles[3], 0, 1, 0)
Temp$curves    <- rotate3d(Temp$curves,    angles[1], 0, 0, 1)
Temp$curves    <- rotate3d(Temp$curves,    angles[2], 1, 0, 0)
Temp$curves    <- rotate3d(Temp$curves,    angles[3], 0, 1, 0)

plot(Temp)

# Cut the face in half to create the right hand side
ind                <- which(Temp$vertices[,1] < 0) 
right.face         <- list(vertices = matrix(Temp$vertices[ind, ], ncol = 3,
                                           dimnames = list(as.character(ind))))
trpls              <- Temp$triangles
indt               <- (trpls[ , 1] %in% ind) & (trpls[ , 2] %in% ind) & (trpls[ , 3] %in% ind)
trpls              <- c(t(trpls[indt, ]))
vec                <- 1:length(ind)
nms                <- format(ind, scientific = FALSE)
nms                <- sub("^ +", "", nms)                   # remove leading zeroes 
names(vec)         <- nms
nms1               <- format(trpls, scientific = FALSE)
nms1               <- sub("^ +", "", nms1)

right.face[["triangles"]]   <- matrix(vec[nms1], ncol = 3, byrow = TRUE)
right.face[["colour"]]    <- Temp$colour[ind]
class(right.face)         <- "face3d"

# Crop the right edge neatly
plot(right.face, new = FALSE)
ind <- grep("brow ridge right", rownames(Temp$curves))
spheres3d(Temp$curves[ind, ])
min.z <- min(Temp$curves[ind, 3])
right.face <- subset(right.face, right.face$vertices[ , 3] >= min.z)
plot(right.face, new = FALSE)

#  Create the left side by reflection
left.face                <- right.face
left.face$vertices[ , 1] <- -left.face$vertices[ , 1]
left.face$triangles      <- matrix(rev(c(t(left.face$triangles))), ncol = 3, byrow = TRUE)
left.face$colour         <- right.face$colour

plot(right.face, display = "mesh")
plot(right.face, display = "normals", add = TRUE)
plot( left.face, display = "mesh",    add = TRUE)
plot( left.face, display = "normals", add = TRUE)

#  Finding the ordered edge points
edge.pts <- edges.face3d(right.face)[[1]]
plot(right.face, display = "mesh", new = FALSE)
plot(left.face,  display = "mesh", add = TRUE)
# for (i in 1:length(edge.r[[1]])) {
#    spheres3d(right.face$vertices[edge.pts[i], ], col = "green")
#    spheres3d( left.face$vertices[edge.pts[i], ], col = "green")
#    scan()
# }
spheres3d(right.face$vertices[edge.pts, ], col = "green", radius = 0.3)
spheres3d( left.face$vertices[edge.pts, ], col = "blue",  radius = 0.3)

# Restrict to points where left and right are close
dst <- apply(right.face$vertices[edge.pts, ] - left.face$vertices[edge.pts, ], 1, function(x) sqrt(sum(x^2)))
ind <- (dst < 10)
edge.pts <- edge.pts[ind]
plot(right.face, display = "mesh", new = FALSE)
plot( left.face, display = "mesh", add = TRUE)
spheres3d(right.face$vertices[edge.pts, ], col = "green", radius = 0.3)
spheres3d( left.face$vertices[edge.pts, ], col = "blue", radius = 0.3)

# Manual adjustment to omit the first and last few, by inspection
edge.pts <- edge.pts[5:(length(edge.pts) - 4)]
plot(right.face, display = "mesh", new = FALSE)
plot( left.face, display = "mesh", add = TRUE)
spheres3d(right.face$vertices[edge.pts, ], col = "green", radius = 0.3)
spheres3d( left.face$vertices[edge.pts, ], col = "blue", radius = 0.3)

# Find the mid-point of each pair of edge points and triangulate the surrounding area
trp.new <- integer(0)
# mid.old  <- (right.face$vertices[edge.pts[1], ] + left.face$vertices[edge.pts[1], ]) / 2
# mid.pts  <- matrix(mid.old, nrow = 1)
mid.pts <- (right.face$vertices[edge.pts, ] + left.face$vertices[edge.pts, ]) / 2
mid.col <- right.face$colour[edge.pts]
mid.ind <- TRUE
ncoords  <- nrow(right.face$vertices)
for (i in 2:length(edge.pts)) {
   # mid.new <- (right.face$vertices[edge.pts[i], ] + left.face$vertices[edge.pts[i], ]) / 2
   projn <- sum((right.face$vertices[edge.pts[i], ] - mid.pts[i - 1, ]) *
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

crds <- rbind(right.face$vertices, left.face$vertices, mid.pts[mid.ind, ])
rownames(crds) <- c(paste("right", 1:nrow(right.face$vertices)),
                    paste("left",  1:nrow(left.face$vertices)),
                    paste("mid",   1:nrow(mid.pts)))
r.triples <- c(t(right.face$triangles))
l.triples <- c(t(left.face$triangles))
template <- as.face3d(list(vertices  = crds,
                           triangles = matrix(c(r.triples, l.triples + ncoords, trp.new), ncol = 3, byrow = TRUE),
                           landmarks = Temp$landmarks, curves = Temp$curves, mesh = Temp$mesh,
                           colour = c(right.face$colour, left.face$colour,  mid.col[mid.ind])))

# Check that the mid-line triangles are correct
# plot(template, display = "mesh", new = FALSE)
# Error here after changing to triangles rather than triples
# trpls <- template$triples[2 * length(right.face$triples) + 1:length(trp.new)]
# trpls <- matrix(trpls, ncol = 3, byrow = TRUE)
# for (i in 1:nrow(trpls)) {
#    spheres3d(t(template$vertices[trpls[i, 1], ]), col = "blue", radius = 0.2)
#    scan()
#    spheres3d(t(template$vertices[trpls[i, 2], ]), col = "blue", radius = 0.2)
#    scan()
#    spheres3d(t(template$vertices[trpls[i, 3], ]), col = "blue", radius = 0.2)
#    scan()
#    pop3d()
#    pop3d()
#    pop3d()
# }

# plot(template, display = "mesh")
# plot(template)
# plot(template, col = "grey")

if (sex == "male") template_male <- template else template_female <- template

# save(template_male,   file = "../face3d/data/template_male.Rda")
# save(template_female, file = "../face3d/data/template_female.Rda")

