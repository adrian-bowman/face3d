#            Tests for the closest.face3d function

install.packages("~/OneDrive - University of Glasgow/research/face3d_0.1-1/face3d",
                 repos = NULL, type = "source")
library(face3d)
library(rgl)
library(fields)

# Examples used in the help file

# Projection onto vertices, edges and faces
shape   <- as.face3d(list(coords = matrix(c(0, 0, 0, 1, 0, 0, 0.5, 1, 0, 0.5, 0.5, 1),
                              ncol = 3, byrow = TRUE),
                        triples = c(1, 2, 4, 2, 3, 4, 1, 4, 3, 3, 2, 1)))
pts     <- matrix(c(1.2, -0.2, 0, 0.5, 1.2, 0.5, 0.3, 0.5, 0, 0, 0, 0, 0.5, 0, 0, 0.5, 0, 0.5),
                  ncol = 3, byrow = TRUE)
closest <- t(apply(pts, 1, function(x) closest.face3d(x, shape)$point))
plot(shape)
spheres3d(pts,     radius = 0.02, col = "red")
spheres3d(closest, radius = 0.03, col = "green")
segments3d(matrix(c(t(cbind(pts, closest))), ncol = 3, byrow = TRUE))

# Projection onto the nose of template_male
x    <- template_male$lmks["pn", ] + c(0.2, 0.1, 1)
sbst <- subset(template_male, rdist(t(x), template_male$coords) < 10)
cx   <- closest.face3d(x, template_male)
plot(cx$subset, display = "mesh")
spheres3d(x, col = "red", radius = 0.1)
spheres3d(cx$point, col = "blue", radius = 0.1)
spheres3d(cx$subset$coords[as.character(cx$id), ], col = "green", radius = 0.1)

# ---------------------------------
# Other examples

spheres3d(x,  col = "red",   radius = 2)
spheres3d(cx$point, col = "green", radius = 1)
plot(cx$subset, display = "mesh")
spheres3d(cx$point, col = "green", radius = 0.1)

# Closest point is on an edge

load("testing/test-data/closest.Rdata")
pp <- planepath.face3d(shape, pt_1, pt_2)$path
cx <- closest.face3d(pp[2,], shape)
cx$id
plot(cx$subset, display = "mesh")
spheres3d(pp[2, ], radius = 0.5, col = "red")
spheres3d(cx$subset$coords[cx$id.sub, ], radius = 0.5, col = "green")
head(cx$subset$coords)
cx$subset$coords[cx$id, ]

all(as.numeric(rownames(cx$subset$triples)) %in% as.numeric(rownames(cx$subset$coords)))
matrix(cx$subset$triples, ncol = 3, byrow = TRUE)
