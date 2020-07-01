#     Test code for the edges.face3d function

install.packages("~/OneDrive - University of Glasgow/research/face3d_0.1-1/face3d",
                 repos = NULL, type = "source")
library(face3d)
library(fields)
library(rgl)
library(geometry)


load("testing/test-data/edges.RData")
plot(newshape)
edges <- edges.face3d(newshape)
lapply(edges, function(x) rgl::lines3d(newshape$coords[x, ], lwd = 5, col = "blue"))


fls <- list.files("Data-Liberty", full.names = TRUE)

i <- 3
load(fls[i])
face  <- subset(face, edist.face3d(face$coords, face$lmks["pn", ]) < 100)
plot(face)
edges <- edges.face3d(face)
spheres3d(face$coords[c(edges), ], col = "yellow")
pop3d()
segments3d(face$coords[c(t(edges)), ], lwd = 3, col = "blue")

# How to put the edges in order

ne <- nrow(edges)
dst <- outer(1:ne, 1:ne, function(x, y) (edges[x, 1] == edges[y, 1]) | (edges[x, 1] == edges[y, 2]) |
                (edges[x, 2] == edges[y, 1]) | (edges[x, 2] == edges[y, 2]))
diag(dst) <- FALSE
apply(dst, 1, sum)
nxt <- t(apply(edges, 1, function(x) which((x[1] == edges[ , 1] | x[1] == edges[ , 2] | 
                                               x[2] == edges[ , 1] | x[2] == edges[ , 2]) &
                                              !(x[1] == edges[ , 1] & x[2] == edges[ , 2]))
))
ord <- 1
for (i in 1:nrow(edges)) {
   
}
