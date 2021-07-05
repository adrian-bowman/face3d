#     Smooth a triangulated surface

library(face3d)

load("~/Dropbox/Phylogeny-from-shape/R work/Reference-Shape/nose.RData")
names(nose)[1] <- "vertices"
nose$triangles <- matrix(nose$triples, ncol = 3, byrow = TRUE)
nose$triples   <- NULL
nvert          <- nrow(nose$vertices)
ind            <- which(1:nvert %in% unique(c(nose$triangles)))
nose           <- subset(nose, ind)
nvert          <- nrow(nose$vertices)

plot(nose)

nbrs.fn <- function(i) {
   connected <- apply(nose$triangles, 1, function(x) i %in% x)
   connected <- unique(c(nose$triangles[connected, ]))
   connected <- connected[-match(i, connected)]
   connected
}
nbrs <- lapply(1:nvert, nbrs.fn)
D    <- diag(1, nvert)
for (i in 1:nvert) D[i, nbrs[[i]]] <- -1 / length(nbrs[[i]])

lambda  <- 100
B1      <- solve(diag(nvert) + lambda * crossprod(D))
sm.nose <- nose
sm.nose$vertices <- B1 %*% nose$vertices
plot(nose, col = "blue")
plot(sm.nose, add = TRUE, col = "green")

