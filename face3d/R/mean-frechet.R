mean_frechet.face3d <- function(shape, values){

  xis           <- shape$vertices
  p             <- shape$vertices
  weights       <- values
  edges <- edges.face3d(shape)
  xis           <-     xis[-unique(c(edges)), ]
  p             <-       p[-unique(c(edges)), ]
  weights       <- weights[-unique(c(edges))]
 #error when i-1, j 381
       #Error in path[nrow(path), ] : incorrect number of dimensions
       
   
   
  distances_pi  <- rep(NA, length(p[ ,1]))
  distances_ps  <- matrix(NA, ncol= length(p[ ,1]), nrow= length(p[ ,1]))
  final.sums.ps <- rep(NA, length(p[ ,1]))
  
  for(i in 1:length(p[ ,1])) {
     for(j in 1:length(xis[ ,1])) {
     	print(c(i, j))
        distances_pi[j] <- max(planepath(shape, p[i, ],xis[j, ], boundary = NA)$arclength)
     }
     distances_ps[ ,i] <- distances_pi
  }
  
  for (z in 1:length(p[,1]))
     final.sums.ps <- sum(weights* distances_ps[ ,z]^2)


min             <- which.min(final.sums.ps)
mean.frechet    <- shape$vertices[min, ]
invisible(mean)
}
