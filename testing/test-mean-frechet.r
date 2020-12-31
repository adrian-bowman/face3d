# setwd("/Users/libertyvittert/PHD/Faces/R codes/R")
for (file in list.files()) source(file)

library(rgl)
library(fields)
  
load("/Users/libertyvittert/PHD/Faces/all images/Liberty-Controls-dmp/Liberty-Controls-dmp-final/controls-liberty-001.dmp")
# picking out a shape to work with
lmks   <- rbind(face$lmks["chR", ], face$lmks["cphR", ])
rng    <- sqrt(sum((lmks[2, ] - lmks[1, ])^2))
bndry  <- c(.1,.3) * rng
unit   <- (lmks[2, ] - lmks[1, ]) / rng
prjn   <- c(sweep(face$coords, 2, lmks[1, ]) %*% unit)
ind1   <- (prjn > - bndry[1]) & (prjn < rng + bndry[1])
prjn2  <- outer(prjn, unit)
ind2   <- apply((sweep(face$coords, 2, lmks[1, ]) - prjn2)^2, 1, function(x) sqrt(sum(x)) < bndry[2])
shape  <- subset.face3d(face, ind1 & ind2, remove.singles = TRUE)	  
shape  <- index.face3d(shape, distance=5)
values <- pmax(abs(shape$kappa1), abs(shape$kappa2)) 
mean   <- mean_frechet.face3d(shape, values)

