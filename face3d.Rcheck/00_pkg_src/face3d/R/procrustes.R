"procrustes.face3d" <- function(lmks, scale = TRUE, threshold = 0.03,
                         rotation = "coronal", rotation.lmks = c("sn", "n", "exR", "exL"), subset) {

  if (!requireNamespace("rgl", quietly = TRUE)) stop("procrustes.face3d requires the rgl package.")

  LNDMS    <- lmks
  id.lndms <- rotation.lmks
  scaling  <- scale
  
  if (missing(subset))
     subset <- FALSE
  else {
     id.subset <- subset
     if (class(id.subset) == "logical") id.subset <- which(subset)
     subset <- TRUE
  }
  
  d <- dim(LNDMS)[2]
  
  if(d == 2){
  	
    if (!subset) {
  k <- dim(LNDMS)[1]
  d <- dim(LNDMS)[2]
  n <- dim(LNDMS)[3]
  ## affine transformations
  PSC.i.all <- apply(LNDMS,3, function(x){
             CENTROID       <- apply(x,2,mean)
             LNDMS.centered <- cbind(outer(x[,1],CENTROID[1],"-"),
                                     outer(x[,2],CENTROID[2],"-"))
             CS             <- sqrt(sum(LNDMS.centered^2))
             if (scaling == TRUE)  LNDMS.centered.rescaled <- LNDMS.centered / CS
             if (scaling == FALSE) LNDMS.centered.rescaled <- LNDMS.centered
             LNDMS.centered.rescaled
  })
  PSC.i.all <- array(c(PSC.i.all), c(k, d, n), dimnames = list(dimnames(LNDMS)[[1]], NULL, NULL))
  CENTROID.all <- apply(LNDMS,3, function(x){
             apply(x,2,mean)
  })
  CS.all <- apply(LNDMS,3, function(x){
             CENTROID       <- apply(x,2,mean)
             LNDMS.centered <- cbind(outer(x[,1],CENTROID[1],"-"),
                                     outer(x[,2],CENTROID[2],"-"))
             sqrt(sum(LNDMS.centered^2))
  })
## optimal rotations
  SOURCE <- PSC.i.all[ , , 1]
  PSC.i.all.rot <- apply(PSC.i.all,3,  function(x){
            SVD <- svd(t(SOURCE)%*%x)
            MATrot <- SVD$v%*%t(SVD$u)
            TARGET.rot <- x%*%MATrot
            TARGET.rot
    })
  PSC.i.all.rot <- array(c(PSC.i.all.rot), c(k, d, n), dimnames = list(dimnames(LNDMS)[[1]], NULL, NULL))
  MEAN.i <- apply(PSC.i.all.rot,c(1,2),mean)
  DIST <- distance.3d.diff(SOURCE - MEAN.i)
  index <- 1
  index.vec <- 1
  while (DIST > threshold) {
  SOURCE <- MEAN.i
  PSC.i.all.rot <- apply(PSC.i.all.rot,3,  function(x){
            SVD <- svd(t(SOURCE)%*%x)
            MATrot <- SVD$v%*%t(SVD$u)
            TARGET.rot <- x%*%MATrot
            TARGET.rot
    })
  PSC.i.all.rot <- array(c(PSC.i.all.rot), c(k, d, n), dimnames = list(dimnames(LNDMS)[[1]], NULL, NULL))
  MEAN.i <- apply(PSC.i.all.rot,c(1,2),mean)
  DIST <- distance.3d.diff(SOURCE - MEAN.i)
  index <- index + 1
  #print(c(index,DIST))
  index.vec <- append(index.vec,index)
  }
}

## subset == TRUE
    if (subset) {
    	
  k1 <- dim(LNDMS)[1]
  LNDMS.subset <- LNDMS[id.subset, , ]
  k <- dim(LNDMS.subset)[1]
  d <- dim(LNDMS.subset)[2]
  n <- dim(LNDMS.subset)[3]
## affine transformations
  PSC.i.all <- apply(LNDMS,3, function(x){
             CENTROID       <- apply(x[id.subset,], 2, mean)
             LNDMS.centered <- cbind(outer(x[id.subset, 1], CENTROID[1], "-"),
                                     outer(x[id.subset, 2], CENTROID[2], "-"))
             CS             <- sqrt(sum(LNDMS.centered^2))
             LNDMS.centered <- cbind(outer(x[,1], CENTROID[1], "-"),
                                     outer(x[,2], CENTROID[2], "-"))
             if (scaling == TRUE)  LNDMS.centered.rescaled <- LNDMS.centered/CS
             if (scaling == FALSE) LNDMS.centered.rescaled <- LNDMS.centered
             LNDMS.centered.rescaled
  })
  dnms <- list(dimnames(LNDMS[id.subset, , ])[[1]], NULL, NULL)
  PSC.i.all <- array(c(PSC.i.all), c(k1, d, n), dimnames = dnms)
  CS.all <- apply(LNDMS.subset,3, function(x){
             CENTROID       <- apply(x, 2, mean)
             LNDMS.centered <- cbind(outer(x[,1], CENTROID[1], "-"),
                                     outer(x[,2], CENTROID[2], "-"))
             sqrt(sum(LNDMS.centered^2))
  })

 ## optimal rotations
  SOURCE <- PSC.i.all[,,1]
  PSC.i.all.rot <- apply(PSC.i.all,3,  function(x){
            SVD <- svd(t(SOURCE[id.subset,])%*%x[id.subset,])
            #print(sign(SVD$d))
            MATrot <- SVD$v%*%t(SVD$u)
            TARGET.rot <- x%*%MATrot
            TARGET.rot
    })
  PSC.i.all.rot <- array(c(PSC.i.all.rot), c(k1, d, n), dimnames = dnms)
  MEAN.i <- apply(PSC.i.all.rot,c(1,2),mean)
  DIST <- distance.3d.diff(SOURCE[id.subset,] - MEAN.i[id.subset,])
  index <- 1
  index.vec <- 1

  while (DIST > threshold) {
  SOURCE <- MEAN.i
  PSC.i.all.rot <- apply(PSC.i.all.rot,3,  function(x){
            SVD <- svd(t(SOURCE[id.subset,])%*%x[id.subset,])
            #print(sign(SVD$d))
            MATrot <- SVD$v%*%t(SVD$u)
            TARGET.rot <- x%*%MATrot
            TARGET.rot
    })
  PSC.i.all.rot <- array(c(PSC.i.all.rot), c(k1, d, n), dimnames = dnms)
  MEAN.i <- apply(PSC.i.all.rot,c(1,2),mean)
  DIST <- distance.3d.diff(SOURCE[id.subset,] - MEAN.i[id.subset,])
  index <- index + 1
  #print(c(index,DIST))
  index.vec <- append(index.vec,index)
  }
}
}

  if(d == 3){
  	
  if (!subset) {
  	
  k <- dim(LNDMS)[1]
  d <- dim(LNDMS)[2]
  n <- dim(LNDMS)[3]
  ## affine transformations
  PSC.i.all <- apply(LNDMS, 3, function(x){
             CENTROID       <- apply(x, 2, mean)  
             LNDMS.centered <- cbind(outer(x[,1], CENTROID[1], "-"),
                                     outer(x[,2], CENTROID[2], "-"),
                                     outer(x[,3], CENTROID[3], "-"))
             CS             <- sqrt(sum(LNDMS.centered^2))
             if (scaling == TRUE)  LNDMS.centered.rescaled <- LNDMS.centered/CS
             if (scaling == FALSE) LNDMS.centered.rescaled <- LNDMS.centered
             LNDMS.centered.rescaled
  })
  PSC.i.all <- array(c(PSC.i.all), c(k, d, n), dimnames = list(dimnames(LNDMS)[[1]], NULL, NULL))
  CENTROID.all <- apply(LNDMS, 3, function(x){
             apply(x, 2, mean)  
  })
  CS.all <- apply(LNDMS, 3, function(x){
             CENTROID       <- apply(x, 2, mean)  
             LNDMS.centered <- cbind(outer(x[,1], CENTROID[1], "-"),
                                     outer(x[,2], CENTROID[2], "-"),
                                     outer(x[,3], CENTROID[3], "-"))
             sqrt(sum(LNDMS.centered^2))
  })
## optimal rotations
  SOURCE <- PSC.i.all[ , , 1]
  PSC.i.all.rot <- apply(PSC.i.all, 3, function(x){
            SVD <- svd(t(SOURCE)%*%x)
            MATrot <- SVD$v%*%t(SVD$u)
            TARGET.rot <- x%*%MATrot
            TARGET.rot     
    })
  PSC.i.all.rot <- array(c(PSC.i.all.rot), c(k, d, n), dimnames = list(dimnames(LNDMS)[[1]], NULL, NULL))
  MEAN.i <- apply(PSC.i.all.rot, c(1, 2), mean)
  DIST <- distance.3d.diff(SOURCE - MEAN.i)
  index <- 1
  index.vec <- 1
  while (DIST > threshold) { 
  SOURCE <- MEAN.i
  PSC.i.all.rot <- apply(PSC.i.all.rot, 3, function(x){
            SVD <- svd(t(SOURCE)%*%x)
            MATrot <- SVD$v%*%t(SVD$u)
            TARGET.rot <- x%*%MATrot
            TARGET.rot     
    })
  PSC.i.all.rot <- array(c(PSC.i.all.rot), c(k, d, n), dimnames = list(dimnames(LNDMS)[[1]], NULL, NULL))
  MEAN.i <- apply(PSC.i.all.rot, c(1, 2), mean)
  DIST <- distance.3d.diff(SOURCE - MEAN.i)
  index <- index + 1
  #print(c(index,DIST))
  index.vec <- append(index.vec,index)
  }
}

## subset == TRUE
    if (subset) {

  k1 <- dim(LNDMS)[1]
  LNDMS.subset <- LNDMS[id.subset, , ]
  k <- dim(LNDMS.subset)[1]
  d <- dim(LNDMS.subset)[2]
  n <- dim(LNDMS.subset)[3]
  dnms <- list(dimnames(LNDMS)[[1]], NULL, NULL)
## affine transformations
  PSC.i.all <- apply(LNDMS,3, function(x){
             CENTROID       <- apply(x[id.subset,],2,mean)  
             LNDMS.centered <- cbind(outer(x[id.subset,1],CENTROID[1],"-"),
                                     outer(x[id.subset,2],CENTROID[2],"-"),
                                     outer(x[id.subset,3],CENTROID[3],"-"))
             CS             <- sqrt(sum(LNDMS.centered^2))
             LNDMS.centered <- cbind(outer(x[,1],CENTROID[1],"-"),
                                     outer(x[,2],CENTROID[2],"-"),
                                     outer(x[,3],CENTROID[3],"-"))
             if (scaling == TRUE)  LNDMS.centered.rescaled <- LNDMS.centered/CS
             if (scaling == FALSE) LNDMS.centered.rescaled <- LNDMS.centered
             LNDMS.centered.rescaled
  })
  PSC.i.all <- array(c(PSC.i.all), c(k1, d, n), dimnames = dnms)
  CS.all <- apply(LNDMS.subset, 3, function(x){
             CENTROID       <- apply(x,2,mean)  
             LNDMS.centered <- cbind(outer(x[,1],CENTROID[1],"-"),
                                     outer(x[,2],CENTROID[2],"-"),
                                     outer(x[,3],CENTROID[3],"-"))
             sqrt(sum(LNDMS.centered^2))
  }) 
 
 ## optimal rotations
  SOURCE <- PSC.i.all[,,1]
  PSC.i.all.rot <- apply(PSC.i.all,3,  function(x){
            SVD <- svd(t(SOURCE[id.subset,])%*%x[id.subset,])
            #print(sign(SVD$d))
            MATrot <- SVD$v%*%t(SVD$u)
            TARGET.rot <- x%*%MATrot
            TARGET.rot     
    })
  PSC.i.all.rot <- array(c(PSC.i.all.rot), c(k1, d, n), dimnames = dnms)
  MEAN.i <- apply(PSC.i.all.rot, c(1, 2), mean)
  DIST <- distance.3d.diff(SOURCE[id.subset,] - MEAN.i[id.subset,])
  index <- 1
  index.vec <- 1

  while (DIST > threshold) { 
  SOURCE <- MEAN.i
  PSC.i.all.rot <- apply(PSC.i.all.rot,3,  function(x){
            SVD <- svd(t(SOURCE[id.subset,])%*%x[id.subset,])
            #print(sign(SVD$d))
            MATrot <- SVD$v%*%t(SVD$u)
            TARGET.rot <- x%*%MATrot
            TARGET.rot     
    })
  PSC.i.all.rot <- array(c(PSC.i.all.rot), c(k1, d, n), dimnames = dnms)
  MEAN.i <- apply(PSC.i.all.rot, c(1, 2), mean)
  DIST <- distance.3d.diff(SOURCE[id.subset,] - MEAN.i[id.subset,])
  index <- index + 1
  #print(c(index,DIST))
  index.vec <- append(index.vec,index)
  }
}

## rotate to frontal
  if (rotation == "coronal") {
     if (missing(id.lndms)) stop("id.lndms missing.") 

  k      <- dim(LNDMS)[1]
  d      <- dim(LNDMS)[2]
  n      <- dim(LNDMS)[3]
  ROT    <- rotate.face3d(MEAN.i, id.lndms, rotation = "coronal")
  MEAN.i <- ROT$ROT.shape
  angles <- as.numeric(ROT$angles)
  PSC.i.all.rot <- apply(PSC.i.all.rot, 3,  function(x){
                         ROT.shape <- rgl::rotate3d(x,         angles[1], 0, 0, 1)
                         ROT.shape <- rgl::rotate3d(ROT.shape, angles[2], 1, 0, 0) 
                         ROT.shape <- rgl::rotate3d(ROT.shape, angles[3], 0, 1, 0)
                         ROT.shape
    })
  PSC.i.all.rot <- array(c(PSC.i.all.rot), c(k, d, n), dimnames = list(dimnames(PSC.i.all.rot)[[1]], NULL, NULL))
  }
    
}
  RESULTS    <- list(rotated = PSC.i.all.rot, mean = MEAN.i, centroid.sizes = CS.all, index = index, DIST = DIST, index.vec = index.vec) 
  return(RESULTS) 
 }
 
 
  "distance.3d.diff" <- function(Xdiff) sqrt(sum((Xdiff)^2))
