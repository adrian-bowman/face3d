"distances.face3d" <- function(source.shape, target.shape, warp = TRUE, procrustes = TRUE, template.shape){
    if (("curves" %in% names(source.shape)) == FALSE){
       }
    if (warp == TRUE){
        if (("lmks" %in% names(source.shape)) & ("curves" %in% names(source.shape))){
            source.shape$coords  <- warp.face3d(template$coords,rbind(template$lmks,template$curves),rbind(source.shape$lmks,source.shape$curves))
            source.shape$triples <- template$triples
            target.shape$coords  <- warp.face3d(template$coords,rbind(template$lmks,template$curves),rbind(target.shape$lmks,target.shape$curves))
            target.shape$triples <- template$triples
            mean.shape$coords    <- (source.shape$coords + target.shape$coords)/2
            mean.shape$lmks      <- (source.shape$lmks + target.shape$lmks)/2
            }
        if (("curves" %in% names(source.shape)) == FALSE){
            source.shape$coords  <- warp.face3d(template$coords,template$lmks,source.shape$lmks)
            source.shape$triples <- template$triples
            target.shape$coords  <- warp.face3d(template$coords,template$lmks,target.shape$lmks)
            target.shape$triples <- template$triples
            mean.shape$coords    <- (source.shape$coords + target.shape$coords)/2
            mean.shape$lmks      <- (source.shape$lmks + target.shape$lmks)/2
                }
        }
    ## area of a triangle
    n                    <- length(source.shape$triples)/3
    source.triples.array <- array(t(source.shape$coords[source.shape$triples,]),c(3,3,n))
    source.triples.areas <- apply(source.triples.array,3, function(x){
                                  1/2 * sqrt(sum((crossproduct(x[,2] - x[,1],x[,3] - x[,1],scale = FALSE))^2))
                            })
    target.triples.array <- array(t(target.shape$coords[target.shape$triples,]),c(3,3,n))
    target.triples.areas <- apply(target.triples.array,3, function(x){
                                  1/2 * sqrt(sum((crossproduct(x[,2] - x[,1],x[,3] - x[,1],scale = FALSE))^2))
                            })
    area.ratio <- log(target.triples.areas/source.triples.areas)

    if (procrustes == TRUE){
        if (("lmks" %in% names(source.shape)) & ("curves" %in% names(source.shape))){
            k1a                  <- dim(source.shape$lmks)[1]
            k1b                  <- dim(source.shape$curves)[1]
            k1                   <- k1a + k1b
            k2                   <- dim(source.shape$coords)[1]
            k                    <- k1 + k2
            LNDMS                <- array(0,c(k,3,2))
            LNDMS[,,1]           <- rbind(source.shape$lmks,source.shape$curves,source.shape$coords)
            LNDMS[,,2]           <- rbind(target.shape$lmks,target.shape$curves,target.shape$coords)
            id.subset            <- 1:k1
            GPA                  <- procrustes.face3d(LNDMS, scale = TRUE, rotation = "coronal", rotation.lmks = c(4,6,8,7), subset = id.subset)
            PSC                  <- GPA$PSC
            id.subset            <- 1:k1
            id.lmks              <- 1:k1a
            id.curves            <- (k1a + 1):k1
            mean.shape$coords    <- GPA$PMS[-c(id.subset),]
            mean.shape$lmks      <- GPA$PMS[id.lmks,]
            mean.shape$curves    <- GPA$PMS[id.curves,]
            source.shape$coords  <- PSC[-c(id.subset),,1]
            source.shape$lmks    <- PSC[id.lmks,,1]
            source.shape$curves  <- PSC[id.curves,,1]
            target.shape$coords  <- PSC[-c(id.subset),,2]
            target.shape$lmks    <- PSC[id.lmks,,2]
            target.shape$curves  <- PSC[id.curves,,2]
            }
        if (("curves" %in% names(source.shape)) == FALSE){
            k1                   <- dim(source.shape$lmks)[1]
            k2                   <- dim(source.shape$coords)[1]
            k                    <- k1 + k2
            LNDMS                <- array(0,c(k,3,2))
            LNDMS[,,1]           <- rbind(source.shape$lmks,source.shape$coords)
            LNDMS[,,2]           <- rbind(target.shape$lmks,target.shape$coords)
            id.subset            <- 1:k1
            GPA                  <- procrustes.face3d(LNDMS, scale = TRUE, rotation = "coronal", rotation.lmks = c(4,6,8,7), subset = id.subset)
            PSC                  <- GPA$PSC
            mean.shape$coords    <- GPA$PMS[-c(id.subset),]
            mean.shape$lmks      <- GPA$PMS[id.subset,]
            source.shape$coords  <- PSC[-c(id.subset),,1]
            source.shape$lmks    <- PSC[id.subset,,1]
            target.shape$coords  <- PSC[-c(id.subset),,2]
            target.shape$lmks    <- PSC[id.subset,,2]
            }
          }
    if (!procrustes){
        mean.shape            <- source.shape
        mean.shape$coords     <- (source.shape$coords + target.shape$coords)/2
        mean.shape$lmks       <- (source.shape$lmks + target.shape$lmks)/2
    }

    ## Distance in normal direction
    local.axes             <- index.face3d(source.shape)$local.axes
    k2                     <- dim(source.shape$coords)[1]
    local.coords.1         <- local.coords.2 <- numeric(k2)
    for (i in 1:k2) local.coords.1[i] <- as.vector(t(source.shape$coords[i,]) %*% local.axes[,,i])[1]
    for (i in 1:k2) local.coords.2[i] <- as.vector(t(target.shape$coords[i,]) %*% local.axes[,,i])[1]
    dist.norm <- local.coords.2 - local.coords.1
    ## 3d sign distance with respect to x,y,z direction
    dist.mat               <- target.shape$coords - source.shape$coords
    ## 3d sign distance with respect to x direction
    sign.lr.x              <- sign(dist.mat[,1]) # left-right on the screen
    index.dir              <- rep(1,k2)
    index                  <- which(template.shape$coords[,1] <= 0)
    index.dir[index]       <- -1
    sign.lr.x              <- sign.lr.x * index.dir
    dist.x    <- sign.lr.x * abs(dist.mat[,1])
    ## 3d sign distance with respect to y direction
    sign.lr.y              <- sign(dist.mat[,2]) # up-down on the screen
    dist.y    <- sign.lr.y * abs(dist.mat[,2])
    ## 3d sign distance with respect to z direction
    sign.lr.z              <- sign(dist.mat[,3]) # in-out of the screen
    dist.z    <- sign.lr.z * abs(dist.mat[,3])

    result <- list(dist.norm = dist.norm, area.ratio = area.ratio,
                   dist.x = dist.x, dist.y = dist.y, dist.z = dist.z)
    cat("Range of dist normal:  ", round(range(dist.norm),4), "\n")
    cat("Range of dist x-axis:  ", round(range(dist.x),4), "\n")
    cat("Range of dist y-axis:  ", round(range(dist.y),4), "\n")
    cat("Range of dist z-axis:  ", round(range(dist.z),4), "\n")
    cat("Range of triangle area:", round(range(area.ratio),4), "\n")
    return(result)
    }
