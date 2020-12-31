"sliding.face3d" <- function(lmks1,lmks2dense,resample = FALSE,n=100){

    if (!requireNamespace("MASS", quietly = TRUE)) stop("the MASS package is required.")
   
    # template = lmks1
    # if resample = TRUE, n > k.dense  
    k           <- dim(lmks1)[1]
    lmks2 <- resample.face3d(lmks2dense,n=k,threshold=1e-5)
    if (resample == FALSE) lmks2.rsmpl <- lmks2dense
    if (resample == TRUE) {
              k.dense           <- dim(lmks2dense)[1]
              if (k.dense > n) stop("number of resampled points is smaller than \n number of points on lmks2dense.") 
              else lmks2.rsmpl <- resample.face3d(lmks2dense,n=n,threshold=1e-5)
       }
    TPS <- tps.face3d(lmks1,lmks2)
    B <- TPS$B
    Be <- TPS$Be
    ## optimal bending energy calculation
    Be.i <- matrix(0,k,n)
    lmks2.optim <- lmks2
    for(i in 1:k){
        for(j in 1:n){
            lmks2.i     <- lmks2
            lmks2.i[i,] <- lmks2.rsmpl[j,]
            Be.i[i,j]   <- bendingenergy.face3d(B,lmks2.i)
        }
    }
    ## reordering
    optim.ord <- numeric(k)
    for(i in 2:(k-1)){
        optim.ord[i]    <- which.min(abs(Be.i[i,]))
        lmks2.optim[i,] <- lmks2.rsmpl[optim.ord[i],]
        lmks2.optim[1,] <- lmks2[1,]
        lmks2.optim[k,] <- lmks2[k,]
    }
    ## optimal bending energy
    Be.optim <- bendingenergy.face3d(B,lmks2.optim)
    ## resutls
    return(lmks2.optim)
}


"cumchord" <- function(MAT){
x <- cumsum(sqrt(apply((MAT-rbind(MAT[1,],MAT[-(dim(MAT)[1]),]))^2,1,sum)))
return(x)
}

"bendingenergy.face3d" <- function(B,lmks){
    lmks.vec   <- c(lmks[,1],lmks[,2],lmks[,3])
    Be         <- t(lmks.vec)%*%B%*%lmks.vec
    Be         <- as.numeric(Be)
    return(Be)
}

"tps.face3d" <- function(lmks1,lmks2){
    if (!requireNamespace("MASS", quietly = TRUE)) stop(" the MASS package is required.")
    k         <- dim(lmks1)[1]
    Sdx        <- outer(lmks1[,1],lmks1[,1],"-")
    Sdy        <- outer(lmks1[,2],lmks1[,2],"-")
    Sdz        <- outer(lmks1[,3],lmks1[,3],"-")
    Sphi       <- sqrt(Sdx * Sdx + Sdy * Sdy + Sdz * Sdz)
    ONE        <- matrix(1, k, 1)
    ZERO       <- matrix(0, 4, 4)
    Q          <- rbind(Sphi, t(ONE), t(lmks1))
    O          <- cbind(ONE, lmks1)
    O          <- as.matrix(O)
    U          <- rbind(O, ZERO)
    GAMMA      <- cbind(Q, U)
    Gi         <- MASS::ginv(GAMMA)
    B          <- kronecker(diag(1,3),Gi[1:k,1:k])
    Be         <- bendingenergy.face3d(B,lmks2)
    results    <- list(L = GAMMA, B = B, Be = Be)
return(results)
}

