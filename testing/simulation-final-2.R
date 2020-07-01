setwd("/Users/libertyvittert/ownCloud/Shared/Face3D_0.1-1/Face3D/R")
for (file in list.files()) source(file)
library(rgl)
library(fields)
library(MASS)
"edist.face3d"  <- function(x, x0) sqrt(rowSums(sweep(x, 2, x0)^2))

##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################


# # 
# x <- y <- seq(-1,1, length=100)
# grid.x <- rep(x, each=length(x))
# grid.y <- rep(y, length(x))

# grid.z <- (grid.x)^2 +grid.x^2

# points3d(cbind(grid.x, grid.y, grid.z))

# library(sm)

# month <- rep(1:12, 10)
# year <- rep(1:10, each=12)
# ddate <- (year+1900)+month/12

# response <- 0.05*(1:length(ddate))

# plot(ddate, response+sin(2*pi*month/12) + cos(2*pi*month/12))

# rs <- response+sin(2*pi*month/12) + cos(2*pi*month/12)
# rs2 <- rs
# rs2[month==6]<- sin(2*pi*month/12)[1:10]
# test <- sm.regression(cbind(month,ddate), rs2, display="rgl")
# ?sm.regression

# plot(rs[month==6]^2)


# plot(sin(2*pi*month/12)[1:10])



ww <-1
p <-5
t <-1


#NOISE and STEEP -mesh size
lm                      <- 0 
cm                      <- 0
mesh.size               <- c(.05,.06,.07,.08,.09,.10)
nsim                    <- seq(0,.01,length=10)
nsim2                   <- seq(0,.5, length=11)
nsim2                   <- nsim2[-1]

final                   <- array(0, dim =c(length(nsim2),20,length(nsim)))
final.average.distances <- rep(NA, length(nsim2))
total.average.distances <- array(0, dim=c(length(nsim2),1,length(nsim)))

for (ww in 1: length(mesh.size)){
ms            <- mesh.size[ww]	

for (totalsim in 1:5){
	
for (p in 1:length(nsim)){
noise           <- nsim[p]     

for(t in 1:length(nsim2)){
steep           <- nsim2[t]
variance        <- 1
max.grid        <- 1
min.grid        <- -1
size.by         <- ms
distance        <- .5
perp.dist.bound <- 1
x               <- seq(min.grid, max.grid, size.by)
grid.x          <- rep(x, each= length(x))
grid.y          <- rep(x, length(x))
grid.z          <- -steep*grid.x^2  


addNoise <- function(mtx) {
  if (!is.matrix(mtx)) mtx <- matrix(mtx, byrow = TRUE, nrow = 1)
  random.stuff  <- matrix(runif(prod(dim(mtx)), min = -noise, max = noise), 
                 nrow =    dim(mtx)[1])
  random.stuff + mtx
}
grid.z          <- as.vector(addNoise(mtx=grid.z))
grid            <- cbind(grid.x, grid.y, grid.z)
fpt1            <- c(lm,-1,-cm)
fpt2            <- c(lm,1,cm)
fixed.y         <- seq(-1,1,.1) 
fixed.x         <- rep(0, length(fixed.y))
fixed.z         <- rep(NA, length(fixed.y))
for (aa in 1: length(fixed.y)) fixed.z[aa] <- -steep*fixed.x[aa]^2 + cm*fixed.y[aa]^3 
rnames          <- matrix(nrow=length(x), ncol=length(x))
for(i in 1:length(x)) rnames[,i ] <- paste(i, 1:length(x), sep= " ")
rownames(grid)                    <- as.vector(rnames)                        
triples         <- NULL
pts             <- length(x)
indices         <- c(1:dim(grid)[1])
paths           <- length(indices)/pts
trips           <- matrix(indices, ncol=paths)
    for (i in 1:(paths-1)){
	    for (j in 1:(pts-1)){
      	   triples <- c(triples, 
      	                      c(trips[j,i], trips[j,i+1], trips[j+1,i+1] , 
		                      trips[j+1, i+1], trips[j+1,i],trips[j,i]))}
    }  
lmks            <- rbind(fpt1,fpt2)
shape           <- list(coords=grid, triples=triples, lmks=rbind(fpt1,fpt2), 
                fixed.ridge = cbind(fixed.x,fixed.y,fixed.z))
shape$lmks[1,]  <-closest.face3d(shape$lmks[1,],shape)$point
shape$lmks[2,]  <-closest.face3d(shape$lmks[2,],shape)$point
shape                          <- index.face3d(shape, distance = distance)
ind                            <- (sign(shape$shape.index) == 1)
shape$shape.index[!ind]        <- 0
shape$kappa1[!ind]             <- 0
shape$kappa2[!ind]             <- 0
values                         <- pmax(-shape$kappa2, 0)


curve               <- shape$fixed.ridge
ccurve              <- closestcurve.face3d(shape, curve)
lambda              <- cumsum(c(0,sqrt(apply((diff(curve))^2, 1, sum))))
xcoord              <- lambda[ccurve$closest.curvept] 
ycoord              <- ccurve$closest.distance
ind                 <- (abs(ycoord) < 10) & (xcoord  > 0.005 * max(xcoord )) & (xcoord  < 0.995 * max(xcoord ))
yc                  <- c(0, ycoord[ind], 0)
xc                  <- c(0, xcoord[ind], max(xcoord))
values              <- c(0, values[ind], 0)
ridge               <- ridge2d.face3d(xc, yc, values, graphics = FALSE)

                    
X                    <- cbind(xcoord,ycoord)
Xnew                 <- cbind(ridge$x,ridge$y)
interp.x             <- interp.barycentric(X, f = shape$coords[ , 1], Xnew)$fnew
interp.y             <- interp.barycentric(X, f = shape$coords[ , 2], Xnew)$fnew
interp.z             <- interp.barycentric(X, f = shape$coords[ , 3], Xnew)$fnew
smooth.path          <- rbind(fpt1, cbind(interp.x, interp.y, interp.z), fpt2)  
        
        
smooth.path          <- resample.face3d(smooth.path, n= dim(shape$fixed.ridge)[1])
shape$fixed.ridge    <- resample.face3d(shape$fixed.ridge, n= dim(shape$fixed.ridge)[1])
distances.btw.points <- sqrt(apply(((shape$fixed.ridge - smooth.path)^2),1 ,sum))
average.distance     <- sum(distances.btw.points)/length(distances.btw.points)
final.average.distances[t]   <- average.distance
#spheres3d(smooth.path, col="red", radius=.05)
#spheres3d(shape$fixed.ridge, radius=.05)
}
total.average.distances[,,p] <- final.average.distances
}
final[ ,totalsim, ]           <- total.average.distances   
print(totalsim)        
}
save(final, nsim, nsim2, file=paste("ms",mesh.size[ww],"_noise_steep.dmp", sep=""))
print("Done with MS round")
}


##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
#NOISE and STEEP -lmks move

cm                      <- 0
ms                      <- .05
lmk.move                <- c(.01,.02,.04,.06,.08,.10)
nsim                    <- seq(0,.01,length=10)
nsim2                   <- seq(0,.5, length=11)
nsim2                   <- nsim2[-1]

final                   <- array(0, dim =c(length(nsim2),20,length(nsim)))
final.average.distances <- rep(NA, length(nsim2))
total.average.distances <- array(0, dim=c(length(nsim2),1,length(nsim)))

for (ww in 1: length(lmk.move)){
lm            <- lmk.move[ww]	

for (totalsim in 1:5){
	
for (p in 1:length(nsim)){
noise           <- nsim[p]     

for(t in 1:length(nsim2)){
steep           <- nsim2[t]
variance        <- 1
max.grid        <- 1
min.grid        <- -1
size.by         <- ms
distance        <- .5
perp.dist.bound <- 1
x               <- seq(min.grid, max.grid, size.by)
grid.x          <- rep(x, each= length(x))
grid.y          <- rep(x, length(x))
grid.z          <- -steep*grid.x^2  + cm*x^3
addNoise <- function(mtx) {
  if (!is.matrix(mtx)) mtx <- matrix(mtx, byrow = TRUE, nrow = 1)
  random.stuff  <- matrix(runif(prod(dim(mtx)), min = -noise, max = noise), 
                 nrow =    dim(mtx)[1])
  random.stuff + mtx
}
grid.z          <- as.vector(addNoise(mtx=grid.z))
grid            <- cbind(grid.x, grid.y, grid.z)
fpt1            <- c(lm,-1,-cm)
fpt2            <- c(lm,1,cm)
fixed.y         <- seq(-1,1,.1) 
fixed.x         <- rep(0, length(fixed.y))
fixed.z         <- rep(NA, length(fixed.y))
for (aa in 1: length(fixed.y)) fixed.z[aa] <- -steep*fixed.x[aa]^2 + cm*fixed.y[aa]^3 
rnames          <- matrix(nrow=length(x), ncol=length(x))
for(i in 1:length(x)) rnames[,i ] <- paste(i, 1:length(x), sep= " ")
rownames(grid)                    <- as.vector(rnames)                        
triples         <- NULL
pts             <- length(x)
indices         <- c(1:dim(grid)[1])
paths           <- length(indices)/pts
trips           <- matrix(indices, ncol=paths)
    for (i in 1:(paths-1)){
	    for (j in 1:(pts-1)){
      	   triples <- c(triples, 
      	                      c(trips[j,i], trips[j,i+1], trips[j+1,i+1] , 
		                      trips[j+1, i+1], trips[j+1,i],trips[j,i]))}
    }  
lmks            <- rbind(fpt1,fpt2)
shape           <- list(coords=grid, triples=triples, lmks=rbind(fpt1,fpt2), 
                fixed.ridge = cbind(fixed.x,fixed.y,fixed.z))
shape$lmks[1,]  <-closest.face3d(shape$lmks[1,],shape)$point
shape$lmks[2,]  <-closest.face3d(shape$lmks[2,],shape)$point
shape                          <- index.face3d(shape, distance = distance)
ind                            <- (sign(shape$shape.index) == 1)
shape$shape.index[!ind]        <- 0
shape$kappa1[!ind]             <- 0
shape$kappa2[!ind]             <- 0
values                         <- pmax(-shape$kappa2, 0)


curve               <- shape$fixed.ridge
ccurve              <- closestcurve.face3d(shape, curve)
lambda              <- cumsum(c(0,sqrt(apply((diff(curve))^2, 1, sum))))
xcoord              <- lambda[ccurve$closest.curvept] 
ycoord              <- ccurve$closest.distance
ind                 <- (abs(ycoord) < 10) & (xcoord  > 0.005 * max(xcoord )) & (xcoord  < 0.995 * max(xcoord ))
yc                  <- c(0, ycoord[ind], 0)
xc                  <- c(0, xcoord[ind], max(xcoord))
values              <- c(0, values[ind], 0)
ridge               <- ridge2d.face3d(xc, yc, values, graphics = FALSE)

                    
X                    <- cbind(xcoord,ycoord)
Xnew                 <- cbind(ridge$x,ridge$y)
interp.x             <- interp.barycentric(X, f = shape$coords[ , 1], Xnew)$fnew
interp.y             <- interp.barycentric(X, f = shape$coords[ , 2], Xnew)$fnew
interp.z             <- interp.barycentric(X, f = shape$coords[ , 3], Xnew)$fnew
smooth.path          <- rbind(fpt1, cbind(interp.x, interp.y, interp.z), fpt2)  
        
smooth.path          <- resample.face3d(smooth.path, n= dim(shape$fixed.ridge)[1])
shape$fixed.ridge    <- resample.face3d(shape$fixed.ridge, n= dim(shape$fixed.ridge)[1])
distances.btw.points <- sqrt(apply(((shape$fixed.ridge - smooth.path)^2),1 ,sum))
average.distance     <- sum(distances.btw.points)/length(distances.btw.points)
final.average.distances[t]   <- average.distance
}
total.average.distances[,,p] <- final.average.distances
}
final[ ,totalsim, ]           <- total.average.distances   
print(totalsim)        
}
save(final, nsim, nsim2, file=paste("lm",lmk.move[ww],"_noise_steep.dmp", sep=""))
print("Done with LM round")
}



##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
#STEEP and CM -mesh size

lm                      <- 0 
mesh.size               <- c(.05,.06,.07,.08,.09,.10)
nsim                    <- seq(0,.5,length=11)
nsim                    <- nsim[-1]
nsim2                   <- seq(0,1, length=10)

final                   <- array(0, dim =c(length(nsim2),20,length(nsim)))
final.average.distances <- rep(NA, length(nsim2))
total.average.distances <- array(0, dim=c(length(nsim2),1,length(nsim)))

for (ww in 1: length(mesh.size)){
ms            <- mesh.size[ww]	
print(ww)
for (totalsim in 1:5){

for (p in 1:length(nsim)){
steep           <- nsim[p]  
print(p)   

for(t in 1:length(nsim2)){

cm              <- nsim2[t]
print(t)
variance        <- 1
max.grid        <- 1
min.grid        <- -1
size.by         <- ms
distance        <- .5
perp.dist.bound <- 1
x               <- seq(min.grid, max.grid, size.by)
grid.x          <- rep(x, each= length(x))
grid.y          <- rep(x, length(x))
grid.z          <- -steep*grid.x^2  + cm*x^4
grid            <- cbind(grid.x, grid.y, grid.z)
fpt1            <- c(lm,-1,-cm)
fpt2            <- c(lm,1,cm)
fixed.y         <- seq(-1,1,.1) 
fixed.x         <- rep(0, length(fixed.y))
fixed.z         <- rep(NA, length(fixed.y))
for (aa in 1: length(fixed.y)) fixed.z[aa] <- -steep*fixed.x[aa]^2 + cm*fixed.y[aa]^3 
rnames          <- matrix(nrow=length(x), ncol=length(x))
for(i in 1:length(x)) rnames[,i ] <- paste(i, 1:length(x), sep= " ")
rownames(grid)                    <- as.vector(rnames)                        
triples         <- NULL
pts             <- length(x)
indices         <- c(1:dim(grid)[1])
paths           <- length(indices)/pts
trips           <- matrix(indices, ncol=paths)
    for (i in 1:(paths-1)){
	    for (j in 1:(pts-1)){
      	   triples <- c(triples, 
      	                      c(trips[j,i], trips[j,i+1], trips[j+1,i+1] , 
		                      trips[j+1, i+1], trips[j+1,i],trips[j,i]))}
    }  
lmks            <- rbind(fpt1,fpt2)
shape           <- list(coords=grid, triples=triples, lmks=rbind(fpt1,fpt2), 
                fixed.ridge = cbind(fixed.x,fixed.y,fixed.z))
shape$lmks[1,]  <-closest.face3d(shape$lmks[1,],shape)$point
shape$lmks[2,]  <-closest.face3d(shape$lmks[2,],shape)$point
shape                          <- index.face3d(shape, distance = distance)
ind                            <- (sign(shape$shape.index) == 1)
shape$shape.index[!ind]        <- 0
shape$kappa1[!ind]             <- 0
shape$kappa2[!ind]             <- 0
values                         <- pmax(-shape$kappa2, 0)


curve               <- shape$fixed.ridge
ccurve              <- closestcurve.face3d(shape, curve)
lambda              <- cumsum(c(0,sqrt(apply((diff(curve))^2, 1, sum))))
xcoord              <- lambda[ccurve$closest.curvept] 
ycoord              <- ccurve$closest.distance
ind                 <- (abs(ycoord) < 10) & (xcoord  > 0.005 * max(xcoord )) & (xcoord  < 0.995 * max(xcoord ))
yc                  <- c(0, ycoord[ind], 0)
xc                  <- c(0, xcoord[ind], max(xcoord))
values              <- c(0, values[ind], 0)
ridge               <- ridge2d.face3d(xc, yc, values, graphics = FALSE)

                    
X                    <- cbind(xcoord,ycoord)
Xnew                 <- cbind(ridge$x,ridge$y)
interp.x             <- interp.barycentric(X, f = shape$coords[ , 1], Xnew)$fnew
interp.y             <- interp.barycentric(X, f = shape$coords[ , 2], Xnew)$fnew
interp.z             <- interp.barycentric(X, f = shape$coords[ , 3], Xnew)$fnew
smooth.path          <- rbind(fpt1, cbind(interp.x, interp.y, interp.z), fpt2)  
        
smooth.path          <- resample.face3d(smooth.path, n= dim(shape$fixed.ridge)[1])
shape$fixed.ridge    <- resample.face3d(shape$fixed.ridge, n= dim(shape$fixed.ridge)[1])
distances.btw.points <- sqrt(apply(((shape$fixed.ridge - smooth.path)^2),1 ,sum))
average.distance     <- sum(distances.btw.points)/length(distances.btw.points)
final.average.distances[t]   <- average.distance

}
total.average.distances[,,p] <- final.average.distances

}
final[ ,totalsim, ]           <- total.average.distances 
  
print(totalsim)        
}
save(final, nsim, nsim2, file=paste("ms",mesh.size[ww],"_steep_cm.dmp", sep=""))
print("Done with MS round")
}

























##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
#STEEP and CM -lmksmove

ms       <- 0.05
lmk.move <- c(.01,.02,.04,.06,.08,.10)
nsim     <- seq(0,.5,length=11)
nsim     <- nsim[-1]
nsim2    <- seq(0,1, length=10)

final                   <- array(0, dim =c(length(nsim2),20,length(nsim)))
final.average.distances <- rep(NA, length(nsim2))
total.average.distances <- array(0, dim=c(length(nsim2),1,length(nsim)))

for (ww in 1: length(lmk.move)){
lm            <- lmk.move[ww]	

for (totalsim in 1:5){
	
for (p in 1:length(nsim)){
steep           <- nsim[p]     

for(t in 1:length(nsim2)){
cm              <- nsim2[t] 
variance        <- 1
max.grid        <- 1
min.grid        <- -1
size.by         <- ms
distance        <- .5
perp.dist.bound <- 1
x               <- seq(min.grid, max.grid, size.by)
grid.x          <- rep(x, each= length(x))
grid.y          <- rep(x, length(x))
grid.z          <- -steep*grid.x^2  + cm*x^3
grid            <- cbind(grid.x, grid.y, grid.z)
fpt1            <- c(lm,-1,-cm)
fpt2            <- c(lm,1,cm)
fixed.y         <- seq(-1,1,.1) 
fixed.x         <- rep(0, length(fixed.y))
fixed.z         <- rep(NA, length(fixed.y))
for (aa in 1: length(fixed.y)) fixed.z[aa] <- -steep*fixed.x[aa]^2 + cm*fixed.y[aa]^3 
rnames          <- matrix(nrow=length(x), ncol=length(x))
for(i in 1:length(x)) rnames[,i ] <- paste(i, 1:length(x), sep= " ")
rownames(grid)                    <- as.vector(rnames)                        
triples         <- NULL
pts             <- length(x)
indices         <- c(1:dim(grid)[1])
paths           <- length(indices)/pts
trips           <- matrix(indices, ncol=paths)
    for (i in 1:(paths-1)){
	    for (j in 1:(pts-1)){
      	   triples <- c(triples, 
      	                      c(trips[j,i], trips[j,i+1], trips[j+1,i+1] , 
		                      trips[j+1, i+1], trips[j+1,i],trips[j,i]))}
    }  
lmks            <- rbind(fpt1,fpt2)
shape           <- list(coords=grid, triples=triples, lmks=rbind(fpt1,fpt2), 
                fixed.ridge = cbind(fixed.x,fixed.y,fixed.z))
shape$lmks[1,]  <-closest.face3d(shape$lmks[1,],shape)$point
shape$lmks[2,]  <-closest.face3d(shape$lmks[2,],shape)$point
shape                          <- index.face3d(shape, distance = distance)
ind                            <- (sign(shape$shape.index) == 1)
shape$shape.index[!ind]        <- 0
shape$kappa1[!ind]             <- 0
shape$kappa2[!ind]             <- 0
values                         <- pmax(-shape$kappa2, 0)


curve               <- shape$fixed.ridge
ccurve              <- closestcurve.face3d(shape, curve)
lambda              <- cumsum(c(0,sqrt(apply((diff(curve))^2, 1, sum))))
xcoord              <- lambda[ccurve$closest.curvept] 
ycoord              <- ccurve$closest.distance
ind                 <- (abs(ycoord) < 10) & (xcoord  > 0.005 * max(xcoord )) & (xcoord  < 0.995 * max(xcoord ))
yc                  <- c(0, ycoord[ind], 0)
xc                  <- c(0, xcoord[ind], max(xcoord))
values              <- c(0, values[ind], 0)
ridge               <- ridge2d.face3d(xc, yc, values, graphics = FALSE)

                    
X                    <- cbind(xcoord,ycoord)
Xnew                 <- cbind(ridge$x,ridge$y)
interp.x             <- interp.barycentric(X, f = shape$coords[ , 1], Xnew)$fnew
interp.y             <- interp.barycentric(X, f = shape$coords[ , 2], Xnew)$fnew
interp.z             <- interp.barycentric(X, f = shape$coords[ , 3], Xnew)$fnew
smooth.path          <- rbind(fpt1, cbind(interp.x, interp.y, interp.z), fpt2)  
        
smooth.path          <- resample.face3d(smooth.path, n= dim(shape$fixed.ridge)[1])
shape$fixed.ridge    <- resample.face3d(shape$fixed.ridge, n= dim(shape$fixed.ridge)[1])
distances.btw.points <- sqrt(apply(((shape$fixed.ridge - smooth.path)^2),1 ,sum))
average.distance     <- sum(distances.btw.points)/length(distances.btw.points)
final.average.distances[t]   <- average.distance
}
total.average.distances[,,p] <- final.average.distances
}
final[ ,totalsim, ]           <- total.average.distances   
print(totalsim)        
}
save(final, nsim, nsim2, file=paste("lm",lmk.move[ww],"_steep_cm.dmp", sep=""))
print("Done with lmk round")
}



##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
#NOISE and CM -mesh size
lm                      <- 0 
steep                   <-.05
mesh.size               <- c(.05,.06,.07,.08,.09,.10)
nsim                    <- seq(0,.01,length=10)
nsim2                   <- seq(0,1, length=11)
nsim2                   <- nsim2[-1]

final                   <- array(0, dim =c(length(nsim2),20,length(nsim)))
final.average.distances <- rep(NA, length(nsim2))
total.average.distances <- array(0, dim=c(length(nsim2),1,length(nsim)))

for (ww in 1: length(mesh.size)){
ms            <- mesh.size[ww]	

for (totalsim in 1:5){
	
for (p in 1:length(nsim)){
noise           <- nsim[p]     

for(t in 1:length(nsim2)){
cm              <- nsim2[t]
variance        <- 1
max.grid        <- 1
min.grid        <- -1
size.by         <- ms
distance        <- .5
perp.dist.bound <- 1
x               <- seq(min.grid, max.grid, size.by)
grid.x          <- rep(x, each= length(x))
grid.y          <- rep(x, length(x))
grid.z          <- -steep*grid.x^2  + cm*x^3
addNoise <- function(mtx) {
  if (!is.matrix(mtx)) mtx <- matrix(mtx, byrow = TRUE, nrow = 1)
  random.stuff  <- matrix(runif(prod(dim(mtx)), min = -noise, max = noise), 
                 nrow =    dim(mtx)[1])
  random.stuff + mtx
}
grid.z          <- as.vector(addNoise(mtx=grid.z))
grid            <- cbind(grid.x, grid.y, grid.z)
fpt1            <- c(lm,-1,-cm)
fpt2            <- c(lm,1,cm)
fixed.y         <- seq(-1,1,.1) 
fixed.x         <- rep(0, length(fixed.y))
fixed.z         <- rep(NA, length(fixed.y))
for (aa in 1: length(fixed.y)) fixed.z[aa] <- -steep*fixed.x[aa]^2 + cm*fixed.y[aa]^3 
rnames          <- matrix(nrow=length(x), ncol=length(x))
for(i in 1:length(x)) rnames[,i ] <- paste(i, 1:length(x), sep= " ")
rownames(grid)                    <- as.vector(rnames)                        
triples         <- NULL
pts             <- length(x)
indices         <- c(1:dim(grid)[1])
paths           <- length(indices)/pts
trips           <- matrix(indices, ncol=paths)
    for (i in 1:(paths-1)){
	    for (j in 1:(pts-1)){
      	   triples <- c(triples, 
      	                      c(trips[j,i], trips[j,i+1], trips[j+1,i+1] , 
		                      trips[j+1, i+1], trips[j+1,i],trips[j,i]))}
    }  
lmks            <- rbind(fpt1,fpt2)
shape           <- list(coords=grid, triples=triples, lmks=rbind(fpt1,fpt2), 
                fixed.ridge = cbind(fixed.x,fixed.y,fixed.z))
shape$lmks[1,]  <-closest.face3d(shape$lmks[1,],shape)$point
shape$lmks[2,]  <-closest.face3d(shape$lmks[2,],shape)$point
shape                          <- index.face3d(shape, distance = distance)
ind                            <- (sign(shape$shape.index) == 1)
shape$shape.index[!ind]        <- 0
shape$kappa1[!ind]             <- 0
shape$kappa2[!ind]             <- 0
values                         <- pmax(-shape$kappa2, 0)


curve               <- shape$fixed.ridge
ccurve              <- closestcurve.face3d(shape, curve)
lambda              <- cumsum(c(0,sqrt(apply((diff(curve))^2, 1, sum))))
xcoord              <- lambda[ccurve$closest.curvept] 
ycoord              <- ccurve$closest.distance
ind                 <- (abs(ycoord) < 10) & (xcoord  > 0.005 * max(xcoord )) & (xcoord  < 0.995 * max(xcoord ))
yc                  <- c(0, ycoord[ind], 0)
xc                  <- c(0, xcoord[ind], max(xcoord))
values              <- c(0, values[ind], 0)
ridge               <- ridge2d.face3d(xc, yc, values, graphics = FALSE)

                    
X                    <- cbind(xcoord,ycoord)
Xnew                 <- cbind(ridge$x,ridge$y)
interp.x             <- interp.barycentric(X, f = shape$coords[ , 1], Xnew)$fnew
interp.y             <- interp.barycentric(X, f = shape$coords[ , 2], Xnew)$fnew
interp.z             <- interp.barycentric(X, f = shape$coords[ , 3], Xnew)$fnew
smooth.path          <- rbind(fpt1, cbind(interp.x, interp.y, interp.z), fpt2)  
        
 
        
smooth.path          <- resample.face3d(smooth.path, n= dim(shape$fixed.ridge)[1])
shape$fixed.ridge    <- resample.face3d(shape$fixed.ridge, n= dim(shape$fixed.ridge)[1])
distances.btw.points <- sqrt(apply(((shape$fixed.ridge - smooth.path)^2),1 ,sum))
average.distance     <- sum(distances.btw.points)/length(distances.btw.points)
final.average.distances[t]   <- average.distance
#spheres3d(smooth.path, col="red", radius=.05)
#spheres3d(shape$fixed.ridge, radius=.05)
}
total.average.distances[,,p] <- final.average.distances
}
final[ ,totalsim, ]           <- total.average.distances   
print(totalsim)        
}
save(final, nsim, nsim2, file=paste("ms",mesh.size[ww],"_noise_cm.dmp", sep=""))
print("Done with MS round")
}

##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
#NOISE and CM -lmk move
ms                      <- .05
steep                   <- .05
lmk.move                <- c(.01,.02,.04,.06,.08,.10)
nsim                    <- seq(0,.01,length=10)
nsim2                   <- seq(0,1, length=11)
nsim2                   <- nsim2[-1]

final                   <- array(0, dim =c(length(nsim2),20,length(nsim)))
final.average.distances <- rep(NA, length(nsim2))
total.average.distances <- array(0, dim=c(length(nsim2),1,length(nsim)))

for (ww in 1: length(lmk.move)){
lm            <- lmk.move[ww]	

for (totalsim in 1:5){
	
for (p in 1:length(nsim)){
noise           <- nsim[p]     

for(t in 1:length(nsim2)){
cm              <- nsim2[t]
variance        <- 1
max.grid        <- 1
min.grid        <- -1
size.by         <- ms
distance        <- .5
perp.dist.bound <- 1
x               <- seq(min.grid, max.grid, size.by)
grid.x          <- rep(x, each= length(x))
grid.y          <- rep(x, length(x))
grid.z          <- -steep*grid.x^2  + cm*x^3
addNoise <- function(mtx) {
  if (!is.matrix(mtx)) mtx <- matrix(mtx, byrow = TRUE, nrow = 1)
  random.stuff  <- matrix(runif(prod(dim(mtx)), min = -noise, max = noise), 
                 nrow =    dim(mtx)[1])
  random.stuff + mtx
}
grid.z          <- as.vector(addNoise(mtx=grid.z))
grid            <- cbind(grid.x, grid.y, grid.z)
fpt1            <- c(lm,-1,-cm)
fpt2            <- c(lm,1,cm)
fixed.y         <- seq(-1,1,.1) 
fixed.x         <- rep(0, length(fixed.y))
fixed.z         <- rep(NA, length(fixed.y))
for (aa in 1: length(fixed.y)) fixed.z[aa] <- -steep*fixed.x[aa]^2 + cm*fixed.y[aa]^3 
rnames          <- matrix(nrow=length(x), ncol=length(x))
for(i in 1:length(x)) rnames[,i ] <- paste(i, 1:length(x), sep= " ")
rownames(grid)                    <- as.vector(rnames)                        
triples         <- NULL
pts             <- length(x)
indices         <- c(1:dim(grid)[1])
paths           <- length(indices)/pts
trips           <- matrix(indices, ncol=paths)
    for (i in 1:(paths-1)){
	    for (j in 1:(pts-1)){
      	   triples <- c(triples, 
      	                      c(trips[j,i], trips[j,i+1], trips[j+1,i+1] , 
		                      trips[j+1, i+1], trips[j+1,i],trips[j,i]))}
    }  
lmks            <- rbind(fpt1,fpt2)
shape           <- list(coords=grid, triples=triples, lmks=rbind(fpt1,fpt2), 
                fixed.ridge = cbind(fixed.x,fixed.y,fixed.z))
shape$lmks[1,]  <-closest.face3d(shape$lmks[1,],shape)$point
shape$lmks[2,]  <-closest.face3d(shape$lmks[2,],shape)$point
shape                          <- index.face3d(shape, distance = distance)
ind                            <- (sign(shape$shape.index) == 1)
shape$shape.index[!ind]        <- 0
shape$kappa1[!ind]             <- 0
shape$kappa2[!ind]             <- 0
values                         <- pmax(-shape$kappa2, 0)


curve               <- shape$fixed.ridge
ccurve              <- closestcurve.face3d(shape, curve)
lambda              <- cumsum(c(0,sqrt(apply((diff(curve))^2, 1, sum))))
xcoord              <- lambda[ccurve$closest.curvept] 
ycoord              <- ccurve$closest.distance
ind                 <- (abs(ycoord) < 10) & (xcoord  > 0.005 * max(xcoord )) & (xcoord  < 0.995 * max(xcoord ))
yc                  <- c(0, ycoord[ind], 0)
xc                  <- c(0, xcoord[ind], max(xcoord))
values              <- c(0, values[ind], 0)
ridge               <- ridge2d.face3d(xc, yc, values, graphics = FALSE)

                    
X                    <- cbind(xcoord,ycoord)
Xnew                 <- cbind(ridge$x,ridge$y)
interp.x             <- interp.barycentric(X, f = shape$coords[ , 1], Xnew)$fnew
interp.y             <- interp.barycentric(X, f = shape$coords[ , 2], Xnew)$fnew
interp.z             <- interp.barycentric(X, f = shape$coords[ , 3], Xnew)$fnew
smooth.path          <- rbind(fpt1, cbind(interp.x, interp.y, interp.z), fpt2)  
        

        
smooth.path          <- resample.face3d(smooth.path, n= dim(shape$fixed.ridge)[1])
shape$fixed.ridge    <- resample.face3d(shape$fixed.ridge, n= dim(shape$fixed.ridge)[1])
distances.btw.points <- sqrt(apply(((shape$fixed.ridge - smooth.path)^2),1 ,sum))
average.distance     <- sum(distances.btw.points)/length(distances.btw.points)
final.average.distances[t]   <- average.distance
#spheres3d(smooth.path, col="red", radius=.05)
#spheres3d(shape$fixed.ridge, radius=.05)
}
total.average.distances[,,p] <- final.average.distances
}
final[ ,totalsim, ]           <- total.average.distances   
print(totalsim)        
}
save(final, nsim, nsim2, file=paste("lm",lmk.move[ww],"_noise_cm.dmp", sep=""))
print("Done with LM round")
}

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
# # #To Input into Latex
library(xtable)


final.means <- array(0, dim=c(10,1,10))
for (i in 1:10){
final.means[,,i] <- rowMeans(final[,,i])
}

dd            <- t(data.frame(round(final.means,3)))

#nsim          <- round(seq(0,.05,length=10),3)   noise

 nsim          <- round(seq(0,2.5,length=11),3)  #steep
 nsim <- nsim[-1]
 
rownames(dd)  <- nsim
colnames(dd)  <- NULL
print(xtable(dd ,type="latex", file="test", floating="TRUE", digits=3, include.rownames ="FALSE", include.colnames = "FALSE"))




