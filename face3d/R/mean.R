mean.face3d <- function(shape, values){
	
mesh.m    <- shape$vertices - mean(shape$vertices)
mesh.pc   <- prcomp(mesh.m)
mesh.pcf  <- mesh.pc$rotation
mesh.r    <- mesh.m %*% mesh.pcf
mesh.2d   <- mesh.r[,1:2]

nbrk <- 30
# cc <- colorRampPalette(c("blue","red"))(30)
# plot(mesh.2d, col = cc[cut(values, nbrk)], pch = 16)
	
mean.x <-(sum(mesh.2d[,1]*values))  / sum(values)
mean.y <-(sum(mesh.2d[,2]*values))  / sum(values)
	
# points(mean.x, mean.y, pch=3)	

X             <- mesh.2d
Xnew          <- matrix(c(mean.x,mesh.2d[4,1],mean.y,mesh.2d[4,2]), ncol=2)
interp.x      <- interp.barycentric(X, f = shape$vertices[ ,1], Xnew)$fnew
interp.y      <- interp.barycentric(X, f = shape$vertices[ ,2], Xnew)$fnew
interp.z      <- interp.barycentric(X, f = shape$vertices[ ,3], Xnew)$fnew
mean          <- cbind(interp.x[1], interp.y[1], interp.z[1])
            
	invisible(mean)
}
