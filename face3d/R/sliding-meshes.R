"slidingmesh.face3d" <- function(shape1,shape2, n = 100, type = "resample")  {
   # choose the particular curve and whether it is slid or not
   
   # if (class(shape1) == "face3d"){ 
   # 	shape1 <- shape1$mesh
   # 	shape2 <- shape2$mesh
   # } 
 


 mesh.names <-          c(paste("mid-face right",       1:20,""),
                          paste("mid-face left",        1:20, ""),
                          paste("upper mid face right", 1:20, ""),
                          paste("upper mid face left",  1:20, ""),
                          paste("lower face right",     1:28, ""),
                          paste("nose left",            1:7, ""),
                          paste("nose right",           1:7, ""),
                          #paste( "upper face right 2",  1:10, ""),
                          #paste( "upper face left 2",   1:10, ""),
                          paste( "lower face left",     1:28, ""))

 mesh.names.notslide <- c(paste("upper lip left",       1:7, ""),
                          paste("lower lip left",       1:7, ""),
                          paste("upper lip right",      1:7, ""),
                          paste("lower lip right",      1:7, ""),
                          paste("philtrum right",       1:9, ""),
                          paste("philtrum left",        1:9, ""),
                          paste( "upper face right 3",  1:7, ""),
                          paste( "upper face left 3",   1:7, ""),
                          paste( "upper face right 1",  1:12, ""),
                          paste( "upper face left 1",   1:12, ""),
                          paste( "upper face right 2",  1:10, ""),
                          paste( "upper face left 2",   1:10, ""))

  n.notslide          <- c(rep(3,7), rep(3,7), rep(3,7), rep(3,7),
                          rep(6,9), rep(6,9), rep(3,7), rep(3,7),
                          rep(3,12),rep(3,12),rep(1,10), rep(1,10))

 n                   <- c(rep(6,20), rep(6,20), rep(5,20), rep(5,20),
                          rep(5,28), rep(5,7), rep(5,7),rep(5,28))

    



if(type == "resample/slide"){
  n <- n * 20
  curves2.optim           <- matrix(0, ncol = 3, nrow = length(shape1$pre.prep.mesh))
  rownames(curves2.optim) <- shape1$pre.prep.mesh
  
  
   for (i in 1:length(mesh.names)) {
    
       ind.i                 <- grep(mesh.names[i],dimnames(curves2.optim)[[1]])
       ind.shape2            <- grep(mesh.names[i],dimnames(shape2$mesh)[[1]])
       curves2.optim[ind.i,] <- sliding.face3d(shape2[ind.i,],shape2[ind.shape2,],
   							                   n=n[i], resample=TRUE)
       print(i)
       }
   
   for (i in 1:length(mesh.names.notslide)) {
       ind.i                 <- grep(mesh.names.notslide[i],dimnames(curves2.optim)[[1]])
       ind.shape2            <- grep(mesh.names.notslide[i],dimnames(shape2$mesh)[[1]])
       curves2.optim[ind.i,] <- resample.face3d(shape2$mesh[ind.shape2,], n=n.notslide[i])
   }
  
  
  
  
}else if(type == "resample"){
  curves2.optim           <- matrix(0, ncol = 3, nrow = length(shape1$pre.prep.mesh))
  rownames(curves2.optim) <- shape1$pre.prep.mesh
  #pre.prep.mesh FIX HEREREEEEEEEE using 1156 names
  
  for (i in 1:length(mesh.names)) {
    ind.i                 <- grep(mesh.names[i], dimnames(curves2.optim)[[1]])
    ind.shape2            <- grep(mesh.names[i],dimnames(shape2$mesh)[[1]])
    curves2.optim[ind.i,] <- resample.face3d(shape2$mesh[ind.shape2,], n=n[i])
       }
  
  for (i in 1:length(mesh.names.notslide)) {
    ind.i                 <- grep(mesh.names.notslide[i],dimnames(curves2.optim)[[1]])
    ind.shape2            <- grep(mesh.names.notslide[i],dimnames(shape2$mesh)[[1]])
    curves2.optim[ind.i,] <- resample.face3d(shape2$mesh[ind.shape2,], n=n.notslide[i])
       }
  
  
}

   return(curves2.optim)
}
  
  
  
  
  
  
  
   #resampling curves to decided amount for the big ones and decided amount for the mesh
     
  # for (i in 1:length(curve.names)) {
      # ind                            <- grep(curve.names[i], rownames(face$curves))
      # x                              <- face$curves[ind,]
      # arc.length                     <- cumsum(sqrt(apply((x-rbind(x[1,],x[-(dim(x)[1]),]))^2,1,sum)))
      # x.coords                       <- approx(arc.length, x[,1],n=n.vec[i])
      # y.coords                       <- approx(arc.length, x[,2],n=n.vec[i])
      # z.coords                       <- approx(arc.length, x[,3],n=n.vec[i])
      # X                              <- cbind(x.coords$y,y.coords$y,z.coords$y)
      # ind                            <- (substr(rownames(face$curves), 1, nchar(curve.names[i])) == curve.names[i])
      # face$curves                    <- face$curves[!ind, ]
      # rownames(X)                    <- paste(curve.names[i], 1:nrow(X))
      # face$curves                    <- rbind(face$curves, X)
   # }
