



"slidingcurves.face3d" <- function(shape1,shape2, n = 100)  {
   # choose the particular curve
 curve.names     <- c(    "cheek-eye right",
                          "cheek-eye left",
                          "cheek-lip right",
                          "cheek-lip left",
                          "cheek-nose right",
                          "cheek-nose left",
                          "brow ridge right", 
                          "brow ridge left",
                          "nasolabial right", 
                          "nasolabial left",
                          "philtrum ridge right",
                          "philtrum ridge left",
                          "upper lip right",
                          "upper lip left",
                          "lower lip right",
                          "lower lip left",
                          "nasal bridge left",
                          "nasal bridge right", 
                          "upper eye socket left",
                          "upper eye socket right",
                          "lower eye socket right",
                          "lower eye socket left",
                          "nasal root left",
                          "nasal root right",
                          "philtrum lip right",
                          "philtrum lip left",
                          "mid-line chin",
                          "mid-line mentolabial",
                          "mid-line upper lip",
                          "mid-line bottom lip",
                          "mid-line nasal root",
                          "mid-line nasal profile",
                          "nasal boundary right",
                          "nasal boundary left",
                          "mid-line philtral",
                          "nasal base right",
                          "nasal base left",
                          "mid-line columella",
                          "mandible right",
                          "mandible left",
                          "mid-line lip left",
                          "mid-line lip right")
                          
   
 
    n             <- c(12, 12, 20, 20, 20, 20, 28, 28, 6, 6, 6, 6, 8, 
                            8,  9,  9,  7,  7, 10, 10,  9,  9, 7, 7, 2, 2,
                            14, 4, 5, 5, 5, 20, 20, 20, 8,
                            8, 8, 8, 28, 28, 9, 9 )
   n <- n * 20   

   curves2.optim <- shape1$curves
   for (i in 1:42) {
   ind.i <- grep(curve.names[i],dimnames(shape1$curves)[[1]])
   ind.shape2 <- grep(curve.names[i],dimnames(shape2$curves)[[1]])
   curves2.optim[ind.i,] <-sliding.face3d(shape1$curves[ind.i,],shape2$curves[ind.shape2,],
   							n=n[i], resample=T)
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
