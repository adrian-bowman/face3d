mesh.face3d <- function(face, meshes, template, graphics=FALSE, monitor=FALSE) {
	
  if (!("meshes" %in% names(face))) face$meshes <- matrix(nrow = 0, ncol = 3)	
  
  
  
  
  
  #names of columns-- curve1, curve2, npts1 (1st resampling), npts2(resampling on the new paths), startpt (which pts to start on),   endpt( which pts to end on) 
  meshes.df <- matrix( c(  "mid-face right",         "cheek-lip right",          "cheek-nose right",                         20,  40, 1, 20,
                           "mid-face left",          "cheek-lip left",           "cheek-nose left",                          20,  40, 1, 20, 
                           "upper mid face right",   "cheek-nose right",         "cheek-eye right",                          20,  40, 1, 20,  
                           "upper mid face left",    "cheek-nose left",          "cheek-eye left",                           20,  40, 1, 20,  
                           "upper lip left",         "upper lip left",           "mid-line lip left",                         9,  20, 3,  9,
                           "lower lip left",         "lower lip left",           "mid-line lip left",                         9,  20, 3,  9,
                           "philtrum right",         "upper lip right",          "nasal base right",                          9,  30, 1,  9,
                           "philtrum left",          "upper lip left",           "nasal base left",                           9,  30, 1,  9, 
                           "upper lip right",        "upper lip right",          "mid-line lip right",                        9,  20, 3,  9,
                           "lower face right",       "mandible right",           "lower lip right",                          28,  50, 1, 28,
                           "nose left",              "nasal bridge left",        "nasal root left",                           7,  40, 1,  7,       
                           "nose right",             "nasal bridge right",       "nasal root right",                          7,  40, 1,  7,
                           "lower lip right",        "lower lip right",          "mid-line lip right",                        9,  20, 3,  9,  
                           "upper face right 1",     "brow ridge right",         "cheek-eye right",                          12,  10, 1, 12, 
                           "upper face right 2",     "brow ridge right",         "upper eye socket right",                   10,  20, 1, 10,
                           "upper face right 3",     "brow ridge right",         "nasal root right",                          7,  20, 1,  7,  
                           "upper face left 1",      "brow ridge left",          "cheek-eye left",                           12,  10, 1, 12, 
                           "upper face left 2",      "brow ridge left",          "upper eye socket left",                    10,  20, 1, 10,
                           "upper face left 3",      "brow ridge left",          "nasal root left",                           7,  20, 1,  7,        
                           "lower face left",        "mandible left",            "lower lip left",                           28,  50, 1, 28), 
                            ncol= 7, byrow = TRUE)	                                  
                           
  meshes.df                   <- as.data.frame(meshes.df[ , 2:7], row.names = meshes.df[ ,1])        
  if (missing(meshes)) meshes <- rownames(meshes.df)
  meshes.df                   <- meshes.df[meshes, ]
  curve.nms                   <- as.matrix(meshes.df[ , 1:2])
  npts1                       <- as.numeric(as.character(meshes.df[ , 3]))
  npts2                       <- as.numeric(as.character(meshes.df[ , 4]))
  startpt                     <- as.numeric(as.character(meshes.df[ , 5])) 
  endpt                       <- as.numeric(as.character(meshes.df[ , 6]))
  
      for (i in 1:length(meshes)){
 	     	   
 	     if (length(grep("upper mid face right", meshes[i])) > 0){
 	     	  ind1        <- grep(curve.nms[i,1],  rownames(face$curves))
              curve1      <- face$curves[ind1,]
              ind2        <- grep(curve.nms[i,2],  rownames(face$curves))
              curve2      <- face$curves[ind2,] 
              ind3        <- grep("lower eye socket right",  rownames(face$curves))
              curve3      <- face$curves[ind3,]
              curve3      <- curve3[-1, ] #removing lmk
 	     	  curve2      <- rbind(curve2,curve3)
 	     	
 	     } else if (length(grep("upper mid face left", meshes[i])) > 0){
 	     	  ind1        <- grep(curve.nms[i,1],  rownames(face$curves))
              curve1      <- face$curves[ind1,]
              ind2        <- grep(curve.nms[i,2],  rownames(face$curves))
              curve2      <- face$curves[ind2,]
              ind3        <- grep("lower eye socket left",  rownames(face$curves))
              curve3      <- face$curves[ind3,] 
              curve3      <- curve3[-1, ]
 	     	  curve2      <- rbind(curve2,curve3)
 	     	
 	      } else if (length(grep("upper lip left", meshes[i])) > 0){
 	     	  ind1            <- grep(curve.nms[i,1], rownames(face$curves))
              curve1          <- face$curves[ind1,]  #upper lip left
              curve3          <- face$curves[grep("philtrum lip left", rownames(face$curves)), ]
              curve1          <- rbind(curve1, curve3[-1,])
              ind2            <- grep(curve.nms[i,2], rownames(face$curves))
              curve2          <- face$curves[ind2,]    # mid line
         
          } else if (length(grep("philtrum left", meshes[i])) > 0){
 	     	      ind1            <- grep(curve.nms[i,1], rownames(face$curves))
              curve1          <- face$curves[ind1,]  #upper lip left
              curve3          <- face$curves[grep("philtrum lip left", rownames(face$curves)), ]
              curve1          <- rbind(curve1, curve3[-1,])
              ind2            <- grep(curve.nms[i,2], rownames(face$curves))
              curve2          <- face$curves[ind2,]    # nasal base
              curve2          <- resample.face3d(curve2, n = npts1[i])


         } else if (length(grep("upper lip right", meshes[i])) > 0){
 	     	  ind1            <- grep(curve.nms[i,1], rownames(face$curves))
              curve1          <- face$curves[ind1,]  #upper lip right
              curve3          <- face$curves[grep("philtrum lip right", rownames(face$curves)), ]
              curve1          <- rbind(curve1, curve3[-1, ])
              ind2            <- grep(curve.nms[i,2], rownames(face$curves))
              curve2          <- face$curves[ind2,]    # mid line
              
          } else if (length(grep("philtrum right", meshes[i])) > 0){
 	     	  ind1            <- grep(curve.nms[i,1], rownames(face$curves))
              curve1          <- face$curves[ind1,]  #upper lip right
              curve3          <- face$curves[grep("philtrum lip right", rownames(face$curves)), ]
              curve1          <- rbind(curve1, curve3[-1, ])
              ind2            <- grep(curve.nms[i,2], rownames(face$curves))
              curve2          <- face$curves[ind2,]    #nasal base
              curve2          <- resample.face3d(curve2, n = npts1[i])


         } else if (length(grep("lower face right", meshes[i])) > 0){
 	     	  ind1            <- grep(curve.nms[i,1], rownames(face$curves))
              curve1          <- face$curves[ind1,]    #mand
              ind2            <- grep(curve.nms[i,2], rownames(face$curves))
              curve2          <- face$curves[ind2,]     #lower lip r
              curve3          <- face$curves[grep("cheek-lip right", rownames(face$curves)),]
              curve2          <- rbind(curve3, curve2[-1,])
          
         } else if (length(grep("OM lower face right", meshes[i])) > 0){
 	     	  ind1            <- grep(curve.nms[i,1], rownames(face$curves))
              curve1          <- face$curves[ind1,]    #mand
              ind2            <- grep(curve.nms[i,2], rownames(face$curves))
              curve2          <- face$curves[ind2,]     #lower lip r
              curve3          <- face$curves[grep("cheek-lip right", rownames(face$curves)),]
              curve2          <- rbind(curve3, curve2[-1,])
              
         } else if (length(grep("upper face right 1", meshes[i])) > 0){
 	     	      ind1              <- grep(curve.nms[i,1],  rownames(face$curves))
              curve1            <- face$curves[ind1,]    # brow ridge
              ind2              <- grep(curve.nms[i,2],  rownames(face$curves))
              curve2            <- face$curves[ind2,]    # cheek eye
              pt                <- face$lmks["exR", ]
              Distance          <- rep(0)
              for (a in 1:length(curve1[,1]))
	          Distance[a]       <- sqrt(apply((diff(rbind(pt, curve1[a,])))^2, 1, sum))
              names(Distance)   <- paste(rownames(curve1))
              pt.brow           <- which.min(Distance)
              curve1            <- curve1[1:pt.brow, ] 
              curve1            <-  resample.face3d(curve1, n=npts1[i])   
            # curve1            <- sliding.face3d(template, curve1, resample=T, n=npts1[i])  
                          
        } else if (length(grep("upper face right 2", meshes[i])) > 0){
 	     	  ind1              <- grep(curve.nms[i,1],  rownames(face$curves))
              curve1            <- face$curves[ind1,]    
              ind2              <- grep(curve.nms[i,2],  rownames(face$curves))
              curve2            <- face$curves[ind2,]   
              pt                <- face$lmks["exR", ]
              Distance          <- rep(0)
              for (a in 1:length(curve1[,1]))
	          Distance[a]       <- sqrt(apply((diff(rbind(pt, curve1[a,])))^2, 1, sum))
              names(Distance)   <- paste(rownames(curve1))
              pt.brow.right     <- which.min(Distance)     #pt closest to exocanthion
              reference.path.r  <- rbind(face$lmks["se", ], face$lmks["exR", ])
              lambda.r          <- cumsum(c(0, sqrt(apply((diff(reference.path.r))^2, 1, sum))))
              Euc.dist.r        <- .25*(lambda.r[2])
              arc.length.R      <- cumsum(sqrt(apply((curve1 - rbind(curve1[1, ],
                                   curve1[-(dim(curve1)[1]),]))^2, 1, sum)))
              closest           <- arc.length.R[length(arc.length.R)] -  Euc.dist.r           
              pt.brow.left      <- which.min(abs(closest -arc.length.R) ) #pt 25% of reference away from nasion
              curve1            <- curve1[pt.brow.right:pt.brow.left, ] 
              curve1            <-  resample.face3d(curve1, n=npts1[i])            

      
        } else if (length(grep("upper face right 3", meshes[i])) > 0){
 	     	  ind1              <- grep(curve.nms[i,1],  rownames(face$curves))
              curve1            <- face$curves[ind1,]    #brow ridge
              ind2              <- grep(curve.nms[i,2],  rownames(face$curves))
              curve2            <- face$curves[ind2,]   
              reference.path.r  <- rbind(face$lmks["se", ], face$lmks["exR", ])
              lambda.r          <- cumsum(c(0, sqrt(apply((diff(reference.path.r))^2, 1, sum))))
              Euc.dist.r        <- .25*(lambda.r[2])
              arc.length.R      <- cumsum(sqrt(apply((curve1 - rbind(curve1[1, ],
                                   curve1[-(dim(curve1)[1]),]))^2, 1, sum)))
              closest           <- arc.length.R[length(arc.length.R)] -  Euc.dist.r           
              pt.brow.left      <- which.min(abs(closest -arc.length.R) ) #pt 25% of reference away from nasion
              curve1            <- curve1[pt.brow.left:length(curve1[,1]), ]
              curve1            <-  resample.face3d(curve1, n=npts1[i])    
               #curve1            <- sliding.face3d(template, curve1, resample=T, n=npts1[i])              
               #curve2            <- sliding.face3d(template, curve2, resample=T, n=npts1[i])              

         } else if (length(grep("upper face left 1", meshes[i])) > 0){
 	     	  ind1              <- grep(curve.nms[i,1],  rownames(face$curves))
              curve1            <- face$curves[ind1,]    # brow ridge
              ind2              <- grep(curve.nms[i,2],  rownames(face$curves))
              curve2            <- face$curves[ind2,]    # cheek eye
              pt                <- face$lmks["exL", ]
              Distance          <- rep(0)
              for (a in 1:length(curve1[,1]))
	          Distance[a]       <- sqrt(apply((diff(rbind(pt, curve1[a,])))^2, 1, sum))
              names(Distance)   <- paste(rownames(curve1))
              pt.brow           <- which.min(Distance)
              curve1            <- curve1[1:pt.brow, ]  
              curve1            <-  resample.face3d(curve1, n=npts1[i])  
              # curve1            <- sliding.face3d(template, curve1, resample=T, n=npts1[i])              
                    
        } else if (length(grep("upper face left 2", meshes[i])) > 0){
 	     	  ind1              <- grep(curve.nms[i,1],  rownames(face$curves))
              curve1            <- face$curves[rev(ind1),]    # brow ridge
              ind2              <- grep(curve.nms[i,2],  rownames(face$curves))
              curve2            <- face$curves[ind2,]   
              pt                <- face$lmks["exL", ]
              Distance          <- rep(0)
              for (a in 1:length(curve1[,1]))
	          Distance[a]       <- sqrt(apply((diff(rbind(pt, curve1[a,])))^2, 1, sum))
              names(Distance)   <- paste(rownames(curve1))
              pt.brow.left      <- which.min(Distance)     #pt closest to exocanthion
              reference.path.l  <- rbind(face$lmks["se", ], face$lmks["exL", ])
              lambda.l          <- cumsum(c(0, sqrt(apply((diff(reference.path.l))^2, 1, sum))))
              Euc.dist.l        <- .25*(lambda.l[2])
              arc.length.L      <- cumsum(sqrt(apply((curve1 - rbind(curve1[1, ],
                                               curve1[-(dim(curve1)[1]),]))^2, 1, sum)))
              closest           <- Euc.dist.l          
              pt.brow.right     <- which.min(abs(closest -arc.length.L) ) #pt 25% of reference away from nasion
              curve1            <- curve1[pt.brow.left:pt.brow.right, ]
              curve1            <-  resample.face3d(curve1, n=npts1[i])      
              # curve1            <- sliding.face3d(template, curve1, resample=T, n=npts1[i])              

              
        } else if (length(grep("upper face left 3", meshes[i])) > 0){
 	     	  ind1              <- grep(curve.nms[i,1],  rownames(face$curves))
              curve1            <- face$curves[rev(ind1),]    #brow ridge
              ind2              <- grep(curve.nms[i,2],  rownames(face$curves))
              curve2            <- face$curves[ind2,]   #nasalroot left
              reference.path.l  <- rbind(face$lmks["se", ], face$lmks["exL", ])
              lambda.l          <- cumsum(c(0, sqrt(apply((diff(reference.path.l))^2, 1, sum))))
              Euc.dist.l        <- .25*(lambda.l[2])
              arc.length.L      <- cumsum(sqrt(apply((curve1 - rbind(curve1[1, ],
                                   curve1[-(dim(curve1)[1]),]))^2, 1, sum)))
              closest           <-  Euc.dist.l           
              pt.brow.left      <- which.min(abs(closest -arc.length.L) ) #pt 25% of reference away from nasion
              curve1            <- curve1[rev(1:pt.brow.left) , ] 
              curve1            <-  resample.face3d(curve1, n=npts1[i]) 
               #curve1            <- sliding.face3d(template, curve1, resample=T, n=npts1[i])              
               #curve2            <- sliding.face3d(template, curve2, resample=T, n=npts1[i])              

              
         } else if (length(grep("lower face left", meshes[i])) > 0){
 	     	  ind1              <- grep(curve.nms[i,1], rownames(face$curves))
              curve1            <- face$curves[ind1,]    #mand
              ind2              <- grep(curve.nms[i,2], rownames(face$curves))
              curve2            <- face$curves[ind2,]     #lower lip l
              curve3            <- face$curves[grep("cheek-lip left", rownames(face$curves)),]
              curve2            <- rbind(curve3, curve2[-1,])
        
   	     } else {
   	     	  ind1               <- grep(curve.nms[i,1],  rownames(face$curves))
              curve1             <- face$curves[ind1,]
   	     	   #curve1             <- sliding.face3d(template, curve1, resample=T, n=npts1[i])              
              ind2               <- grep(curve.nms[i,2],  rownames(face$curves))
 	          curve2             <- face$curves[ind2,]
               #curve2             <- sliding.face3d(template, curve2, resample=T, n=npts1[i])              
         }


if (length(grep("lip", meshes[i])) > 0){

 paths               <- array(0, dim = c(npts2[i] , 3, npts1[i]))
 for (k in startpt[i]:endpt[i]){
	  x1            <- curve1[k, ]
	  x2            <- curve2[k, ]	
	  path          <- planepath.face3d(face, x1, x2, distance = 5, graphics = FALSE, rotation.range= pi/4, boundary = c(1,1))$path
	  path          <- resample.face3d(path, n=npts2[i])
      #paths[ , ,k]  <- sliding.face3d(template, path, resample=T, n=npts2[i])
      paths[ , ,k]  <- path
    }
}

else {
	paths               <- array(0, dim = c(npts2[i] , 3, npts1[i]))
 for (k in startpt[i]:endpt[i]){
	  x1            <- curve1[k, ]
	  x2            <- curve2[k, ]	
	  path          <- planepath.face3d(face, x1, x2, distance = 10, graphics = FALSE, rotation.range= pi/4)$path
	  path          <- resample.face3d(path, n=npts2[i])
	 # paths[ , ,k]  <- sliding.face3d(template, path, resample=T, n=npts2[i])
	  paths[ , ,k]  <- path
	  print(k)
	

          }
}
#Throwing away the endpts of the paths and the unused beginning and end where other curves are
 paths          <- paths[ , , startpt[i]:endpt[i]]
  total.paths    <- endpt[i]- startpt[i] + 1
 # total.pts      <- npts2[i]-2
 # path           <- array(0, dim = c(total.pts , 3, total.paths))
 # for (j in 1:total.paths) { 	   
      # A         <- paths[-npts2[i], ,j]
      # B         <- A[-1,]
      # path[,,j] <- B
     # }
     
 #changes path below to paths and total pts    
     total.pts <- npts2[i]
 X              <- paths[,,1]
 for (a in 2: total.paths){ 
      X         <- rbind(X, paths[,,a])
     }

#Putting mesh into face object
   ind            <- (substr(rownames(face$meshes), 1, nchar(meshes[i])) == meshes[i])
   face$meshes    <- face$meshes[!ind, ]
   rownames(X)    <- paste(meshes[i], paste(rep(1:total.paths, each = total.pts)), paste(rep(1:total.pts, total.paths)))
   face$meshes    <- rbind(face$meshes, X)  
 }


 invisible(face)

}
