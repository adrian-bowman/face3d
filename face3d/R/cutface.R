
cutface.face3d <- function(face, monitor = FALSE){
		
#calculating proportion reference
	reference.path <- rbind(face$landmarks["enR", ], face$landmarks["enL", ])
	lambda         <- cumsum(c(0, sqrt(apply((diff(reference.path))^2, 1, sum))))
	Euc.dist       <- lambda[2]
 	 
#calculating cut pts cheek-eye and new curve
    path.R        <- face$curves[rev(grep("cheek-eye right", rownames(face$curves))), ]
    arc.length.R  <- cumsum(sqrt(apply((path.R - rbind(path.R[1, ],
                          path.R[-(dim(path.R)[1]),]))^2, 1, sum)))
    Cutpt.R       <- which.min(abs(arc.length.R -(1.4*Euc.dist)) )
    Cutpt.R.eye   <- path.R[Cutpt.R, ]
    X             <- path.R[Cutpt.R:1, ]
    ind           <- (substr(rownames(face$curves), 1, nchar("cheek-eye right")) == "cheek-eye right")
    face$curves   <- face$curves[!ind, ]
    rownames(X)   <- paste("cheek-eye right", 1:nrow(X))
    face$curves   <- rbind(face$curves, X)
    path.L        <- face$curves[grep("cheek-eye left", rownames(face$curves)), ]
    arc.length.L  <- cumsum(sqrt(apply((path.L - rbind(path.L[1, ],
                          path.L[-(dim(path.L)[1]),]))^2, 1, sum)))
    Cutpt.L       <- which.min(abs(arc.length.L -(1.4*Euc.dist)) )
    Cutpt.L.eye   <- path.L[Cutpt.L, ]
    X             <- path.L[1:Cutpt.L, ]
    ind           <- (substr(rownames(face$curves), 1, nchar("cheek-eye left")) == "cheek-eye left")
    face$curves   <- face$curves[!ind, ]
    rownames(X)   <- paste("cheek-eye left", 1:nrow(X))
    face$curves   <- rbind(face$curves, X)
   
 #calculating cutpts cheek-lip and new curve
    path.R        <- face$curves[rev(grep("cheek-lip right", rownames(face$curves))), ]
    arc.length.R  <- cumsum(sqrt(apply((path.R - rbind(path.R[1, ],
                          path.R[-(dim(path.R)[1]),]))^2, 1, sum)))
    Cutpt.R       <- which.min(abs(arc.length.R -(2.2*Euc.dist)) )
    Cutpt.R.lip   <- path.R[Cutpt.R, ]
    X             <- path.R[Cutpt.R:1, ]
    ind           <- (substr(rownames(face$curves), 1, nchar("cheek-lip right")) == "cheek-lip right")
    face$curves   <- face$curves[!ind, ]
    rownames(X)   <- paste("cheek-lip right", 1:nrow(X))
    face$curves   <- rbind(face$curves, X)
    path.L        <- face$curves[grep("cheek-lip left", rownames(face$curves)), ]
    arc.length.L  <- cumsum(sqrt(apply((path.L - rbind(path.L[1, ],
                          path.L[-(dim(path.L)[1]),]))^2, 1, sum)))
    Cutpt.L       <- which.min(abs(arc.length.L -(2.2*Euc.dist)) )
    Cutpt.L.lip   <- path.L[Cutpt.L, ]
    X             <- path.L[1:Cutpt.L, ]
    ind           <- (substr(rownames(face$curves), 1, nchar("cheek-lip left")) == "cheek-lip left")
    face$curves   <- face$curves[!ind, ]
    rownames(X)   <- paste("cheek-lip left", 1:nrow(X))
    face$curves   <- rbind(face$curves, X)

  #calculating cutpts cheek-nose
    path.R        <- face$curves[rev(grep("cheek-nose right", rownames(face$curves))), ]
    arc.length.R  <- cumsum(sqrt(apply((path.R - rbind(path.R[1, ],
                          path.R[-(dim(path.R)[1]),]))^2, 1, sum)))
    Cutpt.R       <- which.min(abs(arc.length.R - (2.5* Euc.dist)) )
    Cutpt.R.nose  <- path.R[Cutpt.R, ]
    X             <- path.R[Cutpt.R:1, ]
    ind           <- (substr(rownames(face$curves), 1, nchar("cheek-nose right")) == "cheek-nose right")
    face$curves   <- face$curves[!ind, ]
    rownames(X)   <- paste("cheek-nose right", 1:nrow(X))
    face$curves   <- rbind(face$curves, X)
    path.L        <- face$curves[grep("cheek-nose left", rownames(face$curves)), ]
    arc.length.L  <- cumsum(sqrt(apply((path.L - rbind(path.L[1, ],
                          path.L[-(dim(path.L)[1]),]))^2, 1, sum)))
    Cutpt.L       <- which.min(abs(arc.length.L -(2.5* Euc.dist)) )
    Cutpt.L.nose  <- path.L[Cutpt.L, ]
    X             <- path.L[1:Cutpt.L, ]
    ind           <- (substr(rownames(face$curves), 1, nchar("cheek-nose left")) == "cheek-nose left")
    face$curves   <- face$curves[!ind, ]
    rownames(X)   <- paste("cheek-nose left", 1:nrow(X))
    face$curves   <- rbind(face$curves, X)


#calculating cutpts mandible
    pt              <-  face$landmarks["gn", ]
    Distance        <- rep(0)
    curve1          <- face$curves[grep("Mandible", rownames(face$curves)), ]
      for (a in 1:length(curve1[,1]))
        Distance[a] <- sqrt(apply((diff(rbind(pt, curve1[a,])))^2, 1, sum))
    names(Distance) <- paste(rownames(curve1))
    pt.mand         <- which.min(Distance)
    path.R          <- curve1[pt.mand:1, ]
    arc.length.R    <- cumsum(sqrt(apply((path.R - rbind(path.R[1, ],
                          path.R[-(dim(path.R)[1]),]))^2, 1, sum)))
    Cutpt.R         <- which.min(abs(arc.length.R -(3*Euc.dist)) )
    Cutpt.R.mand    <- path.R[Cutpt.R, ]
    path.L          <- curve1[pt.mand:length(curve1[,1]),]
    arc.length.L    <- cumsum(sqrt(apply((path.L - rbind(path.L[1, ],
                          path.L[-(dim(path.L)[1]),]))^2, 1, sum)))
    Cutpt.L         <- which.min(abs(arc.length.L -(3*Euc.dist)) )
    Cutpt.L.mand    <- path.L[Cutpt.L, ]
    X               <- rbind(path.R[Cutpt.R:1,], path.L[2:Cutpt.L,])
    ind             <- (substr(rownames(face$curves), 1, nchar("Mandible")) == "Mandible")
    face$curves     <- face$curves[!ind, ]
    rownames(X)     <- paste("Mandible", 1:nrow(X))
    face$curves     <- rbind(face$curves, X)
 
 #calculating cutpts brow ridge
    path.R          <- face$curves[rev(grep("brow ridge right", rownames(face$curves))), ]
    arc.length.R    <- cumsum(sqrt(apply((path.R - rbind(path.R[1, ],
                          path.R[-(dim(path.R)[1]),]))^2, 1, sum)))
    Cutpt.R         <- which.min(abs(arc.length.R - (3*Euc.dist)) )
    Cutpt.R.brow    <- path.R[Cutpt.R, ]
    X               <- path.R[Cutpt.R:1, ]
    ind             <- (substr(rownames(face$curves), 1, nchar("brow ridge right")) == "brow ridge right")
    face$curves     <- face$curves[!ind, ]
    rownames(X)     <- paste("brow ridge right", 1:nrow(X))
    face$curves     <- rbind(face$curves, X)
    path.L          <- face$curves[grep("brow ridge left", rownames(face$curves)), ]
    arc.length.L    <- cumsum(sqrt(apply((path.L - rbind(path.L[1, ],
                          path.L[-(dim(path.L)[1]),]))^2, 1, sum)))
    Cutpt.L         <- which.min(abs(arc.length.L - (3*Euc.dist)) )
    Cutpt.L.brow    <- path.L[Cutpt.L, ]
    X               <- path.L[1:Cutpt.L, ]
    ind             <- (substr(rownames(face$curves), 1, nchar("brow ridge left")) == "brow ridge left")
    face$curves     <- face$curves[!ind, ]
    rownames(X)     <- paste("brow ridge left", 1:nrow(X))
    face$curves     <- rbind(face$curves, X)

 #right side cut
    # a                <- planepath.face3d(face, Cutpt.R.brow, Cutpt.R.eye)$path
    # b                <- planepath.face3d(face, Cutpt.R.eye, Cutpt.R.nose)$path
    # c                <- planepath.face3d(face, Cutpt.R.nose, Cutpt.R.lip)$path
    # d                <- planepath.face3d(face, Cutpt.R.lip, Cutpt.R.mand)$path
    # Right.cut        <- rbind(a,b[-1,], c[-1,], d[-1, ])
    # Test.r           <- closestcurve.face3d(face, Right.cut)
    # perp.dist.r      <- Test.r$closest.distance
    # ind.right        <- perp.dist.r > -30
    # names(ind.right) <- paste(c(1: dim(face$vertices)[1]))
 
 #left side cut
    # a                <- planepath.face3d(face, Cutpt.L.brow, Cutpt.L.eye)$path
    # b                <- planepath.face3d(face, Cutpt.L.eye, Cutpt.L.nose)$path
    # c                <- planepath.face3d(face, Cutpt.L.nose, Cutpt.L.lip)$path
    # d                <- planepath.face3d(face, Cutpt.L.lip, Cutpt.L.mand)$path
    # Left.cut         <- rbind(a,b[-1,], c[-1,], d[-1, ])
    # Test.l           <- closestcurve.face3d(face, Left.cut)
    # perp.dist.l      <- Test.l$closest.distance
    # ind.left         <- perp.dist.l < 30
    # names(ind.left)  <- paste(c(1: dim(face$vertices)[1]))
 
 #1st subset   
     # crvs        <- face$curves
     # landmarks        <- face$landmarks
     # face        <- subset.face3d(face, (ind.left & ind.right))
     # face$landmarks   <- landmarks
     # face$curves <- crvs
  
 #neck cuts
    # mand            <- grep("Mandible", rownames(face$curves))
    # mandible        <- face$curves[mand, ]
    # TEST            <- closestcurve.face3d(face, mandible)
    # perp.dist       <- TEST$closest.distance
    # ind.mand        <- perp.dist > -20
    # names(ind.mand) <- paste(c(1: dim(face$vertices)[1]))
    
 #brow cuts
     # ind.R           <- grep("brow ridge right", rownames(face$curves))
     # brow.right      <- face$curves[ind.R,]
     # ind.L           <- grep("brow ridge left", rownames(face$curves))
     # brow.left       <- face$curves[ind.L,] 
     # brow.ridge      <- rbind(brow.right, brow.left[-1, ])
     # TEST            <- closestcurve.face3d(face, brow.ridge)
     # perp.dist       <- TEST$closest.distance
     # ind.brow        <- perp.dist < 10
     # names(ind.brow) <-  paste(c(1: dim(face$vertices)[1]))
  
  #2nd subset   
     # crvs        <- face$curves
     # landmarks        <- face$landmarks
     # face        <- subset.face3d(face, (ind.mand & ind.brow))
     # face$landmarks   <- landmarks
     # face$curves <- crvs

  #fixing gnathion curve
     pt              <-  face$landmarks["gn", ]
     Distance        <- rep(0)
     curve1          <- face$curves[grep("Mandible", rownames(face$curves)), ]
       for (a in 1:length(curve1[,1]))
	        Distance[a]      <- sqrt(apply((diff(rbind(pt, curve1[a,])))^2, 1, sum))
     names(Distance) <- paste(rownames(curve1))
     pt.mand         <- which.min(Distance)
     pt.mand         <- curve1[pt.mand, ]
    #   Subset the object by a tube around x1 and the relevant direction
      lmks     <- rbind(pt, pt.mand)
      boundary <- c(20,20)
      rng      <- sqrt(sum((lmks[2, ] - lmks[1, ])^2))
      bndry    <- boundary * rng
      unit     <- (lmks[2, ] - lmks[1, ]) / rng
      prjn     <- c(sweep(face$vertices, 2, lmks[1, ]) %*% unit)
      ind1     <- (prjn > - bndry[1]) & (prjn < rng + bndry[1])
      prjn2    <- outer(prjn, unit)
      ind2     <- apply((sweep(face$vertices, 2, lmks[1, ]) - prjn2)^2, 1, function(x) sqrt(sum(x)) < bndry[2])
      shape    <- subset.face3d(face, ind1 & ind2, remove.singles = TRUE)	  
    #is pt above or below the new mand curve(pt.mand included)?
	#which point on the shape vertices it is closest to
	  Distance        <- rep(0)
        for (a in 1:length(shape$vertices[,1]))
	        Distance[a]  <- sqrt(apply((diff(rbind(pt, shape$vertices[a,])))^2, 1, sum))
	    closest.pt.pt    <-   rownames(shape$vertices )[which.min(Distance)]
	    cutt             <- closestcurve.face3d(shape, curve1)
	    sign             <- sign(cutt$closest.distance[closest.pt.pt])
	 
      if  (sign == -1){
      midline           <- face$curves[grep("mid-line chin", rownames(face$curves)), ]
      Distance        <- rep(0)
        for (a in 1:length(midline[,1]))
	        Distance[a]  <- sqrt(apply((diff(rbind(pt.mand, midline[a,])))^2, 1, sum))
      newpt.mand        <- midline[which.min(Distance),]
      face$landmarks["gn", ] <-  newpt.mand
      X                 <- midline[c(1:which.min(Distance)),]
      ind               <- (substr(rownames(face$curves), 1, nchar("mid-line chin")) == "mid-line chin")
      face$curves       <- face$curves[!ind, ]
      rownames(X)       <- paste("mid-line chin", 1:nrow(X))
      face$curves       <- rbind(face$curves, X)
     	
     }

     else  {
      face$landmarks["gn", ] <-  pt.mand
      midline           <- face$curves[grep("mid-line chin", rownames(face$curves)), ]
  	  X                 <- rbind(midline, pt.mand)
      ind               <- (substr(rownames(face$curves), 1, nchar("mid-line chin")) == "mid-line chin")
      face$curves       <- face$curves[!ind, ]
      rownames(X)       <- paste("mid-line chin", 1:nrow(X))
      face$curves       <- rbind(face$curves, X)
     }
     
  
  
  
  
          
  	     
  	     
  	     
  	     
  	     
  	     
  	     
  	     
  	     
  	     
   #removing mL, mR      
     ind         <- grep("mandible left", rownames(face$curves))
     face$curves <- face$curves[-ind, ]
     ind         <- grep("mandible right", rownames(face$curves))
     face$curves <- face$curves[-ind, ]
      
    
       #splitting mid-line 
             mid.line        <- face$curves[grep("mid-line lip", rownames(face$curves)),]	
             pt              <- face$landmarks["st", ]
             Distance        <- rep(0)
               for (a in 1:length(mid.line[,1])) Distance[a]  <- sqrt(apply((diff(rbind(pt, mid.line[a,])))^2, 1, sum))
             names(Distance) <- paste(rownames(mid.line))
             new.st          <- which.min(Distance)           	     	   
	         ml.l            <- mid.line[new.st: length(mid.line[,1]), ]
	         X               <- ml.l
             ind             <- (substr(rownames(face$curves), 1, nchar("mid-line lip left")) ==    								"mid-line lip left")
             face$curves     <- face$curves[!ind, ]
             rownames(X)     <- paste("mid-line lip left", 1:nrow(X))
             face$curves     <- rbind(face$curves, X)
             ml.r            <- mid.line[1: new.st, ]
             X               <- ml.r
             ind             <- (substr(rownames(face$curves), 1, nchar("mid-line lip right")) 									== "mid-line lip right")
             face$curves     <- face$curves[!ind, ]
             rownames(X)     <- paste("mid-line lip right", 1:nrow(X))
             face$curves     <- rbind(face$curves, X)

#splitting mandible
             Mandible  <- face$curves[grep("Mandible", rownames(face$curves)),]	
             pt              <- face$landmarks["gn", ]
             Distance        <- rep(0)
               for (a in 1:length(Mandible[,1])) Distance[a]  <- sqrt(apply((diff(rbind(pt, Mandible[a,])))^2, 1, sum))
             names(Distance) <- paste(rownames(Mandible))
             new.gn          <- which.min(Distance)          
             mandible.left   <- Mandible[new.gn: length(Mandible[,1]), ]
	         X               <- mandible.left
             ind             <- (substr(rownames(face$curves), 1, nchar("mandible left")) == 								"mandible left")
             face$curves     <- face$curves[!ind, ]
             rownames(X)     <- paste("mandible left", 1:nrow(X))
             face$curves     <- rbind(face$curves, X)
             mandible.right  <- Mandible[1: new.gn, ] 
             X               <- mandible.right
             ind             <- (substr(rownames(face$curves), 1, nchar("mandible right")) == 								"mandible right")
             face$curves     <- face$curves[!ind, ]
             rownames(X)     <- paste("mandible right", 1:nrow(X))
             face$curves     <- rbind(face$curves, X)

#removing mand, mid     
     ind         <- grep("Mandible", rownames(face$curves))
     face$curves <- face$curves[-ind, ]
     ind         <- grep("mid-line lip", rownames(face$curves))
     face$curves <- face$curves[-ind[1:24], ]
      


  curve.names        <- c("nasal root left", "nasal boundary left", "nasal bridge left", 
  						  "nasal base left", "mid-line lip left", "lower lip left", 
 						  "upper lip left", "philtrum ridge left", "nasolabial left", 
                          "cheek-lip left",  "cheek-eye left", "cheek-nose left", 
                          "philtrum lip left", "mandible left", 
                          "lower eye socket left", "brow ridge left")  
                          
                                                                
  for (i in 1:length(curve.names)) {
  	  X             <- face$curves[rev(grep(curve.names[i], rownames(face$curves))),] 
      ind           <- (substr(rownames(face$curves), 1, nchar(curve.names[i])) == curve.names[i])
      face$curves   <- face$curves[!ind, ]
      rownames(X)   <- paste(curve.names[i], 1:nrow(X))
      face$curves   <- rbind(face$curves, X)
     
   }
   
 
 

   invisible(face)
}
