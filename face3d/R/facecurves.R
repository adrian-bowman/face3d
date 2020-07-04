
# plot(face)
# spheres3d(face$curves)



facecurves.face3d <- function(face, curves, monitor = FALSE, pcrv.path = FALSE, start.pc = NA) {

   crvs.missing <- missing(curves)
   if (!("curves" %in% names(face))) face$curves <- matrix(nrow = 0, ncol = 3)
   paths        <- list()
# Col Names of matrix - curve name, lmk1, lmk2, df, si.target, resampled points/npts , boundary1, boundary2, perp.dist.bound,   distance, rotation range
   curves.df    <- matrix(c("nasal root left",            "se",  "enL",  0,   NA, 50,  0.1,    0.3,    6,  5, pi/2,
                            "nasal root right",          "enR",   "se",  0,   NA, 50,  0.1,    0.3,    6,  5, pi/2,
                            "nasal boundary left",       "acL",  "enL",  0,   NA, 50,  0.2,    0.5,    3,  5, pi/4,
                            "nasal boundary right",      "enR",  "acR",  0,   NA, 50,  0.2,    0.5,    3,  5, pi/4,
                            "nasal bridge left",          "pn",  "acL",  5,  0.5, 40,  0.2,    0.5,    10,  5, pi/2,
                            "nasal bridge right",        "acR",   "pn",  5,  0.5, 40,  0.2,    0.5,    10,  5, pi/2,
                            "mid-line philtral",          "sn",   "ls",  0,   NA, 20,  0.6,    1.0,    5,  5, pi/4,
                            "nasal base left",            "sn",  "acL",  0,   NA, 40,  0.2,    0.4,    4,  5, pi/2,
                            "nasal base right",          "acR",   "sn",  0,   NA, 40,  0.2,    0.4,    4,  5, pi/2,
                            "mid-line lip",              "chR",  "chL",  5, -0.5, 80,  0.2,    0.5,    6,  5, pi/2,
                            "lower lip left",             "li",  "chL",  5,  0.5, 40,  0.1,    0.3,    6,  5, pi/2,
                            "lower lip right",           "chR",   "li",  5,  0.5, 40,  0.1,    0.3,    6,  5, pi/2,
                            "upper lip left",           "cphL",  "chL",  5,  0.5, 40,  0.1,    0.3,   10,  5, pi/2,
                            "upper lip right",           "chR", "cphR",  5,  0.5, 40,  0.1,    0.3,   10,  5, pi/2,
                            "philtrum ridge right",     "cphR",   "sn",  0,   NA, 20,  0.1,    0.3,    6,  5, pi/2,
                            "philtrum ridge left",        "sn", "cphL",  0,   NA, 20,  0.1,    0.3,    6,  5, pi/2,
                            "philtrum lip left",          "ls", "cphL",  0,  0.5, 10,  0.4,    0.8,    3,  5, pi/4,
                            "philtrum lip right",       "cphR",   "ls",  0,  0.5, 10,  0.4,    0.8,    3,  5, pi/4,
                            "mid-line nasal profile",     "se",   "pn",  0,   NA, 40,  0.2,    0.5,   15,  5, pi/2,
                            "mid-line nasal root",         "n",   "se",  0,   NA, 20,  0.2,    0.5,    6,  5, pi/2,
                            "mid-line columella",         "pn",   "sn",  0,   NA, 20,  0.5,    1.0,    2,  5, pi/2,
                            "mid-line upper-lip",         "ls",   "st",  0,   NA, 15,  0.5,    1.0,    5,  5, pi/2,
                            "mid-line bottom lip",        "st",   "li",  0,   NA, 15,  0.5,    1.0,    6,  5, pi/4,
                            "mid-line mentolabial",       "li",   "sl",  0,   NA, 20,  0.2,    0.5,    6,  5, pi/2,
                            "mid-line chin",              "sl",   "gn",  0,   NA, 30,  0.2,    0.5,    3,  5, pi/2,
                            "nasolabial left",           "acL",  "chL",  0,   NA, 50,  0.2,    0.5,   10,  5, pi/4,
                            "nasolabial right",          "chR",  "acR",  0,   NA, 50,  0.2,    0.5,   10,  5, pi/4,
                            "cheek-lip left",            "chL",  "oiL",  0,   NA,  0,  0.1,    0.3,    7,  5, pi/2,
                            "cheek-lip right",           "oiR",  "chR",  0,   NA,  0,  0.1,    0.3,    7,  5, pi/2,
                            "cheek-nose left",           "acL",   "tL",  0,   NA,  0,  0.1,    0.3,    6,  5, pi/2,
                            "cheek-nose right",           "tR",  "acR",  0,   NA,  0,  0.1,    0.3,    6,  5, pi/2,
                            "cheek-eye left",            "exL",   "tL",  0,   NA,  0,  0.1,    0.3,    6,  5, pi/2,
                            "cheek-eye right",            "tR",  "exR",  0,   NA,  0,  0.1,    0.3,    6,  5, pi/2,
                            "ear left",                   "tL",  "oiL",  0,   NA,  8,  0.2,    0.5,    6,  5, pi/2,
                            "ear right",                  "tR",  "oiR",  0,   NA,  8,  0.2,    0.5,    6,  5, pi/2,
                            "mandible right",            "oiR",   "gn",  0,  NA, 250, 0.1,    0.3,   12, 10, pi/2,
                            "mandible left",              "gn",  "oiL",  0,  NA, 250, 0.1,    0.3,   12, 10, pi/2,
                            "Mandible",                  "oiR",  "oiL",  5,  0.5, 500, 0.2,    1.0,   10, 15, pi/2,
                            "upper eye socket right",    "exR",  "enR",  0, -0.5, 10,  0.4,    1.0,   20,  5, pi/2,
                            "upper eye socket left",     "enL",  "exL",  0, -0.5, 10,  0.4,    1.0,   20,  5, pi/2,
                            "lower eye socket right",    "exR",  "enR",  5,  0.5, 30,  0.2,    0.5,   15,  5, pi/2,
                            "les right",                 "exR",  "enR",  5, -0.5, 16,  0.4,    1.0,   15,  5, pi/2,
                            "lower eye socket left",     "enL",  "exL",  5,  0.5, 30,  0.2,    0.5,   15,  5, pi/2,
                            "les left",                  "enL",  "exL",  5, -0.5, 16,  0.4,    1.0,   15,  5, pi/2,
                            "brow ridge right",           "tR",    "n",  5,  0.5, 250, 0.1,    0.35,  10, 10, pi/2,
                            "brow ridge left",             "n",   "tL",  5,  0.5, 250, 0.1,    0.35,  10, 10, pi/2,
                            "OM upper lip right",        "chR", "cphR",  5,  0.5,  40, 0.2,    0.5,    8,  5, pi/2,
                            "OM upper lip left",        "cphL", "chL",   5,  0.5,  40, 0.2,    0.5,    8,  5, pi/2,
                            "OM lower lip right",        "chR", "li",    5,  0.5,  40,  0.2,    0.5,    8,  5, pi/2,
                            "OM lower lip left",          "li", "chL",   5,  0.5,  40,  0.2,    0.5,    8,  5, pi/2,
                            "OM mid-line lip",            "chR", "chL",  0,   NA, 80,  0.2,    0.5,    8,  5, pi/2,
                            "Brow ridge",                 "tR",   "tL",  5,  0.5, 80,  0.2,    1.0,   15, 15, pi/2),
                          ncol = 11, byrow = TRUE)
   curves.df       <- as.data.frame(curves.df[ , 2:11], row.names = curves.df[ ,1])
   if (crvs.missing) curves <- rownames(curves.df)
   if (!all(curves %in% rownames(curves.df))) stop("curve not recognised.")
   curves.df <- curves.df[curves, ]
   # if (df.missing) 
   df              <- as.numeric(as.character(curves.df[ , 3]))
   lmks.nms        <- as.matrix(curves.df[ , 1:2])
   si.target       <- as.numeric(as.character(curves.df[ , 4]))
   npts            <- as.numeric(as.character(curves.df[ , 5]))
   boundary        <- cbind(as.numeric(as.character(curves.df[ , 6])), as.numeric(as.character(curves.df[ , 7])))
   perp.dist.bound <- as.numeric(as.character(curves.df[ , 8]))
   distance        <- as.numeric(as.character(curves.df[ , 9]))
   rotation.range  <- as.numeric(as.character(curves.df[ , 10]))
   ind             <- !(c(lmks.nms) %in% rownames(face$landmarks))
   if (any(ind))
      stop(paste("landmarks not found:", unique(c(lmks.nms)[ind]), collapse = " "))
   nms.crvs   <- character(0)
   curves.mat <- matrix(nrow = 0, ncol = 3)
   if (monitor) {
      cat(paste(c("Target:   ", rep(".", length(curves))), collapse = ""), "\n")
      cat("Progress: ")
   }
   # for (k in 1:nrow(face$landmarks)) face$landmarks[k, ] <- closest.face3d(face$landmarks[k, ], face)$points
   
   for (i in 1:length(curves)) {
      lmks <- face$landmarks[lmks.nms[i, ], ]
      reference.path <- NA
      shape.smooth   <- NA
      #  Subset the object by a tube around x1 and the relevant direction
      rng   <- sqrt(sum((lmks[2, ] - lmks[1, ])^2))
      bndry <- boundary[i, ] * rng
      unit  <- (lmks[2, ] - lmks[1, ]) / rng
      prjn  <- c(sweep(face$vertices, 2, lmks[1, ]) %*% unit)
      near  <- function(x, pts, distance)
                    sqrt(.rowSums((sweep(pts, 2, x))^2, nrow(pts), 3)) < distance
      ind1  <- (prjn > - bndry[1]) & (prjn < rng + bndry[1])
      prjn2 <- outer(prjn, unit)
      ind2  <- apply((sweep(face$vertices, 2, lmks[1, ]) - prjn2)^2, 1, function(x) sqrt(sum(x)) < bndry[2])
      shape <- subset.face3d(face, ind1 & ind2, remove.singles = TRUE)	  
 	  if (length(grep("OM upper lip", curves[i])) > 0) {
     	    path.ml                        <- planepath.face3d(shape, lmks[1,], lmks[2,])$path   
   	  	 	cutt                           <- closestcurve.face3d(shape, path.ml)
   	        shape                          <- index.face3d(shape, distance=distance[i]) 
   	        ind1                           <- sign(shape$shape.index)!=sign(si.target[i])
   	        ind11                          <- which(cutt$closest.distance < 0) 
            shape$shape.index[ind11]       <- 0
            shape$kappa1[ind11]            <- 0
            shape$kappa2[ind11]            <- 0   
            shape$shape.index[ind1]        <- 0
            shape$kappa1[ind1]             <- 0
            shape$kappa2[ind1]             <- 0 
            reference.path                 <- planepath.face3d(shape, lmks[1,], lmks[2,])$path
            

                   }
     if(length(grep("nasal bridge", curves[i])) > 0) {
     	    path <- planepath.face3d(shape, lmks[1,], lmks[2,])$path
     	    cutt <- closestcurve.face3d(shape, path)
            shape <- index.face3d(shape, distance=distance[i]) 
   	        ind1 <- sign(shape$shape.index)!=sign(si.target[i])
   	        ind11 <- which(cutt$closest.distance < 3) 
     	    shape$shape.index[ind11]       <- 0
            shape$kappa1[ind11]            <- 0
            shape$kappa2[ind11]            <- 0   
            shape$shape.index[ind1]        <- 0
            shape$kappa1[ind1]             <- 0
            shape$kappa2[ind1]             <- 0 
            reference.path                 <- planepath.face3d(shape, lmks[1,], lmks[2,], si.target=si.target[i])$path
            values                         <- pmax(-shape$kappa2, 0)
            if (monitor)                  plot(shape, colour=values*100)
       
     	
     	            }              
     if (length(grep("OM lower lip", curves[i])) > 0) {
     	    path.ml                        <- planepath.face3d(shape, lmks[1,], lmks[2,])$path   
   	  	 	cutt                           <- closestcurve.face3d(shape, path.ml)
   	        shape                          <- index.face3d(shape, distance=distance[i]) 
   	        ind1                           <- sign(shape$shape.index)!=sign(si.target[i])
   	        ind11                          <- which(cutt$closest.distance > 0) 
            shape$shape.index[ind11]       <- 0
            shape$kappa1[ind11]            <- 0
            shape$kappa2[ind11]            <- 0   
            shape$shape.index[ind1]        <- 0
            shape$kappa1[ind1]             <- 0
            shape$kappa2[ind1]             <- 0 
            reference.path                 <- planepath.face3d(shape, lmks[1,], lmks[2,])$path
                   }

     if (length(grep("Mandible", curves[i])) > 0) {
     	    indL                           <- grep("mandible left",  rownames(face$curves))
 	 	    mL                             <- face$curves[indL,]
 	     	indR                           <- grep("mandible right",  rownames(face$curves))
 	   	    mR                             <- face$curves[indR,]
   	   	  	reference.path                 <- rbind(mR, mL)
   	   	  	xy                             <- closestcurve.face3d(shape, reference.path)
            lambda                         <- cumsum(c(0, sqrt(apply((diff(reference.path))^2, 1, sum))))
            xcoord                         <- lambda[xy$closest.curvept] # distance along the curve
            ycoord                         <- xy$closest.distance
            indd                           <- (ycoord < 20 & ycoord > -20)
            shape                          <- subset.face3d(shape, indd)
            xy1                             <- closestcurve.face3d(shape, reference.path)
            lambda1                         <- cumsum(c(0, sqrt(apply((diff(reference.path))^2, 1, sum))))
            xcoord1                         <- lambda1[xy1$closest.curvept] # distance along the curve
            shape                          <- index.face3d(shape, distance=15)
            ind2                           <- (xcoord1 < (0.2 * rng)) & (xcoord1 > (0.8 * rng))
            shape$shape.index[ind2]         <- 0
            shape$kappa1[ind2]              <- 0
            shape$kappa2[ind2]              <- 0 
                  }    


     if (length(grep("upper lip", curves[i])) > 0) {   
     	    path.ml                        <- face$curves[grep("mid-line lip",  rownames(face$curves)), ]               
   	  	 	cutt                           <- closestcurve.face3d(shape, path.ml)
   	        ind                            <- which(cutt$closest.distance >= -2) 
   	        shape                          <- index.face3d(shape, distance=distance[i]) 
   	        ind1                           <- sign(shape$shape.index)!=sign(si.target[i])
   	        ind11                          <- which(cutt$closest.distance < -1) 
            shape$shape.index[ind11]       <- 0
            shape$kappa1[ind11]            <- 0
            shape$kappa2[ind11]            <- 0   
            shape$shape.index[ind1]        <- 0
            shape$kappa1[ind1]             <- 0
            shape$kappa2[ind1]             <- 0   
            #shape.smooth                   <- subset.face3d(shape, ind, remove.singles = FALSE)
            shape.smooth                   <- shape
            ind2                           <- sign(shape.smooth$shape.index)!=sign(si.target[i])
            shape.smooth$shape.index[ind2] <- 0
            shape.smooth$kappa1[ind2]      <- 0
            shape.smooth$kappa2[ind2]      <- 0  
            reference.path <- planepath.face3d(shape.smooth, lmks[1,], lmks[2,], si.target=si.target[i], monitor = monitor)$path   	   }  	  

     if (length(grep("mid-line lip", curves[i])) > 0) {    
     	    shape <- index.face3d(shape)
   	        ind1                           <- sign(shape$shape.index)!=sign(si.target[i]) 
            shape$shape.index[ind1]        <- 0
            shape$kappa1[ind1]             <- 0
            shape$kappa2[ind1]             <- 0   
           # plot(shape.smooth, colour=values*100)
            reference.path <- planepath.face3d(shape, lmks[1,], lmks[2,])$path
  	   }  
  	   
     if (length(grep("lower lip", curves[i])) > 0) {
     	    path.ml                        <- face$curves[grep("mid-line lip",  rownames(face$curves)), ]               
   	  	 	cutt                           <- closestcurve.face3d(shape, path.ml)
   	        ind                            <- which(cutt$closest.distance <= 1) 
   	        shape                          <- index.face3d(shape, distance=distance[i]) 
   	        ind1                           <- sign(shape$shape.index)!=sign(si.target[i])
   	        ind11                          <- which(cutt$closest.distance > -1) 
            shape$shape.index[ind11]       <- 0
            shape$kappa1[ind11]            <- 0
            shape$kappa2[ind11]            <- 0   
            shape$shape.index[ind1]        <- 0
            shape$kappa1[ind1]             <- 0
            shape$kappa2[ind1]             <- 0   
            #shape.smooth                   <- subset.face3d(shape, ind, remove.singles = FALSE)
            shape.smooth                   <- shape
            ind2                           <- sign(shape.smooth$shape.index)!=sign(si.target[i])
            shape.smooth$shape.index[ind2] <- 0
            shape.smooth$kappa1[ind2]      <- 0
            shape.smooth$kappa2[ind2]      <- 0  
            reference.path <- planepath.face3d(shape.smooth, lmks[1,], lmks[2,], si.target=si.target[i], monitor = monitor)$path 
       }  	 

   	if (length(grep("brow ridge left", curves[i])) > 0) {
   	 		lmksa                          <- rbind(face$landmarks["exL" , ], face$landmarks[ "enL", ])
   	 		boundaryy                      <- c(.2, .7)
   	 		rng2                            <- sqrt(sum((lmksa[2, ] - lmksa[1, ])^2))
            bndry                          <- boundaryy * rng2
            unit                           <- (lmksa[2, ] - lmksa[1, ]) / rng2
            prjn                           <- c(sweep(face$vertices, 2, lmksa[1, ]) %*% unit)
            ind1                           <- (prjn > - bndry[1]) & (prjn < rng2 + bndry[1])
            prjn                           <- outer(prjn, unit)
            ind2                           <- apply((sweep(face$vertices, 2, lmksa[1, ]) - prjn)^2, 1, 
            	                                  function(x) sqrt(sum(x)) < bndry[2])
            shape.eye                      <- subset.face3d(face, ind1 & ind2, remove.singles = FALSE) 
            shape                          <- index.face3d(shape, distance=distance[i])
            ind3                           <- sign(shape$shape.index)!=sign(si.target[i])
            shape$shape.index[ind3]        <- 0
            shape$kappa1[ind3]             <- 0
            shape$kappa2[ind3]             <- 0
            Match                          <- match(rownames(shape.eye$vertices), rownames(shape$vertices))
            shape$shape.index[Match]       <- 0
            shape$kappa1[Match]            <- 0
            shape$kappa2[Match]            <- 0  
            path.eye                       <- face$curves[grep("cheek-eye left",  rownames(face$curves)), ] 
            cutt                           <- closestcurve.face3d(shape, path.eye)
   	        ind                            <- which(cutt$closest.distance <= 0) 
   	        shape$shape.index[ind]         <- 0
            shape$kappa1[ind]              <- 0
            shape$kappa2[ind]              <- 0
            ref.path                       <- planepath.face3d(shape, lmks[1, ], lmks[2, ], boundary = boundary[i, ], 
                                              distance = distance[i], rotation.range=rotation.range[i], si.target =    
                                              si.target[i], monitor = monitor)$path
            cutt1                          <- closestcurve.face3d(shape, ref.path)
            lambda                         <- cumsum(c(0, sqrt(apply((diff(ref.path))^2, 1, sum))))
            xcoord                         <- lambda[cutt1$closest.curvept] # distance along the curve
            ind2                           <- xcoord > (0.7 * rng) 
            shape$shape.index[ind2]         <- 0
            shape$kappa1[ind2]              <- 0
            shape$kappa2[ind2]              <- 0 
            shape.smooth <- shape                     
            values                         <- pmax(-shape$kappa2, 0)  
            if (monitor)	 	           plot(shape, colour = 6 + values)
       }


     if (length(grep("brow ridge right", curves[i])) > 0) {
   	 		lmksa                          <- rbind(face$landmarks["enR" , ], face$landmarks[ "exR", ])
   	 		boundaryy                      <- c(.20, .7)
   	 		rng2                            <- sqrt(sum((lmksa[2, ] - lmksa[1, ])^2))
            bndry                          <- boundaryy * rng2
            unit                           <- (lmksa[2, ] - lmksa[1, ]) / rng2
            prjn                           <- c(sweep(face$vertices, 2, lmksa[1, ]) %*% unit)
            ind1                           <- (prjn > - bndry[1]) & (prjn < rng2 + bndry[1])
            prjn                           <- outer(prjn, unit)
            ind2                           <- apply((sweep(face$vertices, 2, lmksa[1, ]) - prjn)^2, 1, 
            	                                  function(x) sqrt(sum(x)) < bndry[2])
            shape.eye                      <- subset.face3d(face, ind1 & ind2, remove.singles = FALSE) 
            shape                          <- index.face3d(shape, distance=distance[i])
            ind3                           <- sign(shape$shape.index)!=sign(si.target[i])
            shape$shape.index[ind3]        <- 0
            shape$kappa1[ind3]             <- 0
            shape$kappa2[ind3]             <- 0
            Match                          <- match(rownames(shape.eye$vertices), rownames(shape$vertices))
            shape$shape.index[Match]       <- 0
            shape$kappa1[Match]            <- 0
            shape$kappa2[Match]            <- 0  
            path.eye                       <- face$curves[grep("cheek-eye right",  rownames(face$curves)), ] 
            cutt                           <- closestcurve.face3d(shape, path.eye)
   	        ind                            <- which(cutt$closest.distance <= 0) 
   	        shape$shape.index[ind]         <- 0
            shape$kappa1[ind]              <- 0
            shape$kappa2[ind]              <- 0
            ycoord                         <- cutt$closest.distance
            indd                           <- ycoord >25
            shape$shape.index[ind3]        <- 0
            shape$kappa1[ind3]             <- 0
            shape$kappa2[ind3]             <- 0 
            ref.path                       <- planepath.face3d(shape, lmks[1, ], lmks[2, ], boundary = boundary[i, ], 
                                              distance = distance[i], rotation.range=rotation.range[i], si.target =    
                                              si.target[i], monitor = monitor)$path
            cutt1                          <- closestcurve.face3d(shape, ref.path)
            lambda                         <- cumsum(c(0, sqrt(apply((diff(ref.path))^2, 1, sum))))
            xcoord                         <- lambda[cutt1$closest.curvept] # distance along the curve
            ind2                           <- xcoord < 0.3 * rng 
            shape$shape.index[ind2]         <- 0
            shape$kappa1[ind2]              <- 0
            shape$kappa2[ind2]              <- 0 
            shape$shape.index[ind3]         <- 0
            shape$kappa1[ind3]              <- 0
            shape$kappa2[ind3]              <- 0 
            shape.smooth <- shape                     
            values                         <- pmax(-shape.smooth$kappa2, 0)  
            if (monitor)	 	           plot(shape.smooth, colour = 6 + values)
          
   	  }
   
     if (length(grep("lower eye socket", curves[i])) > 0) {
     	   # reference.path                 <- apply(lmks, 2, function(x) seq(x[1], x[2], length = 50))
     	     reference.path                 <- planepath.face3d(shape, lmks[1, ], lmks[2, ])$path
     	    cut1                            <- closestcurve.face3d(shape, reference.path)
     	    ind1                            <- which(cut1$closest.distance >= -4)  
     	    shape                           <- index.face3d(shape, distance=distance)
     	    shape$shape.index[ind1]         <- 0
     	    shape$kappa1[ind1]              <- 0
     	    shape$kappa2[ind1]              <- 0
     	    ind                             <-  sign(shape$shape.index)!=sign(si.target[i])
     	    shape$shape.index[ind]          <- 0
     	    shape$kappa1[ind]               <- 0
     	    shape$kappa2[ind]               <- 0
     	    shape.smooth                    <- shape
 
     	 
      }
   	    if (length(grep("les left", curves[i])) > 0) {
     	    path.lid                       <- face$curves[grep("lower eye socket left",  rownames(face$curves)), ]              
   	  	 	cutt                           <- closestcurve.face3d(shape, path.lid)
   	        shape                          <- index.face3d(shape, distance=distance[i])
   	        ind                            <- sign(shape$shape.index)!=sign(si.target[i])
   	        ind1                           <- which(cutt$closest.distance < 2) 
            shape$shape.index[ind]         <- 0
            shape$kappa1[ind]              <- 0
            shape$kappa2[ind]              <- 0   
            shape$shape.index[ind1]        <- 0
            shape$kappa1[ind1]             <- 0
            shape$kappa2[ind1]             <- 0   
    	        values                         <- pmax(-shape.smooth$kappa2, 0)  
            pcurve3d                       <- pcurve.face3d(shape, weights = values, fixed=lmks, initial = path.lid, stretch=0)
            pcurve                         <- pcurve3d$s
            pcurve                         <- pcurve[order(pcurve[,1]),]
            for (k in 1:nrow(pcurve)) pcurve[k, ] <- closest.face3d(pcurve[k, ], face)$points
            reference.path                 <- pcurve
             if (monitor) { plot(shape, colour= 2 + values)
   	                         rgl::spheres3d(lmks, col = "black", radius = 1)
                             rgl::spheres3d(reference.path, col = "black", radius = 0.5)}

   	  }
   	  
   	  if (length(grep("les right", curves[i])) > 0) {
     	    path.lid                       <- face$curves[grep("lower eye socket right",  rownames(face$curves)), ]              
   	  	 	cutt                           <- closestcurve.face3d(shape, path.lid)
   	        shape                          <- index.face3d(shape, distance=distance[i])
   	        ind                            <- sign(shape$shape.index)!=sign(si.target[i])
   	        ind1                           <- which(cutt$closest.distance < 2) 
            shape$shape.index[ind]         <- 0
            shape$kappa1[ind]              <- 0
            shape$kappa2[ind]              <- 0   
            shape$shape.index[ind1]        <- 0
            shape$kappa1[ind1]             <- 0
            shape$kappa2[ind1]             <- 0   
    	    values                         <- pmax(-shape.smooth$kappa2, 0)  
            pcurve3d                       <- pcurve.face3d(shape, weights = values, fixed=lmks, initial = path.lid, stretch=0)
            pcurve                         <- pcurve3d$s
            pcurve                         <- pcurve[order(pcurve[,1]),]
            save(shape,pcurve,values, file="~/Desktop/data.dmp")
            stop()
            for (k in 1:nrow(pcurve)) pcurve[k, ] <- closest.face3d(pcurve[k, ], face)$points
            reference.path                 <- pcurve
             if (monitor) { plot(shape, colour= 2 + values)
   	                         rgl::spheres3d(lmks, col = "black", radius = 1)
                             rgl::spheres3d(reference.path, col = "black", radius = 0.5)}

   	  }


     if (length(grep("upper eye socket", curves[i])) > 0) {
            reference.path                 <- planepath.face3d(shape, lmks[1, ], lmks[2, ], boundary = boundary[i, ], 
                                                               distance = distance[i], ngrid=50, monitor = FALSE)$path    
   	  	 	cutt                           <- closestcurve.face3d(shape, reference.path)
   	        ind                            <- which(cutt$closest.distance > 0) 
  	        # shape                          <- subset.face3d(shape, ind, remove.singles = FALSE)
  	        reference.path                 <- apply(lmks, 2, function(x) seq(x[1], x[2], length = 50))
            shape                          <- index.face3d(shape, distance = distance)
            
            xy                             <- closestcurve.face3d(shape, reference.path)
            lambda                         <- cumsum(c(0, sqrt(apply((diff(reference.path))^2, 1, sum))))
            xcoord                         <- lambda[xy$closest.curvept] # distance along the curve
            ycoord                         <- xy$closest.distance
            indd                           <- (sign(shape$shape.index) != sign(si.target[i]))   | ycoord < 10 |ycoord > 17
            shape.smooth                   <- subset.face3d(shape, !indd)
   	        values                         <- pmax(shape.smooth$kappa1, 0)    
   	        if (monitor)                  plot(shape.smooth, colour= 2 + values)   	   
   	  }

      print("here")

      if ((df[i] == 0) && !is.na(si.target[i]))
         path <- planepath.face3d(shape, lmks[1, ], lmks[2, ], boundary = boundary[i, ], 
                                  distance = distance[i], rotation.range=rotation.range[i], si.target = si.target[i], 
                                  monitor = monitor)$path
      else if ((df[i] == 0) && is.na(si.target[i]))
         path <- planepath.face3d(shape, lmks[1, ], lmks[2, ], boundary = boundary[i, ],
                                  distance = distance[i], rotation.range=rotation.range[i], monitor = monitor)$path
      else
         path <- smoothpath.ridge.face3d(shape, shape.smooth, lmks[1, ], lmks[2, ], df = df[i], reference.path = reference.path, 
                                   si.target = si.target[i], rotation.range=rotation.range[i], curve.name = curves[i], boundary = boundary[i, ],
                                    npts = npts[i], distance = distance[i], perp.dist.bound = perp.dist.bound[i], monitor = monitor, 
                                    pcrv.path = pcrv.path, start.pc = start.pc)$path

      # Resampling path
        # if (npts == 0){
        	# path <- path
        # }
        # else path <- resample.face3d(path,n=npts,threshold=1e-15)

        ind            <- (substr(rownames(face$curves), 1, nchar(curves[i])) == curves[i])
        face$curves    <- face$curves[!ind, ]
        rownames(path) <- paste(curves[i], 1:nrow(path))
        face$curves    <- rbind(face$curves, path)

   	 # if (monitor) cat(curves[i], "... ")
   	   if (monitor) cat(".")
   }
   if (monitor) cat("\n")
   invisible(face)
}
