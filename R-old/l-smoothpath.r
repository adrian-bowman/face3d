smoothpath.face3d <- function(shape, shape.smooth, lmk1, lmk2, df = 5, distance = 5, curve.name = NA,
                              reference.path = NA, rotation.range=rotation.range, si.target, boundary = c(0.2, 0.5),
                              npts, graphics = FALSE,        
                              perp.dist.bound = 10, pcrv.path = FALSE, start.pc = NA, smooth.method = NULL) {  
	
   if (!require(fields))   stop("the fields package is required.")
   if (!require(rgl))      stop("the rgl package is required.")
   if (!require(Matrix))   stop("the Matrix package is required.")
   if (!require(geometry)) stop("the geometry package is required.")
   
   
   if(pcrv.path) { 
   	if(missing(start.pc)) stop("Need to create a start for principal curve")
		ind <- !(shape.smooth$kappa1 == 0 & shape.smooth$kappa2 ==0)
		shape.smooth2 <- shape.smooth$coords[ind, ]
		shape.smooth2 <- subset.face3d(shape.smooth, ind, remove.singles = FALSE)

		chR <- lmk1 
		chL <- lmk2
		u <- (chL-chR)/norm(chL-chR)
		crds <- sweep(shape.smooth2$coords, 2, chL)
		v <- c(crds%*%u)
		shape.smooth2 <- subset.face3d(shape.smooth2, v < 1, remove.singles = FALSE)
		crds <- crds[v<1,]
		for( i in 1:nrow(crds)) crds[i,] <- crds[i,]+chL 
		u <- (chR-chL)/norm(chR-chL)
		crds <- sweep(crds, 2, chR)
		v <- c(crds%*%u)
		shape.smooth2 <- subset.face3d(shape.smooth2, v < 1, remove.singles = FALSE)
		crds <- crds[v<1,]
		for( i in 1:nrow(crds)) crds[i,] <- crds[i,]+chR 
		            values       <- abs(shape.smooth2$kappa1)
		fit <- pcurve(shape.smooth2, df = c(4, 4, 4), weights = values,
		              fixed = rbind(lmk1, lmk2), initial = start.pc)
		for (i in 1:nrow(fit$s)) fit$s[i, ] <- closest.face3d(fit$s[i, ], shape.smooth)$point
		values       <- pmax(abs(shape$kappa1), abs(shape$kappa2))
		curve <- fit$s[fit$tag,]
		curve <- unique(curve)	
		shape.smooth <- shape 
	}   else  { 
  
# Planepath curve 
   if (!any(is.na(shape.smooth)) & !any(is.na(reference.path))){
   	           curve        <- reference.path
   	           values       <- pmax(abs(shape.smooth$kappa1), abs(shape.smooth$kappa2))
    } else  if (!any(is.na(reference.path)) & ( "shape.index" %in% names(shape))){
   	           curve        <- reference.path
   	           values       <- pmax(abs(shape$kappa1), abs(shape$kappa2))
   	           shape.smooth <- shape  	          
    } else if (!any(is.na(shape.smooth))) { 
    	           values       <- pmax(abs(shape.smooth$kappa1), abs(shape.smooth$kappa2))
    	           curve        <- planepath.face3d(shape, lmk1, lmk2, distance = distance, si.target = 
                                         si.target, boundary = boundary, graphics = graphics)$path                                          
    } else if ( "shape.index" %in% names(shape)) {
               values       <- pmax(abs(shape$kappa1), abs(shape$kappa2))
               curve        <- planepath.face3d(shape, lmk1, lmk2, distance = distance, si.target = 
                                         si.target,  boundary = boundary, graphics = graphics)$path
               shape.smooth <- shape   
                            
    } else {
               pp           <- planepath.face3d(shape, lmk1, lmk2, distance = distance, rotation.range=rotation.range, si.target =  
                                      si.target, boundary = boundary, graphics = graphics)
               curve        <- pp$path
               values       <- pp$values
               shape        <- pp$shape
               shape.smooth <- shape
               
    }
 }
 
 
 

   
    if(graphics) 	display.face3d(shape.smooth, colour = shape.smooth$shape.index)
    if(graphics)    display.face3d(shape.smooth, colour = values)   
    if(graphics){
       display.face3d(shape.smooth, colour=values)
       spheres3d(curve)
     }
     
#Project onto 2d surface and readying for smoothing procedure
   ccurve  <- closestcurve.face3d(shape.smooth, curve)
   lambda  <- cumsum(c(0, sqrt(apply((diff(curve))^2, 1, sum))))
   xcoord  <- lambda[ccurve$closest.curvept] # distance along the curve
   ycoord  <- ccurve$closest.distance
     
   if(graphics){
   	   spheres3d(rbind(lmk1, lmk2), col = "red", radius = 1)
       spheres3d(shape.smooth$coords, col = topo.colors(20)[cut(values, 20, labels = FALSE)], 
                                 radius = 0.3, alpha = 0.6)
       spheres3d(curve, col="black", radius = 1 )
       }
      
   #Removing excess perp.dist pts
   ind     <- abs(ycoord) < perp.dist.bound
   #Removing excess side pts
   lowend  <- (xcoord > 1e-8 )
   highend <- (xcoord < max(xcoord) - 1e-8)
   ind     <- ind & lowend & highend
   xc      <- xcoord[ind]
   yc      <- ycoord[ind]
   va      <- values[ind]
  
   #Making shape only the coordinates that are smoothed for tps
   shape.smooth   <- subset.face3d(shape.smooth, ind, remove.singles=FALSE)
        
   #adding back on landmarks
   xc  <- c(0, xc, max(lambda))
   yc  <- c(0, yc, 0)
   va1 <- c(0, va, 0)
   if (NA %in% va1){ 
   	ind  <- which(va1 != "NA")
   	va1  <- va1[ind]
   	xc   <- xc[ind]
   	yc   <- yc[ind]}
   	
#Smoothing step 
 if(!missing(smooth.method)) {	
   if (smooth.method =="quad"){
     result.smquad <- smquad.face3d(xc, yc, va1, df = df,
                                  fixed = matrix(c(0, max(lambda), 0, 0), nrow = 2),
                                  eval.points = seq(0, max(lambda), length = npts),
                                  display="none")    
                                  
     if(graphics){ plot(result$x, result$pdist, col = topo.colors(20)[cut(result$y, 20, labels = FALSE)])
                  abline(h=0)
                  lines(result$eval.points, result$estimate, col = "red")}
     new.pdist    <- result$estimate
     arclength    <- result$eval.points
     }                             
 
    else if (smooth.method =="pcurve"){  
  	 result.pcurve  <- pcurve.face3d(shape = cbind(xc,yc), df=c(df,df), weights = va1,
                                   fixed = matrix(c(0, 0, max(xc), 0), ncol=2, byrow=T),
                                    initial = cbind(seq(min(xc), max(xc),
                                    length = (length(xc)+ 2)), rep(0, (length(xc)+ 2))))  					  
     curve.pc    <- result.pcurve$s
     curve.pc    <- curve.pc[order(curve.pc[,1]),]

     if(graphics){ plot(xc, yc, col = topo.colors(20)[cut(va1, 20, labels = FALSE)], 
                   ylab= "Perpindicular Distance to the Reference Curve", 
                   xlab="Arclength of Reference Curve", main="2-Dimensional subspace")}
     abline(h=0)
     lines(curve.pc, col = "red")
     new.pdist     <- curve.pc[,2]
     arclength     <- curve.pc[,1]
                                        
     }} 
     
 #uses sm.psp as a default for smoothing method
      best.lambda <- sm.map(xc, yc, weights = va1, fixed = matrix(c(0, max(xc), 0, 0), nrow = 2))$best.lambda
      result.psp     <- sm.psp(xc, yc, weights = va1, df=df, lambda = best.lambda,
   						fixed = matrix(c(0, max(lambda), 0, 0), nrow = 2), 
   						eval.points = seq(0, max(lambda), length = npts), display= "none") 	
   				
      if(graphics){ plot(xc, yc, col = topo.colors(20)[cut(va1, 20, labels = FALSE)],
                    ylab= "Perpindicular Distance to the Reference Curve", xlab="Arclength of   
                    Reference Curve", main="2-Dimensional subspace")
                    abline(h=0)
                    lines(result.psp$eval.points, result.psp$estimate, col = "red")}
      new.pdist     <- result.psp$estimate
      arclength     <- result.psp$eval.points

     




   						
     
# Project back onto the surface full pts (interp bary)
   full.ccurve   <- closestcurve.face3d(shape, curve)
   lambda        <- cumsum(c(0, sqrt(apply((diff(curve))^2, 1, sum))))
   full.xc       <- lambda[full.ccurve$closest.curvept] # distance along the curve
       full.xc   <- c(0, full.xc, max(lambda))
   full.yc       <- full.ccurve$closest.distance
       full.yc   <- c(0, full.yc, 0)
   full.va1      <- c(0, pmax(abs(shape$kappa1), abs(shape$kappa2)),0)
   
   X             <- cbind(full.xc, full.yc)
   Xnew          <- cbind(arclength, new.pdist)
   interp.x      <- interp.barycentric(X, f = c(lmk1[1], shape$coords[ , 1], lmk2[1]), Xnew)$fnew
   interp.y      <- interp.barycentric(X, f = c(lmk1[2], shape$coords[ , 2], lmk2[2]), Xnew)$fnew
   interp.z      <- interp.barycentric(X, f = c(lmk1[3], shape$coords[ , 3], lmk2[3]), Xnew)$fnew
   interp.values <- interp.barycentric(X, f = full.va1, Xnew)$fnew
   smooth.path   <- cbind(interp.x, interp.y, interp.z)
        
  if(graphics){
	  display.face3d(face)
	  spheres3d(smooth.path, radius=1)
  }

  if(graphics){
	   display.face3d(shape, colour="grey")
   	   spheres3d(rbind(lmk1, lmk2), col = "red", radius = 1)
       spheres3d(shape.smooth$coords, col = topo.colors(20)[cut(values, 20, labels = FALSE)], 
                                 radius = 0.3, alpha = 0.6)
       spheres3d(smooth.path, col="black", radius = .5 )
       }

    invisible(list(path = smooth.path, values = interp.values))
    	}
   	
   	
   	
   	
