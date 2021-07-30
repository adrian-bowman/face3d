smoothpath.ridge.face3d <- function(shape, shape.smooth, lmk1, lmk2, df = 5, distance = 5, curve.name = NA, penalty = .002,
                              reference.path = NA, rotation.range=rotation.range, si.target, boundary = c(0.2, 0.5), npts, monitor = FALSE,        
                              perp.dist.bound = 10, pcrv.path = FALSE, start.pc = NA, smooth.method = NULL, average.alpha = average.alpha) {  
	
   if (!requireNamespace("fields"))   stop("the fields package is required.")
   if (!requireNamespace("rgl"))      stop("the rgl package is required.")
   if (!requireNamespace("Matrix"))   stop("the Matrix package is required.")
   if (!requireNamespace("geometry")) stop("the geometry package is required.")
   
   if(pcrv.path) { 
   	if(missing(start.pc)) stop("Need to create a start for principal curve")
		ind <- !(shape.smooth$kappa1 == 0 & shape.smooth$kappa2 ==0)
		shape.smooth2 <- shape.smooth$vertices[ind, ]
		shape.smooth2 <- subset.face3d(shape.smooth, ind, remove.singles = FALSE)

		chR <- lmk1 
		chL <- lmk2
		u <- (chL-chR)/norm(chL-chR)
		crds <- sweep(shape.smooth2$vertices, 2, chL)
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
		for (i in 1:nrow(fit$s)) fit$s[i, ] <- closest.face3d(fit$s[i, ], shape.smooth)$points
		values       <- pmax(abs(shape$kappa1), abs(shape$kappa2))
		curve <- fit$s[fit$tag,]
		curve <- unique(curve)	
		shape.smooth <- shape 
	}   else  { 
  
# Planepath curve (!any(is.na(shape.smooth))

   if (!any(is.na(shape.smooth)) & !any(is.na(reference.path))){
   	           curve        <- reference.path
   	           if(sign(si.target) == -1){
                  values <- pmax(shape.smooth$kappa1, 0) #valley
               } else {
                	values <- pmax(-shape.smooth$kappa2, 0) #ridge
               }

    } else  if (!any(is.na(reference.path)) & ( "shape.index" %in% names(shape))){
   	           curve        <- reference.path
   	          if(sign(si.target) == -1){
                 values <- pmax(shape$kappa1, 0) #valley
              } else {
      	         values <- pmax(-shape$kappa2, 0) #ridge
              }

   	           shape.smooth <- shape  
   	          
   	             	            	          
    } else if ( ("shape.index" %in% names(shape.smooth)) & (is.na(reference.path))) { 
              if(sign(si.target) == -1){
                 values <- pmax(shape.smooth$kappa1, 0) #valley
              } else {
      	         values <- pmax(-shape.smooth$kappa2, 0) #ridge
              }    	          
               curve        <- planepath(shape, lmk1, lmk2, distance = distance, si.target = 
                                         si.target, boundary = boundary, monitor = monitor)$path
    	          
                                        
    } else if ( "shape.index" %in% names(shape)) {
                if(sign(si.target) == -1){
                 values <- pmax(shape$kappa1, 0) #valley
              } else {
      	         values <- pmax(-shape$kappa2, 0) #ridge
              }
               curve        <- planepath(shape, lmk1, lmk2, distance = distance, si.target = 
                                         si.target,  boundary = boundary, monitor = monitor)$path
               shape.smooth <- shape 
    } else {
               pp           <- planepath(shape, lmk1, lmk2, distance = distance, rotation.range=rotation.range, si.target =  
                                      si.target, boundary = boundary, monitor = monitor)
               curve        <- pp$path
               values       <- pp$values
               shape        <- pp$shape
               shape.smooth <- shape
    }
 }
 
   
    if(monitor) 	plot(shape.smooth, colour = shape.smooth$shape.index)
    if(monitor)  plot(shape.smooth, colour = values)   
    if(monitor){
       plot(shape.smooth, colour=values)
       rgl::spheres3d(curve)
     }
   
 #save(shape, lmk1, lmk2, shape.smooth, curve, values file="/Users/liberty/Desktop/data.dmp")         
#Project onto 2d surface and readying for smoothing procedure
   ccurve  <- closestcurve.face3d(shape.smooth, curve)
   lambda  <- cumsum(c(0, sqrt(apply((diff(curve))^2, 1, sum))))
   xcoord  <- lambda[ccurve$closest.curvept] # distance along the curve
   ycoord  <- ccurve$closest.distance
   si      <- shape.smooth$shape.index
   values[sign(si)!= sign(si.target)] <- 0 
   ind      <- (abs(ycoord) < perp.dist.bound) & (xcoord > 0.005 * max(xcoord)) & (xcoord < 0.995 * max(xcoord))
   ycoord   <- c(0, ycoord[ind], 0)
   xcoord   <- c(0, xcoord[ind], max(xcoord))
   values   <- c(0, values[ind], 0)
   
   xcoord <- xcoord
   ycoord <- ycoord
   if (monitor){clr       <- topo.colors(20)[cut(values, 20, labels = FALSE)]
   plot(xcoord, ycoord, col = clr)
   abline(h = 0)}
   
   # save(xcoord,ycoord,values,file="ridge_issue.dmp")
   
   ridge           <- ridge2d.face3d(xcoord, ycoord, values, lambda = penalty, monitor = monitor)
   shape.smooth    <- subset.face3d(shape.smooth, ind, remove.singles=FALSE)
   arclength       <- ridge$x 
   new.pdist       <- ridge$y

   # Project back onto the surface full pts (interp bary)
   full.ccurve   <- closestcurve.face3d(shape, curve)
   lambda        <- cumsum(c(0, sqrt(apply((diff(curve))^2, 1, sum))))
   full.xc       <- lambda[full.ccurve$closest.curvept] # distance along the curve
       full.xc   <- c(0, full.xc, max(lambda))
   full.yc       <- full.ccurve$closest.distance
       full.yc   <- c(0, full.yc, 0)
   full.va1      <- c(0, pmax(abs(shape$kappa1), abs(shape$kappa2)),0)
   
   X             <- cbind(full.xc, full.yc)
   # X             <- cbind(xc,yc)
   Xnew          <- cbind(arclength, new.pdist)
   interp.x      <- interp.barycentric(X, f = c(lmk1[1], shape$vertices[ , 1], lmk2[1]), Xnew)$fnew
   interp.y      <- interp.barycentric(X, f = c(lmk1[2], shape$vertices[ , 2], lmk2[2]), Xnew)$fnew
   interp.z      <- interp.barycentric(X, f = c(lmk1[3], shape$vertices[ , 3], lmk2[3]), Xnew)$fnew
   interp.values <- interp.barycentric(X, f = full.va1, Xnew)$fnew
   smooth.path   <- cbind(interp.x, interp.y, interp.z)
   aver.path     <- smooth.path
   #if(average==TRUE){ 
  	#if(!any(is.na(average.alpha))){
  
   #Weighting of quadratic path and optimized planepath   3d     
     # qp        <- smooth.path
     # pp        <- resample.face3d(curve, n = dim(qp)[1])  
     # lambda    <- cumsum(c(0, sqrt(apply((diff(pp))^2, 1, sum))))
     # a         <- average.alpha
     # aver.path <- matrix(nrow=dim(qp)[1], ncol=3)
     # for (i in 1:dim(qp)[1]){
        # d         <- lambda[i]/max(lambda) 
        # wt        <- exp((-.5)*(d^2/a^2)) + exp((-.5)*((1-d)^2/a^2))  
        # x         <- wt*pp[i,1] + (1-wt)*qp[i,1] 
        # y         <- wt*pp[i,2] + (1-wt)*qp[i,2]
        # z         <- wt*pp[i,3] + (1-wt)*qp[i,3]
        # aver.path[i,] <-c(x,y,z) 
      # }
    # }    else {aver.path <- smooth.path }  

   if(monitor){
	   plot(face)
	   rgl::spheres3d(aver.path, radius=1)
   }

   if(monitor){
	   plot(shape, colour="grey")
   	 rgl::spheres3d(rbind(lmk1, lmk2), col = "red", radius = 1)
       rgl::spheres3d(shape.smooth$vertices, col = topo.colors(20)[cut(values, 20, labels = FALSE)], 
                                 radius = 0.3, alpha = 0.6)
       rgl::spheres3d(aver.path, col="black", radius = .5 )
       }

   invisible(list(path = aver.path, values = interp.values))
}
