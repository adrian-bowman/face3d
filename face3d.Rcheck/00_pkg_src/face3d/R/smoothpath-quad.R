



smoothpath.quad.face3d <- function(shape, shape.smooth, lmk1, lmk2, df = 5, distance = 5, curve.name = NA,
                              reference.path = NA, rotation.range=rotation.range, si.target, boundary = c(0.2, 0.5), npts, graphics = FALSE,        
                              perp.dist.bound = 10, pcrv.path = FALSE, start.pc = NA, smooth.method = NULL, average.alpha = average.alpha) {  
	
   if (!requireNamespace("fields",   quietly = TRUE)) stop("the fields package is required.")
   if (!requireNamespace("rgl",      quietly = TRUE)) stop("the rgl package is required.")
   if (!requireNamespace("Matrix",   quietly = TRUE)) stop("the Matrix package is required.")
   if (!requireNamespace("geometry", quietly = TRUE)) stop("the geometry package is required.")
   
 
   
   
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
               curve        <- planepath.face3d(shape, lmk1, lmk2, distance = distance, si.target = 
                                         si.target, boundary = boundary, graphics = graphics)$path
    	          
                                        
    } else if ( "shape.index" %in% names(shape)) {
                if(sign(si.target) == -1){
                 values <- pmax(shape$kappa1, 0) #valley
              } else {
      	         values <- pmax(-shape$kappa2, 0) #ridge
              }
               curve        <- planepath.face3d(shape, lmk1, lmk2, distance = distance, si.target = 
                                         si.target,  boundary = boundary, graphics = graphics)$path
               shape.smooth <- shape 
               print("here1")
                            
    } else {
               pp           <- planepath.face3d(shape, lmk1, lmk2, distance = distance, rotation.range=rotation.range, si.target =  
                                      si.target, boundary = boundary, graphics = graphics)
               curve        <- pp$path
               values       <- pp$values
               shape        <- pp$shape
               shape.smooth <- shape
               print("here")
               
    }
 }
 
   
    if(graphics) 	display.face3d(shape.smooth, colour = shape.smooth$shape.index)
    if(graphics)    display.face3d(shape.smooth, colour = values)   
    if(graphics){
       display.face3d(shape.smooth, colour=values)
       rgl::spheres3d(curve)
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
   
  # shape.smooth <- index.face3d(shape.smooth, distance=5)
   # values <- pmax(shape.smooth$kappa1, 0)     valley
  # values <- pmax(-shape.smooth$kappa2, 0)     #ridge
   
 # display.face3d(shape.smooth, colour="grey")
   #spheres3d(shape.smooth$coords, col = topo.colors(20)[cut(values, 20, labels = FALSE)], 
   #         radius = 0.3, alpha = 0.6)
   
   
   #Removing excess perp.dist pts
   ind     <- abs(ycoord) < perp.dist.bound
   #Removing excess side pts
   lowend  <- (xcoord > 1e-8 )
   highend <- (xcoord < max(xcoord) - 1e-8)
   ind     <- ind & lowend & highend & (values >= 0) 
   #ind     <- ind & (values > 0) 
   xc      <- xcoord[ind]
   yc      <- ycoord[ind]
   va      <- values[ind]
  
   #Making shape only the coordinates that are smoothed for tps
   shape.smooth   <- subset.face3d(shape.smooth, ind, remove.singles=FALSE)
        
   #adding back on landmarks
   #va[which(va<=0)]<- .001

   xc  <- c(0, xc, max(lambda))
   yc  <- c(0, yc, 0)
  # va1 <- c(values[1], va, values[length(values)])
     va1 <- c(.001, va, .001)
   # if (length(grep("lower lip right", curve.name)) > 0    | length(grep("upper lip right", curve.name)) > 0 ){
   # #trying to fix move away - just weighted the beginning

  
  
    # xc  <- c(seq(0,3,by=.25), xc)
    # yc  <- c(rep(0,13), yc)
    # va1 <- c(rep(.3,13), va1)
# }
   # if (length(grep("lower lip left", curve.name)) > 0    | length(grep("upper lip left", curve.name)) > 0 ){
   # #trying to fix move away - just weighted the beginning

  
  
    # xc  <- c(xc,seq(max(lambda)-6, max(lambda), by=.25))
    # yc  <- c(yc, rep(0,25))
    # va1 <- c(va1,rep(.3,25))
# }

   if (NA %in% va1){ 
   	ind  <- which(va1 != "NA")
   	va1  <- va1[ind]
   	xc   <- xc[ind]
   	yc   <- yc[ind]}
 
#Smoothing step 
#save(shape.smooth, xc,yc,va1, df, lambda, file="/Users/libertyvittert/Desktop/data.dmp")  

# display.face3d(shape.smooth, colour="grey")
# spheres3d(shape.smooth$coords, col = topo.colors(20)[cut(va1, 20, labels = FALSE)], 
 #          radius = 0.3, alpha = 0.6)
# spheres3d(aver.path, col="black", radius = .5 )
# ssss<- index.face3d(shape.smooth, distance=5)
# display.face3d(ssss, colour="grey")
# spheres3d(ssss$coords, col = topo.colors(20)[cut(pmax(abs(ssss$kappa1), abs(ssss$kappa2)), 20, labels = FALSE)], 
#           radius = 0.3, alpha = 0.6)
# 
# ccc<- pmax(abs(ssss$kappa1),0)
# display.face3d(ssss, colour=ccc*100)

#load("/Users/libertyvittert/Desktop/data.dmp")
     result <- smquad.face3d(xc, yc, va1, df = df,
                                  fixed = matrix(c(0, max(lambda), 0, 0), nrow = 2),
                                  eval.points = seq(0, max(lambda), length = 103),
                                  display="none")    
                                  
     if(graphics){ plot(result$x, result$pdist, col = topo.colors(20)[cut(result$y, 20, labels = FALSE)])
                  abline(h=0)
                  lines(result$eval.points, result$estimate, col = "red")}
     
     #log(va1[order(result$x)])
     #abs(result$curvature[order(result$x)])
     
# #Weighting of quadratic path and optimized planepath
       
     # pdist.est    <- result$estimate
     # arclength    <- result$eval.points
     # qp           <- cbind(arclength, pdist.est)
     # curvature.qp <- result$curvature
     # pp           <- cbind(result$eval.points, rep(0, length=length(result$eval.points)))
     # cs           <- max(-curvature.qp, 0)
     # ws           <- -curvature.qp/max(cs)
     # new.pdist    <- ws * qp[,2]
     # 
     # 
     # new.curvature           <- va1/max(cs)
     curv <- c(result$curvature)
     #plot(curv ~ xc)
     # # # # # NOTE curvature is always negative coming out- why??
    curv[sign(curv) >0] <- 0
    curv <- curv^4
     # curv[is.infinite(curv)]<- 0
    curvature <- curv/ mean(curv)
    # plot(curvature ~ xc)
    # plot(result$curvature ~ xc)
    # m <- length(result$beta) / 3
    # ind <- 2 * m + 1:m
    # plot(result$beta[ind])
    # cs <- result$B[ , ind] %*%result$beta[ind]
    # plot(curv ~ xc)
    
    if(graphics){ plot(result$x, result$pdist, col = topo.colors(20)[cut(curvature, 20, labels = FALSE)])
                  abline(h=0)}
    
    result <- smquad.face3d(xc, yc, va1, df = df, weights = curvature,
                                  fixed = matrix(c(0, max(lambda), 0, 0), nrow = 2),
                                  eval.points = seq(0, max(lambda), length = 103),
                                  display="none")    
    
          
   arclength    <- result$eval.points 
   new.pdist    <- result$estimate
     if(graphics){ plot(result$x, result$pdist, col = topo.colors(20)[cut(result$y, 20, labels = FALSE)])
                  abline(h=0)
                  lines(arclength, new.pdist, col = "red")}
    
             # arclength    <- result$eval.points
             # new.pdist    <- result$estimate

 


   						
     
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
   interp.x      <- interp.barycentric(X, f = c(lmk1[1], shape$coords[ , 1], lmk2[1]), Xnew)$fnew
   interp.y      <- interp.barycentric(X, f = c(lmk1[2], shape$coords[ , 2], lmk2[2]), Xnew)$fnew
   interp.z      <- interp.barycentric(X, f = c(lmk1[3], shape$coords[ , 3], lmk2[3]), Xnew)$fnew
   interp.values <- interp.barycentric(X, f = full.va1, Xnew)$fnew
   smooth.path     <- cbind(interp.x, interp.y, interp.z)
   aver.path       <- smooth.path
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
        
        
        
        
        
        
  if(graphics){
	  display.face3d(face)
	  spheres3d(aver.path, radius=1)
  }

  if(graphics){
	   display.face3d(shape, colour="grey")
   	   spheres3d(rbind(lmk1, lmk2), col = "red", radius = 1)
       spheres3d(shape.smooth$coords, col = topo.colors(20)[cut(values, 20, labels = FALSE)], 
                                 radius = 0.3, alpha = 0.6)
       spheres3d(aver.path, col="black", radius = .5 )
       }

    invisible(list(path = aver.path, values = interp.values))
    	}
   	
   	
   	
   	
