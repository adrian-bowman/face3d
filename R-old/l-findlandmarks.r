findlandmarks.face3d <- function(face, lmks, orient = TRUE, monitor = FALSE, graphics = FALSE,
                                 overwrite = FALSE, niter = 1, monitor.extra = FALSE) {

   if (missing(lmks)) lmks <- c("pn",  "acL",  "acR",  "sn",  "se",  "n", "exL", "exR", "enL", "enR",
                                "ls", "cphL", "cphR", "chL", "chR", "st",  "li",  "sl", "gn")

   if ("lmks" %in% names(face)) {
      ind <- which(is.na(match(lmks, rownames(face$lmks))))
      if (length(ind) > 0)
         face$lmks <- rbind(face$lmks, matrix(nrow = length(ind), ncol = 3, dimnames = list(lmks[ind])))
   }
   else
      face$lmks <- matrix(nrow = length(lmks), ncol = 3, dimnames = list(lmks))

   # lmks4 <- matrix(nrow = 4, ncol = 3, dimnames = list(c("pn", "enL", "enR", "se"), NULL))

   if (orient) {
      if (monitor) cat("orienting ... ")
      face <- orient.face3d(face)
      orient <- FALSE
   }
   else if (!("nearest" %in% names(face)))
      stop("orient is set to FALSE but face$nearest is not present.")

   quadlocate <- function(face, constraint = NA, si, si.within, pair = FALSE, graphics = FALSE) {
   	suppressWarnings(if(is.na(constraint)!= TRUE){
      sbst   <- subset.face3d(face, constraint)
      sbst   <- index.face3d(sbst)}  )
      if(graphics){display.face3d(sbst, colour = 2 + sbst$kappa1 * sbst$kappa2, new=FALSE)
      	            print("prior of interest curvature values on face")
      	            scan()}
                    
      ind1 <- which(abs(sbst$shape.index - si) < si.within)
      sbst  <- subset.face3d(sbst, ind1, remove.singles=FALSE)
      # crv   <- sbst$kappa1 * sbst$kappa2
      # sbst  <- subset.face3d(sbst, crv > quantile(crv, 0.9))
      crv   <- sbst$kappa1 * sbst$kappa2
      if (graphics) {display.face3d(sbst, colour = (2 + crv), new = FALSE)
                    print("prior of interest curvature values after threshold")
                    scan()}
                   
      parts <- connected.face3d(sbst)
      # ord   <- order(tapply(crv, parts, mean), decreasing = TRUE)
      # parts <- ord[parts]
      ########################
      for (i in 1:(1 + as.numeric(pair))) {
      	ind2   <- which(parts==i)
         sbst1 <- subset.face3d(sbst, parts == i) 
         # sbst1 <- index.face3d(sbst1, overwrite = TRUE)
         # sbst1  <- subset.face3d(sbst1, abs(sbst1$shape.index - si) < si.within)
         dst   <- as.matrix(dist(sbst1$coords))
         crv   <- sbst1$kappa1
         crv   <- sbst1$kappa1 + sbst1$kappa2
         crv   <- sbst1$kappa1 * sbst1$kappa2
         rss   <- apply(dst, 1, function(x) sum(lsfit(x^2, crv)$residuals^2))
         # suppressWarnings(if(is.na(constraint)!= TRUE) {
         # rss   <- apply(dst, 1, function(x) sum(lsfit(x^2, crv)$residuals^2))}
         # else {sigma <- .1
               # rss   <- (-1/(2*(sigma^2))) * apply(dst, 1, function(x) sum(lsfit(x^2, crv)$residuals^2))})
         beta  <- apply(dst, 1, function(x) lsfit(x^2, crv)$coeff[2])
         neg   <- (abs(si) > si.within)
         ind   <- if (neg) beta > 0 else beta < 0
         rss[ind] <- max(rss)
         lmk.i <- sbst1$coords[which.min(rss), ]
         lmk   <- if (i == 1) lmk.i else rbind(lmk, lmk.i)
         suppressWarnings(if (i==1) sbst.pair <- sbst1$coords
         else if (is.na(constraint)!= FALSE){
         	    sbst.pair.f <- array(0, dim=c(100,3,2))
         	    sbst.pair.f[c(1:dim(sbst.pair)[1]),,1] <- sbst.pair
         	    sbst.pair.f[c(1:dim(sbst1$coords)[1]),,2]<- sbst1$coords
         	    sbst.pair <- sbst.pair.f}
         else sbst.pair.f <- NA)
         if (graphics) {
            # display.face3d(sbst1, new = FALSE, colour = "shape index", add = (i == 2))
            # scan()
            display.face3d(sbst1, colour = 2 + crv, new = FALSE, add = (i == 2))
            spheres3d(lmk, radius = 1, col = "black")
            print("final subset with new landmark black")
            scan()}
             
            # ind <- which.min(rss)
            # plot(dst[ind, ], crv)
        
      }
      if (pair && (lmk[1, 1] < lmk[2, 1])) lmk <- lmk[2:1, ]
      invisible(list(lmk = lmk, rss = rss, sbst = sbst, indices1 =ind1, indices2 = ind2,sbst.pair=sbst.pair))
   }
   

   #---------------------------------------------------------------------------------
   #                       locate pn, enL, enR, se
   #---------------------------------------------------------------------------------

      monitor.extra <- FALSE
   
      if (any(is.na(face$lmks["pn", ]))) {
         if (monitor) cat("Locating pn ... ")
         nrst <- face$nearest
         nrst <- nrst[!is.na(nrst)] 
         nrst <- face$coords[nrst, ]
         ord  <- order(nrst[ , 3], decreasing = TRUE)
         nrst <- nrst[ord, ]
         flag <- FALSE
         i    <- 0
         while (!flag & (i < nrow(nrst))) {
            i    <- i + 1
            ind  <- apply(sweep(face$coords, 2, nrst[i, ]), 1, function(x) sqrt(sum(x^2))) < 20
            face <- index.face3d(face, subset = ind)
            crv  <- pmin(-face$kappa1[ind], -face$kappa2[ind])
            crv[face$shape.index[ind] < 0.8] <- 0
            if (length(which(crv > 0.1)) > 10) flag <- TRUE
            # display.face3d(face, colour = "shape index", new = FALSE)
         }
         if (!flag) stop("pn could not be identified.")                        
         face$lmks["pn", ] <- quadlocate(face, edist.face3d(face$coords, nrst[i, ]) < 30,
                                         1, 0.25, graphics = graphics)$lmk
         # face$lmks["pn", ] <- quadlocate(face, edist.face3d(face$coords, face$lmks["pn", ]) < 10,
                                         # 1, 0.25, graphics = graphics)$lmk
      }
      
      if (any(is.na(c(face$lmks[c("enL", "enR"), ])))) {
         if (monitor) cat("enL/R ... ")
         face <- index.face3d(face, subset = 
           (edist.face3d(face$coords, face$lmks["pn", ]) < 70) & (face$coords[ , 2] > face$lmks["pn", 2] + 10))
         face$lmks[c("enL", "enR"), ] <- quadlocate(face,
                 (edist.face3d(face$coords, face$lmks["pn", ]) < 70) & (face$coords[ , 2] > face$lmks["pn", 2] + 10),
                 -1, 0.3, pair = TRUE, graphics = graphics)$lmk
         # face$lmks["enL", ] <- quadlocate(face, edist.face3d(face$coords, face$lmks["enL", ]) < 20,
                 # -1, 0.25, graphics = graphics)$lmk
         # face$lmks["enR", ] <- quadlocate(face, edist.face3d(face$coords, face$lmks["enR", ]) < 20,
                 # -1, 0.25, graphics = graphics)$lmk
      }
      
      if (any(is.na(face$lmks["se", ]))) {
         if (monitor) cat("se ... ")
         curve <- planepath.face3d(face, face$lmks["enL", ], face$lmks["enR", ], boundary = c(0.2, 2))$path
         gcrv  <- gcurvature.face3d(curve, 4)
         ind   <- which.max(gcrv$gcurvature)
         face$lmks["se", ] <- gcrv$resampled.curve[ind, ]
         if (monitor.extra) {
            display.face3d(face, new = FALSE)
            spheres3d(gcrv$resampled.curve, col = "green")
            spheres3d(face$lmks["se", ], radius = 1.5, col = "blue")
            scan()
         }
         face <- index.face3d(face, subset = 
           edist.face3d(face$coords, face$lmks["se", ]) < 10)
         
         # se is a saddlepoint but the negative curvature will be stronger so go from -0.25 to ridge (0.5).
         face$lmks["se", ] <- quadlocate(face, edist.face3d(face$coords, face$lmks["se", ]) < 10,
                                                    0.25, 0.35, graphics = graphics)$lmk
         # face$lmks["se", ] <- quadlocate(face, edist.face3d(face$coords, face$lmks["se", ]) < 10,
                                                    # 0.125, 0.375, graphics = graphics)$lmk
      }
      
      if (monitor) cat("completed.\n")
      
   #return(invisible(face))
#s}

   #---------------------------------------------------------------------------------
   #                       locate acL, acR
   #---------------------------------------------------------------------------------
   
         # ind   <- apply(face$coords, 1, function(x) sqrt(sum((x - face$lmks["pn", ])^2))) < 70
         # face  <- index.face3d(face, distance = 10, subset = ind)
         # ind   <- (face$shape.index < -0.25) & abs(face$coords[ , 2] - face$lmks["pn", 2]) < 20 &
                       # apply(face$coords, 1, function(x) sqrt(sum((x - face$lmks["pn", ])^2))) < 50
         # sbst  <- subset.face3d(face, ind)
         # parts <- connected.face3d(sbst)
         # ind   <- as.numeric(rownames(sbst$coords)[parts %in% 1:2])
         # sbst  <- subset.face3d(sbst, parts %in% 1:2)
         # crv   <- pmin(sbst$kappa1, sbst$kappa2)
         # brks  <- seq(0.0, 0.2, length = 21)
         # brks  <- c(min(crv) - 1, quantile(crv, seq(0.05, 0.95, 0.05)), max(crv) + 1)
         # parts <- parts[parts %in% 1:2]
         # sbst1 <- subset.face3d(sbst, parts == 1)
         # crv1  <- crv[parts == 1]
         # alL   <- sbst1$coords[which.max(crv1), ]
         # sbst2 <- subset.face3d(sbst, parts == 2)
         # crv2  <- crv[parts == 2]
         # alR   <- sbst2$coords[which.max(crv2), ]
         # al    <- if (alL[1] > alR[1]) cbind(alL, alR) else cbind(alR, alL)
         
         # ind   <- apply(face$coords, 1, function(x) sqrt(sum((x - face$lmks["pn", ])^2))) < 40
         # face  <- index.face3d(face, distance = 10, subset = ind)
         # ind   <- ind & (face$shape.index > 0.25)
         # sbst  <- subset.face3d(face, ind)
         # sbst  <- subset(sbst, connected.face3d(sbst) == 1)

         # if ("acL" %in% lmks) face$lmks["acL", ] <- al[ , 1]
         # if ("acR" %in% lmks) face$lmks["acR", ] <- al[ , 2]


   #---------------------------------------------------------------------------------
   #                       locate the others
   #---------------------------------------------------------------------------------
   #make sure prior.dmp is loaded
   #load("/Volumes/LACIE SHARE/to work on april19/R codes/bayes/bayes lmks/prior.dmp")

   # if (graphics) {for (i in 1:n) spheres3d(glmks[ , , i], radius = 1, col = 1:k)
                 # spheres3d(gpa$mean, radius = 1)}
   #reorder lmks to be same order that choose in::

   reorder         <-  c("pn", "se", "enL", "enR", "sn", "acL",  "acR", "n", "exL", "exR", 
                   "ls", "cphL", "cphR", "chL", "chR", "st",  "li",  "sl", "gn")   
   reorder.num     <- c(1, 5, 9, 10, 4, 2, 3, 6, 7, 8, 11:19 )                 
   glmks           <- gpa$rotated
   rownames(glmks) <- c("pn",  "acL",  "acR",  "sn",  "se",  "n", "exL", "exR", "enL", "enR",
                                "ls", "cphL", "cphR", "chL", "chR", "st",  "li",  "sl", "gn")
   glmks           <- glmks[reorder.num, ,] 
   gmean           <- gpa$mean
   gmean           <- gmean[reorder.num, ]
   face$lmks       <- face$lmks[reorder.num, ]
   gmean.mean <- apply(gmean, 2, mean)
   k          <- dim(glmks)[1]
   n          <- dim(glmks)[3]
   X          <- matrix(c(aperm(sweep(glmks, 1:2, gmean), c(2, 1, 3))), nrow = n, byrow = TRUE)
   covmat     <- cov(X)
   if (!require(shapes)) stop("the shapes package is required.")
   g4          <- gmean[c("pn","se", "enL", "enR" ), ]
   g4mean      <- apply(g4, 2, mean)
   #covmat4     <- cov(matrix(c(aperm(sweep(glmks[c(1,5,9,10),,], 1:2, 
   #                        gmean[c(1,5,9,10),]), c(2, 1, 3))), nrow = n, byrow = TRUE))
    covmat4     <- cov(matrix(c(aperm(sweep(glmks[c(1:4),,], 1:2, 
                           gmean[c(1:4),]), c(2, 1, 3))), nrow = n, byrow = TRUE))
                        
     
      # if (graphics){
       #lmks4             <- face$lmks[c("pn", "se", "enL", "enR"), ]   
       #lmean             <- apply(lmks4, 2, mean)
       #opa               <- procOPA(g4, lmks4, scale = FALSE)
       #face1$lmks        <- sweep(sweep(face$lmks,   2, lmean) %*% opa$R, 2, g4mean, "+")
       #face1$coords      <- sweep(sweep(face$coords, 2, lmean) %*% opa$R, 2, g4mean, "+")
      # display.face3d(face1, new = FALSE)
                  # for (i in 1:k) {
                        # ind      <- (i - 1) * 3 + 1:3
                         # covmat.i <- covmat[ind, ind]
                         # ell      <- ellipse3d(covmat.i, centre = gmean[i, ], level = 0.99)
                         # plot3d(ell, add = TRUE, col = i, alpha = 0.5)
                      # }
                  # }
                # scan()
                  
# Now bring in prior information on the four known landmarks 
   #for pn/se
   si        <- c(1,.25, -1, -1)
   si.within <- c(.25,.35,.3, .3)
   pair      <- c(FALSE, FALSE, FALSE, FALSE)
   stddev    <- c(10,10,10, 10)
   
        
for ( i in 1:4) {
   for (iter in 1:niter) { 
    	#i <- 2   #just doing sellion at the moment
    	  # doing opa on "new placement"
      lmks4             <- face$lmks[c("pn", "se", "enL", "enR"), ]   
      id                <- rownames(lmks4)[i]
      lmean             <- apply(lmks4, 2, mean)
      opa               <- procOPA(g4, lmks4, scale = FALSE)
      face$lmks         <- sweep(sweep(face$lmks,   2, lmean) %*% opa$R, 2, g4mean, "+")
      face$coords       <- sweep(sweep(face$coords, 2, lmean) %*% opa$R, 2, g4mean, "+")
      
      #picking out 4 landmarks
      lmks.cond         <- face$lmks[c("pn", "se", "enL", "enR"), ]   
      if(graphics){display.face3d(face, new=FALSE)
      	           spheres3d(face$lmks)
      	           print("found landmarks")
      	           scan()}
      	           
      lmk               <- lmks.cond[i, ]
      ind               <- (i - 1) * 3 + 1:3
      covmat.i          <- covmat4[ind, ind]                         
      cond.cov          <- covmat.i - covmat4[ind, -ind] %*% solve(covmat4[-ind,-ind]) %*%covmat4[-ind, ind]   
      cond.mean         <-  c(t(g4))[ind] +  covmat4[ind, -ind] %*% solve(covmat4[-ind,-ind]) %*% 
                            (c(t(lmks.cond[-i,]))  - c(t(g4))[-ind] ) 
      crds              <- sweep(face$coords, 2, cond.mean)
      M.dst             <- rowSums((crds %*% solve(cond.cov)) * crds)
      crds              <- sweep(crds, 2, cond.mean, "+")
      ind               <- (M.dst <= stddev[i])
      logprior          <- -M.dst[ind]  
      face              <- index.face3d(face, subset = ind)  
      sbst              <- subset.face3d(face, ind)
     if(graphics){display.face3d(face, new=FALSE)
      	           spheres3d(sbst$coords)
      	           print("prior of interest")
      	           scan()}
     if(graphics){display.face3d(sbst,colour="shape index", new=FALSE)
     	         print("shape index of prior")
     	         scan() }
      	           	           
      ql                <- quadlocate(face=sbst, si = si[i], si.within = si.within[i], pair = pair[i],
                           constraint=NA, graphics = FALSE)
      loglik            <- ql$rss
      if (graphics){display.face3d(face, new=FALSE)
      	           spheres3d(ql$sbst.pair, col = loglik)
      	           print("loglik")
      	           scan()}
      	           
      logprior          <- logprior[ql$indices2]  
            if (graphics){display.face3d(face, new=FALSE)
      	           spheres3d(ql$sbst.pair, col = -logprior)
      	                 	           print("logprior")

      	           scan()}
      	                                                     
      logpost           <- loglik + logprior
            if (graphics){display.face3d(face, new=FALSE)
      	           spheres3d(ql$sbst.pair, col = -logpost)
      	                 	           print("logpost")
      	           scan()}
      	           
      new.cond.lmk      <- ql$sbst.pair[as.numeric(names(which.min(logpost))), ]        
      if (graphics){       display.face3d(face, new = FALSE)
                           spheres3d(face$lmks)  
                           spheres3d(new.cond.lmk, radius=.5,col="red")
                           print("red is new landmark")
                           scan() 
                    }
                     
        } 
        
        face$lmks[id,] <- new.cond.lmk
     
     }
      
 ## Now that first four landmarks have stabilized, begin the next landmark sequentially  

    #IF k =6 need it to be si.within of .25 or it goes off. 
    #if K=7 , if .25 then the parts are unconnected and the true landmark gets dumped. 
    
      si          <- c(-.25,    -1,   -1,   1,  -.8,  -.8,   .8,   .8,   .8,   -.8,  -.8,   -.8,   .8,  -.7,   .7)
      si.within   <- c(  .3,   .25,  .3,  .3,   .3,   .3,   .3,   .2,   .2,    .2,   .2,    .3,   .3,   .3,   .3)
      stdd        <- c(  20,    20,   20,  10,   10,   10,   20,   20,   20,    70,   70,    50,   20,   20,   20)
      #stddev      <- stdd*15
      tt<-1
      
     for (k in 5: dim(gmean)[1]) { 
    
     lmk               <- gmean[k, ]
     ind               <- (k - 1) * 3 + 1:3 
     covmat.i          <- covmat[ind, ind]                
     crds              <- sweep(face$coords, 2, lmk)
     dst               <- rowSums((crds %*% solve(covmat.i)) * crds)
     crds              <- sweep(crds, 2, lmk, "+")
     ind               <- (dst <= stdd[tt])
     sbst              <- subset.face3d(face, ind)
     sbst              <- index.face3d(sbst)
     ql                <- quadlocate(sbst, rep(TRUE, nrow(sbst$coords)), 
     					  si= si[tt], si.within= si.within[tt], pair = FALSE, graphics = FALSE)
     logprior          <- -dst
     logprior          <- logprior[ql$indices2]
     loglik            <- ql$rss
     logpost           <- loglik + logprior
     id                <- rownames(gmean)[k] 
     face$lmks[id, ]   <- ql$sbst.pair[as.numeric(which.min(logpost)), ]  
     
     
      lmks.new             <- face$lmks[c(1:k), ]   
      lmean             <- apply(lmks.new, 2, mean)
      opa               <- procOPA(gmean[c(1:k), ], lmks.new, scale = FALSE)
      face$lmks         <- sweep(sweep(face$lmks,   2, lmean) %*% opa$R, 2, (apply(gmean[c(1:k), ], 2, mean)), "+")
      face$coords       <- sweep(sweep(face$coords, 2, lmean) %*% opa$R, 2, (apply(gmean[c(1:k), ], 2, mean)), "+")
tt<- tt+1
     }
     
     
     
     
     
  
  
  
  
  
  
  
      invisible(face)
      }
  
     
     
     
     
     
   # tt <- 2  
  # #for (k in lmks.undone) { 
   # for (k in 6: dim(gmean)[1]) {
   	# g.new        <- gmean[c(1:(k-1)), ]
     # g.new.mean   <- apply(g.new, 2, mean)
     # covmat.new   <- cov(matrix(c(aperm(sweep(glmks[c(1:(k-1)),,], 1:2, 
                           # gmean[c(1:(k-1)),]), c(2, 1, 3))), nrow = n, byrow = TRUE))
     # lmks.new           <- rbind(face$lmks[c(1:(k-1)), ])
     # id                 <- rownames(lmks.new)[(k-1)]
     # lmean              <- apply(lmks.new, 2, mean)
     # opa                <- procOPA(g.new, lmks.new, scale = FALSE)
     # face$lmks          <- sweep(sweep(face$lmks,   2, lmean) %*% opa$R, 2, g.new.mean, "+")
     # face$coords        <- sweep(sweep(face$coords, 2, lmean) %*% opa$R, 2, g.new.mean, "+")

     # lmk               <- gmean[k, ]
     # ind               <- (k - 1) * 3 + 1:3 
     # covmat.i          <- covmat[ind, ind]                
     # crds              <- sweep(face$coords, 2, lmk)
     # dst               <- rowSums((crds %*% solve(covmat.i)) * crds)
     # crds              <- sweep(crds, 2, lmk, "+")
     # ind               <- (dst <= stdd[tt])
     # sbst              <- subset.face3d(face, ind)
     # sbst              <- index.face3d(sbst)
     # ql                <- quadlocate(sbst, rep(TRUE, nrow(sbst$coords)), 
     					  # si= si[tt], si.within= si.within[tt], pair = FALSE, graphics = FALSE)
     # logprior          <- -dst
     # logprior          <- logprior[ql$indices2]
     # loglik            <- ql$rss
     # logpost           <- loglik + logprior
     # id                <- rownames(gmean)[k] 
     # face$lmks[id, ]   <- ql$sbst.pair[as.numeric(which.min(logpost)), ]  
      

      
   #generalize from here
  
     # g.new        <- gmean[c(1:k), ]
     # g.new.mean   <- apply(g.new, 2, mean)
     # covmat.new   <- cov(matrix(c(aperm(sweep(glmks[c(1:k),,], 1:2, 
                           # gmean[c(1:k),]), c(2, 1, 3))), nrow = n, byrow = TRUE))
     # lmks.new           <- rbind(face$lmks[c(1:k), ])
     # id                 <- rownames(lmks.new)[k]
     # lmean              <- apply(lmks.new, 2, mean)
     # opa                <- procOPA(g.new, lmks.new, scale = FALSE)
     # face$lmks          <- sweep(sweep(face$lmks,   2, lmean) %*% opa$R, 2, g.new.mean, "+")
     # face$coords        <- sweep(sweep(face$coords, 2, lmean) %*% opa$R, 2, g.new.mean, "+")
      # lmks.cond         <- face$lmks[c(1:k), ]    	           
      # lmk               <- lmks.cond[k, ]
      # ind               <- (k - 1) * 3 + 1:3
      # covmat.i          <- covmat.new[ind, ind]                         
      # cond.cov          <- covmat.i - covmat.new[ind, -ind] %*% solve(covmat.new[-ind,-ind]) %*%covmat.new[-ind, ind]   
      # cond.mean         <-  c(t(g.new))[ind] +  covmat.new[ind, -ind] %*% solve(covmat.new[-ind,-ind]) %*% 
                            # (c(t(lmks.cond[-i,]))  - c(t(g.new))[-ind] ) 
      # crds              <- sweep(face$coords, 2, cond.mean)
       # M.dst             <- rowSums((crds %*% solve(cond.cov)) * crds)
       # crds              <- sweep(crds, 2, cond.mean, "+")
       # ind               <- (M.dst <= stddev[tt])
      # logprior          <- -M.dst[ind]  
      # face              <- index.face3d(face, subset = ind)  
      # sbst              <- subset.face3d(face, ind)           	           
      # ql                <- quadlocate(face=sbst, si = si[tt], si.within = si.within[tt], pair = FALSE,
                           # constraint=NA, graphics = FALSE)
      # loglik            <- ql$rss    	           
      # logprior          <- logprior[ql$indices2]                                
      # logpost           <- loglik + logprior
      # new.cond.lmk      <- ql$sbst.pair[as.numeric(which.min(logpost)), ] 
      # id                <- rownames(gmean)[k]                     
      # face$lmks[id,]    <- new.cond.lmk
       # tt                 <- tt + 1
 # print(k)
      # }


  
  
  
  
  
     
      
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   

      # # if (monitor) cat("Locating sn ... ")
        # # i                 <- 4
        # # lmk               <- gmean[i, ]
        # # ind               <- (i - 1) * 3 + 1:3
        # # covmat.i          <- covmat[ind, ind]                
        # # crds              <- sweep(face$coords, 2, lmk)
        # # dst               <- rowSums((crds %*% solve(covmat.i)) * crds)
        # # crds              <- sweep(crds, 2, lmk, "+")
        # # ind               <- (dst <= 20)
        # # sbst              <- subset.face3d(face, ind)
        # # sbst              <- index.face3d(sbst)
        # # display.face3d(sbst, colour = "shape index", new = FALSE)
        # # spheres3d(lmk, col = "green")
        # # scan()
        # # display.face3d(sbst, colour = 2 + sbst$kappa1 * sbst$kappa2, new = FALSE)
        # # spheres3d(lmk, col = "green")        
        # # scan()
        # # ql <- quadlocate(sbst, rep(TRUE, nrow(sbst$coords)), 0, 0.25, pair = FALSE, graphics = FALSE)
        # # print(ql)
        # # print(nrow(sbst$coords))
        # # display.face3d(ql$sbst, colour = 2 + ql$rss, new = FALSE)
        # # spheres3d(lmk, col = "green")        
        # # spheres3d(ql$lmk, col = "red")        
        # # scan()
        # # # plot3d(ellipse3d(covmat.i, lmk, level = 0.99), add = TRUE, col = "blue", alpha = 0.5)
        # # logprior          <- -dst
        # # gc                <- face$kappa1 * face$kappa2
        # # loglik            <- log(gc[which.max(gc)])
        # # logpost.i         <- loglik + logprior
        # # face$lmks["sn", ] <- crds[which.max(logpost.i), ]                
      
   





    
# #      for (iter in 1:niter){ 
   # # #do pair  - find each prior separately, colour subset only by mahal distance. 
      # covmat4           <- cov(matrix(c(aperm(sweep(glmks[c(1,5,9,10),,], 1:2, 
                           # gmean[c(1,5,9,10),]), c(2, 1, 3))), nrow = n, byrow = TRUE))
     # lmks.cond         <- face$lmks[c("pn", "se", "enL", "enR"), ]   
      # lmk               <- lmks.cond[3, ]
      # ind               <- (3 - 1) * 3 + 1:3
      # covmat.i          <- covmat4[ind, ind]                         
      # cond.cov          <- covmat.i - covmat4[ind, -ind] %*% solve(covmat4[-ind,-ind]) %*%covmat4[-ind, ind]   
      # cond.mean         <-  c(t(g4))[ind] +  covmat4[ind, -ind] %*% solve(covmat4[-ind,-ind])%*%   
                            # (c(t(lmks.cond[-3,]))  - c(t(g4))[-ind] )   
      # crds              <- sweep(face$coords, 2, cond.mean)
      # M.dst             <- rowSums((crds %*% solve(cond.cov)) * crds)
      # crds              <- sweep(crds, 2, cond.mean, "+")
      # ind               <- (M.dst <= 20)
      # logprior          <- -M.dst[ind]
      # M.dst.L           <- M.dst
      # ind.L             <- ind
      # logprior.L        <- logprior
      # lmk               <- lmks.cond[4, ]
      # ind               <- (4 - 1) * 3 + 1:3
      # covmat.i          <- covmat4[ind, ind]                         
      # cond.cov          <- covmat.i - covmat4[ind, -ind] %*% solve(covmat4[-ind,-ind]) %*%covmat4[-ind, ind]   
      # cond.mean         <-  c(t(g4))[ind] +  covmat4[ind, -ind] %*% solve(covmat4[-ind,-ind]) %*% 
                            # (c(t(lmks.cond[-4,]))  - c(t(g4))[-ind] )   
      # crds              <- sweep(face$coords, 2, cond.mean)
      # M.dst             <- rowSums((crds %*% solve(cond.cov)) * crds)
      # crds              <- sweep(crds, 2, cond.mean, "+")
      # ind               <- (M.dst <= 20)
      # logprior          <- -M.dst[ind]
      # M.dst.R           <- M.dst
      # ind.R             <- ind
      # logprior.R        <- logprior

      # sbst                                <- subset.face3d(face,(edist.face3d(face$coords, 
                                             # face$lmks["pn", ]) < 70) & (face$coords[ , 2] > face$lmks["pn", 2] + 10) ) #initial subset on area
      # sbst                                <- index.face3d(sbst)    
      # ind.final                           <- c( which(ind.L== "TRUE"), which(ind.R==TRUE)) 
      # ind.good                            <- match( ind.final,rownames(sbst$coords))
      # sbst$shape.index[-ind.good]         <- 1
      # sbst$kappa1[-ind.good]              <- -1       
      # sbst$kappa2[-ind.good]              <- -1 
      # ql                                  <- quadlocate(face=sbst, si = -1, si.within = 0.3, 
                                              # pair = TRUE, constraint= NA, graphics = FALSE)
      # loglik                              <- ql$rss
      # logpri.L                            <- logprior.L[ql$indices2]
      # logpri.R                            <- logprior.R[ql$indices2]
      # logpost.L                           <- logpri.L +loglik
      # logpost.R                           <- logpri.R +loglik 
      # new.cond.lmk.R                      <- ql$sbst.pair[,,1][which.max(logpost.R), ]   
      # new.cond.lmk.L                      <- ql$sbst.pair[,,2][which.max(logpost.L), ] 
      # if (graphics){                       display.face3d(face, new=FALSE)
                                           # spheres3d(face$lmks)  
                                           # spheres3d(lmks.cond, col="red") 	
                   # } 
# face$lmks[c("enL", "enR"),] <- rbind(new.cond.lmk.L, new.cond.lmk.R)

  # } 
    
         
               
                
 # # Now that first four landmarks have stabilized, begin the next landmark sequentially  


   # # #locating the sn
   # # monitor.extra <- FALSE
   
      # # if (monitor) cat("Locating sn ... ")
        # # i                 <- 4
        # # lmk               <- gpa$mean[i, ]
        # # ind               <- (i - 1) * 3 + 1:3
        # # covmat.i          <- covmat[ind, ind]                
        # # crds              <- sweep(face$coords, 2, lmk)
        # # dst               <- rowSums((crds %*% solve(covmat.i)) * crds)
        # # crds              <- sweep(crds, 2, lmk, "+")
        # # ind               <- (dst <= 20)
        # # sbst              <- subset.face3d(face, ind)
        # # sbst              <- index.face3d(sbst)
        # # display.face3d(sbst, colour = "shape index", new = FALSE)
        # # spheres3d(lmk, col = "green")
        # # scan()
        # # display.face3d(sbst, colour = 2 + sbst$kappa1 * sbst$kappa2, new = FALSE)
        # # spheres3d(lmk, col = "green")        
        # # scan()
        # # ql <- quadlocate(sbst, rep(TRUE, nrow(sbst$coords)), 0, 0.25, pair = FALSE, graphics = FALSE)
        # # print(ql)
        # # print(nrow(sbst$coords))
        # # display.face3d(ql$sbst, colour = 2 + ql$rss, new = FALSE)
        # # spheres3d(lmk, col = "green")        
        # # spheres3d(ql$lmk, col = "red")        
        # # scan()
        # # # plot3d(ellipse3d(covmat.i, lmk, level = 0.99), add = TRUE, col = "blue", alpha = 0.5)
        # # logprior          <- -dst
        # # gc                <- face$kappa1 * face$kappa2
        # # loglik            <- log(gc[which.max(gc)])
        # # logpost.i         <- loglik + logprior
        # # face$lmks["sn", ] <- crds[which.max(logpost.i), ]                