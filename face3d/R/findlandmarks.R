# findlandmarks.face3d <- function(face, lmks, orient = TRUE, monitor = FALSE, monitor = FALSE,
                                 # overwrite = FALSE, niter = 1, monitor.extra = FALSE) {

   # if (missing(lmks)) lmks <- c("pn",  "acL",  "acR",  "sn",  "se",  "n", "exL", "exR", "enL", "enR",
                                # "ls", "cphL", "cphR", "chL", "chR", "st",  "li",  "sl", "gn")

   # if ("landmarks" %in% names(face)) {
      # ind <- which(is.na(match(lmks, rownames(face$landmarks))))
      # if (length(ind) > 0)
         # face$landmarks <- rbind(face$landmarks, matrix(nrow = length(ind), ncol = 3, dimnames = list(lmks[ind])))
   # }
   # else
      # face$landmarks <- matrix(nrow = length(lmks), ncol = 3, dimnames = list(lmks))

   # # lmks4 <- matrix(nrow = 4, ncol = 3, dimnames = list(c("pn", "enL", "enR", "se"), NULL))

   # if (orient) {
      # if (monitor) cat("orienting ... ")
      # face <- orient.face3d(face)
      # orient <- FALSE
   # }
   # else if (!("nearest" %in% names(face)))
      # stop("orient is set to FALSE but face$nearest is not present.")

   # # quadlocate <- function(face, constraint = NA, si, si.within, pair = FALSE, monitor = FALSE) {
   	# # suppressWarnings(if(is.na(constraint)!= TRUE){
      # # sbst   <- subset.face3d(face, constraint)
      # # sbst   <- curvatures(sbst)}  )
      # # if(monitor){plot(sbst, colour = 2 + sbst$kappa1 * sbst$kappa2, new=FALSE)
      	            # # print("prior of interest curvature values on face")
      	            # # scan()}
                    
      # # ind1 <- which(abs(sbst$shape.index - si) < si.within)
      # # sbst  <- subset.face3d(sbst, ind1, remove.singles=FALSE)
      # # # crv   <- sbst$kappa1 * sbst$kappa2
      # # # sbst  <- subset.face3d(sbst, crv > quantile(crv, 0.9))
      # # crv   <- sbst$kappa1 * sbst$kappa2
      # # if (monitor) {plot(sbst, colour = (2 + crv), new = FALSE)
                    # # print("prior of interest curvature values after threshold")
                    # # scan()}
                   
      # # parts <- connected.face3d(sbst)
      # # # ord   <- order(tapply(crv, parts, mean), decreasing = TRUE)
      # # # parts <- ord[parts]
      # # ########################
      # # for (i in 1:(1 + as.numeric(pair))) {
      	# # ind2   <- which(parts==i)
         # # sbst1 <- subset.face3d(sbst, parts == i) 
         # # # sbst1 <- curvatures(sbst1, overwrite = TRUE)
         # # # sbst1  <- subset.face3d(sbst1, abs(sbst1$shape.index - si) < si.within)
         # # dst   <- as.matrix(dist(sbst1$vertices))
         # # crv   <- sbst1$kappa1
         # # crv   <- sbst1$kappa1 + sbst1$kappa2
         # # crv   <- sbst1$kappa1 * sbst1$kappa2
         # # rss   <- apply(dst, 1, function(x) sum(lsfit(x^2, crv)$residuals^2))
         # # # suppressWarnings(if(is.na(constraint)!= TRUE) {
         # # # rss   <- apply(dst, 1, function(x) sum(lsfit(x^2, crv)$residuals^2))}
         # # # else {sigma <- .1
               # # # rss   <- (-1/(2*(sigma^2))) * apply(dst, 1, function(x) sum(lsfit(x^2, crv)$residuals^2))})
         # # beta  <- apply(dst, 1, function(x) lsfit(x^2, crv)$coeff[2])
         # # neg   <- (abs(si) > si.within)
         # # ind   <- if (neg) beta > 0 else beta < 0
         # # rss[ind] <- max(rss)
         # # lmk.i <- sbst1$vertices[which.min(rss), ]
         # # lmk   <- if (i == 1) lmk.i else rbind(lmk, lmk.i)
         # # suppressWarnings(if (i==1) sbst.pair <- sbst1$vertices
         # # else if (is.na(constraint)!= FALSE){
         	    # # sbst.pair.f <- array(0, dim=c(100,3,2))
         	    # # sbst.pair.f[c(1:dim(sbst.pair)[1]),,1] <- sbst.pair
         	    # # sbst.pair.f[c(1:dim(sbst1$vertices)[1]),,2]<- sbst1$vertices
         	    # # sbst.pair <- sbst.pair.f}
         # # else sbst.pair.f <- NA)
         # # if (monitor) {
            # # # plot(sbst1, new = FALSE, colour = "shape index", add = (i == 2))
            # # # scan()
            # # plot(sbst1, colour = 2 + crv, new = FALSE, add = (i == 2))
            # # spheres3d(lmk, radius = 1, col = "black")
            # # print("final subset with new landmark black")
            # # scan()}
             
            # # # ind <- which.min(rss)
            # # # plot(dst[ind, ], crv)
        
      # # }
      # # if (pair && (lmk[1, 1] < lmk[2, 1])) lmk <- lmk[2:1, ]
      # # invisible(list(lmk = lmk, rss = rss, sbst = sbst, indices1 =ind1, indices2 = ind2,sbst.pair=sbst.pair))
   # # }
   

   # #---------------------------------------------------------------------------------
   # #                       locate pn, enL, enR, se
   # #---------------------------------------------------------------------------------

      # monitor.extra <- FALSE
   
      # if (any(is.na(face$landmarks["pn", ]))) {
         # if (monitor) cat("Locating pn ... ")
         # nrst <- face$nearest
         # nrst <- nrst[!is.na(nrst)] 
         # nrst <- face$vertices[nrst, ]
         # ord  <- order(nrst[ , 3], decreasing = TRUE)
         # nrst <- nrst[ord, ]
         # flag <- FALSE
         # i    <- 0
         # while (!flag & (i < nrow(nrst))) {
            # i    <- i + 1
            # ind  <- apply(sweep(face$vertices, 2, nrst[i, ]), 1, function(x) sqrt(sum(x^2))) < 20
            # face <- curvatures(face, subset = ind)
            # crv  <- pmin(-face$kappa1[ind], -face$kappa2[ind])
            # crv[face$shape.index[ind] < 0.8] <- 0
            # if (length(which(crv > 0.1)) > 10) flag <- TRUE
            # # plot(face, colour = "shape index", new = FALSE)
         # }
         # if (!flag) stop("pn could not be identified.")                        
         # face$landmarks["pn", ] <- quadlocate(face, edist.face3d(face$vertices, nrst[i, ]) < 30,
                                         # 1, 0.25, monitor = monitor)$lmk
         # # face$landmarks["pn", ] <- quadlocate(face, edist.face3d(face$vertices, face$landmarks["pn", ]) < 10,
                                         # # 1, 0.25, monitor = monitor)$lmk
      # }
      
      # if (any(is.na(c(face$landmarks[c("enL", "enR"), ])))) {
         # if (monitor) cat("enL/R ... ")
         # face <- curvatures(face, subset = 
           # (edist.face3d(face$vertices, face$landmarks["pn", ]) < 70) & (face$vertices[ , 2] > face$landmarks["pn", 2] + 10))
         # face$landmarks[c("enL", "enR"), ] <- quadlocate(face,
                 # (edist.face3d(face$vertices, face$landmarks["pn", ]) < 70) & (face$vertices[ , 2] > face$landmarks["pn", 2] + 10),
                 # -1, 0.3, pair = TRUE, monitor = monitor)$lmk
         # # face$landmarks["enL", ] <- quadlocate(face, edist.face3d(face$vertices, face$landmarks["enL", ]) < 20,
                 # # -1, 0.25, monitor = monitor)$lmk
         # # face$landmarks["enR", ] <- quadlocate(face, edist.face3d(face$vertices, face$landmarks["enR", ]) < 20,
                 # # -1, 0.25, monitor = monitor)$lmk
      # }
      
      # if (any(is.na(face$landmarks["se", ]))) {
         # if (monitor) cat("se ... ")
         # curve <- planepath(face, face$landmarks["enL", ], face$landmarks["enR", ], boundary = c(0.2, 2))$path
         # gcrv  <- gcurvature.face3d(curve, 4)
         # ind   <- which.max(gcrv$gcurvature)
         # face$landmarks["se", ] <- gcrv$resampled.curve[ind, ]
         # if (monitor.extra) {
            # plot(face, new = FALSE)
            # spheres3d(gcrv$resampled.curve, col = "green")
            # spheres3d(face$landmarks["se", ], radius = 1.5, col = "blue")
            # scan()
         # }
         # face <- curvatures(face, subset = 
           # edist.face3d(face$vertices, face$landmarks["se", ]) < 10)
         
         # # se is a saddlepoint but the negative curvature will be stronger so go from -0.25 to ridge (0.5).
         # face$landmarks["se", ] <- quadlocate(face, edist.face3d(face$vertices, face$landmarks["se", ]) < 10,
                                                    # 0.25, 0.35, monitor = monitor)$lmk
         # # face$landmarks["se", ] <- quadlocate(face, edist.face3d(face$vertices, face$landmarks["se", ]) < 10,
                                                    # # 0.125, 0.375, monitor = monitor)$lmk
      # }
      
      # if (monitor) cat("completed.\n")
      
   # #return(invisible(face))
# #s}

   # #---------------------------------------------------------------------------------
   # #                       locate acL, acR
   # #---------------------------------------------------------------------------------
   
         # # ind   <- apply(face$vertices, 1, function(x) sqrt(sum((x - face$landmarks["pn", ])^2))) < 70
         # # face  <- curvatures(face, distance = 10, subset = ind)
         # # ind   <- (face$shape.index < -0.25) & abs(face$vertices[ , 2] - face$landmarks["pn", 2]) < 20 &
                       # # apply(face$vertices, 1, function(x) sqrt(sum((x - face$landmarks["pn", ])^2))) < 50
         # # sbst  <- subset.face3d(face, ind)
         # # parts <- connected.face3d(sbst)
         # # ind   <- as.numeric(rownames(sbst$vertices)[parts %in% 1:2])
         # # sbst  <- subset.face3d(sbst, parts %in% 1:2)
         # # crv   <- pmin(sbst$kappa1, sbst$kappa2)
         # # brks  <- seq(0.0, 0.2, length = 21)
         # # brks  <- c(min(crv) - 1, quantile(crv, seq(0.05, 0.95, 0.05)), max(crv) + 1)
         # # parts <- parts[parts %in% 1:2]
         # # sbst1 <- subset.face3d(sbst, parts == 1)
         # # crv1  <- crv[parts == 1]
         # # alL   <- sbst1$vertices[which.max(crv1), ]
         # # sbst2 <- subset.face3d(sbst, parts == 2)
         # # crv2  <- crv[parts == 2]
         # # alR   <- sbst2$vertices[which.max(crv2), ]
         # # al    <- if (alL[1] > alR[1]) cbind(alL, alR) else cbind(alR, alL)
         
         # # ind   <- apply(face$vertices, 1, function(x) sqrt(sum((x - face$landmarks["pn", ])^2))) < 40
         # # face  <- curvatures(face, distance = 10, subset = ind)
         # # ind   <- ind & (face$shape.index > 0.25)
         # # sbst  <- subset.face3d(face, ind)
         # # sbst  <- subset(sbst, connected.face3d(sbst) == 1)

         # # if ("acL" %in% lmks) face$landmarks["acL", ] <- al[ , 1]
         # # if ("acR" %in% lmks) face$landmarks["acR", ] <- al[ , 2]


   # #---------------------------------------------------------------------------------
   # #                       locate the others
   # #---------------------------------------------------------------------------------
   # #make sure prior.dmp is loaded
   # #load("/Volumes/LACIE SHARE/to work on april19/R codes/bayes/bayes lmks/prior.dmp")

   # # if (monitor) {for (i in 1:n) spheres3d(glmks[ , , i], radius = 1, col = 1:k)
                 # # spheres3d(gpa$mean, radius = 1)}
   # #reorder lmks to be same order that choose in::

   # reorder         <-  c("pn", "se", "enL", "enR", "sn", "acL",  "acR", "n", "exL", "exR", 
                   # "ls", "cphL", "cphR", "chL", "chR", "st",  "li",  "sl", "gn")   
   # reorder.num     <- c(1, 5, 9, 10, 4, 2, 3, 6, 7, 8, 11:19 )                 
   # glmks           <- gpa$rotated
   # rownames(glmks) <- c("pn",  "acL",  "acR",  "sn",  "se",  "n", "exL", "exR", "enL", "enR",
                                # "ls", "cphL", "cphR", "chL", "chR", "st",  "li",  "sl", "gn")
   # glmks           <- glmks[reorder.num, ,] 
   # gmean           <- gpa$mean
   # gmean           <- gmean[reorder.num, ]
   # face$landmarks       <- face$landmarks[reorder.num, ]
   # gmean.mean <- apply(gmean, 2, mean)
   # k          <- dim(glmks)[1]
   # n          <- dim(glmks)[3]
   # X          <- matrix(c(aperm(sweep(glmks, 1:2, gmean), c(2, 1, 3))), nrow = n, byrow = TRUE)
   # covmat     <- cov(X)
   # if (!require(shapes)) stop("the shapes package is required.")
   # g4          <- gmean[c("pn","se", "enL", "enR" ), ]
   # g4mean      <- apply(g4, 2, mean)
   # #covmat4     <- cov(matrix(c(aperm(sweep(glmks[c(1,5,9,10),,], 1:2, 
   # #                        gmean[c(1,5,9,10),]), c(2, 1, 3))), nrow = n, byrow = TRUE))
    # covmat4     <- cov(matrix(c(aperm(sweep(glmks[c(1:4),,], 1:2, 
                           # gmean[c(1:4),]), c(2, 1, 3))), nrow = n, byrow = TRUE))
                        
     
      # # if (monitor){
       # #lmks4             <- face$landmarks[c("pn", "se", "enL", "enR"), ]   
       # #lmean             <- apply(lmks4, 2, mean)
       # #opa               <- procOPA(g4, lmks4, scale = FALSE)
       # #face1$landmarks        <- sweep(sweep(face$landmarks,   2, lmean) %*% opa$R, 2, g4mean, "+")
       # #face1$vertices      <- sweep(sweep(face$vertices, 2, lmean) %*% opa$R, 2, g4mean, "+")
      # # plot(face1, new = FALSE)
                  # # for (i in 1:k) {
                        # # ind      <- (i - 1) * 3 + 1:3
                         # # covmat.i <- covmat[ind, ind]
                         # # ell      <- ellipse3d(covmat.i, centre = gmean[i, ], level = 0.99)
                         # # plot3d(ell, add = TRUE, col = i, alpha = 0.5)
                      # # }
                  # # }
                # # scan()
                  
# # Now bring in prior information on the four known landmarks 
   # #for pn/se
   # si        <- c(1,.25, -1, -1)
   # si.within <- c(.25,.35,.3, .3)
   # pair      <- c(FALSE, FALSE, FALSE, FALSE)
   # stddev    <- c(10,10,10, 10)
   
        
# for ( i in 1:4) {
   # for (iter in 1:niter) { 
    	# #i <- 2   #just doing sellion at the moment
    	  # # doing opa on "new placement"
      # lmks4             <- face$landmarks[c("pn", "se", "enL", "enR"), ]   
      # id                <- rownames(lmks4)[i]
      # lmean             <- apply(lmks4, 2, mean)
      # opa               <- procOPA(g4, lmks4, scale = FALSE)
      # face$landmarks         <- sweep(sweep(face$landmarks,   2, lmean) %*% opa$R, 2, g4mean, "+")
      # face$vertices       <- sweep(sweep(face$vertices, 2, lmean) %*% opa$R, 2, g4mean, "+")
      
      # #picking out 4 landmarks
      # lmks.cond         <- face$landmarks[c("pn", "se", "enL", "enR"), ]   
      # if(monitor){plot(face, new=FALSE)
      	           # spheres3d(face$landmarks)
      	           # print("found landmarks")
      	           # scan()}
      	           
      # lmk               <- lmks.cond[i, ]
      # ind               <- (i - 1) * 3 + 1:3
      # covmat.i          <- covmat4[ind, ind]                         
      # cond.cov          <- covmat.i - covmat4[ind, -ind] %*% solve(covmat4[-ind,-ind]) %*%covmat4[-ind, ind]   
      # cond.mean         <-  c(t(g4))[ind] +  covmat4[ind, -ind] %*% solve(covmat4[-ind,-ind]) %*% 
                            # (c(t(lmks.cond[-i,]))  - c(t(g4))[-ind] ) 
      # crds              <- sweep(face$vertices, 2, cond.mean)
      # M.dst             <- rowSums((crds %*% solve(cond.cov)) * crds)
      # crds              <- sweep(crds, 2, cond.mean, "+")
      # ind               <- (M.dst <= stddev[i])
      # logprior          <- -M.dst[ind]  
      # #face              <- curvatures(face, subset = ind)  
      # sbst              <- subset.face3d(face, ind)
     # if(monitor){plot(face, new=FALSE)
      	           # spheres3d(sbst$vertices)
      	           # print("prior of interest")
      	           # scan()}
     # if(monitor){plot(sbst,colour="shape index", new=FALSE)
     	         # print("shape index of prior")
     	         # scan() }
      	           	           
      # ql                <- quadlocate(face=sbst, si = si[i], si.within = si.within[i], pair = pair[i],
                           # constraint=NA, monitor = FALSE)
      # loglik            <- ql$rss
      # if (monitor){plot(face, new=FALSE)
      	           # spheres3d(ql$sbst.pair, col = loglik)
      	           # print("loglik")
      	           # scan()}
      	           
      # logprior          <- logprior[ql$indices2]  
            # if (monitor){plot(face, new=FALSE)
      	           # spheres3d(ql$sbst.pair, col = -logprior)
      	                 	           # print("logprior")

      	           # scan()}
      	                                                     
      # logpost           <- loglik + logprior
            # if (monitor){plot(face, new=FALSE)
      	           # spheres3d(ql$sbst.pair, col = -logpost)
      	                 	           # print("logpost")
      	           # scan()}
      	           
      # new.cond.lmk      <- ql$sbst.pair[as.numeric(names(which.min(logpost))), ]        
      # if (monitor){       plot(face, new = FALSE)
                           # spheres3d(face$landmarks)  
                           # spheres3d(new.cond.lmk, radius=.5,col="red")
                           # print("red is new landmark")
                           # scan() 
                    # }
                     
        # } 
        
        # face$landmarks[id,] <- new.cond.lmk
     
     # }
      
 # ## Now that first four landmarks have stabilized, begin the next landmark sequentially  

    # #IF k =5 show log prior and log lik then show them added
    # #if K=7 , can be driven too much by the log prior
    # #if k=9/10 , gets confused by the eyelid
      # si          <- c(-.25,  -.8,   -.8,   .8,  -.8,  -.8,   .6,   .7,   .7,   -.5,  -.5,   -.8,   .8,  -.7,   .7)
      # si.within   <- c(  .3,   .4,    .4,   .3,   .3,   .3,   .2,   .2,   .2,    .2,   .2,    .3,   .3,   .3,   .3)
      # stdd        <- c(  20,   20,    20,   10,    7,    7,   20,   20,   20,    30,   30,    20,   20,   20,   20)
      
      # tt<-1
      
     # for (k in 5: dim(gmean)[1]) { 
   
     # lmk                <- gmean[k, ]
     # ind                <- (k - 1) * 3 + 1:3 
     # covmat.i           <- covmat[ind, ind]                
     # crds               <- sweep(face$vertices, 2, lmk)
     # dst                <- rowSums((crds %*% solve(covmat.i)) * crds)
     # crds               <- sweep(crds, 2, lmk, "+")
     # ind                <- (dst <= stdd[tt])
     # sbst               <- subset.face3d(face, ind)
     # sbst               <- curvatures(sbst, distance = 5)
     # ql                 <- quadlocate(sbst, rep(TRUE, nrow(sbst$vertices)), 
     					  # si= si[tt], si.within= si.within[tt], pair = FALSE, monitor = FALSE)
     # logprior           <- -dst
     # logprior           <- logprior[ql$indices2]
     # loglik             <- -ql$rss
     # logpost            <- loglik + logprior
     # id                 <- rownames(gmean)[k] 
     # face$landmarks[id, ]    <- ql$sbst.pair[as.numeric(which.max(logpost)), ] 
     # if(monitor){ 
     # plot(face)
     # spheres3d(ql$sbst.pair[as.numeric(which.max(loglik)), ], col="orange" )
     # spheres3d(ql$sbst.pair[as.numeric(which.max(logprior)), ], col="blue" )
     # spheres3d(ql$sbst.pair[as.numeric(which.max(logpost)), ] , radius=1.1) }
     
      # lmks.new          <- face$landmarks[c(1:k), ]   
      # lmean             <- apply(lmks.new, 2, mean)
      # opa               <- procOPA(gmean[c(1:k), ], lmks.new, scale = FALSE)
      # face$landmarks         <- sweep(sweep(face$landmarks,   2, lmean) %*% opa$R, 2, (apply(gmean[c(1:k), ], 2, mean)), "+")
      # face$vertices       <- sweep(sweep(face$vertices, 2, lmean) %*% opa$R, 2, (apply(gmean[c(1:k), ], 2, mean)), "+")
      # tt                <- tt+1
     
     # }
     # if (monitor){
     # spheres3d(sbst$vertices)
     # plot(face)
     # spheres3d(face$landmarks) 
 # }
  
  
  
  
  
      # invisible(face)
      # }
  
  
  
  
  
  
# # #   
    # # # #to do in pairs
    
     # k<- 6
     # lmk                <- gmean[k, ]
     # ind                <- (k - 1) * 3 + 1:3
     # covmat.i           <- covmat[ind, ind]
     # crds               <- sweep(face$vertices, 2, lmk)
     # dst                <- rowSums((crds %*% solve(covmat.i)) * crds)
     # crds               <- sweep(crds, 2, lmk, "+")
     # ind                <- (dst <= stdd[tt])
     # sbst1              <- subset.face3d(face, ind)
     # sbst1              <- curvatures(sbst1, distance = 5)
     # logprior           <- -dst
     # ind1               <- which(abs(sbst1$shape.index - si[tt]) < si.within[tt])
     # sub1               <- subset.face3d(sbst1, ind1, remove.singles=FALSE)
     # parts              <- connected.face3d(sub1)
     # ind2               <- which(parts==1)
     # sub1                <- subset.face3d(sub1, parts == 1)
     # #logprior1          <- logprior[ind2]
 
  # k<-7
     # lmk                <- gmean[k, ]
     # ind                <- (k - 1) * 3 + 1:3
     # covmat.i           <- covmat[ind, ind]
     # crds               <- sweep(face$vertices, 2, lmk)
     # dst                <- rowSums((crds %*% solve(covmat.i)) * crds)
     # crds               <- sweep(crds, 2, lmk, "+")
     # ind                <- (dst <= stdd[tt])
     # sbst               <- subset.face3d(face, ind)
     # sbst               <- curvatures(sbst, distance = 5)
     # logprior           <- -dst
     # ind1               <- which(abs(sbst$shape.index - si[tt]) < si.within[tt])
     # sub                <- subset.face3d(sbst, ind1, remove.singles=FALSE)
     # parts              <- connected.face3d(sub)
     # ind2               <- which(parts==1)
     # sub                <- subset.face3d(sub, parts == 1)
     # #logprior           <- logprior[ind2]




  # g7                <- gmean[c("pn","se", "enL", "enR" , "sn", "acL", "acR"), ]    
  # g7mean            <- apply(g7, 2, mean)     
  # covmat7           <- cov(matrix(c(aperm(sweep(glmks[c(1:7),,], 1:2, 
                           # gmean[c(1:7),]), c(2, 1, 3))), nrow = n, byrow = TRUE))          
 
  
  
  
  
  
                
  # prod     <- rep(NA, dim(sub1$vertices)[1])
  # allprod  <- matrix(NA, ncol= dim(sub1$vertices)[1], nrow= dim(sub$vertices)[1])

  
 # for(i in 1:dim(sub$vertices)[1]) {
 # for(j in 1:dim(sub1$vertices)[1]) {
  
  # lmks.cond         <- rbind(face$landmarks[c("pn", "se", "enL", "enR", "sn"),], sub1$vertices[j,], sub$vertices[i,])  
  # #lmk               <- face$landmarks[c("acL", "acR"),]         
  
  # ind               <- (7 - 2) * 3 + 1:6
  # covmat.i          <- covmat7[ind, ind]                         
  # cond.cov          <- covmat.i - covmat7[ind, -ind] %*% solve(covmat7[-ind,-ind]) %*%covmat7[-ind, ind]   
  # cond.mean         <-  c(t(g7))[ind] +  covmat7[ind, -ind] %*% solve(covmat7[-ind,-ind]) %*% 
                            # (c(t(lmks.cond[c(1:5),]))  - c(t(g7))[-ind] ) 
                            
   # crds              <- sweep(face$vertices, 2, cond.mean)
      # M.dst             <- rowSums((crds %*% solve(cond.cov)) * crds)
      # crds              <- sweep(crds, 2, cond.mean, "+")
      # ind               <- (M.dst <= 20)
      # lgp               <- -M.dst[ind]
      # aa                <- which.max(lgp) 
      # logprior          <- lgp[aa]
      # if (length(logprior) == 0){
         # prod[j] <-   0 
         # } else{ 
      # prod[j] <- logprior
      # }
    
    
     # }
     # allprod[i, ] <- prod
  
  # }
  

  
   # # # # a <- which(allprod == max(allprod), arr.ind=TRUE)
   # # # # logprior  <- logprior[a[1]]
   # # # # logprior1 <- logprior1[a[2]]



  # #prior value for the 1st in sbst max
  # #logprior  <- apply(allprod, 1, max)
  # logprior1 <- apply(allprod, 2, max)

# #redo for loop the other way around for ind to get logprior i.e. acR


 # g7                <- gmean[c("pn","se", "enL", "enR" , "sn", "acR", "acL"), ]    
  # g7mean            <- apply(g7, 2, mean)     
  # covmat7           <- cov(matrix(c(aperm(sweep(glmks[c(1:5, 7, 6),,], 1:2, 
                           # gmean[c(1:5, 7, 6),]), c(2, 1, 3))), nrow = n, byrow = TRUE))          
 
 
# prod     <- rep(NA, dim(sub1$vertices)[1])
  # allprod  <- matrix(NA, ncol= dim(sub1$vertices)[1], nrow= dim(sub$vertices)[1])

  
 # for(i in 1:dim(sub$vertices)[1]) {
 # for(j in 1:dim(sub1$vertices)[1]) {
  
  
  # lmks.cond         <- rbind(face$landmarks[c("pn", "se", "enL", "enR", "sn"),], sub$vertices[i,], sub1$vertices[j,])  
  # #lmk               <- face$landmarks[c("acL", "acR"),]         
  
  # ind               <- (7 - 1) * 3 + 1:3
  # covmat.i          <- covmat7[ind, ind]                         
  # cond.cov          <- covmat.i - covmat7[ind, -ind] %*% solve(covmat7[-ind,-ind]) %*%covmat7[-ind, ind]   
  # cond.mean         <-  c(t(g7))[ind] +  covmat7[ind, -ind] %*% solve(covmat7[-ind,-ind]) %*% 
                            # (c(t(lmks.cond[-7,]))  - c(t(g7))[-ind] ) 
                            
   # crds              <- sweep(face$vertices, 2, cond.mean)
      # M.dst             <- rowSums((crds %*% solve(cond.cov)) * crds)
      # crds              <- sweep(crds, 2, cond.mean, "+")
      # ind               <- (M.dst <= 20)
      # lgp               <- -M.dst[ind]
      # aa                <- which.max(lgp) 
      # logprior          <- lgp[aa]
      # if (length(logprior) == 0){
         # prod[j] <-   0 
         # } else{ 
      # prod[j] <- logprior
      # }
    
    
     # }
     # allprod[i, ] <- prod
  # }
  

 # #prior value for the 1st in sbst max
  # logprior  <- apply(allprod, 1, max)
  # #logprior1 <- apply(allprod, 2, max)

















# ql                 <- quadlocate(sbst, rep(TRUE, nrow(sbst$vertices)),
     					  # si= si[tt], si.within= si.within[tt], pair = FALSE, monitor = FALSE)

# ql1                 <- quadlocate(sbst1, rep(TRUE, nrow(sbst1$vertices)),
     					  # si= si[tt], si.within= si.within[tt], pair = FALSE, monitor = FALSE)


    # loglik             <- ql$rss
    # loglik1            <- ql1$rss

# #joint loglik


   # prod     <- rep(NA, length(loglik1))
  # allprod  <- matrix(NA, ncol= length(loglik1), nrow= length(loglik))

    # for(i in 1:length(loglik)) {
     # for(j in 1:length(loglik1)) {
        # prod[j] <-loglik[i]* loglik1[j]
     # }
     # allprod[i, ] <- prod
  # }

  # loglik  <- apply(allprod, 1, max)
  # loglik1 <- apply(allprod, 2, max)





   # logpost  <- loglik + logprior
   # logpost1 <- loglik1 + logprior1

     # k<-6
     # id                 <- rownames(gmean)[k]
     # face$landmarks[id, ]    <- sub$vertices[as.numeric(which.max(logpost)), ]

     # k<-7
     # id                 <- rownames(gmean)[k]
     # face$landmarks[id, ]    <- sub1$vertices[as.numeric(which.max(logpost1)), ]







     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
   # # tt <- 2  
  # # #for (k in lmks.undone) { 
   # # for (k in 6: dim(gmean)[1]) {
   	# # g.new        <- gmean[c(1:(k-1)), ]
     # # g.new.mean   <- apply(g.new, 2, mean)
     # # covmat.new   <- cov(matrix(c(aperm(sweep(glmks[c(1:(k-1)),,], 1:2, 
                           # # gmean[c(1:(k-1)),]), c(2, 1, 3))), nrow = n, byrow = TRUE))
     # # lmks.new           <- rbind(face$landmarks[c(1:(k-1)), ])
     # # id                 <- rownames(lmks.new)[(k-1)]
     # # lmean              <- apply(lmks.new, 2, mean)
     # # opa                <- procOPA(g.new, lmks.new, scale = FALSE)
     # # face$landmarks          <- sweep(sweep(face$landmarks,   2, lmean) %*% opa$R, 2, g.new.mean, "+")
     # # face$vertices        <- sweep(sweep(face$vertices, 2, lmean) %*% opa$R, 2, g.new.mean, "+")

     # # lmk               <- gmean[k, ]
     # # ind               <- (k - 1) * 3 + 1:3 
     # # covmat.i          <- covmat[ind, ind]                
     # # crds              <- sweep(face$vertices, 2, lmk)
     # # dst               <- rowSums((crds %*% solve(covmat.i)) * crds)
     # # crds              <- sweep(crds, 2, lmk, "+")
     # # ind               <- (dst <= stdd[tt])
     # # sbst              <- subset.face3d(face, ind)
     # # sbst              <- curvatures(sbst)
     # # ql                <- quadlocate(sbst, rep(TRUE, nrow(sbst$vertices)), 
     					  # # si= si[tt], si.within= si.within[tt], pair = FALSE, monitor = FALSE)
     # # logprior          <- -dst
     # # logprior          <- logprior[ql$indices2]
     # # loglik            <- ql$rss
     # # logpost           <- loglik + logprior
     # # id                <- rownames(gmean)[k] 
     # # face$landmarks[id, ]   <- ql$sbst.pair[as.numeric(which.min(logpost)), ]  
      

      
   # #generalize from here
  
     # # g.new        <- gmean[c(1:k), ]
     # # g.new.mean   <- apply(g.new, 2, mean)
     # # covmat.new   <- cov(matrix(c(aperm(sweep(glmks[c(1:k),,], 1:2, 
                           # # gmean[c(1:k),]), c(2, 1, 3))), nrow = n, byrow = TRUE))
     # # lmks.new           <- rbind(face$landmarks[c(1:k), ])
     # # id                 <- rownames(lmks.new)[k]
     # # lmean              <- apply(lmks.new, 2, mean)
     # # opa                <- procOPA(g.new, lmks.new, scale = FALSE)
     # # face$landmarks          <- sweep(sweep(face$landmarks,   2, lmean) %*% opa$R, 2, g.new.mean, "+")
     # # face$vertices        <- sweep(sweep(face$vertices, 2, lmean) %*% opa$R, 2, g.new.mean, "+")
      # # lmks.cond         <- face$landmarks[c(1:k), ]    	           
      # # lmk               <- lmks.cond[k, ]
      # # ind               <- (k - 1) * 3 + 1:3
      # # covmat.i          <- covmat.new[ind, ind]                         
      # # cond.cov          <- covmat.i - covmat.new[ind, -ind] %*% solve(covmat.new[-ind,-ind]) %*%covmat.new[-ind, ind]   
      # # cond.mean         <-  c(t(g.new))[ind] +  covmat.new[ind, -ind] %*% solve(covmat.new[-ind,-ind]) %*% 
                            # # (c(t(lmks.cond[-i,]))  - c(t(g.new))[-ind] ) 
      # # crds              <- sweep(face$vertices, 2, cond.mean)
       # # M.dst             <- rowSums((crds %*% solve(cond.cov)) * crds)
       # # crds              <- sweep(crds, 2, cond.mean, "+")
       # # ind               <- (M.dst <= stddev[tt])
      # # logprior          <- -M.dst[ind]  
      # # face              <- curvatures(face, subset = ind)  
      # # sbst              <- subset.face3d(face, ind)           	           
      # # ql                <- quadlocate(face=sbst, si = si[tt], si.within = si.within[tt], pair = FALSE,
                           # # constraint=NA, monitor = FALSE)
      # # loglik            <- ql$rss    	           
      # # logprior          <- logprior[ql$indices2]                                
      # # logpost           <- loglik + logprior
      # # new.cond.lmk      <- ql$sbst.pair[as.numeric(which.min(logpost)), ] 
      # # id                <- rownames(gmean)[k]                     
      # # face$landmarks[id,]    <- new.cond.lmk
       # # tt                 <- tt + 1
 # # print(k)
      # # }


  
  
  
  
  
     
      
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   

      # # # if (monitor) cat("Locating sn ... ")
        # # # i                 <- 4
        # # # lmk               <- gmean[i, ]
        # # # ind               <- (i - 1) * 3 + 1:3
        # # # covmat.i          <- covmat[ind, ind]                
        # # # crds              <- sweep(face$vertices, 2, lmk)
        # # # dst               <- rowSums((crds %*% solve(covmat.i)) * crds)
        # # # crds              <- sweep(crds, 2, lmk, "+")
        # # # ind               <- (dst <= 20)
        # # # sbst              <- subset.face3d(face, ind)
        # # # sbst              <- curvatures(sbst)
        # # # plot(sbst, colour = "shape index", new = FALSE)
        # # # spheres3d(lmk, col = "green")
        # # # scan()
        # # # plot(sbst, colour = 2 + sbst$kappa1 * sbst$kappa2, new = FALSE)
        # # # spheres3d(lmk, col = "green")        
        # # # scan()
        # # # ql <- quadlocate(sbst, rep(TRUE, nrow(sbst$vertices)), 0, 0.25, pair = FALSE, monitor = FALSE)
        # # # print(ql)
        # # # print(nrow(sbst$vertices))
        # # # plot(ql$sbst, colour = 2 + ql$rss, new = FALSE)
        # # # spheres3d(lmk, col = "green")        
        # # # spheres3d(ql$lmk, col = "red")        
        # # # scan()
        # # # # plot3d(ellipse3d(covmat.i, lmk, level = 0.99), add = TRUE, col = "blue", alpha = 0.5)
        # # # logprior          <- -dst
        # # # gc                <- face$kappa1 * face$kappa2
        # # # loglik            <- log(gc[which.max(gc)])
        # # # logpost.i         <- loglik + logprior
        # # # face$landmarks["sn", ] <- crds[which.max(logpost.i), ]                
      
   





    
# # #      for (iter in 1:niter){ 
   # # # #do pair  - find each prior separately, colour subset only by mahal distance. 
      # # covmat4           <- cov(matrix(c(aperm(sweep(glmks[c(1,5,9,10),,], 1:2, 
                           # # gmean[c(1,5,9,10),]), c(2, 1, 3))), nrow = n, byrow = TRUE))
     # # lmks.cond         <- face$landmarks[c("pn", "se", "enL", "enR"), ]   
      # # lmk               <- lmks.cond[3, ]
      # # ind               <- (3 - 1) * 3 + 1:3
      # # covmat.i          <- covmat4[ind, ind]                         
      # # cond.cov          <- covmat.i - covmat4[ind, -ind] %*% solve(covmat4[-ind,-ind]) %*%covmat4[-ind, ind]   
      # # cond.mean         <-  c(t(g4))[ind] +  covmat4[ind, -ind] %*% solve(covmat4[-ind,-ind])%*%   
                            # # (c(t(lmks.cond[-3,]))  - c(t(g4))[-ind] )   
      # # crds              <- sweep(face$vertices, 2, cond.mean)
      # # M.dst             <- rowSums((crds %*% solve(cond.cov)) * crds)
      # # crds              <- sweep(crds, 2, cond.mean, "+")
      # # ind               <- (M.dst <= 20)
      # # logprior          <- -M.dst[ind]
      # # M.dst.L           <- M.dst
      # # ind.L             <- ind
      # # logprior.L        <- logprior
      # # lmk               <- lmks.cond[4, ]
      # # ind               <- (4 - 1) * 3 + 1:3
      # # covmat.i          <- covmat4[ind, ind]                         
      # # cond.cov          <- covmat.i - covmat4[ind, -ind] %*% solve(covmat4[-ind,-ind]) %*%covmat4[-ind, ind]   
      # # cond.mean         <-  c(t(g4))[ind] +  covmat4[ind, -ind] %*% solve(covmat4[-ind,-ind]) %*% 
                            # # (c(t(lmks.cond[-4,]))  - c(t(g4))[-ind] )   
      # # crds              <- sweep(face$vertices, 2, cond.mean)
      # # M.dst             <- rowSums((crds %*% solve(cond.cov)) * crds)
      # # crds              <- sweep(crds, 2, cond.mean, "+")
      # # ind               <- (M.dst <= 20)
      # # logprior          <- -M.dst[ind]
      # # M.dst.R           <- M.dst
      # # ind.R             <- ind
      # # logprior.R        <- logprior

      # # sbst                                <- subset.face3d(face,(edist.face3d(face$vertices, 
                                             # # face$landmarks["pn", ]) < 70) & (face$vertices[ , 2] > face$landmarks["pn", 2] + 10) ) #initial subset on area
      # # sbst                                <- curvatures(sbst)    
      # # ind.final                           <- c( which(ind.L== "TRUE"), which(ind.R==TRUE)) 
      # # ind.good                            <- match( ind.final,rownames(sbst$vertices))
      # # sbst$shape.index[-ind.good]         <- 1
      # # sbst$kappa1[-ind.good]              <- -1       
      # # sbst$kappa2[-ind.good]              <- -1 
      # # ql                                  <- quadlocate(face=sbst, si = -1, si.within = 0.3, 
                                              # # pair = TRUE, constraint= NA, monitor = FALSE)
      # # loglik                              <- ql$rss
      # # logpri.L                            <- logprior.L[ql$indices2]
      # # logpri.R                            <- logprior.R[ql$indices2]
      # # logpost.L                           <- logpri.L +loglik
      # # logpost.R                           <- logpri.R +loglik 
      # # new.cond.lmk.R                      <- ql$sbst.pair[,,1][which.max(logpost.R), ]   
      # # new.cond.lmk.L                      <- ql$sbst.pair[,,2][which.max(logpost.L), ] 
      # # if (monitor){                       plot(face, new=FALSE)
                                           # # spheres3d(face$landmarks)  
                                           # # spheres3d(lmks.cond, col="red") 	
                   # # } 
# # face$landmarks[c("enL", "enR"),] <- rbind(new.cond.lmk.L, new.cond.lmk.R)

  # # } 
    
         
               
                
 # # # Now that first four landmarks have stabilized, begin the next landmark sequentially  


   # # # #locating the sn
   # # # monitor.extra <- FALSE
   
      # # # if (monitor) cat("Locating sn ... ")
        # # # i                 <- 4
        # # # lmk               <- gpa$mean[i, ]
        # # # ind               <- (i - 1) * 3 + 1:3
        # # # covmat.i          <- covmat[ind, ind]                
        # # # crds              <- sweep(face$vertices, 2, lmk)
        # # # dst               <- rowSums((crds %*% solve(covmat.i)) * crds)
        # # # crds              <- sweep(crds, 2, lmk, "+")
        # # # ind               <- (dst <= 20)
        # # # sbst              <- subset.face3d(face, ind)
        # # # sbst              <- curvatures(sbst)
        # # # plot(sbst, colour = "shape index", new = FALSE)
        # # # spheres3d(lmk, col = "green")
        # # # scan()
        # # # plot(sbst, colour = 2 + sbst$kappa1 * sbst$kappa2, new = FALSE)
        # # # spheres3d(lmk, col = "green")        
        # # # scan()
        # # # ql <- quadlocate(sbst, rep(TRUE, nrow(sbst$vertices)), 0, 0.25, pair = FALSE, monitor = FALSE)
        # # # print(ql)
        # # # print(nrow(sbst$vertices))
        # # # plot(ql$sbst, colour = 2 + ql$rss, new = FALSE)
        # # # spheres3d(lmk, col = "green")        
        # # # spheres3d(ql$lmk, col = "red")        
        # # # scan()
        # # # # plot3d(ellipse3d(covmat.i, lmk, level = 0.99), add = TRUE, col = "blue", alpha = 0.5)
        # # # logprior          <- -dst
        # # # gc                <- face$kappa1 * face$kappa2
        # # # loglik            <- log(gc[which.max(gc)])
        # # # logpost.i         <- loglik + logprior
        # # # face$landmarks["sn", ] <- crds[which.max(logpost.i), ]               
