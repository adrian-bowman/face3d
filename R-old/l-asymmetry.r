
asymmetry.face3d <- function(landmarks, coding =  NULL, 
                      landmark.names = rownames(landmarks), match.landmarks =     
                      unlist(coding), subset.string = NULL, Centroid=NULL, display = FALSE) {

   if (!require(shapes)) stop("this function requires the shapes package.")
 
   
    if (coding =="lmks"){
    	coding     <- list(single.lmks = c("pn", "se", "sn", "n", "ls", "st", "li", "sl", "gn"),
                           paired.lmks = matrix(c("acL", "exL", "enL", "tL", "cphL", "chL", "oiL",
                                          "acR", "exR", "enR", "tR", "cphR", "chR", "oiR"),
                                          ncol = 2))
    }
   
    else if (coding =="curves"){
    	       # if (grep("upper eye socket", rownames(landmarks)) == 1) stop("this function does not account for upper eye sockets. ")
               	ind           <- grep("mid-line", rownames(landmarks))
                midline.1     <- landmarks[ind, ]
                ind.rm        <- grep("mid-line lip", rownames(midline.1))
                midline       <- midline.1[-ind.rm, ]
 paired.curves <-  matrix(c(rownames(landmarks[grep("lower eye socket right", rownames(landmarks)),]), 
                            rownames(landmarks[grep("nasal root right", rownames(landmarks)),]),   
                            rownames(landmarks[grep("nasal boundary right", rownames(landmarks)),]),
                            rownames(landmarks[grep("nasal bridge right", rownames(landmarks)),]),   
                            rownames(landmarks[grep("nasal base right", rownames(landmarks)),]),                                            
                            rownames(landmarks[grep("nasolabial right", rownames(landmarks)),]),
                            rownames(landmarks[grep("mid-line lip right", rownames(landmarks)),]),
                            rownames(landmarks[grep("lower lip right", rownames(landmarks)),]),
                            rownames(landmarks[grep("upper lip right", rownames(landmarks)),]),
                            rownames(landmarks[grep("philtrum ridge right", rownames(landmarks)),]),
                            rownames(landmarks[grep("cheek-lip right",  rownames(landmarks)),]),
                            rownames(landmarks[grep("cheek-nose right", rownames(landmarks)),]),
                            rownames(landmarks[grep("cheek-eye right",  rownames(landmarks)),]),
                            rownames(landmarks[grep("brow ridge right", rownames(landmarks)),]),
                            rownames(landmarks[grep("mandible right", rownames(landmarks)),]),
                            rownames(landmarks[grep("philtrum lip right",  rownames(landmarks)),]),
                            rownames(landmarks[grep("lower eye socket left",  rownames(landmarks)),]),
                            rownames(landmarks[grep("nasal root left",  rownames(landmarks)),]),
                            rownames(landmarks[grep("nasal boundary left",  rownames(landmarks)),]),
                            rownames(landmarks[grep("nasal bridge left", rownames(landmarks)),]),
                            rownames(landmarks[grep("nasal base left",  rownames(landmarks)),]),
                            rownames(landmarks[grep("nasolabial left",  rownames(landmarks)),]),
                            rownames(landmarks[grep("mid-line lip left",  rownames(landmarks)),]),
                            rownames(landmarks[grep("lower lip left", rownames(landmarks)),]),
                            rownames(landmarks[grep("upper lip left", rownames(landmarks)),]),
                            rownames(landmarks[grep("philtrum ridge left",  rownames(landmarks)),]),
                            rownames(landmarks[grep("cheek-lip left", rownames(landmarks)),]),
                            rownames(landmarks[grep("cheek-nose left", rownames(landmarks)),]),
                            rownames(landmarks[grep("cheek-eye left",  rownames(landmarks)),]),
                            rownames(landmarks[grep("brow ridge left", rownames(landmarks)),]),
                            rownames(landmarks[grep("mandible left",  rownames(landmarks)),]),
                            rownames(landmarks[grep("philtrum lip left", rownames(landmarks)),])), ncol=2)
                                            


 coding        <- list(single.lmks = rownames(midline),     
                       paired.lmks = paired.curves)              
    }
   
   
   
   
   else if (coding =="meshes"){
   	coding <- list(single.lmks= c(rownames(landmarks[grep("mid-line columella", rownames(landmarks)),]),
  	                              rownames(landmarks[grep("mid-line chin", rownames(landmarks)),]),
  	                              rownames(landmarks[grep("mid-line nasal root", rownames(landmarks)),]),
  	                              rownames(landmarks[grep("mid-line nasal profile", rownames(landmarks)),]),
  	                              rownames(landmarks[grep("mid-line philtral", rownames(landmarks)),]),
  	                              rownames(landmarks[grep("mid-line upper lip", rownames(landmarks)),]),
  	                              rownames(landmarks[grep("mid-line bottom lip", rownames(landmarks)),])),
  	                             
  	paired.lmks = matrix(c(rownames(landmarks[grep("ULR", rownames(landmarks)),]),
                           rownames(landmarks[grep("LLR", rownames(landmarks)),]),
                           rownames(landmarks[grep("MLR", rownames(landmarks)),]),
  	                       rownames(landmarks[grep("lower face right", rownames(landmarks)),]),
  	                       rownames(landmarks[grep("mid-face right", rownames(landmarks)),]),
                           rownames(landmarks[grep("upper mid face right", rownames(landmarks)),]),
                           rownames(landmarks[grep("upper face right 3", rownames(landmarks)),]),
                           rownames(landmarks[grep("upper face right 1", rownames(landmarks)),]),
                           rownames(landmarks[grep("upper face right 2", rownames(landmarks)),]),
                           rownames(landmarks[grep("nose right", rownames(landmarks)),]),
                           rownames(landmarks[grep("philtrum right", rownames(landmarks)),]),
                           rownames(landmarks[grep("upper lip right", rownames(landmarks)),]),
                           rownames(landmarks[grep("lower lip right", rownames(landmarks)),]), 
                           rownames(landmarks[grep("ULL", rownames(landmarks)),]),
                           rownames(landmarks[grep("LLL", rownames(landmarks)),]),
                           rownames(landmarks[grep("MLL", rownames(landmarks)),]),
                           rownames(landmarks[grep("lower face left", rownames(landmarks)),]),
                           rownames(landmarks[grep("mid-face left", rownames(landmarks)),]),
                           rownames(landmarks[grep("upper mid face left", rownames(landmarks)),]),
                           rownames(landmarks[grep("upper face left 3", rownames(landmarks)),]),
                           rownames(landmarks[grep("upper face left 1", rownames(landmarks)),]),
                           rownames(landmarks[grep("upper face left 2", rownames(landmarks)),]),
                           rownames(landmarks[grep("nose left", rownames(landmarks)),]),
                           rownames(landmarks[grep("philtrum left", rownames(landmarks)),]),
                           rownames(landmarks[grep("upper lip left", rownames(landmarks)),]),
                           rownames(landmarks[grep("lower lip left", rownames(landmarks)),])),ncol=2))
  }  
  
  else if (coding == NULL){
  	stop("this function requires a coding.")
  }
    

#to extract specific features
     
   if (!missing(subset.string)) {
   	
   	   if (subset.string == "lower face"){
   			coding <- list(single.lmks= rownames(landmarks[grep("mid-line chin", rownames(landmarks)),]),
  	                             
  	                       paired.lmks = matrix(c(rownames(landmarks[grep("lower face right", rownames(landmarks)),]),
                                                  rownames(landmarks[grep("lower face left", rownames(landmarks)),])),ncol=2))
                          
            ind            <- match(coding$single.lmks, rownames(landmarks))
            indd           <- match(coding$paired.lmks, rownames(landmarks))
            landmarks      <- landmarks[c(ind,indd), ]           
            landmark.names <- rownames(landmarks)
    
         }  
   		
      	else if (subset.string == "nose"){
   			coding <- list(single.lmks= c(rownames(landmarks[grep("mid-line columella", rownames(landmarks)),]),
  	                                      rownames(landmarks[grep("mid-line nasal profile", rownames(landmarks)),])),
  	                             
  	                       paired.lmks = matrix(c(rownames(landmarks[grep("nose right", rownames(landmarks)),]),
                                                  rownames(landmarks[grep("nose left", rownames(landmarks)),])),ncol=2))
                          
            ind            <- match(coding$single.lmks, rownames(landmarks))
            indd           <- match(coding$paired.lmks, rownames(landmarks))
            landmarks      <- landmarks[c(ind,indd), ]           
            landmark.names <- rownames(landmarks)
    
         }  
       else if (subset.string == "upper lips"){
   			coding <- list(single.lmks= rownames(landmarks[grep("mid-line upper lip", rownames(landmarks)),]),
  	                             
  	                       paired.lmks = matrix(c(rownames(landmarks[grep("MLR", rownames(landmarks)),]),
  	                                              rownames(landmarks[grep("ULR", rownames(landmarks)),]),
  	                                              rownames(landmarks[grep("upper lip right", rownames(landmarks)),]),
  	                                              rownames(landmarks[grep("MLL", rownames(landmarks)),]),
  	                                              rownames(landmarks[grep("ULL", rownames(landmarks)),]),
                                                  rownames(landmarks[grep("upper lip left", rownames(landmarks)),])),ncol=2))
                          
            ind            <- match(coding$single.lmks, rownames(landmarks))
            indd           <- match(coding$paired.lmks, rownames(landmarks))
            landmarks      <- landmarks[c(ind,indd), ]           
            landmark.names <- rownames(landmarks)
    
         }  

        else if (subset.string == "lower lips"){
   			coding <- list(single.lmks= rownames(landmarks[grep("mid-line bottom lip", rownames(landmarks)),]),
  	                             
  	                       paired.lmks = matrix(c(rownames(landmarks[grep("MLR", rownames(landmarks)),]),
  	                                              rownames(landmarks[grep("LLR", rownames(landmarks)),]),
  	                                              rownames(landmarks[grep("lower lip right",rownames(landmarks)),]),
  	                                              rownames(landmarks[grep("MLL", rownames(landmarks)),]),
  	                                              rownames(landmarks[grep("LLL", rownames(landmarks)),]),
                                                  rownames(landmarks[grep("lower lip left", rownames(landmarks)),])),ncol=2))
                          
            ind            <- match(coding$single.lmks, rownames(landmarks))
            indd           <- match(coding$paired.lmks, rownames(landmarks))
            landmarks      <- landmarks[c(ind,indd), ]           
            landmark.names <- rownames(landmarks)
    
         }  


   		else if (subset.string == "lips"){
   			coding <- list(single.lmks= c(rownames(landmarks[grep("mid-line bottom lip", rownames(landmarks)),]),
  	                                      rownames(landmarks[grep("mid-line upper lip", rownames(landmarks)),])),
  	                             
  	                       paired.lmks = matrix(c(rownames(landmarks[grep("ULR", rownames(landmarks)),]),
  	                                              rownames(landmarks[grep("LLR", rownames(landmarks)),]),
  	                                              rownames(landmarks[grep("MLR", rownames(landmarks)),]),
  	                                              rownames(landmarks[grep("lower lip right", rownames(landmarks)),]),
  	                                              rownames(landmarks[grep("upper lip right", rownames(landmarks)),]),
  	                                              rownames(landmarks[grep("ULL", rownames(landmarks)),]),
  	                                              rownames(landmarks[grep("LLL", rownames(landmarks)),]),
  	                                              rownames(landmarks[grep("MLL", rownames(landmarks)),]),
  	                                              rownames(landmarks[grep("lower lip left", rownames(landmarks)),]),
  	                                              rownames(landmarks[grep("upper lip left", rownames(landmarks)),])),
                                                  ncol=2))
                          
            ind            <- match(coding$single.lmks, rownames(landmarks))
            indd           <- match(coding$paired.lmks, rownames(landmarks))
            landmarks      <- landmarks[c(ind,indd), ]           
            landmark.names <- rownames(landmarks)
    
         }  

   		
   		
}
   		
 
   
   	# if (subset.string == "lower face"){
   	     # subset.string   <- c("lower face", "lower lip", "LLL", "LLR", "MLL", "MLR", "mid-line chin", "mid-line bottom lip")
   	     # coding1         <- list()
   	     # landmark.names1 <- rep(NA)
   	# for (i in 1:length(subset.string)) {
   	     # ind               <- grep(subset.string[i], coding$single.lmks)
   	     # coding1$single.lmks <- rbind(coding1$single.lmks,  coding$single.lmks[ind])
   	     # ind                <- grep(subset.string[i], coding$paired.lmks[ , 1])
   	     # coding1$paired.lmks <- rbind(coding1$paired.lmks, coding$paired.lmks[ind, ])
   	     # ind1               <- landmarks[grep(subset.string[i],rownames(landmarks)),]
   	     # landmarks          <- as.matrix(ind1)
   	     # landmark.names1     <- c(landmark.names1, rownames(landmarks))   
   	    # } 
   	       # }  
   	  	  
   



#landmarks:       a q x m x n array defining landmark configurations ((k+l) lmks in m dimensions)
#currently needs to be a matrix but can change subset.string and coding so that works if it is an array as well. 

#coding:          a list with components paired.lmks and single.lmks
#landmark.names:  a character vector defining the complete list of landmark labels i.e. the rows of the data matrix
#match.landmarks: a character vector defining the landmarks to be used in matching. (unlisted coding) - must be a subset of paired and single
#paired.lmks:     (q/2 x 2) character matrix with the indices of paired lmks in each config. in landmarks
#single.lmks:     l character vector with the indices of the single lmks in each config. in landmarks

   q         <- dim(landmarks)[1]
   m         <- dim(landmarks)[2]
   n         <- dim(landmarks)[3]

# instead of doing for loop can do a full array analysis at once	 
# - if only one set of lmks: then make as a basic matrix
   if (is.na(n)) {
      n         <- 1
      landmarks <- array(landmarks, dim = c(dim(landmarks), 1))
      }
      object    <- coding
 
#Function to calculate asymmetry scores for configurations of landmarks (see below)
   result.G <- asymmetry.prepare(landmarks = landmarks, object$paired.lmks, object$single.lmks,
                                 landmark.names = landmark.names, match.landmarks = match.landmarks)
 #sweep mean of both 
  
  landmarks   <- result.G$landmarks
  ref.rel.opa <- result.G$ref.rel.opa
  ref.rel.opa <- sweep(ref.rel.opa, 2, apply(ref.rel.opa, 2, mean))
  landmarks   <- sweep(landmarks, 2, apply(landmarks, 2, mean))
  if(!missing(Centroid)){
  ref.rel.opa <- ref.rel.opa * (Centroid/centroid.size(ref.rel.opa))
  landmarks   <- landmarks * (Centroid/centroid.size(landmarks))
  }
# funtion so match elements of first argument to the second
   ind      <- pmatch(c(object$single.lmks, object$paired.lmks),
   			          result.G$landmark.names)
# Calculating the actual score   			          
   GA <- rep(0, n)
   for (i in 1:n) {
      object.G            <- landmarks[ind, , i] #for each lmk the sym
      object.G.reflected  <- ref.rel.opa[ind, , i]   #the reflected version
      k.object            <- nrow(object.G)   #the number of lmks
      GA[i]               <- sqrt(sum((object.G - object.G.reflected)^2) / k.object)  #taking the symm version sum
      }
   
# could put the scaling here to make the centroid size 1? and sweep the mean      
# Final output 

    symmetric.data <- (landmarks[ind,,] + ref.rel.opa[ind,,])/2 
        
   results <- list(global = GA, landmarks = landmarks[ind, , ],
                                ref.rel.opa = ref.rel.opa[ind, , ], symmetric = symmetric.data )# , landmark.names.final = lmksnames.final)

   class(results) <- "asymmetry"
   invisible(results)
   
   
  
   }

############################################################################################################################################

#everything below is used for the function above

#to calculate the scaling constant if desired
size <- function(config)
   {
   m <- nrow(config)
   centroid <- apply(config, 2, mean)
   config.centred <- config - matrix(rep(centroid, m), nrow = m, byrow = TRUE)
   sqrt(sum(diag(t(config.centred) %*% config.centred)))
   }


#sizes: [optional for asymmetry.prepare] n-vector containing the scaling constant to use for each config. in landmarks
# If "sizes" is present, the function scales by these amounts to
# 'preshape' the data, rather than by the size of each configuration
# This is used when the results based on different subsets of lmks are
# to be compared, in which case the "sizes" usually correspond to the
# maximum subset of landmarks to be considered.

#Function to calculate asymmetry scores for configurations of landmarks
asymmetry.prepare <- function(landmarks, paired.lmks, single.lmks = NA, landmark.names = NA,
		match.landmarks = c(single.lmks, paired.lmks), sizes = NA) {

   n <- dim(landmarks)[3]
   if (is.na(n)) {
      n <- 1
      landmarks <- array(landmarks,dim=c(dim(landmarks),1))
      }
# alllmks and match.landmarks are the same thing but match.landmarks has the single/paired coding names
   alllmks   <- c(single.lmks, paired.lmks)
   all.ind   <- pmatch(alllmks, landmark.names)
   match.ind <- pmatch(match.landmarks, landmark.names)
   if (any(is.na(all.ind))) 
      print(paste('Error:',alllmks[is.na(all.ind)],
		          ' does not appear in landmark.names',sep=' '))
   if (any(is.na(match.ind)))
      print(paste('Error:',match.landmarks[is.na(match.ind)],
		          ' does not appear in landmark.names',sep=' '))
   if (any(is.na(match(match.ind, all.ind))))
      print('Some matching landmarks do not appear in specified lmks')

   # landmarks to be used and order them
   landmark.names  <- landmark.names[all.ind]
   landmarks       <- landmarks[all.ind,,]
   
   lmksnames.final <- landmark.names
   #print(landmark.names)
   #save(landmark.names, file="~/Desktop/landmarknames.dmp")
   

   # Define indices of single and paired landmarks
   single.ind      <- pmatch(single.lmks, landmark.names)
   if (length(paired.lmks) > 0) {
      paired.ind1  <- pmatch(paired.lmks[,1], landmark.names)
      paired.ind2  <- pmatch(paired.lmks[,2], landmark.names)
      paired.ind   <- as.matrix(cbind(paired.ind1, paired.ind2))
      }
   else
      paired.ind <- paired.lmks


   #######################################
   # Setting up/checking the data

   results <- list(landmark.names = NA, landmarks = NA, sizes = NA, ref.rel.opa = NA)

   m                     <- dim(landmarks)[2]           # Number of dimensions
   k                     <- 2*dim(paired.ind)[1]        # Number of paired landmarks
   if (length(k) == 0) k <- 0       # To handle case of no paired landmarks
   l                     <- length(single.ind)       	# Number of unpaired landmarks
   
   if (length(dim(landmarks)) < 3) {
      landmarks <- array(landmarks, dim = c(dim(landmarks), 1))
      }
   results$landmarks <- landmarks
   Xi.s              <- landmarks
   Yi.s              <- landmarks

   #######################################
   # Reflect and relabel the data...

   if (k == 0)
      for (i in 1:n) Yi.s[,1,i] <- -Yi.s[,1,i]
   else
      for (i in 1:n) Yi.s[,,i] <- reflect.relabel(Xi.s[,,i], paired.ind)

   results$ref.rel.opa <- array(NA, dim = dim(landmarks))

   n.match.pts <- length(match.ind)
   for (i in 1:n) {
      if (!any(is.na(rbind(Xi.s[,,i], Yi.s[,,i])))) {
         diff.centroid <- 
               (matrix(1, ncol= n.match.pts, nrow=k+l) / n.match.pts) %*%
	           (Xi.s[match.ind,,i] - Yi.s[match.ind,,i])
     match.result             <- procOPA(Xi.s[match.ind, , i], Yi.s[match.ind, , i], scale = FALSE)
    	 results$ref.rel.opa[,,i] <- Yi.s[,,i] %*% match.result$R + diff.centroid
    	 	     }
      }
    
   results$landmark.names <- landmark.names
   invisible(results)

   }
 
reflect.relabel <- function(mu, paired.lmks) {
	# Inputs: mu          - (k+l) x m configuration matrix
	#         paired.lmks - (k/2 x 2) matrix with the indices of paired lmks in mu
	# Output: reflected and relabeled version of mu via Q(mu)A 
	n                        <- dim(mu)[1]
	m                        <- dim(mu)[2]
	k                        <- dim(paired.lmks)[1]*2
	l                        <- n-k
	A                        <- diag(rep(1, m))
	A[1, 1]                  <- -1  # This satisfies MBM's criteria that A be orthogonal with det(A)=-1
	output                   <- mu%*%A
	temp                     <- output
	output[paired.lmks[,1],] <- temp[paired.lmks[,2],]
	output[paired.lmks[,2],] <- temp[paired.lmks[,1],]
	invisible(output)
    }
