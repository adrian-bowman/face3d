as.face3d <- function(object, model.mesh = FALSE) {

if (!model.mesh) {
   class(object) <- "face3d"
   if (is.face3d(object, report = TRUE))
      object
   else
      NULL
}
else {
   if (!is.matrix(object)) stop("'object' is a not a matrix.")
   if (dim(object)[2] != 3) stop("'object' must have three columns.")
   if (is.null(rownames(object))) stop("'object' must have rownames.")
   
   #All Patches  
   mesh.names    <- c("nose right",        "mid-face left",      "upper mid face left",  
                      "lower face left",   "upper lip right",    "lower lip left",
                      "upper face right 1","upper face right 3", "philtrum left")
   pts.names     <- paste(mesh.names, "1 ")                               
   triples.right <- NULL
   for (k in 1:length(mesh.names)){
      pts          <- length(grep(pts.names[k], rownames(object)))
      indices      <- grep(mesh.names[k], rownames(object))
      paths        <- length(indices)/pts
      trips        <- matrix(indices, ncol = paths)
      for (i in 1:(paths - 1)){
         for (j in 1:(pts - 1)){
            triples.right <- c(triples.right, 
                             c(trips[j + 1, i + 1], trips[j,     i + 1], trips[j, i], 
                               trips[j,     i],     trips[j + 1, i],     trips[j + 1, i + 1]))
         }
      }
   }
   
   mesh.names       <- c("nose left",         "mid-face right",    "upper mid face right",  
                         "lower face right",  "upper lip left",    "lower lip right",
                         "upper face left 1", "upper face left 3", "philtrum right")
   pts.names        <- paste(mesh.names, "1 ")                               
   triples.left     <- NULL
   for (k in 1:length(mesh.names)){
      pts          <- length(grep(pts.names[k], rownames(object)))
      indices      <- grep(mesh.names[k], rownames(object))
      paths        <- length(indices) / pts
      trips        <- matrix(indices, ncol = paths)
      for (i in 1:(paths-1)){
         for (j in 1:(pts-1)){
            triples.left <- c(triples.left, 
                              c(trips[j,     i],     trips[j,     i + 1], trips[j + 1, i + 1] , 
                                trips[j + 1, i + 1], trips[j + 1, i],     trips[j,     i]))
         }
      }
   }
   
   
   #Stitch Together Gaps
   
   # FOR j,i
   pts.names.1 <- matrix( c(paste("ULR", 1:2), paste("upper lip right 1 1"), rep(NA, 16), 
                            paste("LLL", 1:2), paste("lower lip left 1 1"), rep(NA, 16),
                            paste("lower face left", 19:22, "4"), rep(NA, 15),
                            paste("ULL 1"), paste("upper lip left", 1:6, "1"), rep(NA, 12),
                            paste("philtrum left 7", 6:2), paste("upper lip left", "6", 1:3), paste("lower lip left","6", 2:1), rep(NA, 9), 
                            paste("upper face right 3", 1:6, "2"), rep(NA, 13),
                            paste("upper face left 3", "6", 1:2), paste("nose left 6", 5:1), rep(NA, 12),
                            paste("mid-line chin", 1:4), paste("mid-line bottom lip 3"), rep(NA, 14),
                            paste("mid-face left",1:19, "5"),
                            paste("lower face left",1:19, "4"), 
                            paste("upper mid face left",1:12, "5"), rep(NA,7),
                            paste("mid-face right","19", 1:5),paste("upper mid face right 19 1"),  rep(NA,13),
                            paste("lower lip right",1:6, "1"), rep(NA, 13),
                            paste("lower lip left",1:6, "2"), rep(NA, 13),
                            paste("upper mid face right","19",1:5), rep(NA, 14)),
                          ncol=15)
   
   
   
   pts.names.2 <- matrix(c(paste("ULR 1"), paste("MLR 2"), paste("upper lip right 1 3"), rep(NA, 16),
                           paste("ULL 1"), paste("MLL 2"), paste("upper lip left 1 3"), rep(NA, 16),
                           paste("mid-face left 19 1"), paste("ULL 1"), paste("LLL 2"), paste("lower lip left 1 1"), rep(NA, 15),
                           paste("philtrum left",1:7, "2"), rep(NA, 12),
                           paste("mid-line philtral", 1:5), paste("mid-line upper lip", 1:3), paste("mid-line bottom lip",2:3), rep(NA, 9),
                           paste("nose right", 1:6, "5"),  rep(NA, 13),
                           paste("mid-line nasal root", 1:2), paste("mid-line nasal profile", 1:5), rep(NA, 12),
                           paste("lower face left", "27", 1:4), paste("lower lip left 6 1"), rep(NA,14),
                           paste("upper mid face left", 1:19, "1"),
                           paste("mid-face left", 1:19, "1"),
                           paste("upper face left 1", 1:12, "2"), rep(NA,7),
                           paste("ULR 1"), paste("philtrum right","1", 2:6), rep(NA,13),
                           paste("lower face right", 22:27, "4"), rep(NA, 13),
                           paste("upper lip left", 1:6, "3"), rep(NA,13),
                           paste("nose right", "1", 1:5), rep(NA,14)),
                         ncol=15)
   
   points       <- c(3, 3, 4, 7, 10, 6, 7, 5, 19, 19, 12, 6, 6, 6, 5)
   triples.ij   <- NULL
   
   for(k in 1:15) {
      pts     <- points[k]
      indices <- c(match(pts.names.1[c(1:pts), k], rownames(object)),
                   match(pts.names.2[c(1:pts), k], rownames(object)))
      paths   <- length(indices) / pts
      trips   <- matrix(indices, ncol = paths)
      for (i in 1:(paths-1)){
         for (j in 1:(pts-1)){
            triples.ij <- c(triples.ij, 
                            c( trips[j,i], trips[j,i+1],trips[j+1,i+1], 
                               trips[j+1, i+1] , trips[j+1,i],trips[j,i]))
         }
      }
   }  
   
   
   
   # FOR j+1,i+1
   pts.names.1 <- matrix( c(paste("ULL", 1:2), paste("upper lip left 1 1"), rep(NA, 16), 
                            paste("LLR", 1:2), paste("lower lip right 1 1"), rep(NA, 16),
                            paste("lower face right", 19:22, "4"), rep(NA, 15),
                            paste("ULR 1"), paste("upper lip right", 1:6, "1"), rep(NA, 12),
                            paste("philtrum right 7", 6:2), paste("upper lip right", "6", 1:3), paste("lower lip right","6", 2:1), rep(NA, 9), 
                            paste("upper face left 3", 1:6, "2"), rep(NA, 13),
                            paste("upper face right 3", "6", 1:2), paste("nose right 6", 5:1), rep(NA, 12),
                            paste("mid-line chin", 1:4), paste("mid-line bottom lip 3"), rep(NA, 14),
                            paste("mid-face right",1:19, "5"),
                            paste("lower face right",1:19, "4"), 
                            paste("upper mid face right",1:12, "5"), rep(NA,7),
                            paste("mid-face left","19", 1:5),paste("upper mid face left 19 1"),  rep(NA,13),
                            paste("lower lip left",1:6, "1"), rep(NA, 13),
                            paste("lower lip right",1:6, "2"), rep(NA, 13),
                            paste("upper mid face left","19",1:5), rep(NA, 14)),
                          ncol=15)
   
   
   
   pts.names.2 <- matrix(c(paste("ULL 1"), paste("MLL 2"), paste("upper lip left 1 3"), rep(NA, 16),
                           paste("ULR 1"), paste("MLR 2"), paste("upper lip right 1 3"), rep(NA, 16),
                           paste("mid-face right 19 1"), paste("ULR 1"), paste("LLR 2"), paste("lower lip right 1 1"), rep(NA, 15),
                           paste("philtrum right",1:7, "2"), rep(NA, 12),
                           paste("mid-line philtral", 1:5), paste("mid-line upper lip", 1:3), paste("mid-line bottom lip",2:3), rep(NA, 9),
                           paste("nose left", 1:6, "5"),  rep(NA, 13),
                           paste("mid-line nasal root", 1:2), paste("mid-line nasal profile", 1:5), rep(NA, 12),
                           paste("lower face right", "27", 1:4), paste("lower lip right 6 1"), rep(NA,14),
                           paste("upper mid face right", 1:19, "1"),
                           paste("mid-face right", 1:19, "1"),
                           paste("upper face right 1", 1:12, "2"), rep(NA,7),
                           paste("ULL 1"), paste("philtrum left","1", 2:6), rep(NA,13),
                           paste("lower face left", 22:27, "4"), rep(NA, 13),
                           paste("upper lip right", 1:6, "3"), rep(NA,13),
                           paste("nose left", "1", 1:5), rep(NA,14)),
                         ncol=15)
   
   points         <- c(3, 3, 4, 7, 10, 6, 7, 5, 19, 19, 12, 6, 6, 6, 5)
   triples.ijij   <- NULL
   
   for (k in 1:15){
      pts     <- points[k]
      indices <- c(match(pts.names.1[c(1:pts),k], rownames(object)),
                   match(pts.names.2[c(1:pts),k], rownames(object)))
      paths   <- length(indices)/pts
      trips   <- matrix(indices, ncol=paths)
      for (i in 1:(paths-1)){
         for (j in 1:(pts-1)){
            triples.ijij <- c(triples.ijij, 
                              c(trips[j + 1, i + 1], trips[j,     i + 1], trips[j,     i], 
                                trips[j,     i],     trips[j + 1, i],     trips[j + 1, i + 1]))	    
         }
      }
      
   }  
   
   result        <- list(coords = object,  triples = c(triples.left, triples.right, triples.ij, triples.ijij))                       
   class(result) <- "face3d"
   invisible(result)
}
   
}
