sample.face3d <- function(x, spacing) {

   #     Create a sample of vertices with minimum spacing from a shape or set of vertices
   
   if (is.face3d(x))
      x <- x$vertices
   else
      if (!(is.matrix(x) && ncol(x) == 3))
         stop("x must be a face3d object or a three-column matrix of vertices.")
   
   nvert   <- nrow(x)
   sampled <- as.integer(1)
   mindist <- edist(x[sampled, ], x)
   while(max(mindist) > spacing & length(sampled) < nvert) {
      isampled <- which.max(mindist)
      sampled  <- c(sampled, isampled)
      idist    <- edist(x[isampled, ], x)
      mindist  <- pmin(idist, mindist)
   }
   
   invisible(sampled)
   
}
