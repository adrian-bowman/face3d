normals.face3d <- function(shape) {
   
   if (!("face3d" %in% class(shape))) stop("shape is not a face3d object.")

   clist <- function(list1, list2) {
      nms1         <- sapply(list1, length)
      nms2         <- sapply(list2, length)
      nms1         <- rep(names(nms1), nms1)
      nms2         <- rep(names(nms2), nms2)
      list1        <- unlist(list1)
      list2        <- unlist(list2)
      names(list1) <- nms1
      names(list2) <- nms2
      list2        <- c(list1, list2)
      nms          <- names(list2)
      names(list2) <- NULL
      list2        <- tapply(list2, nms, c)
      list2        <- list2[order(as.numeric(names(list2)))]
      list2        <- lapply(list2, unique)
      list2
   }
  
   local.axes <- function(trngs) {
   	  d21    <- shape$vertices[triples[trngs, 2], ] - shape$vertices[triples[trngs, 1], ]
   	  d31    <- shape$vertices[triples[trngs, 3], ] - shape$vertices[triples[trngs, 1], ]
   	  normal <- apply(crossproduct(d21, d31), 2, mean)
   	  normal <- normal / sqrt(sum(normal^2))
   	  x1     <- crossproduct(normal, c(t(d21))[1:3])
   	  x2     <- crossproduct(normal, x1)
   	  cbind(normal, c(x1), c(x2))
   }

   triples       <- shape$triangles
   triangles     <- 1:nrow(triples)
   trngs         <-              tapply(triangles, triples[ , 1], c)
   trngs         <- clist(trngs, tapply(triangles, triples[ , 2], c))
   trngs         <- clist(trngs, tapply(triangles, triples[ , 3], c))
   axes          <- lapply(trngs, local.axes)
   shape$axes    <- axes
   axes          <- array(unlist(axes), dim = c(3, 3, length(axes)),
                          dimnames = list(c("x", "y", "z"), c("normal", "axis 1", "axis 2"), NULL))
   axes          <- aperm(axes, c(3, 1, 2))
   shape$normals <- axes[ , , 1]

   invisible(shape)
}
