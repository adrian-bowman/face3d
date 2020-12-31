normals.face3d <- function(shape) {

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
   	  d21    <- shape$coords[triples[trngs, 2], ] - shape$coords[triples[trngs, 1], ]
   	  d31    <- shape$coords[triples[trngs, 3], ] - shape$coords[triples[trngs, 1], ]
   	  normal <- apply(crossproduct(d21, d31), 2, mean)
   	  normal <- normal / sqrt(sum(normal^2))
   	  x1     <- crossproduct(normal, c(t(d21))[1:3])
   	  x2     <- crossproduct(normal, x1)
   	  cbind(normal, c(x1), c(x2))
   }

   triples   <- matrix(shape$triples, ncol = 3, byrow = TRUE)
   triangles <- 1:nrow(triples)
   trngs     <-              tapply(triangles, triples[ , 1], c)
   trngs     <- clist(trngs, tapply(triangles, triples[ , 2], c))
   trngs     <- clist(trngs, tapply(triangles, triples[ , 3], c))
   trngs1    <- trngs

   axes          <- lapply(trngs1, local.axes)
   axes          <- array(unlist(axes), dim = c(3, 3, nrow(shape$coords)))
   shape$normals <- t(axes[ , 1, ])

   class(shape) <- "face3d"
   invisible(shape)
}
