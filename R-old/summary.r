summary.face3d <- function(object, ...) {
   mtch <- match(c("coords", "triples"), names(object))
   if (any(is.na(mtch)))
      stop("This object does not contain both coords and triples.\n")
   n.vertices <- nrow(object$coords)
   n.faces    <- length(object$triples) / 3
   ranges <- apply(object$coords, 2, range)
   dimnames(ranges) <- list(c("min", "max"), c("x", "y", "z"))
   mtch.si <- match(c("shape.index"), names(object))
   mtch <- c(mtch, mtch.si)
   mtch <- mtch[!is.na(mtch)]
   lst <- list(...)
   prnt <- if ("print" %in% names(lst)) lst$print else TRUE
   if (prnt) {
      cat("Number of vertices:", n.vertices, "\n")
      cat("Number of faces:", n.faces, "\n")
      cat("Range of x:", ranges[ , 1], "\n")
      cat("Range of y:", ranges[ , 2], "\n")
      cat("Range of z:", ranges[ , 3], "\n")
      if (!is.na(mtch.si))
         cat("shape index: distance =", object$si.distance, "\n")
      if (length(names(object)) > length(mtch))
         cat("Other information available: ", names(object)[-mtch], "\n")
   }
   invisible(list(n.vertices = n.vertices , n.faces = n.faces, ranges = ranges))
}
