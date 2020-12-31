summary.face3d <- function(object, checks = FALSE, fix = FALSE, ...) {

   if (!is.face3d(object)) stop("this is not a face3d object.")
   checks <- checks | fix

   n.vertices <- nrow(object$vertices)
   n.faces <- nrow(object$triangles)
   ranges  <- apply(object$vertices, 2, range)
   dimnames(ranges) <- list(c("min", "max"), c("x", "y", "z"))

   lst     <- list(...)
   prnt    <- if ("print" %in% names(lst)) lst$print else TRUE
   if (prnt) {
      cat("Number of vertices: ", n.vertices, "\n")
      cat("Number of triangles:", n.faces, "\n")
      cat("Range of x:", ranges[ , 1], "\n")
      cat("Range of y:", ranges[ , 2], "\n")
      cat("Range of z:", ranges[ , 3], "\n")
      if ("shape index" %in% names(object))
         cat("shape index: distance =", object$si.distance, "\n")
      more <- which(!(names(object) %in% c("vertices", "triangles", "shape index")))
      if (length(more) > 0)
         cat("Other information available: ", names(object)[more], "\n")
   }
   result <- list(n.vertices = n.vertices , n.faces = n.faces, ranges = ranges)
   
   if (any(is.na(object$triangles)))
      cat("Warning: the triangles contain missing values.")
   if (min(c(object$triangles), na.rm = TRUE) < 1)
      cat("Warning: the smallest index in triangles is < 1.\n")
   if (max(c(object$triangles), na.rm = TRUE) > n.vertices)
      cat("Warning: the largest index in triangles is greater than the number of vertices.\n")
   
   if (!checks) return(invisible(result))
   
   problem <- FALSE
   
   # Check for isolated vertices which are not in the triangulation
   isolated <- (1:n.vertices)[-unique(c(t(object$triangles)))]
   if (length(isolated) > 0) {
      result$isolated <- isolated
      if (fix)
         result$object <- subset(object, -isolated)
      else {
         if (prnt) {
            cat("Warning: this object contains isolated vertices.\n")
            cat("  The indices are returned in component 'isolated.\n")
         }
         problem <- TRUE
      }
   }
   
   # Check for duplicated co-ordinates
   ind <- which(duplicated(object$vertices) | duplicated(object$vertices, fromLast = TRUE))
   if (length(ind) > 0) {
      ind1      <- duplicated(object$coord[ind, ])
      retained  <- ind[which(!ind1)]
      discarded <- ind[which(ind1)]
      map <- matrix(nrow = 0, ncol = 2)
      for (i in retained) {
         ind1 <- which(duplicated(object$vertices[c(i, discarded), ])) - 1
         # ind1 <- which(edist.face3d(object$vertices[discarded, ], object$vertices[i, ]) < 100 * .Machine$double.eps)
         map  <- rbind(map, cbind(discarded[ind1], i))
      }
      result$duplicated <- map[ , 1]
      result$matched    <- map[ , 2]
      if (fix) {
         triples <- c(t(object$triangles))
         for (i in 1:length(result$duplicated))
            triples[triples == result$duplicated[i]] <- result$matched[i]
         trpls          <- matrix(triples, ncol = 3, byrow = TRUE)
         ind            <- apply(trpls, 1, function(x) length(unique(x)) < 3)
         object$triangles <- matrix(c(t(trpls[!ind, ])), ncol = 3, byrow = TRUE)
         object         <- subset(object, -result$duplicated)
         result$object  <- object
      }
      else {
         if (prnt) {
            cat("Warning: this object contains duplicated co-ordinates.\n")
            cat("  The indices are returned in component 'duplicated'.\n")
         }
         problem <- TRUE
         }
   }
   
   trpls <- object$triangles
   # This test of collinearity misses some - alternative below
   # fn <- function(a) {
   #    ln <- c(sqrt((object$vertices[a[1], ] - object$vertices[a[2], ])^2),
   #            sqrt((object$vertices[a[1], ] - object$vertices[a[3], ])^2),
   #            sqrt((object$vertices[a[2], ] - object$vertices[a[3], ])^2))
   #    mx <- which.max(ln)
   #    ln[mx] == sum(ln[-mx])
   # }
   fn <- function(a) {
      d21    <- object$vertices[a[2], ] - object$vertices[a[1], ]
      d31    <- object$vertices[a[3], ] - object$vertices[a[1], ]
      any(is.nan(crossproduct(d21, d31)))
   }
   ind <- apply(trpls, 1, fn)
   ind <- which(ind)
   if (length(ind) > 0) {
      result$collinear <- trpls[ind, ]
      dimnames(result$collinear) <- list(
         paste("triangle", 1:length(ind)), paste("vertex", 1:3))
      if (fix) {
         object$triangles <- trpls[-ind, ]
         result$object  <- object
      }
      else {
         if (prnt) {
            cat("Warning: this object contains collinear triangles.\n")
            cat("  These are returned in component 'collinear'.\n")
         }
         problem <- TRUE
      }
   }
   
   if (prnt) {
      if (problem) cat("Rerun the 'summary' function with 'fix = TRUE' to address these issues.\n")
      else         cat("All checks passed.\n")
   }
   if (fix) result <- result$object
   
   invisible(result)
}
