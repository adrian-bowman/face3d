"edist" <- function(x1, x2, type = "pairwise", index = FALSE, x2curve = FALSE) {

   if (missing(x2))   x2 <- rep(0, 3)
   if (is.vector(x1)) x1 <- t(as.matrix(x1))
   if (is.vector(x2)) x2 <- t(as.matrix(x2))
   x1f3d <- is.face3d(x1)
   if (x1f3d) {
      x1.face3d <- x1
      x1        <- x1$vertices
   }
   if (is.face3d(x2)) x2 <- x2$vertices

   if (ncol(x1) != 3) stop("the dimension of x1 do not correspond to 3d co-ordinates.")
   if (ncol(x2) != 3) stop("the dimension of x2 do not correspond to 3d co-ordinates.")
   if (type != "min") {
      index   <- FALSE
      x2curve <- FALSE
   }
   if (x2curve & !x1f3d)
      stop("when x2 is a curve, x1 must be a face3d object.")
   
   dst <- fields::rdist(x1, x2)
   if (type == "pairwise") {
      if (any(dim(dst) == 1)) dst <- c(dst)
   }
   else if (type == "min") {
      if (index | x2curve) x2ind  <- apply(dst, 1, which.min)
      if (x2curve) {
         if (!("normals" %in% names(x1.face3d))) x1.face3d <- normals.face3d(x1.face3d)
         meshpt <- apply(dst, 2, which.min)
      }
      dst <- apply(dst, 1, min)
      if (x2curve) {
         # Remove trailing digits
         if (is.null(rownames(x2))) rownames(x2) <- paste("point", 1:nrow(x2))
         cnms <- rownames(x2)
         fn   <- function(x) {
                    while (nchar(x) > 0 && grepl("[0-9]", substr(x, nchar(x), nchar(x))))
                       x <- substr(x, 1, nchar(x) - 1)
                    x
         }
         cnms <- sapply(cnms, fn)
         if (any(nchar(cnms) == 0)) stop("the rownames of x2 cannot be all integer.")
         # Create a differencing matrix for vectors along the curve
         cnms <- unique(cnms)
         dff.fn <- function(x) {
            nms <- grep(x, rownames(x2))
            lng <- length(nms)
            cbind(nms[c(2:lng, lng)], nms[c(1, 1:(lng - 1))])
         }
         diffmat <- matrix(ncol = 2, nrow = 0)
         for (nm in cnms) diffmat <- rbind(diffmat, dff.fn(nm))
         # Project the distances onto the curve vectors
         normals <- x1.face3d$normal[meshpt, ]
         x2vec   <- x2[diffmat[ , 1], ] - x2[diffmat[ , 2], ]
         x1vec   <- x1 - x2[x2ind, ]
         axis1   <- crossproduct(normals[x2ind, ], x2vec[x2ind, ])
         dst     <- apply(x1vec * axis1, 1, sum)
      }
      if (index) dst <- list(distance = dst, index = x2ind)
   }
   else if (type == "minsum") {
      dst <- sum(apply(dst, 1, min))
   }
   else
      stop("invalid value of the type argument.")
   invisible(dst)
}
