crossproduct <- function(a, b, scale = TRUE) {
   if (is.vector(a)) a  <- t(as.matrix(a))
   if (is.vector(b)) b  <- t(as.matrix(b))
   if (nrow(a) != nrow(b))
      stop("the rows of the input do not match in crossproduct.")
   if (ncol(a) != 3 | ncol(b) != 3)
      stop("the input columns should be 3 in crossproduct.")
   cp <- cbind(a[ , 2]*b[ , 3] - a[ , 3]*b[ , 2],
               a[ , 3]*b[ , 1] - a[ , 1]*b[ , 3],
               a[ , 1]*b[ , 2] - a[ , 2]*b[ , 1])
   if (scale == TRUE) cp <- cp / sqrt(apply(cp^2, 1, sum))
#   if (length(c(cp)) == 3) cp <- c(cp)
   cp
}
