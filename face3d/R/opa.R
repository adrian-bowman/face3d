opa.face3d <- function(from, to, carry = from, scale = TRUE, weights, triangles,
                       return.parameters = FALSE, exclude = NULL) {
   
   if (!all(dim(from) == dim(to)))
      stop("the dimensions of 'from' and 'to' do not match.")
   
   if (missing(weights)) {
      if (!missing(triangles))
         weights <- area.face3d(as.face3d(list(vertices = to, triangles = triangles)))$points
      else
         weights <- rep(1, nrow(from))
   }

   mn.from <- apply(from,   2, function(x) sum(weights  * x) / sum(weights ))
   mn.to   <- apply(to,     2, function(x) sum(weights  * x) / sum(weights ))
   x.from  <- sweep(from,   2, mn.from,   "-")
   x.to    <- sweep(to,     2, mn.to,     "-")
   x.from  <- sweep(x.from, 1, sqrt(weights), "*")
   x.to    <- sweep(x.to,   1, sqrt(weights), "*")
   A       <- crossprod(x.to, x.from) /
                 (sqrt(sum(diag(crossprod(x.to)))) * sqrt(sum(diag(crossprod(x.from)))))
   decomp  <- svd(A)
   Gamma   <- decomp$v %*% t(decomp$u)
   if (scale)
      beta <- sum(diag(t(x.to) %*% x.from %*% Gamma)) / sum(diag(t(x.from) %*% x.from))
   else
      beta <- 1
   # x.from  <- sweep(x.from, 1, sqrt(weights), "/")
   # opa     <- beta * x.from %*% Gamma
   # opa     <- sweep(opa, 2, mn.to, "+")

   # Apply the transformation to the carry object
   opa.fn <- function(x) {
      x <- sweep(x, 2, mn.from)
      if (scale) x <- beta * x
      x <- x %*% Gamma
      x <- sweep(x, 2, mn.to, "+")
   }
   if (is.face3d(carry)) {
      ind     <- sapply(carry, function(x) is.matrix(x) && ncol(x) == 3)
      ind     <- which(ind)
      exclude <- c("triangles", exclude)
      ind     <- ind[-match(exclude, names(ind))]
      opa      <- carry
      opa[ind] <- lapply(opa[ind], opa.fn) 
   }
   else
      opa  <- opa.fn(carry)

   # Return
   if (return.parameters)
      result <- list(opa = opa, mean.from = mn.from, mean.to = mn.to,
                     rotation = Gamma, scale = beta)
   else
      result <- opa
   invisible(result)
}
