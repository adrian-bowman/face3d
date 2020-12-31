opaold.face3d <- function(from, to, carry = from, scale = TRUE, weights, model.mesh = FALSE,
                       return.parameters = FALSE) {
   
   if (!all(dim(from) == dim(to)))
      stop("the dimensions of 'from' and 'to' do not match.")
   
   if (missing(weights)) {
      if (model.mesh) weights <- area.face3d(as.face3d(to, model.mesh = TRUE))$points
      else            weights <- rep(1, nrow(from))
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

   result  <- sweep(carry, 2, mn.from)
   result  <- result %*% Gamma
   if (scale) result  <- beta * result
   opa     <- sweep(result, 2, mn.to, "+")
   rownames(opa) <- rownames(carry)
   
   if (return.parameters)
      result <- list(opa = opa, mean.from = mn.from, mean.to = mn.to,
                     rotation = Gamma, scale = beta)
   else
      result <- opa
   invisible(result)
}
