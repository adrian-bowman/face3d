opa.face3d <- function(from, to, scale = TRUE, weights, match.subset, model.mesh = TRUE) {
   
   if (!all(dim(from) == dim(to)))
      stop("the dimensions of 'from' and 'to' do not match.")
   
   weights.missing <- missing(weights)
   if (!weights.missing) model.mesh <- TRUE
   if (model.mesh | missing(match.subset))
      match.subset <- 1:nrow(from)
   else
      if (is.logical(match.subset)) match.subset <- which(match.subset)
   else
      if (!is.integer(match.subset))
         stop("'match.subset' is neither logical nor integer.")
   else
      if (all(match.subset < 0)) match.subset <- (1:nrow(from))[match.subset]

   if (length(match.subset) > nrow(from))
      stop("the length of 'match.subset' is greater than the first dimension of 'from'.")
   if (min(match.subset) < 1 | max(match.subset) > nrow(from))
      stop("some values in match.subset are out of range.")
   
   if (model.mesh) {
      wts     <- if (weights.missing) area.face3d(as.face3d(to, model.mesh = TRUE))$points else weights
      mn.from <- apply(from,   2, function(x) sum(wts  * x) / sum(wts ))
      mn.to   <- apply(to,     2, function(x) sum(wts  * x) / sum(wts ))
      x.from  <- sweep(from,   2, mn.from,   "-")
      x.to    <- sweep(to,     2, mn.to,     "-")
      x.from  <- sweep(x.from, 1, sqrt(wts), "*")
      x.to    <- sweep(x.to,   1, sqrt(wts), "*")
      # result  <- procOPA(x.to, x.from, scale = scale)$Bhat
      A       <- crossprod(x.to, x.from) / (sqrt(sum(diag(crossprod(x.to)))) * sqrt(sum(diag(crossprod(x.from)))))
      decomp  <- svd(A)
      Gamma   <- decomp$v %*% t(decomp$u)
      beta    <- sum(diag(t(x.to) %*% x.from %*% Gamma)) / sum(diag(t(x.from) %*% x.from))
      x.from  <- sweep(x.from, 1, sqrt(wts), "/")
      opa     <- beta * x.from %*% Gamma
      opa     <- sweep(opa, 2, mn.to, "+")
   }
   else {
      opa     <- procOPA(to[match.subset, ], from[match.subset, ], scale = scale)
      Gamma   <- opa$R
      beta    <- opa$s
      mn.from <- apply(from[match.subset, ], 2, mean)
      mn.to   <- apply(  to[match.subset, ], 2, mean)
      result  <- sweep(from, 2, mn.from)
      result  <- result %*% Gamma
      if (scale) result  <- beta * result
      opa     <- sweep(result, 2, mn.to, "+")
   }
   rownames(opa) <- rownames(from)
      
   invisible(opa)
}
