interpbary <- function(X, f, triangles, Xnew) {
   # 2D barycentric interpolation at points Xnew of values f measured at locations X.
   # For each point in Xnew, tsearch identifies the surrounding triangle.
   # The interpolation is constructed as a sparse matrix operation for speed.
   if (!requireNamespace("Matrix", quietly = TRUE)) stop("the Matrix package is required.")
   tri    <- tsearch(X[ , 1], X[ , 2], triangles, Xnew[ , 1], Xnew[ , 2], bary = TRUE)
   active <- triangles[tri$idx, ]
   ind    <- !is.na(rowSums(active))
   active <- active[ind, ]
   Xnew   <- Xnew[ind, ]
   tri$p  <- tri$p[ind, ]
   M      <- Matrix::sparseMatrix(i = rep(1:nrow(Xnew), each = 3), j = as.numeric(t(active)),
                                  x = as.numeric(t(tri$p)), dims = c(nrow(Xnew), length(f)))
   list(Xnew = Xnew, fnew = as.numeric(M %*% f))
}

