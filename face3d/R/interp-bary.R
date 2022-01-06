interp.barycentric <- function(X, f, Xnew) {
   # 2D barycentric interpolation at points Xnew of values f measured at locations X.
   # For each point in Xnew, geometry::tsearch identifies the surrounding triangle.
   # The interpolation is constructed as a sparse matrix operation for speed.
   if (!requireNamespace("geometry", quietly = TRUE)) stop("the geometry package is required.")
   if (!requireNamespace("Matrix", quietly = TRUE)) stop("the Matrix package is required.")
   dn     <- geometry::delaunayn(X)
   tri    <- geometry::tsearch(X[ , 1], X[ , 2], dn, Xnew[ , 1], Xnew[ , 2], bary = TRUE)
   active <- dn[tri$idx, ]
   ind    <- !is.na(rowSums(active))
   active <- active[ind, ]
   Xnew   <- Xnew[ind, ]
   tri$p  <- tri$p[ind, ]
   M      <- Matrix::sparseMatrix(i = rep(1:nrow(Xnew), each = 3), j = as.numeric(t(active)),
                          x = as.numeric(t(tri$p)), dims = c(nrow(Xnew), length(f)))
   list(Xnew = Xnew, fnew = as.numeric(M %*% f))
}

# X     <- cbind(rnorm(250), rnorm(250))
# f     <- X[ , 1]^2 + X[ , 2]^2
# xgrid <- runif(250, -1.5, 1.5)
# ygrid <- runif(250, -1.5, 1.5)
# Xnew  <- cbind(xgrid, ygrid)
# intrp <- interp.barycentric(X, f, Xnew)

# library(rpanel)
# rp.plot3d(c(X[,1], intrp$Xnew[,1]), c(f, intrp$fnew), c(X[,2], intrp$Xnew[,2]),
        # col = rep(1:2, c(nrow(X), nrow(intrp$Xnew))))
