warp.face3d <- function(from, to, carry = from) {

  if (!requireNamespace("MASS", quietly = TRUE))
    stop("the MASS package is not available.")
  if (!is.matrix(from)) 
    stop("from must be a matrix.")
  if (!is.matrix(to) & !is.vector(to)) 
    stop("to must be a matrix.")
  if (is.vector(to))
    to <- as.matrix(to)
  if (!is.matrix(carry) && !is.list(carry)) 
     stop("carry is not a matrix or a list.")
  if (nrow(from) != nrow(to)) 
    stop("from and to must have the same number of rows.")
  # if (ncol(from) != ncol(to)) 
  #    stop("from and to must have the same number of columns.")
   
  if (is.list(carry)) ID_carry_mat <- FALSE
  if (is.matrix(carry)) {
    ID_carry_mat <- TRUE
    carry <- list(carry)
  }
  
  clr <- NULL
  triples.present <- ("triples" %in% names(carry))
  if (is.list(carry) && (ncol(to) == 3) && ("coords" %in% names(carry)) && (triples.present)) {
    trpls  <- carry$triples
    clr    <- if ("colour" %in% names(carry)) carry$colour else NULL
    carry1 <- list()
    if ("lmks" %in% names(carry))   carry1$lmks   <- carry$lmks
    if ("curves" %in% names(carry)) carry1$curves <- carry$curves
    if ("mesh" %in% names(carry))   carry1$mesh   <- carry$mesh
    if ("coords" %in% names(carry)) carry1$coords <- carry$coords
    carry <- carry1
    to_object <- list()
    if ("lmks" %in% names(carry))   to_object$lmks   <- carry$lmks
    if ("curves" %in% names(carry)) to_object$curves <- carry$curves
    if ("mesh" %in% names(carry))   to_object$mesh   <- carry$mesh
    if ("coords" %in% names(carry)) to_object$coords <- carry$coords
  }
  
  if (is.list(carry)) to_object <- carry
  
  Gi <- ginvtps.face3d(from)
  tps.coeffs <- Gi %*% rbind(to, matrix(0, 4, ncol(to))) 
  k  <- length(carry)
  k1 <- dim(from)[1]
  for (i in 1:k) {
    if (is.matrix(carry[[i]]) && (ncol(from) == ncol(carry[[i]]))) {
      k2 <- nrow(carry[[i]]) 
      affine <- rep(1, k2) %o% tps.coeffs[k1 + 1, ] + 
                carry[[i]][, 1] %o% tps.coeffs[k1 + 2, ] + 
                carry[[i]][, 2] %o% tps.coeffs[k1 + 3, ] + 
                carry[[i]][, 3] %o% tps.coeffs[k1 + 4, ]
      dx <- outer(from[, 1], carry[[i]][, 1], "-")
      dy <- outer(from[, 2], carry[[i]][, 2], "-")
      dz <- outer(from[, 3], carry[[i]][, 3], "-")
      phi <- sqrt(dx * dx + dy * dy + dz * dz)
      nonaff <- t(phi) %*% tps.coeffs[1:k1, ]
      to_object[[i]] <- affine + nonaff 
    }
    if (ncol(from) != ncol(carry[[i]])) {
      to_object[[i]] <- NA
      print(paste(i,"-th element of the list doesn't have ",ncol(from), " columns.", sep = ""))
    }
    if (!is.matrix(carry[[i]])) {
      to_object[[i]] <- NA
      print(paste(i,"-th element of the list is not a matrix.", sep = ""))
    }
    
  }
  
  if ("coords" %in% names(carry) & triples.present) {
     to_object$triples <- trpls
     class(to_object) <- "face3d"
  }
  if (!is.null(clr)) to_object$colour <- clr
  
  if (ID_carry_mat) to_object <- to_object[[1]]
  if (is.matrix(to_object) && ncol(to_object) == 1) to_object <- c(to_object)
  return(to_object)
}
