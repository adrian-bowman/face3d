index.face3d <- function(shape, extent = 2, distance = 10, type = "euclidean",
                         subset = 1:nrow(shape$vertices),
                         extension = "NA", overwrite = FALSE, directions = FALSE) {

   sbst <- subset
   if (all(is.logical(subset)) & (length(sbst) == nrow(shape$vertices))) sbst <- which(sbst)
   if (length(sbst) == 1 && sbst <= 1) sbst <- round(nrow(shape$vertices) * sbst)

   if (!("shape.index" %in% names(shape))) shape$shape.index <- rep(NA, nrow(shape$vertices))
   if (!overwrite) sbst <- sbst[is.na(shape$shape.index[sbst])]
   if (length(sbst) == 0) return(invisible(shape))

   extension.missing <- missing(extension)
  
   clist <- function(list1, list2) {
      nms1         <- sapply(list1, length)
      nms2         <- sapply(list2, length)
      nms1         <- rep(names(nms1), nms1)
      nms2         <- rep(names(nms2), nms2)
      list1        <- unlist(list1)
      list2        <- unlist(list2)
      names(list1) <- nms1
      names(list2) <- nms2
      list2        <- c(list1, list2)
      nms          <- names(list2)
      names(list2) <- NULL
      list2        <- tapply(list2, nms, c)
      list2        <- list2[order(as.numeric(names(list2)))]
      list2        <- lapply(list2, unique)
      list2
   }
   
   index.fn <- function(i, axes, directions) {
   	  pt.i  <- shape$vertices[sbst[i], ]
   	  if (type == "mesh") {
   	     nbrs <- unique(c(triples[trngs1[[i]], ]))
           if (extent > 1) {
              for (i in 2:extent) {
                 trngs <- lapply(nbrs,  function(x) unique(trngs1[[x]]))
                 nbrs  <- unique(unlist(trngs))
                 nbrs  <- unique(c(triples[nbrs, ]))
              }
           }
   	  }
   	  else
           nbrs  <- which(edist.face3d(pts.ext, pt.i) < distance)
   	  if (length(nbrs) < 7) {
   	     shape.index        <- NA
   	     kappa              <- rep(NA, 2)
   	     rss                <- NA
   	     mean.curvature     <- NA
   	     gaussian.curvature <- NA
   	     if (directions) drns <- c(matrix(NA, 2, 2))
   	  }
   	  else {
   	     y                    <- c(sweep(pts.ext[nbrs, ], 2, pt.i) %*% axes[ , 1])
   	     x1                   <- c(sweep(pts.ext[nbrs, ], 2, pt.i) %*% axes[ , 2])
   	     x2                   <- c(sweep(pts.ext[nbrs, ], 2, pt.i) %*% axes[ , 3])
   	     fit                  <- lsfit(cbind(x1^2/2, x1*x2, x2^2/2, x1^3, x1^2*x2, x1*x2^2, x2^3), y,
   	                                      intercept = FALSE)
   	     beta                 <- fit$coef
   	     rss                  <- sum((fit$residuals)^2)
           W                    <- matrix(beta[c(1, 2, 2, 3)], 2, 2)
           eig                  <- eigen(W, only.values = !directions)
           kappa                <- eig$values
   	     if (directions) drns <- c(eig$vectors)
           # kappa                <- rev(eig$values)
   	     # if (directions) drns <- eig$vectors[ , 2:1]
           shape.index          <- 2 / pi * (atan((kappa[2] + kappa[1]) / (kappa[2] - kappa[1])))
           # mean.curvature       <- det(W)
           # gaussian.curvature   <- 1/2 * (sum(diag(W)))
   	  }
      results <- list(shape.index = shape.index, kappa1 = kappa[1], kappa2 = kappa[2], rss = rss)
      if (directions) results$drns <- drns
      invisible(results)
   }
   
   if (overwrite | !("axes" %in% names(shape))) shape <- normals.face3d(shape)
   axes <- shape$axes
#    vec      <- apply(shape$normals, 1, function(x) which.min(x)[1])
#    mat      <- matrix(0, nrow = length(vec), ncol = 3)
#    ind      <- cbind(1:length(vec), vec)
#    # id <- which(sapply(vec, length) != 1)
#    # print(id)
#    # print(vec[[id]])
#    mat[ind] <- 1
#    a2       <- crossproduct(shape$normals, mat)
#    a3       <- crossproduct(shape$normals, a2)
#    axes     <- lapply(1:length(vec), function(i) cbind(shape$normals[i, ], a2[i, ], a3[i, ]))

   # sbst <- subset
   # if (all(is.logical(subset)) & (length(sbst) == nrow(shape$vertices))) sbst <- which(sbst)
   # if (length(sbst) == 1 && sbst <= 1) sbst <- round(nrow(shape$vertices) * sbst)
   pts <- shape$vertices[sbst, ]
   
   # Deal with extension
   if (extension.missing) {
      if (length(sbst) < 10) extensn <- TRUE
      else {
        ind.ext <- sbst
        extensn <- FALSE
      }
   }
   else
      extensn <- extension
   if (is.logical(extensn) && extensn) {
      dst     <- rdist(shape$vertices[sbst, ], shape$vertices)
      ind.ext <- apply(dst, 1, function(x) which(x <= distance))
      ind.ext <- unique(c(unlist(ind.ext)))
   }
   else if (is.integer(extensn)) ind.ext <- c(sbst, extensn)
   pts.ext <- shape$vertices[ind.ext, ]
   
   if (type == "mesh") {
   	triples   <- shape$triangles
      triangles <- 1:nrow(triples)
      trngs1    <-               tapply(triangles, triples[ , 1], c)
      trngs1    <- clist(trngs1, tapply(triangles, triples[ , 2], c))
      trngs1    <- clist(trngs1, tapply(triangles, triples[ , 3], c))
   }
   
   if (!("kappa1"      %in% names(shape))) shape$kappa1      <- rep(NA, nrow(shape$vertices))
   if (!("kappa2"      %in% names(shape))) shape$kappa2      <- rep(NA, nrow(shape$vertices))
   if (!("rss"         %in% names(shape))) shape$rss         <- rep(NA, nrow(shape$vertices))
   index                   <- mapply(index.fn, 1:length(sbst), axes[sbst], directions)
   index                   <- matrix(unlist(index), nrow = length(sbst), byrow = TRUE)
   shape$shape.index[sbst] <- index[ , 1]
   shape$kappa1[sbst]      <- index[ , 2]
   shape$kappa2[sbst]      <- index[ , 3]
   shape$rss[sbst]         <- index[ , 4]
   shape$si.distance       <- distance
   if (directions & (overwrite | !("directions" %in% names(shape)))) {
      axes       <- array(unlist(axes), dim = c(3, 3, length(axes)))
      axes       <- axes[ , , sbst]
      # shape$directions <- array(NA, dim = c(2, 2, nrow(shape$vertices)))
      # shape$directions[ , , sbst] <- array(c(t(index[ , 5:8])), dim = c(2, 2, length(sbst)))
      # a <- shape$directions[1, 1, sbst]
      # b <- shape$directions[2, 1, sbst]
      dirns <- array(NA, dim = c(2, 2, nrow(shape$vertices)))
      dirns[ , , sbst] <- array(c(t(index[ , 5:8])), dim = c(2, 2, length(sbst)))
      a <- dirns[1, 1, sbst]
      b <- dirns[2, 1, sbst]
      if (!("directions" %in% names(shape)))
         shape$directions <- array(NA, dim = c(3, 2, nrow(shape$vertices)))
      shape$directions[ , 1, sbst] <- t(a * t(axes[ , 2, ]) + b * t(axes[ , 3, ]))
      shape$directions[ , 2, sbst] <- t(b * t(axes[ , 2, ]) - a * t(axes[ , 3, ]))
      # print(a * t(shape$axes[ , 2, ]) + b * t(shape$axes[ , 3, ]))
      # print(shape$directions[ , 1, ])
      # print(a^2 + b^2)
      # print(apply(shape$axes[ , 2, ], 2, function(x) sqrt(sum(x^2))))
      # print(apply(shape$axes[ , 3, ], 2, function(x) sqrt(sum(x^2))))
      # print(apply(shape$directions[, 1, ], 2, function(x) sqrt(sum(x^2))))
   }
   
   class(shape) <- "face3d"
   invisible(shape)
}
