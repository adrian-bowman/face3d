subset.face3d <- function(shape, subset, remove.singles = TRUE, retain.indices = FALSE, ...) {

   ind   <- subset
   nvert <- nrow(shape$coords)
   if (is.logical(ind)) ind <- which(ind)
   if (any(ind < 0)) {
      if (any(ind >= 0))
         stop("the subset definition contains both positive and non-negative indices.")
      ind <- (1:nvert)[ind]
   }
   if (length(ind) == 0) stop("the subset is empty.")
   
   dnms        <- dimnames(shape$coords)
   rnms        <- if (is.null(dnms[[1]])) as.character(1:nvert) else dnms[[1]]
   sbst        <- shape
   sbst$coords <- matrix(c(shape$coords[ind, ]), ncol = 3, dimnames = list(rnms[ind], dnms[[2]]))
   
   if ("triples" %in% names(shape)) {
      trpls        <- matrix(shape$triples, ncol = 3, byrow = TRUE)
      indt         <- (trpls[ , 1] %in% ind) & (trpls[ , 2] %in% ind) & (trpls[ , 3] %in% ind)
      trpls        <- c(t(trpls[indt, ]))
      vec          <- 1:length(ind)
      nms          <- format(ind, scientific = FALSE)
      nms          <- sub("^ +", "", nms)                   # remove leading zeroes 
      names(vec)   <- nms
      nms1         <- format(trpls, scientific = FALSE)
      nms1         <- sub("^ +", "", nms1)
      sbst$triples <- vec[nms1]
   }

   if ("colour" %in% names(shape))
      sbst$colour <- shape$colour[ind]
      
   if (remove.singles & ("triples" %in% names(shape))) {
      ind1         <- (vec %in% unique(sbst$triples))
      vec1         <- 1:length(ind1)
      nms          <- format((1:length(ind))[ind1], scientific = FALSE)
      nms          <- sub("^ +", "", nms)
      names(vec1)  <- nms
      nms1         <- format(sbst$triples, scientific = FALSE)
      nms1         <- sub("^ +", "", nms1)
      sbst$triples <- vec1[nms1]
      sbst$coords  <- sbst$coords[ind1, ]
      if ("colour" %in% names(shape))
         sbst$colour <- sbst$colour[ind1]
      ind <- ind[ind1]
   }

   # Subset information which has a dimension matching length of ind
   inms <- 1:length(shape)
   inms <- inms[-match(c("coords", "triples"), names(shape))]
   if ("colour" %in% names(shape)) inms <- inms[inms != match("colour", names(shape))]
   for (i in inms) {
      if (is.vector(shape[[i]]) & length(shape[[i]]) == nvert)
         sbst[[i]] <- shape[[i]][ind]
      else if (is.array(shape[[i]]) & (sum(dim(shape[[i]]) == nvert) == 1)) {
         dm  <- dim(shape[[i]])
         idm <- which(dm == nvert)
         lst <- list()
         for (j in 1:length(dm)) lst[[j]] <- 1:dm[j]
         lst[[idm]] <- ind
         vals       <- shape[[i]][as.matrix(expand.grid(lst))]
         dm[idm]    <- length(ind)
         dimnms <- dimnames(shape[[i]])
         if (!is.null(dimnms) && !is.null(dimnms[[idm]])) dimnms[[idm]] <- dimnms[[idm]][ind]
         sbst[[i]]  <- array(vals, dim = dm, dimnames = dimnms)
      }
      else if (is.list(shape[[i]] & (length(shape[[i]]) == nvert)))
         sbst[[i]] <- shape[[i]][ind]
      else
         sbst[[i]] <- shape[[i]]
   }
   
   if (retain.indices) sbst$subset <- ind
   sbst
}
