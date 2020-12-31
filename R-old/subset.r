subset.face3d <- function(x, ind, remove.singles = TRUE, ...) {

   shape <- x
   if (is.logical(ind)) ind <- which(ind)
   if (any(ind < 0)) {
      if (any(ind >= 0))
         stop("the subset definition contains both positive and non-negative indices.")
      ind <- (1:nrow(shape$coords))[ind]
   }
   if (length(ind) == 0) stop("the subset is empty.")
   sbst <- list(coords = matrix(shape$coords[ind, ], ncol = 3,
                           dimnames = list(as.character(ind))))
   
   if ("triples" %in% names(shape)) {
      trpls      <- matrix(shape$triples, ncol = 3, byrow = TRUE)
      indt       <- (trpls[ , 1] %in% ind) & (trpls[ , 2] %in% ind) & (trpls[ , 3] %in% ind)
      trpls      <- c(t(trpls[indt, ]))
      vec        <- 1:length(ind)
      nms        <- format(ind, scientific = FALSE)
      nms        <- sub("^ +", "", nms)                   # remove leading zeroes 
      names(vec) <- nms
      nms1       <- format(trpls, scientific = FALSE)
      nms1       <- sub("^ +", "", nms1)
      sbst[["triples"]] <- vec[nms1]
   }

   if ("colour" %in% names(shape))
      sbst[["colour"]] <- shape$colour[ind]
      
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

   ind <- as.numeric(rownames(sbst$coords))
   for (comp in c("shape.index", "kappa1", "kappa2"))
      if (comp %in% names(shape)) sbst[[comp]] <- shape[[comp]][ind]
   if ("normals" %in% names(shape))
      sbst$normals <- shape$normals[ind, ]
   if ("local.axes" %in% names(shape))
      sbst$local.axes <- shape$local.axes[ , , ind]

   if ("lndms" %in% names(shape))
      sbst$lndms <- shape$lndms
   if ("lmks" %in% names(shape))
      sbst$lmks <- shape$lmks
      
   nms <- which(!(names(shape) %in% c("coords", "triples", "colour", "shape.index", "kappa1", "kappa2",
                                "normals", "local.axes")))
   nms <- names(shape)[nms]
   if (length(ind) > 0) {
      for (i in nms) sbst[[i]] <- shape[[i]]
   }
   
   class(sbst) <- "face3d"
   sbst
}
