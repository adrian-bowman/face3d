trimimage.face3d <- function(shape, lmks.vertical = c("sn", "n"),
                                    lmks.horizontal = c("exR", "exL")) {

   # Reduce an image using landmarks as a guide

   if (!all(c("lmks", "mesh") %in% names(shape)))
      stop("both lmks and a mesh are required in shape.")
   
   shape <- rotate.face3d(shape, id.lndms = c(lmks.vertical, lmks.horizontal),
                         rotation = "coronal")
   ind <- shape$vertices[ , 3] > min(shape$mesh[ , 3])
   ind <- ind & (shape$vertices[ , 1] >= min(shape$mesh[ , 1]))
   ind <- ind & (shape$vertices[ , 1] <= max(shape$mesh[ , 1]))
   ind <- ind & (shape$vertices[ , 2] >= min(shape$mesh[ , 2]) - 0.1 * diff(range(shape$mesh[ , 2])))
   ind <- ind & (shape$vertices[ , 2] <= max(shape$mesh[ , 2]) + 0.1 * diff(range(shape$mesh[ , 2])))
   shape <- subset.face3d(shape, ind)
   shape
}
