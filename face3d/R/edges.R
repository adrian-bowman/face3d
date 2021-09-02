edges <- function(shape) {

   # Find the pairs of points which define the edge pieces
   edges <- shape$triangles
   edges <- rbind(edges[ , 1:2], edges[ , 2:3], edges[ , c(1, 3)])
   edges <- cbind(pmin(edges[ , 1], edges[ , 2]),
                  pmax(edges[ , 1], edges[ , 2]))
   ecode <- paste(edges[ , 1], edges[ , 2], sep = "-")
   tbl   <- table(ecode)
   ind   <- which(tbl == 1)
   edges <- as.numeric(unlist(strsplit(names(tbl)[ind], "-")))
   edges <- matrix(edges, ncol = 2, byrow = TRUE)
   
   # Find the order of the points along each edge
   lst   <- list()
   while (nrow(edges) > 0) {
      strt  <- edges[1, 1]
      tgt   <- edges[1, 2]
      nxt   <- 1
      edg   <- strt
      lngth <- 1
      while (tgt != strt & lngth > 0) {
         edges <- edges[-nxt, , drop = FALSE]
         nxt   <- which(edges[ , 1] == tgt | edges[ , 2] == tgt)
         lngth <- length(nxt)
         if (lngth > 0) {
            edg <- c(edg, tgt)
            nxt <- nxt[1]
            tgt <- if (edges[nxt, 1] == tgt) edges[nxt, 2] else edges[nxt, 1]
         }
      }
      edg   <- c(edg, tgt)
      edges <- edges[-nxt, , drop = FALSE]
      lst[[length(lst) + 1]] <- edg
   }
   
   # Reorder the list by arclength of each edge (longest first)
   fn  <- function(ids) max(arclengths(shape$vertices[ids, ]))
   sz  <- sapply(lst, fn)
   lst <- lst[order(sz, decreasing = TRUE)]

   invisible(lst)
}
