edges.face3d <- function(shape) {
   edges <- matrix(shape$triples, ncol = 3, byrow = TRUE)
   edges <- rbind(edges[ , 1:2], edges[ , 2:3], edges[ , c(1, 3)])
   edges <- cbind(pmin(edges[ , 1], edges[ , 2]),
                  pmax(edges[ , 1], edges[ , 2]))
   ecode <- paste(edges[ , 1], edges[ , 2], sep = "-")
   tbl   <- table(ecode)
   ind   <- which(tbl < 2)
   edges <- as.numeric(unlist(strsplit(names(tbl)[ind], "-")))
   matrix(edges, ncol = 2, byrow = TRUE)
}
