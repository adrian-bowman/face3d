connected.face3d <- function(shape) {
   
   summary(shape, print = FALSE)
   
   triples       <- matrix(shape$triples, ncol = 3, byrow = TRUE)
   current       <- triples[1, ]
   prt           <- 1
   part          <- rep(0, nrow(shape$coords))
   part[current] <- prt
   
   while (length(current) > 0) {
      ind           <- (triples[ , 1] %in% current) | (triples[ , 2] %in% current) |
                       (triples[ , 3] %in% current)
      current       <- unique(c(triples[ind, ]))
      ind           <- which(part[current] == 0)
      current       <- current[ind]
      if (length(current) > 0)
         part[current] <- prt
      else if (any(part == 0)) {
      	 flag <- TRUE
         while (flag) {
            ind0 <- which(part == 0)[1]
            ind1 <- (triples[ , 1] == ind0) | (triples[ , 2] == ind0) |
                    (triples[ , 3] == ind0)
            ind1 <- which(ind1)
            if (length(ind1) == 0) part[ind0] <- -1
            flag <- (length(ind1) == 0) & any(part == 0)
         }
         if (length(ind1) > 0) {
            current <- c(triples[ind1[1], ])
            prt     <- prt + 1
         }
      }
   }
   
   part[part == -1] <- 0
   part0            <- part[part > 0]
   ord              <- order(table(part0), decreasing = TRUE)
   code             <- 1:length(ord)
   names(code)      <- ord
   part[part > 0]   <- code[as.character(part0)]
   invisible(part)
}

# Code for triangle areas.
# n                    <- length(source.shape$triples)/3
# source.triples.array <- array(t(source.shape$coords[source.shape$triples,]),c(3,3,n))
# source.triples.areas <- apply(source.triples.array,3, function(x){
#    1/2 * sqrt(sum((crossproduct(x[,2] - x[,1],x[,3] - x[,1],scale = FALSE))^2))
# })
# 
