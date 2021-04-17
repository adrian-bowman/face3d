connected.face3d <- function(shape) {
   
   if (!is.face3d(shape)) stop("this is not a face3d object.")
   
   triples       <- shape$triangles
   current       <- triples[1, ]
   prt           <- 1
   part          <- rep(0, nrow(shape$vertices))
   part[current] <- prt
   
   while (length(current) > 0) {
      indt       <- (triples[ , 1] %in% current) | (triples[ , 2] %in% current) |
                    (triples[ , 3] %in% current)
      if ("lines" %in% names(shape)) {
         indl    <- (shape$lines[ , 1] %in% current) | (shape$lines[ , 2] %in% current)
         current <- unique(c(c(triples[indt, ]), c(shape$lines[indl, ])))
      }
      else
         current    <- unique(c(triples[indt, ]))
      ind           <- which(part[current] == 0)
      current       <- current[ind]
      if (length(current) > 0)
         part[current] <- prt
      else if (any(part == 0)) {
      	flag <- TRUE
         while (flag) {
            ind0  <- which(part == 0)[1]
            ind1t <- (triples[ , 1] == ind0) | (triples[ , 2] == ind0) |
                     (triples[ , 3] == ind0)
            ind1t <- which(ind1t)
            ind1l <- numeric(0)
            if ("lines" %in% names(shape)) {
               ind1l <- (shape$lines[ , 1] == ind0) | (shape$lines[ , 2] == ind0)
               ind1l <- which(ind1l)
            }
            if (length(ind1t) + length(ind1l) == 0) part[ind0] <- -1
            flag <- (length(ind1t) + length(ind1l) == 0) & any(part == 0)
         }
         if (length(ind1t) + length(ind1l) > 0) {
            if (length(ind1t) > 0)
               current <- c(triples[ind1t[1], ])
            else
               current <- c(shape$lines[ind1l[1], ])
            prt <- prt + 1
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
