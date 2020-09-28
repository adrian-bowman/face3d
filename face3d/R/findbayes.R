findbayes <- function(face, mn.popn, covinv) {
   
   # Prior parameters
   mu <- 0
   nu <- 1000
   V  <- nu * diag(2)
   a  <- 0.001
   b  <- 0.001
   
   lmks  <- face$landmarks[rownames(mn.popn), ]

   ivec <- match("se", rownames(lmks))

   for (i in ivec) {
      
      nm  <- rownames(lmks)[i]
      ind <- face[[paste(nm, ".ind", sep = "")]]
      y   <- face[[paste(nm, ".crv", sep = "")]]
      n   <- length(ind)
      
      lpost <- rep(0, length(ind))
      prnt  <- matrix(nrow = 0, ncol = 4)
      
      for (j in 1:length(ind)) {
         
         lmks[i, ] <- face$vertices[ind[j], ]
         face      <- opa.face3d(lmks, mn.popn, face, scale = FALSE)
         
         logpost <- 0
         
         dst <- c(rdist(face$vertices[ind, ], t(lmks[i, ])))
         B   <- cbind(1, dst^2)
      
         V1  <- solve(crossprod(B) + diag(2) / nu)
         mu1 <- V1 %*% (t(B) %*% y)
         a1  <- a + n / 2
         b1  <- b + (sum(y^2) - c(t(mu1) %*% (crossprod(B) + diag(2) / nu) %*% mu1)) / 2
                  
         logpost <- logpost + lgamma(a1) + 0.5 * determinant(V1)$modulus - a1 * log(b1)
   
         vec     <- c(lmks) - c(mn.popn)
         logpost <- logpost - c(0.5 * t(vec) %*% covinv %*% vec)
         attr(logpost, "logarithm") <- NULL
         lpost[j] <- logpost
         
         # print(logpost)
         dt <- 0.5 * determinant(V1)$modulus
         attr(dt, "logarithm") <- NULL
         prnt <- rbind(prnt, c(lgamma(a1), -a1 * log(b1), dt,
                       -c(0.5 * t(vec) %*% covinv %*% vec)))
      }
      print(prnt)
      
   }
   
   lpost
}