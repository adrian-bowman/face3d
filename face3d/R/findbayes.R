findbayes <- function(face, mn.popn, covinv) {
   
   # Prior parameters
   mu <- 0
   nu <- 1000
   V  <- nu * diag(2)
   a  <- 0.001
   b  <- 0.001
   
   lmks  <- face$landmarks[rownames(mn.popn), ]
   lpost <- rep(0, length(face$se.ind))
   
   for (j in 1:length(face$se.ind)) {
      lmks["se", ] <- face$vertices[face$se.ind[j], ]
   
   face <- opa.face3d(lmks, mn.popn, face, scale = FALSE)
   
   logpost <- 0
   
   for (i in 1:nrow(lmks)) {
      
      ind <- face[[paste(rownames(lmks)[i], ".ind", sep = "")]]
      n   <- length(ind)
      dst <- c(rdist(face$vertices[ind, ], t(lmks[i, ])))
      B   <- cbind(1, dst^2)
      y   <- face[[paste(rownames(lmks)[i], ".crv", sep = "")]]
      
      V1  <- solve(crossprod(B) + diag(2) / nu)
      mu1 <- V1 %*% (t(B) %*% y)
      a1  <- a + n / 2
      b1  <- b + (sum(y^2) - c(t(mu1) %*% (crossprod(B) + diag(2) / nu) %*% mu1)) / 2
                  
      logpost <- logpost + lgamma(a1) + 0.5 * determinant(V1)$modulus - a1 * log(b1)
      
   }
   
   vec     <- c(lmks) - c(mn.popn)
   logpost <- logpost - c(0.5 * t(vec) %*% covinv %*% vec)
   attr(logpost, "logarithm") <- NULL
   print(logpost)
   lpost[j] <- logpost
   
   }
   
   lpost
}