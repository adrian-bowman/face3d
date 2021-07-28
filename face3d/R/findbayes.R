findbayes <- function(face, mn.popn, covinv, monitor = 0) {
   
   # Prior parameters
   mu <- 0
   nu <- 1000
   V  <- nu * diag(2)
   a  <- 0.001
   b  <- 0.001
   
   lmks            <- face$landmarks[rownames(mn.popn), ]
   face$lmks.bayes <- lmks
   lmks.previous   <- matrix(0, nrow = nrow(lmks), ncol = 3)
   delta           <- mean(apply(lmks - lmks.previous, 1, function(x) sqrt(sum(x^2))))
   iter            <- 0
   
   while(delta > 1e-5) {

   if (iter > 0) lmks.previous <- face$lmks.bayes
   
   for (i in 1:nrow(lmks)) {
      
      nm    <- rownames(lmks)[i]
      ind   <- face[[paste(nm, ".ind", sep = "")]]
      n     <- length(ind)
      lpost <- rep(0, length(ind))
      prnt  <- matrix(nrow = 0, ncol = 4)
      
      for (j in 1:length(ind)) {
         
         face$lmks.bayes[i, ] <- face$vertices[ind[j], ]
         face      <- opa.face3d(face$lmks.bayes, mn.popn, face, scale = FALSE)
         lmks      <- face$lmks.bayes
         logpost   <- 0

         for (k in 1:nrow(lmks)) {
            nm      <- rownames(lmks)[k]
            indk    <- face[[paste(nm, ".ind", sep = "")]]
            y       <- face[[paste(nm, ".crv", sep = "")]]
            dst     <- c(rdist(face$vertices[indk, ], t(lmks[k, ])))
            B       <- cbind(1, dst^2)
            V1      <- solve(crossprod(B) + diag(2) / nu)
            mu1     <- V1 %*% (t(B) %*% y)
            a1      <- a + n / 2
            b1      <- b + (sum(y^2) - c(t(mu1) %*% (crossprod(B) + diag(2) / nu) %*% mu1)) / 2
            logpost <- logpost + lgamma(a1) + 0.5 * determinant(V1)$modulus - a1 * log(b1)
         }
         
         vec     <- c(lmks) - c(mn.popn)
         logpost <- logpost - c(0.5 * t(vec) %*% covinv %*% vec)
         attr(logpost, "logarithm") <- NULL
         lpost[j] <- logpost
         
         # plot(face)
         # spheres3d(mn.popn, col = "blue", radius = 1.5)
         # for (k in 1:nrow(lmks)) {
         #    indk   <- c(k, k + nrow(lmks), k + 2 * nrow(lmks))
         #    print(indk)
         #    plot3d(ellipse3d(covt[indk, indk], centre = mn.popn[k, ]),
         #           col = "lightblue", alpha = 0.5, add = TRUE)
         # }
         # spheres3d(lmks,    col = "red",  radius = 1.5)
         # print(logpost)
         # scan()
         
         dt <- 0.5 * determinant(V1)$modulus
         attr(dt, "logarithm") <- NULL
         prnt <- rbind(prnt, c(lgamma(a1), -a1 * log(b1), dt,
                       -c(0.5 * t(vec) %*% covinv %*% vec)))
      }
      
      face$lmks.bayes[i, ] <- face$vertices[ind[which.max(lpost)], ]
      face <- opa.face3d(face$lmks.bayes, mn.popn, face, scale = FALSE)
      
      if (monitor > 0) {
         sbst <- subset(face, ind, remove.singles = FALSE)
         plot(face)
         plot(sbst, col = lpost,  display = "spheres", palette = topo.colors(20), add = TRUE)
         spheres3d(face$lmks.bayes[i, ], col = "red", radius = 1.5)
         if (monitor > 2) invisible(readline(prompt = "Press [enter] to continue"))
         spherespresent <- TRUE
      }

      # print(prnt)
   }
   
   delta <- mean(apply(face$lmks.bayes - lmks.previous, 1, function(x) sqrt(sum(x^2))))
   iter  <- iter + 1
   if (monitor > 0) cat("iteration", iter, ": delta =", delta, "\n")
   }

   face
}
