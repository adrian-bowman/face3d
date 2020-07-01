
 # animate.face3d <- function(shape1, shape2, rot.shape1, rot.shape2, n.animations = 10, n.cycles = 5, new = TRUE) {
   # if (new) display.face3d(shape1)
   # shape.temp  <- shape1
   # shape.temp2 <- rot.shape1
   # id <- rgl.ids()[ , 1]
   # i  <- 1
   # sgn <- 1
   # for (j in 1:(n.cycles * 2 * n.animations)) {
      # i <- i + 1 * sgn
      # if (i == n.animations | i == 0) sgn <- -sgn
      # shape.temp$coords <- (shape1$coords * (n.animations - i) + shape2$coords * i) / n.animations
      # shape.temp2$coords <- (rot.shape1$coords * (n.animations - i) + rot.shape2$coords * i) / n.animations
      # # normals           <- normals.face3d(shape.temp)$normal
      # # difference        <- shape.temp$coords - shape2$coords
      # # euc.dist          <- sqrt(apply(difference^2, 1, sum))
      # # ind               <- apply(difference * normals, 1, sum)
      # # pts               <- euc.dist * sign(ind)
      # # nbrk              <- 50
      # # cc                <- topo.colors(nbrk)
      # # colours           <- cc[cut(pts, nbrk)]
      # # shape.temp$colour <- colours
      # display.face3d(shape.temp, add = TRUE)
      # display.face3d(shape.temp2, add=TRUE)
      # lines3d(shape.temp$coords[grep("upper face right 2", rownames(shape.temp$coords)), ])
      # lines3d(shape.temp$coords[grep("upper face left 2", rownames(shape.temp$coords)), ])
      # lines3d(shape.temp$coords[grep("mid-line columella", rownames(shape.temp$coords)), ])
      # lines3d(shape.temp2$coords[grep("upper face right 2", rownames(shape.temp2$coords)), ])
      # lines3d(shape.temp2$coords[grep("upper face left 2", rownames(shape.temp2$coords)), ])
      # lines3d(shape.temp2$coords[grep("mid-line columella", rownames(shape.temp2$coords)), ])
      # rgl.snapshot(paste(id, ".png",sep=""), "png")
      # pop3d(id = id)
      # id <- rgl.ids()[ , 1]
   # }
# }



# animate.face3d(rendering$rend.lmks, rendering$rend.rro)



animate.face3d <- function(shape1, shape2, n.animations = 10, n.cycles = 5, new = TRUE) {
   if (new) display.face3d(shape1)
   shape.temp <- shape1
   id <- rgl.ids()[ , 1]
   i  <- 1
   sgn <- 1
   for (j in 1:(n.cycles * 2 * n.animations)) {
      i <- i + 1 * sgn
      if (i == n.animations | i == 0) sgn <- -sgn
      shape.temp$coords <- (shape1$coords * (n.animations - i) + shape2$coords * i) / n.animations
      print(i)
      #also need to calculate new colour ie. how far from symmetric is it?
      normals           <- normals.face3d(shape.temp)$normal
      difference        <- shape.temp$coords - shape2$coords
      euc.dist          <- sqrt(apply(difference^2, 1, sum))
      ind               <- apply(difference * normals, 1, sum)
      pts               <- euc.dist * sign(ind)
      nbrk              <- c(seq(-3,-1,1), seq(1,3,1))
      # nbrk              <- c(-40, seq(-10,-1,1), seq(1,10,1),40)
      cc                <- topo.colors(length(nbrk))
      colours           <- cc[cut(pts, nbrk)]
      
      if (any(is.na(colours))==TRUE){
          corners       <- which(is.na(colours)==TRUE)
        # print(length(corners))
          colours[corners]      <- colours[corners + 1 ]
       }

      shape.temp$colour <- colours
      display.face3d(shape.temp, add = TRUE)
      
      lines3d(shape.temp$coords[grep("upper face right 2", rownames(shape.temp$coords)), ])
      lines3d(shape.temp$coords[grep("upper face left 2", rownames(shape.temp$coords)), ])
      lines3d(shape.temp$coords[grep("mid-line columella", rownames(shape.temp$coords)), ])
      rgl.snapshot(paste(id, ".png",sep=""), "png")
      pop3d(id = id)
      id <- rgl.ids()[ , 1]
   }
}
