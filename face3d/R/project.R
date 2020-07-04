project.face3d <- function(x, normals, shape) {
   trpls <- shape$triangles
   proj  <- 0 * x
   for (i in 1:nrow(x)) {
      ind <- which(edist.face3d(shape$vertices, x[i, ]) <= 10)
      # ind <- which(order(edist.face3d(shape$vertices, x[i, ])) <= 100)
      ind <- which(apply(trpls, 1, function(x) all(x %in% ind)))
      po  <- sweep(shape$vertices[trpls[ind, 1], ], 2, x[i, ])
      n   <- crossproduct(shape$vertices[trpls[ind, 2], ] - shape$vertices[trpls[ind, 1], ],
                          shape$vertices[trpls[ind, 3], ] - shape$vertices[trpls[ind, 1], ])
      tt  <- apply(po * n, 1, sum) / c(n %*% normals[i, ])
      xp  <- rep(1, length(tt)) %o% x[i, ] + tt %o% normals[i, ]
      # spheres3d(xp)
      # spheres3d(shape$vertices[c(trpls[ind, ]), ])
      # pop3d()
      l1 <- crossproduct(shape$vertices[trpls[ind, 2], ] - shape$vertices[trpls[ind, 1], ], xp - shape$vertices[trpls[ind, 1], ])
      l2 <- crossproduct(shape$vertices[trpls[ind, 3], ] - shape$vertices[trpls[ind, 2], ], xp - shape$vertices[trpls[ind, 2], ])
      l3 <- crossproduct(shape$vertices[trpls[ind, 1], ] - shape$vertices[trpls[ind, 3], ], xp - shape$vertices[trpls[ind, 3], ])
      l1 <- apply(l1 * n, 1, sum)
      l2 <- apply(l2 * n, 1, sum)
      l3 <- apply(l3 * n, 1, sum)
      ind <- (l1 >= 0) & (l2 >= 0) & (l3 >= 0)
      proj[i, ] <- xp[which(ind)[1], ]
      # spheres3d(x[i, ], col = "yellow", radius = 1.1)
      # spheres3d(proj[i, ], col = "green", radius = 1.3)
      # print(i)
      # scan()
   }
   invisible(proj)
}
