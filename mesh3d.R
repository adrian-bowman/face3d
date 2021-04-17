vertices <- matrix(c(0, 0, 0, 1, 0, 0, 1, 1, 0), ncol = 3, byrow = TRUE)
indices  <- matrix(c(1, 2, 3), nrow = 1)
object <- tmesh3d(vertices, indices, homogeneous = FALSE)

temp <- cube3d()
