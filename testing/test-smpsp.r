#     Test code for subset.face3d

setwd("/Volumes/adrian/research/shape/Face3D/Face3D-package/Face3D_0.1-1/Face3D/R")
# setwd("~/Desktop/Face3D-package/Face3D_0.1-1/Face3D/R")

n  <- 50
x1 <- runif(n)
x2 <- cbind(x1, runif(n))
x3 <- cbind(x2, runif(n))
y  <- rnorm(n)

source("sm-psp.r")
model <- sm.psp(x1, y)
model <- sm.psp(x1, y, period = 1)
model <- sm.psp(x2, y)
model <- sm.psp(x2, y, period = c(1, NA))
model <- sm.psp(x3, y)
model <- sm.psp(x3, y, period = c(1, NA, NA))

# This doesn't handle periodic cases properly
model <- sm.psp(x1, y, period = 1)
model$beta
