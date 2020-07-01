a <- 0.0001
b <- 0.0001
a <- 1
b <- 1
curve(x^(-a-1) * exp(-b/x), 0.01, 10)

xgrid <- seq(0.01, 10, length = 100)
dgrid <- xgrid^(-a-1) * exp(-b/xgrid)
plot(xgrid, cumsum(dgrid) / sum(dgrid))
