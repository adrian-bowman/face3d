

x <- seq(0, 4*pi, len = 100)
y <- sin(x) + rnorm(length(x))
plot(x, y)
lines(x, sin(x), lty=2, col="blue")
solution <- sm.map(x, y)
lines(solution$eval.points, solution$estimate, col = "green")
