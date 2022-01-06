#     Ridgepath algorithm

library(face3d)
if (reinstall) devtools::install("face3d")

library(rgl)

test_label("ridgepath", test.prompt)
nose <- subset(template_male, edist(x1, template_male) < 30 | edist(x2, template_male) < 30)
nose <- curvatures(nose)
plot(nose, col = "shape index")
plot(nose, col = nose$kappa1)
plot(nose, col = nose$kappa2)
x1   <- template_male$landmarks["se", ]
x2   <- template_male$landmarks["pn", ]
ppath <- planepath(nose, x1, x2, si.target = 1, rotation = 0, monitor = 2)
ppath <- planepath(nose, x1, x2, si.target = 1, monitor = 2)
plot(template_male)
spheres3d(rbind(x1, x2), col = "yellow")
spheres3d(ppath$path, col = "blue", radius = 0.9)
nose_ridge <- ridgepath(nose, x1, x2, ppath, penalty = 0.002, si.target = 1, rotation = 0, monitor = 1)
spheres3d(nose_ridge, col = "green", radius = 0.95)
cat("\n*** Why is this not straight when the template is symmetric? ***\n\n")
