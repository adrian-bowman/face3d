#     Ridgepath algorithm

library(face3d)
if (reinstall) devtools::install("face3d")

library(rgl)

test_label("ridgepath", test.prompt)
plot(template_male)
spheres3d(template_male$landmarks[c("se", "pn"), ], col = "yellow")
ppath <- planepath(template_male, template_male$landmarks["se", ], template_male$landmarks["pn", ],
                   si.target = 1, monitor = 2)
spheres3d(nose_ridge, col = "blue", radius = 0.9)
nose_ridge <- ridgepath(template_male, template_male$landmarks["se", ], template_male$landmarks["pn", ],
                        ppath, penalty = 1, si.target = 1, monitor = 1)
spheres3d(nose_ridge, col = "green", radius = 0.95)
cat("\n*** Why is this not straight when the template is symmetric? ***\n\n")
