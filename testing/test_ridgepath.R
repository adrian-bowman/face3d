#     Ridgepath algorithm

library(face3d)
if (reinstall) devtools::install("face3d")

library(rgl)

test_label("ridgepath", test.prompt)
plot(template_male)
spheres3d(template_male$landmarks[c("se", "pn"), ], col = "yellow")
nose_ridge <- ridgepath(template_male, template_male$landmarks["se", ], template_male$landmarks["pn", ],
                        penalty = 0.002, si.target = 1, monitor = 1)
