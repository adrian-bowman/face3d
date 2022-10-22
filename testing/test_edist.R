#     Tests for the edist function

library(face3d)
if (reinstall) devtools::install("face3d")

library(rgl)

test_label("edist: pairwise", test.prompt)
cat("  number of vertices:", nrow(template_male$vertices), "\n")
cat("  number of curve points:", nrow(template_male$curves), "\n")
cat("  matrix of distances:", dim(edist(template_male, template_male$curves)), "\n")

test_label("edist: pairwise with single points", test.prompt)
cat("  ", length(edist(template_male, template_male$curves[1, ])), "\n")
cat("  ", length(edist(template_male$curves[1, ], template_male)), "\n")
cat("  ", edist(template_male$curves[1, ], template_male$curves[10, ]), "\n")

test_label("edist: only one set of points", test.prompt)
cat("  ", length(edist(template_male$curves)), "\n")

test_label("edist: min", test.prompt)
dst <- edist(template_male, template_male$curves, "min")
plot(template_male, col = dst)
spheres3d(template_male$curves, col = "red")
cat("  the colour shows the shortest distance of each vertex from the curve.\n")

test_label("edist: min with a curve", test.prompt)
dst <- edist(template_male, template_male$curves, "min", x2curve = TRUE)
plot(template_male, col = dst)
spheres3d(template_male$curves, col = "yellow")
cat("  the colour shows the signed distance of each vertex from the curve.\n")

test_label("edist: minsum", test.prompt)
cat("  ", edist(template_male$curves, template_male, "minsum"), "\n")
cat("  This distance is not close to 0 because the curves extend beyond the template.\n")

