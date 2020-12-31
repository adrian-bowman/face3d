# Test code for the procrustes.face3d function

library(Face3D)
library(rgl)
library(shapes)

load("~/Desktop/lmks.dmp")

for (i in 1:dim(lmks)[3]) spheres3d(lmks[ , , i], radius = 1, col = i)

gpa <- procrustes.face3d(lmks)
gpa <- procrustes.face3d(lmks, scale = FALSE)
gpa <- procrustes.face3d(lmks, subset = 1:3)
gpa <- procGPA(lmks)
gpa <- procGPA(lmks, scale = FALSE)

centroid.size(lmks)
centroid.size(gpa$rotated)

clear3d()
for (i in 1:dim(lmks)[3]) spheres3d(gpa$rotated[ , , i], radius = 1, col = 1:dim(lmks)[1])
spheres3d(gpa$mean, radius = 0.01)

gpa <- procrustes.face3d(lmks)
gpa <- procrustes.face3d(lmks, subset = 1:3)

plot(cbind(c(gpa$rotated[ , 1, ]), c(gpa$rotated[ , 1, ])), type = "n")
for (i in 1:dim(lmks)[3]) points(gpa$rotated[ , , i])