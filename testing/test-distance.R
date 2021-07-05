#     Test code for the edges.face3d function

install.packages("~/OneDrive - University of Glasgow/research/face3d_0.1-1/face3d",
                 repos = NULL, type = "source")
library(face3d)
library(fields)
library(rgl)
library(geometry)

distance.face3d(1:3, 2:4)
distance.face3d(1:3, matrix(1:12, ncol = 3))
distance.face3d(matrix(1:12, ncol = 3), matrix(2:13, ncol = 3))
boxplot(distance.face3d(1:3, template_male))
temp <- template_male
temp$vertices[ , 1] <- -template_male$vertices[ , 1]
temp$normals[ , 1]  <- -template_male$normals[ , 1]
dst <- distance.face3d(temp, template_male)
plot(template_male, col = dst$x)
plot(template_male, col = dst$y)
plot(template_male, col = dst$z)
plot(template_male, col = dst$xyz)
plot(template_male, col = dst$normal)

# Calls which should produce errors
distance.face3d(1:3, 2:3)
distance.face3d(1:3, matrix(1:12, ncol = 4))

