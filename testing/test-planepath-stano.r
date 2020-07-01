################################################################################
"distance.3d.diff" <- function(Xdiff) sqrt(sum((Xdiff)^2))
################################################################################
"closest.point" <- function(MESH,POINT){
k <- dim(MESH)[1]
POINTS <- rep(1,k)%*%t(as.vector(POINT))
DIST <- apply(MESH - POINTS,1,distance.3d.diff)
index <- which.min(DIST)
RESULTS <- list(POINT = MESH[index,],index = index)
}


setwd("/Volumes/Face3D/Glasgow-controls/Kathryn/kathryn-dmp")
files <- list.files()

################################################################################
################################################################################
## Notes:
## FAce3D server: Glasgow-controls\Kathryn\kathryn-dmp
## from "kathryn001.dmp" to "kathryn010.dmp"   []
## landmarks 5 (sellion), 10 (endocathion R),9 (endocathion L)
## landmarks 4 (subnasale), 3 (alare crest R), 2 (alare crest L)
## pseudo-geodesics: 5-10; 5-9; 4-3; 4-2
################################################################################

################################################################################
### for(i in 1:10){
 ## i = 1
 setwd("I:/Face.project/Face3dDATA/SK01/glasgow-controls-dmp-complete/") 

i <- 4
   dmpfile <-files[i] 
   load(dmpfile)
## calculate closest point to the lndm on the triangulation
CLOSE.1 <- closest.point(face$coords,face$lndms[5,])$POINT
CLOSE.2 <- closest.point(face$coords,face$lndms[10,])$POINT
CLOSE.3 <- closest.point(face$coords,face$lndms[9,])$POINT
CLOSE.5 <- closest.point(face$coords,face$lndms[4,])$POINT
CLOSE.6 <- closest.point(face$coords,face$lndms[3,])$POINT
CLOSE.7 <- closest.point(face$coords,face$lndms[2,])$POINT
################################################################################
display.face3d(face, type = "mesh")
display.face3d(face, type = "mesh",
       subset = abs(face$coords[ , 2] - 0) < 15, new = FALSE)

shape <- subset.face3d(face, abs(face$coords[ , 2] - 0) < 15)
display.face3d(shape, type = "mesh")
spheres3d(rbind(CLOSE.1, CLOSE.2), radius = 0.5, col = "blue")

min(apply(sweep(shape$coords, 2, CLOSE.1), 1, function(x) sum(abs(x))))


pop3d()
source("/Volumes/adrian/research/shape/Face3D/Face3D-package/Face3D_0.1-1/Face3D/R/planepath.r")
SPATH1 <- planepath.face3d(face, x1 = CLOSE.1, x2 = CLOSE.2)
spheres3d(SPATH1$path, col = "green", radius = 0.2)
#Error: 2,4,7
#Warning message: 1,5,9
SPATH2 <- planepath.face3d(face, x1 = CLOSE.1, x2 = CLOSE.3)
#Error:
#Warning message:  2,4,5,6,8,9
SPATH3 <- planepath.face3d(face, x1 = CLOSE.5, x2 = CLOSE.6)
#Error: 5
#Warning message: 3,4
SPATH4 <- planepath.face3d(face, x1 = CLOSE.5, x2 = CLOSE.7)
#Error:
#Warning message: 1,8,9,10
################################################################################
#Error in if (idx1 <= nind1) idx1 + nind1 else idx1 - nind1 : 
#argument is of length zero
################################################################################
#Warning message:
#In ind1 < ind2 :
#  longer object length is not a multiple of shorter object length
################################################################################
open3d()
rgl.viewpoint(theta = 0,phi = 0,fov=30,zoom=0.7)
bg3d("white")
par3d(windowRect=c(100,100,600,600))
triangles3d(face$coords[face$triples, ], col = "green")
spheres3d(face$lndms[c(2:5,9:10),],radius = 1.5,col="red",add=TRUE)
points3d(SPATH1$path, radius = 1, col = "red",add=TRUE)
points3d(SPATH2$path, radius = 1, col = "red",add=TRUE)
points3d(SPATH3$path, radius = 1, col = "red",add=TRUE)
points3d(SPATH4$path, radius = 1, col = "red",add=TRUE)
spheres3d(rbind(CLOSE.1,CLOSE.2,CLOSE.3,CLOSE.5,CLOSE.6,CLOSE.7), radius = 2, col = "blue",add=TRUE)
### }
################################################################################

  