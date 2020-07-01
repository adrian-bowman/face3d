#     Test code to orient faces

library(Face3D)

display.face3d(face)
face <- orient.face3d(face)
display.face3d(face, new = FALSE)

files <- list.files(
            "/Volumes/Face3D/Glasgow-controls/Kathryn/kathryn-dmp", 
            full.names = TRUE)
ind   <- grep("kathryn0", files)
files <- files[ind]
load(files[1])
display.face3d(face)

i     <- 1
for (fnm in files) {
   load(fnm)
   angles       <- runif(3, -pi/8, pi/8)
   face$coords  <- rotate3d(face$coords,  angles[3], 0, 0, 1)
   face$coords  <- rotate3d(face$coords,  angles[2], 0, 1, 0)
   face$coords  <- rotate3d(face$coords,  angles[1], 1, 0, 0)
   face$normals <- rotate3d(face$normals, angles[3], 0, 0, 1)
   face$normals <- rotate3d(face$normals, angles[2], 0, 1, 0)
   face$normals <- rotate3d(face$normals, angles[1], 1, 0, 0)
   display.face3d(face, new = FALSE)
   snapshot3d(paste("temp/old", i, ".png", sep = ""))
   face <- orient.face3d(face)
   display.face3d(face, new = FALSE)
   spheres3d(face$coords[face$nearest, ], col = "green", radius = 2)
   snapshot3d(paste("temp/new", i, ".png", sep = ""))
   print(i)
   print(angles)
   print(face$rotations)
   i <- i + 1
}

