

library(fields)
library(MASS)
library(rgl)
"edist.face3d"  <- function(x, x0) sqrt(rowSums(sweep(x, 2, x0)^2))

setwd("/Users/libertyvittert/ownCloud/Shared/Face3D_0.1-1/Face3D/R")
for (file in list.files())source(file)

source("/Users/libertyvittert/ownCloud/Shared/Face3D_0.1-1/Face3D/R/l-smoothpath-ridge.r")



for (i in 1:length(dmp)){
	load(dmp[i])
display.face3d(face)
spheres3d(face$curves)
scan()
}




setwd("/Users/libertyvittert/ownCloud/Shared/Face3D_0.1-1/testing")
dmp <- list.files()
dmp <- dmp[-46]
for (zz in 2:length(dmp)){
load(dmp[zz])
	face <- facecurves.face3d(face, "nasal bridge left", graphics = FALSE)
	face <- facecurves.face3d(face, "nasal bridge right", graphics = FALSE)
	face <- facecurves.face3d(face, "mid-line lip", graphics = FALSE)
	face <- facecurves.face3d(face, "lower lip right", graphics = FALSE)
	face <- facecurves.face3d(face, "lower lip left", graphics = FALSE)
	face <- facecurves.face3d(face, "upper lip right", graphics = FALSE)
	face <- facecurves.face3d(face, "upper lip left", graphics = FALSE)
	

     # face <- facecurves.face3d(face, "mandible right", graphics = FALSE)
     # face <- facecurves.face3d(face, "mandible left", graphics = FALSE)
     # face <- facecurves.face3d(face, "Mandible", graphics = FALSE)

     # face <- facecurves.face3d(face, "brow ridge left", graphics = FALSE)
     # face <- facecurves.face3d(face, "brow ridge right", graphics = FALSE)

save(face, file=paste(dmp[zz]))
print(zz)

}

