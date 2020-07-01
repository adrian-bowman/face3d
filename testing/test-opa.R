#     Test code to orient faces

install.packages("~/OneDrive - University of Glasgow/research/face3d_0.1-1/face3d",
                 repos = NULL, type = "source")
library(face3d)
library(rgl)
library(shapes)

fls <- list.files("~/Desktop/Glasgow-controls", recursive = TRUE, full.names = TRUE)
fls <- fls[-grep("details", fls)]
load(fls[81])
face1 <- face
load(fls[82])
lmks <- scale(face$lmks, scale = FALSE)

opa  <- opa.face3d(face1$lmks, lmks)
opa1 <- procOPA(lmks, face1$lmks)$Bhat
max(abs(opa - opa1))

opa  <- opa.face3d(face1$lmks, lmks, scale = FALSE)
opa1 <- procOPA(lmks, face1$lmks, scale = FALSE)$Bhat
max(abs(opa - opa1))
