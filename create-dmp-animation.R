setwd("/Volumes/Face3D/Face4D/NewFaces")

library(Face3D)

resizejpg.face3d()

flo <- list.files(recursive = TRUE, pattern = ".obj")
fld <- sub(".obj", ".dmp", flo)
for (i in 1:length(flo)) {
  face <- read.face3d(flo[i])
  save(face, file = fld[i])
}

setwd("/Users/libertyvittert/Desktop/Esco_2174/")

face<- read.face3d("150427094228.obj")
filename<- "150427094228.obj"


jpgfile <- paste(substr(filename, 1, nchar(filename) - 4), jpgfile.addition, ".jpg", sep = "")
            image.mat <- readJPEG(jpgfile)
            m <- dim(image.mat)[1]
            n <- dim(image.mat)[2]
            IC1 <- 1 + round((1 - coordst2[, 2]) * (m - 1))
            IC2 <- 1 + round((coordst2[, 1]) * (n - 1))
            r <- image.mat[cbind(IC1, IC2, 1)]
            g <- image.mat[cbind(IC1, IC2, 2)]
            b <- image.mat[cbind(IC1, IC2, 3)]
            ind <- which(is.na(r + g + b))
            r[ind] <- 0.5
            g[ind] <- 0.5
            b[ind] <- 0.5
            clr <- rgb(r, g, b)
            result$colour <- clr
            class(result) <- "face3d"

result$colour <- coordst2


for (fname in fld) {
  load(fname)
  display.face3d(face)
  view3d(phi = 0)
  fnm <- sub(".dmp", "", fname)
  # id <- strsplit(fnm, "-")[[1]][3]
  sq <- seq(0, 360, length = 100)
  for (i in 1:length(sq)) {
    view3d(theta = sq[i], phi = 0)
    snapshot3d(paste("temp/temp-", i, ".png", sep = ""))
  }
  system(paste("ffmpeg -i temp/", 
               "%d.png  -r 20 -b:v 1024k ", "gun", ".mpg", sep = ""))
  rgl.close()
}
setwd("/Users/liberty/Documents/consultancy/talks/esri/guns/")
