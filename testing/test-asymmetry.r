# Tests for the asymmetry.face3d function

setwd("/Volumes/adrian/research/shape/Face3D/Face3D-package/Face3D_0.1-1/Face3D/R")
# setwd("~/Desktop/Face3D-package/Face3D_0.1-1/Face3D/R")

library(Face3D)

load("/Volumes/Face3D/Glasgow-controls/glasgow-controls-liberty/lmks.dmp")
coding <- list(paired.lmks = matrix(c("acL", "acR", "nbL", "nbR",
                    "exL", "exR", "enL", "enR", "tL", "tR", 
                    "cphL", "cphR", "chL", "chR", "aL", "aR"),
                    ncol = 2, byrow = TRUE),
               single.lmks = c("pm", "sn", "se", "n", "ls", "sto",
                    "li", "sl", "pg")
                    )

source("asymmetry.r")
score <- asymmetry.face3d(lmks, coding)
