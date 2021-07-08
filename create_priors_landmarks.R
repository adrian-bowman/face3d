#     Create prior information for landmarks

library(Face3D)
library(fields)
library(rgl)
library(shapes)
library(geometry)
library(colorspace)

source("Face3D/R/a-midline.R")
source("Face3D/R/a-initiallandmarks.R")

liberty.data <- "~/Desktop/Data-Liberty/"
# liberty.data <- "/Volumes/Face3D/Glasgow-controls/Liberty/Liberty-Controls-dmp/Liberty-Controls-dmp-final/Data-Liberty/"
fls <- list.files(liberty.data, full.names = TRUE)

load("~/Desktop/lmks-liberty.rda")
load(fls[1])
dimnames(lmks.liberty)[[1]] <- rownames(face$lmks)

lmk.nms <- c("pn", "se", "enL", "enR")
lmk.nms <- "se"
priors  <- array(dim = c(length(fls), 3, length(lmk.nms)),
                 dimnames = list(NULL, NULL, lmk.nms))
nhd.dst <- 10
for (i in 1:length(fls)) {
   load(fls[i])
   for (lmk.nm in lmk.nms) {
      dst   <- edist.face3d(face$coords, face$lmks[lmk.nm, ])
      ind   <- which(dst < nhd.dst)
      dst   <- dst[ind]
      sbst  <- subset(face, ind)
      if (lmk.nm == "pn") {
         crv <- face$shape.index[ind]
      }
      else if (lmk.nm == "se") {
         pp  <- planepath.face3d(sbst, face$lmks[lmk.nm, ],
                                 direction = face$lmks[lmk.nm, ] - face$lmks["pn", ],
                                 bothways = TRUE, rotation = 0)
         gc  <- gcurvature.face3d(pp$path, 3)
         crv <- gc$gcurvature
         dst <- gc$arclength - pp$x1.arclength
         se  <- gc$resampled.curve[which.max(gc$gcurvature), ]
      }
      else
         face$kappa1[ind] * face$kappa2[ind]
      model <- lm(crv ~ I(dst^2))
      plot(crv ~ dst)
      points(dst, fitted(model), pch = 16, col = "red")
      plot(sbst, type = "mesh", new = FALSE)
      spheres3d(face$lmks[lmk.nm, ], radius = 0.4, col = "yellow")
      spheres3d(pp$path, radius = 0.2)
      spheres3d(se, radius = 0.4, col = "red")
      
      # plot(sbst, col = "shape index")
      # plot(sbst, col = sbst$kappa1 * sbst$kappa2)
      # spheres3d(face$lmks[lmk.nm, ])
      # plot(sbst$shape.index ~ dst^2)
      scan()
      # priors[i, , lmk.nm] <- c(model$coefficients,
      #                          sqrt(sum(model$residuals^2) / model$df.residual))
   }
   cat(i, "")
}

# apply(priors, 2:3, mean)
# covs <- apply(priors, 3, function(x) cov(x[ , 1:2]))
# covs <- array(c(covs), dim = c(2, 2, dim(priors)[3]),
#               dimnames = list(NULL, NULL, dimnames(priors)[[3]]))
# 
# plot(priors[ , 1:2, "pn"])
# plot(priors[ , 1:2, "se"])
# plot(priors[ , 1:2, "enL"])
# plot(priors[ , 1:2, "enR"])
# hist(priors[ , 3, "pn"])
# hist(priors[ , 3, "se"])
# hist(priors[ , 3, "enL"])
# hist(priors[ , 3, "enR"])

