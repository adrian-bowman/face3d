# Fit a model to an image by pca

setwd("~/ownCloud/Face3D_0.1-1")

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
fls          <- list.files(liberty.data, full.names = TRUE)

i <- 32
load(fls[i])
plot(face, new = FALSE)
face <- subset(face, c(rdist(face$coords, t(face$lmks["pn", ]))) < 40)
plot(face, new = FALSE, col = "shape index")
face <- index.face3d(face, distance = 10, overwrite = TRUE)
plot(face, new = FALSE, col = "shape index", key = TRUE)
plot(face, new = FALSE, col = face$kappa1, key = TRUE)
plot(face, new = FALSE, col = face$kappa2)

# Create pc's
load("lmks-liberty.rda")
load(fls[1])
dimnames(lmks.liberty)[[1]] <- rownames(face$lmks)
lmk.names <- c("pn", "enL", "enR", "se", "acL", "acR", "sn", "exL", "exR",
               "chL", "chR", "sl", "gn")
gpa   <- procGPA(lmks.liberty[lmk.names, , ])
mnt   <- apply(gpa$tan, 1, mean)
covt  <- cov(t(gpa$tan))
indy  <- 1:4
mny   <- gpa$mshape[indy, ]
mnmny <- apply(mny, 2, mean)
indy  <- rep(1:4, 3) + rep((0:2) * length(lmk.names), each = length(indy))
# ind   <- c(rbind(ind * 3 - 2, ind * 3 - 1, ind * 3))
Sxx   <- covt[-indy, -indy]
Sxy   <- covt[-indy,  indy]
Syx   <- covt[ indy, -indy]
Syy1  <- solve(covt[ indy,  indy])
Sc    <- Sxx - Sxy %*% Syy1 %*% Syx

clear3d()
for (i in 1:gpa$n) spheres3d(lmks.liberty[lmk.names, , i], col = i)
clear3d()
for (i in 1:gpa$n) spheres3d(gpa$rotated[ , , i], col = i)
for (i in 1:13) {
   ind <- c(i, i + 13, i + 26)
   plot3d(ellipse3d(covt[ind, ind], centre = gpa$mshape[i, ]), col = "yellow",
          alpha = 0.7, add = TRUE)
}
spheres3d(mny, col = "blue", radius = 4)

# Find initial landmarks and store these for later use
# for (i in 1:length(fls)) {
#    load(fls[i])
#    lmks <- face$lmks
#    face$lmks <- NULL
#    face <- initiallandmarks.face3d(face, orient = TRUE, monitor = FALSE, graphics = FALSE, overwrite = TRUE)
#    face$initial.lmks <- face$lmks
#    face$lmks <- lmks
#    save(face, file = fls[i])
#    cat(i, "")
# }
# cat("\n")

# Review the initial landmarks
for (i in 1:length(fls)) {
   load(fls[i])
   plot(face, new = FALSE)
   spheres3d(face$initial.lmks, col = "yellow", radius = 2)
   cat(i)
   scan()
}

# Find the columella
# for (i in 1:length(fls)) {
   i <- 11
   load(fls[i])
   face$lmks <- face$initial.lmks
   face <- subset(face, c(rdist(face$coords, t(face$lmks["pn", ]))) < 40)
   plot(face, new = FALSE)
   plot(face, col = face$kappa1, new = FALSE)
   # for (angle in seq(-pi/8, pi/8, length = 12)) {
   pp <- planepath.face3d(face, face$lmks["se", ], face$lmks["pn", ])
   spheres3d(pp$path)
      # pp <- planepath.face3d(face, face$lmks["pn", ], direction = face$lmks["pn", ] - face$lmks["se", ],
      #                        rotation = angle)
      # spheres3d(face$lmks[c("pn", "se"), ], col = "green", radius = 1.1)
      # segments3d(rbind(face$lmks["pn", ], face$lmks["pn", ] + 5 * pp$normal))
      # spheres3d(pp$path)
#    cat(i, "")
#    scan()
# }


# Find points on the nasal base and the mouth valley
for (i in 1:length(fls)) {
load(fls[i])
face$lmks <- face$initial.lmks
face <- subset(face, c(rdist(face$coords, t(face$lmks["pn", ]))) < 100)
plot(face, new = FALSE)
spheres3d(face$lmks[c("pn", "se"), ], col = "green", radius = 1.1)
for (angle in seq(-pi/8, pi/8, length = 12)) {
   pp <- planepath.face3d(face, face$lmks["pn", ], direction = face$lmks["pn", ] - face$lmks["se", ],
                          rotation = angle)
   # spheres3d(pp$path)
   gcrv <- gcurvature.face3d(pp$path, 12)
   # plot(gcrv$gcurvature ~ gcrv$arc.length)
   # ind  <- which.max(gcrv$gcurvature)
   turn <- c(0, diff(sign(diff(gcrv$d2.z))), 0)
   inds <- which(turn < 0)[1]
   indm <- which(turn < 0)[2]
   subn <- gcrv$resampled.curve[inds, ]
   mth  <- gcrv$resampled.curve[indm, ]
   spheres3d(subn, col = "green",  radius = 2)
   spheres3d(mth,  col = "yellow", radius = 2)
   plot(gcrv$arc.length, gcrv$d2.z)
   points(gcrv$arc.length[inds], gcrv$d2.z[inds], pch = 16, col = "green")
   points(gcrv$arc.length[indm], gcrv$d2.z[indm], pch = 16, col = "yellow")
   scan()
}
cat(i, "")
scan()
}


# Use the prior to find other landmarks
plot(face, new = FALSE)
for (i in 1:length(fls)) {
   load(fls[i])
   face$lmks <- face$initial.lmks
   plot(face, type = "mesh", new = FALSE)
   spheres3d(face$lmks, col = "yellow", radius = 2)
   ind      <- 5:13
   mnx      <- if (length(ind) == 1) matrix(gpa$mshape[ind, ], ncol = 3) else gpa$mshape[ind, ]
   opa      <- procOPA(mny, face$lmks)
   y        <- c(sweep(opa$Bhat, 2, mnmny, "+"))
   indS     <- c(ind - 4, ind - 4 + (length(lmk.names) - 4),
                          ind - 4 + 2 * (length(lmk.names) - 4))
   pred      <- mnx + matrix(Sxy[indS, ] %*% Syy1 %*% (y - c(mny)), ncol = 3, byrow = TRUE)
   # for (j in (ind - 4)) {
   #    Sp <- Sc[c(j, j + 9, j + 18), c(j, j + 9, j + 18)]
   #    plot3d(ellipse3d(Sp, centre = pred[j, ]), col = "yellow", add = TRUE)
   # }
   # spheres3d(pred,  col = "red", radius = 5)
   # spheres3d(mnmny, col = "green", radius = 5)
   # spheres3d(rep(0, 3), col = "lightblue", radius = 5)
   pred     <- sweep(pred, 2, mnmny) %*% t(opa$R) / opa$s
   pred     <- sweep(pred, 2, apply(face$lmks, 2, mean), "+")
   # pop3d()
   # spheres3d(pred, col = "lightgreen", radius = 5)
   for (j in 1:3) {
      indS     <- c(j, j + 9, j + 18)
      Spred    <- opa$R %*% Sc[indS, indS] %*% t(opa$R)
      mhd      <- apply(face$coords, 1, function(x) t(x - pred[j, ]) %*% Spred %*% (x - pred[j, ]))
      spheres3d(face$coords[order(mhd)[1], ], col = "blue")
      # plot3d(ellipse3d(Spred, centre = pred[j, ]), add = TRUE)
   }
   cat(i, "")
   scan()
}

# Find the mouth
for (i in 1:length(fls)) {
   load(fls[i])
   face$lmks <- NULL
   plot(face, new = FALSE)
   face   <- initiallandmarks.face3d(face, orient = FALSE, monitor = TRUE, graphics = "add")
   sbst   <- subset(face, c(rdist(face$coords, t(face$lmks["pn", ]))) < 70)
   plot(sbst, new = FALSE, col = sbst$kappa1)
   lmk.se <- face$lmks["se", ]
   lmk.pn <- face$lmks["pn", ]
   mp     <- (lmk.se + lmk.pn) / 2
   proj   <- sweep(sbst$coords, 2, mp) %*% (lmk.se - mp)
   sbst   <- subset(sbst, proj < 0)
   plot(sbst, new = FALSE, col = sbst$kappa1)

   p      <- 0.8
   parts  <- 1:2
   while (length(table(parts)) < 3) {
      ind    <- sbst$kappa1 > quantile(sbst$kappa1, p, na.rm = TRUE)
      sbst0 <- subset(sbst, ind)
      parts  <- connected.face3d(sbst0)
      p      <- p - 0.05
   }
   # spheres3d(sbst0$coords[parts == 1, ], col = "yellow")
   # spheres3d(sbst0$coords[parts == 2, ], col = "green")
   # spheres3d(sbst0$coords[parts == 3, ], col = "grey")
   sbst1  <- subset(sbst0, parts == 1)
   sbst2  <- subset(sbst0, parts == 2)
   sbst3  <- subset(sbst0, parts == 3)
   dist1  <- mean(rdist(sbst1$coords, t(lmk.pn)))
   dist2  <- mean(rdist(sbst2$coords, t(lmk.pn)))
   dist3  <- mean(rdist(sbst3$coords, t(lmk.pn)))
   ord    <- order(c(dist1, dist2, dist3))
   nose   <- switch(match(1, ord), sbst1, sbst2, sbst3)
   mouth  <- switch(match(2, ord), sbst1, sbst2, sbst3)
   sulcus <- switch(match(3, ord), sbst1, sbst2, sbst3)
   plot(nose,   type = "spheres", col = "yellow", add= TRUE)
   plot(mouth,  type = "spheres", col = "green",  add= TRUE)
   plot(sulcus, type = "spheres", col = "grey",   add= TRUE)
   cat(i, "\n")
   scan()
}
   
plot(sbst, new = FALSE)
# sbst <- index.face3d(sbst, distance = 5, overwrite = TRUE)
plot(sbst, new = FALSE, col = "shape index")
plot(sbst, new = FALSE, col = sbst$kappa1)

p <- 0.7
k1   <- sbst$kappa1
proj <- sweep(sbst$coords, 2, mp) %*% (lmk.se - mp)
k1[proj > 0] <- NA
ind    <- k1 > quantile(k1, p, na.rm = TRUE)
plot(sbst, new = FALSE, col = k1)
spheres3d(sbst$coords[ind, ])

plot(face, new = FALSE, col = face$kappa1)
plot(face, new = FALSE)
plot(mouth, type = "spheres", col = "green", add= TRUE)
plot(sbst, type = "spheres", col = "green", add= TRUE)

cat(i, "\n")
scan()
}

plot(face, new = FALSE, col = face$kappa1)

nose   <- index.face3d(nose,   distance = 3, directions = TRUE, overwrite = TRUE)
mouth  <- index.face3d(mouth,  distance = 3, directions = TRUE, overwrite = TRUE)
sulcus <- index.face3d(sulcus, distance = 3, directions = TRUE, overwrite = TRUE)

# Warp the template

warp.template <- function(lmks, template) {
   ind   <- match(rownames(lmks), rownames(template$lmks))
   mat1  <- rbind(template$lmks[ind, ], template$lmks, template$curves, template$coords)
   wp    <- warp.face3d(mat1, lmks, subset = (nrow(lmks) + 1):nrow(mat1),
                        general = TRUE, project = FALSE)
   model <- list(coords = wp[nrow(template$lmks) + nrow(template$curves) + 1:nrow(template$coords), ],
                 triples = template$triples, lmks = wp[1:nrow(template$lmks), ],
                 curves = wp[nrow(template$lmks) + 1:nrow(template$curves), ])
   rownames(model$lmks)   <- rownames(template$lmks)
   rownames(model$curves) <- rownames(template$curves)
   class(model) <- "face3d"
   model
}

pcmodel <- function(scores, lmk.names, gpa, face, template) {
   lmks <- gpa$mshape + sweep(matrix(gpa$pcar[ , 1], ncol = 3), 1, scores * gpa$pcasd[1], "*") 
   rownames(lmks) <- lmk.names
   warp.template(lmks[lmk.names, ], template)
}

template <- template.male
clr <- rainbow_hcl(2)

i <- 85
load(fls[i])
face$lmks <- NULL
face <- initiallandmarks.face3d(face, orient = TRUE, monitor = FALSE, graphics = FALSE, overwrite = TRUE)
face <- subset(face, edist.face3d(face$coords, face$lmks["pn", ]) < 100)
plot(face, type = "mesh", col = clr[1], new = FALSE)
spheres3d(face$lmks, col = "yellow")

model <- warp.template(face$lmks[c("pn", "se", "enL", "enR"), ], template)
model <- subset(model, edist.face3d(model$coords, model$lmks["pn", ]) < 100)
plot(model, type = "mesh", colour = clr[2], add= TRUE)
pcm <- pcmodel(c(-3, rep(0, 12)), lmk.names, gpa, face, model)
pcm <- warp.template(face$lmks[c("pn", "se", "enL", "enR"), ], pcm)
pop3d()
plot(pcm, type = "mesh", colour = clr[2], add= TRUE)


open3d()
spheres3d(gpa$rotated)
spheres3d(gpa$mshape)
shp <- gpa$mshape - 3 * gpa$pcasd[1] * matrix(gpa$pcar[ , 1], ncol = 3)

pop3d()
spheres3d(shp, col = "red", radius = 1.1)


