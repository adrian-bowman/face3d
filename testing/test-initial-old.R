#     Test initiallandmarks.face3d

setwd("~/ownCloud/Face3D_0.1-1")

library(Face3D)
library(fields)
library(rgl)
library(geometry)

source("Face3D/R/a-orient.R")
source("Face3D/R/a-midline.R")
source("Face3D/R/a-ridge2d.R")
source("Face3D/R/a-display.R")
source("Face3D/R/a-initiallandmarks.R")
source("Face3D/R/a-initialcurves.R")
source("Face3D/R/a-utilities.R")
source("Face3D/R/a-planepath.R")
source("Face3D/R/a-gcurvature.R")

fls <- list.files("Data-Liberty", full.names = TRUE)

# Store lmks and curves
# for (i in 1:length(fls)) {
#    load(fls[i])
#    face$manual.lmks   <- face$lmks
#    face$manual.curves  <- face$curves
#    face$lmks          <- NULL
#    face$curves        <- NULL
#    save(face, file = fls[i])
#    cat(i, "")
# }

# Orient
for (i in 1:length(fls)) {
   load(fls[i])
   face   <- orient.face3d(face)
   save(face, file = fls[i])
   cat(i, "")
}

# Inspect
for (i in 1:length(fls)) {
   load(fls[i])
   view3d(0, 0)
   plot(face, new = FALSE)
   face$normals <- face$normals + face$coords
   face$coords  <- face$orient(face$coords,  face$centroid, face$rotation.angles)
   face$normals <- face$orient(face$normals, face$centroid, face$rotation.angles)
   face$normals <- face$normals - face$coords
   scan()
   plot(face, new = FALSE)
   scan()
   spheres3d(face$coords[face$nearest, ], col = "blue")
   scan()
   segments3d(matrix(c(t(cbind(face$coords, face$coords + face$normals))), ncol = 3, byrow = TRUE))
   scan()
}

# Find initial landmarks
for (i in (1:length(fls))) {
   load(fls[i])
   face$lmks <- NULL
   face <- initiallandmarks.face3d(face, orient = FALSE)
   save(face, file = fls[i])
   cat(i, "")
}

# Inspect
for (i in (1:length(fls))) {
   load(fls[i])
   clear3d()
   mfrow3d(1, 2, sharedMouse = TRUE)
   plot(face, new = FALSE)
   spheres3d(face$lmks, radius = 3, col = "yellow")
   view3d(180 * face$rotation.angles[2] / pi, 180 * face$rotation.angles[1] / pi, zoom = 0.7)
   next3d()
   plot(face, new = FALSE)
   spheres3d(face$lmks, radius = 3, col = "yellow")
   view3d(180 * face$rotation.angles[2] / pi - 90, 180 * face$rotation.angles[1] / pi, zoom = 0.7)
   cat(i, "")
   scan()
}

# Compute shape index
for (i in 1:length(fls)) {
   load(fls[i])
   ind  <- (edist.face3d(face$coords, face$lmks["pn", ]) < 100)
   face <- index.face3d(face, distance = 10, subset = ind)
   save(face, file = fls[i])
   cat(i, "")
}

# Inspect
for (i in 1:length(fls)) {
   load(fls[i])
   clear3d()
   mfrow3d(1, 2)
   plot(face, new = FALSE, colour = face$kappa1)
   view3d(180 * face$rotation.angles[2] / pi, 180 * face$rotation.angles[1] / pi, zoom = 0.7)
   next3d()
   plot(face, new = FALSE, colour = face$kappa2)
   view3d(180 * face$rotation.angles[2] / pi, 180 * face$rotation.angles[1] / pi, zoom = 0.7)
   scan()
}

# Find initial curves
for (i in (1:length(fls))) {
   load(fls[i])
   face$curves <- NULL
   face <- initialcurves.face3d(face)
   save(face, file = fls[i])
   cat(i, "")
}

i <- 3
load(fls[i])
face$curves <- NULL
face <- initialcurves.face3d(face, graphics = TRUE, new.window = FALSE)
plot(face, new = FALSE)
spheres3d(face$lmks, radius = 2, col = "yellow")
sapply(face$curves, spheres3d, col = "red")


# Inspect
for (i in (1:length(fls))) {
   load(fls[i])
   plot(face, new = FALSE)
   spheres3d(face$lmks, radius = 2, col = "yellow")
   sapply(face$curves, spheres3d, col = "red")
   # snapshot3d(paste("temp/image-", i, ".png", sep = ""))
   cat(i, "")
   scan()
}

i <- 21
load(fls[i])
face$lmks <- NULL
face <- initiallandmarks.face3d(face, orient = FALSE, overwrite = TRUE, graphics = TRUE)
plot(face, new = FALSE)
spheres3d(face$lmks, radius = 2, col = "yellow")
face <- initialcurves.face3d(face, graphics = TRUE, new.window = FALSE)





i <- 1
load(fls[i])
# face$lmks.orig     <- face$lmks
face$curves.manual <- face$curves
face$lmks          <- NULL
face$curves        <- NULL
face <- initiallandmarks.face3d(face, orient = FALSE, graphics = TRUE)
plot(face, new = FALSE)
spheres3d(face$lmks, radius = 2, col = "yellow")
face <- initialcurves.face3d(face, graphics = TRUE)
plot(face, new = FALSE)
spheres3d(face$lmks, radius = 2, col = "yellow")
lapply(face$curves, spheres3d, col = "red")

face1 <- subset(face, edist.face3d(face$coords, face$lmks["pn", ]) < 100)
plot(face1, new = FALSE, col = "shape index")
plot(face1, new = FALSE, col = face1$kappa1)
plot(face1, new = FALSE, col = face1$kappa2)
plot(face1, new = FALSE, col = face1$kappa1 * face1$kappa2)


for (i in (1:length(fls))) {
   load(fls[i])
   face <- subset(face, edist.face3d(face$coords, face$lmks["pn", ]) < 40)
   mfrow3d(1, 2)
   plot(face, new = FALSE, col = face$kappa1)
   next3d()
   plot(face, new = FALSE, col = face$kappa2)
   scan()
}

i <- 5
load(fls[i])
se.old <- face$lmks["se", ]
face <- initialcurves.face3d(face, graphics = TRUE)
plot(face, new = FALSE)
spheres3d(face$lmks, radius = 2, col = "blue")
spheres3d(se.old, radius = 2, col = "yellow")
pp <- planepath.face3d(face, face$lmks["acL", ], face$lmks["acR", ])
spheres3d(pp$path)
gc <- gcurvature.face3d(pp$path, 6)
plot(gc$gcurvature ~ gc$arc.length)

plot(initial$sbst, new = FALSE, col = initial$sbst$kappa2)
lapply(initial$curves, spheres3d, col = "red")


sbst <- index.face3d(sbst, distance = 3, directions = TRUE, overwrite = TRUE)
plot(sbst, type = "mesh", new = FALSE, col = sbst$kappa1)
spheres3d(x)
plot(sbst, type = "points", new = FALSE, col = sbst$kappa1)
crds <- cbind(sbst$coords, sbst$coords + t(sbst$directions[ , 1, ]))
crds <- matrix(c(t(crds)), ncol = 3, byrow = TRUE)
segments3d(crds)
drn  <- c(crossproduct(face$lmks["pn", ] - face$lmks["se", ], nrm))



face1 <- subset(face, edist.face3d(face$coords, x) >= 10)
plot(face1, new = FALSE)
plot(sbst, add = TRUE, col = sbst$kappa2)

plot(face, new = FALSE, col = face$kappa1)
plot(face, col = face$kappa2)



# Find the midline ridge
plot(sbst2, new = FALSE, col = sbst2$kappa2)
spheres3d(face$lmks[c("se", "pn"), ], col = "red")
ppv <- planepath.face3d(sbst2, face$lmks["se", ], face$lmks["pn", ])
spheres3d(ppv$path)
drn <- c(crossproduct(face$lmks["se", ] - face$lmks["pn", ], ppv$normal))
pph <- planepath.face3d(sbst2, face$lmks["pn", ], si.target = 1, graphics = FALSE,
                        direction = drn, rotation.range = pi/8, bothways = TRUE)
spheres3d(pph$path)





# Check the template
# load("~/ownCloud/Face3D_0.1-1/Face3D/data/templateMale.rda")
# template <- template.male
# plot(template, type = "mesh", new = FALSE)
# plot(template, type = "normals", add = TRUE)
# spheres3d(template$curves, col = "red", radius = 0.2)
# lns <- matrix(c(t(cbind(template$curves, template$curves + 5 * cnormals))), ncol = 3, byrow = TRUE)
# segments3d(lns, col = "red")

# Locate sn by tracing minimum asymmetry down the columella
i <- 85
load(fls[i])
face$lmks <- NULL
face <- initiallandmarks.face3d(face, orient = TRUE, monitor = TRUE, graphics = FALSE, overwrite = TRUE)
face <- subset(face, edist.face3d(face$coords, face$lmks["pn", ]) < 100)
plot(face, new = FALSE)
spheres3d(face$lmks, col = "yellow")
drn  <- face$lmks["pn", ] - face$lmks["se", ]
reference <- planepath.face3d(face, face$lmks["pn", ], direction = face$lmks["pn", ] - face$lmks["se", ],
                              rotation = 0, boundary = 30, directions = TRUE)
spheres3d(reference$path)
mline <- midline.face3d(face, reference = reference, lambda = 0.005, df = 24, d.asym = 25, nv = 31,
                         quantile = 0.25, extension = c(0, 0), monitor = TRUE)


# Locate sn by the shape of planepaths heading south from pn
i <- 30
load(fls[i])
face <- initiallandmarks.face3d(face, orient = TRUE, monitor = TRUE, graphics = FALSE, overwrite = TRUE)
face <- subset(face, edist.face3d(face$coords, face$lmks["pn", ]) < 100)
drn  <- face$lmks["pn", ] - face$lmks["se", ]
pp   <- planepath.face3d(face, face$lmks["pn", ], rotation = 0, direction = drn, bothways = FALSE, boundary = c(100, 100))
n    <- pp$normal

plot(face, new = FALSE)
spheres3d(face$lmks[c("pn", "se"), ], col = "yellow")
segments3d(rbind(face$lmks["pn", ], face$lmks["pn", ] + 4 * n), col = j)

surf   <- matrix(nrow = 0, ncol = 200)
alen   <- matrix(nrow = 0, ncol = 200)
angles <- seq(-pi/5, pi/5, length = 21)
for (angle in angles) {
   dr  <- rgl::rotate3d(drn, angle, n[1], n[2], n[3])
   pp  <- planepath.face3d(face, face$lmks["pn", ], rotation = 0, direction = dr, bothways = FALSE, boundary = c(100, 100))
   ind <- (pp$arclength <= 1 * sqrt(sum(drn^2)))
   al  <- pp$arclength[ind]
   pp  <- pp$path[ind, ]
   dr  <- dr / sqrt(sum(dr^2))
   axs <- cbind(c(crossproduct(dr, -n)), -n, dr)
   ppt <- sweep(pp, 2, face$lmks["pn", ]) %*% axs
   gc  <- gcurvature.face3d(ppt, 6)
   if (nrow(rgl.ids()) == 6) for (j in 1:3) pop3d()
   for (j in c(1, 3))
      segments3d(rbind(face$lmks["pn", ], face$lmks["pn", ] + 4 * axs[ , j]), col = j)
   spheres3d(pp, col = "blue")
   # plot(gc$arc.length, gc$gcurvature,
   #      ylim = c(0, 0.35), main = round(angle, 2))
   # plot(gc$arc.length, gc$gcurvature * sign(gc$d2.y),
   #      ylim = c(-0.35, 0.35), main = round(angle, 2))
   plot(gc$arc.length, gc$d2.z, xlim = c(0, max(al)),
        ylim = c(-0.35, 0.35), main = round(angle, 2))
   surf <- rbind(surf, gc$d2.z)
   alen <- rbind(alen, gc$arc.length)
}
image(x = angles, z = surf, xlab = "angle")
contour(x = angles, z = surf, add = TRUE)
intrp <- interp(rep(angles, each = 200), c(t(alen)), c(t(surf)))
image(intrp)
contour(intrp$x, intrp$y, intrp$z, add = TRUE)
ind <- which.max(intrp$z)
rw  <- 1 + (ind - 1) %/% nrow(intrp$z)
cl  <- 1 + (ind - 1)  %% nrow(intrp$z)
points(intrp$x[cl], intrp$y[rw], pch = 16, col = "blue")
persp3d(intrp$x, intrp$y, intrp$z, col = "darkgreen")

# Inspect the kappa1 and kappa2 values across the face
i <- 1
load(fls[i])
plot(face, new = FALSE)
quants <- quantile(face$kappa1, seq(0.7, 1, by = 0.05), na.rm = TRUE)
for (i in 1:length(quants)) {
   spheres3d(face$coords[face$kappa1 > rev(quants)[i], ], col = "yellow")
   scan()
}

plot(face, col = "shape index",   new = FALSE)
plot(face, new = FALSE)
spheres3d(face$coords[-face$kappa2 > 0.4, ], col = "red")
spheres3d(face$coords[ face$kappa1 > 0.4, ], col = "red")
k1 <- face$kappa1
k1[k1 < 0.3] <- 0
plot(face, col = k1 + 1, new = FALSE)
plot(face, col = face$kappa2 + 1, new = FALSE)


# Look for 6 positions with cup shapes
for (i in 1:length(fls)) {
   load(fls[i])
   plot(face, new = FALSE)
   ed     <- edist.face3d(face$coords, face$lmks["pn", ])
   ind    <- (ed < 70)
   face   <- index.face3d(face, distance = 10, subset = ind)
   ind    <- (face$shape.index < -0.65) &
                  apply(face$coords, 1, function(x) sqrt(sum((x - face$lmks["pn", ])^2))) > 30
   sbst   <- subset.face3d(face, ind)
   spheres3d(sbst$coords, col = "blue")
   parts <- connected.face3d(sbst)
   sbst  <- subset(sbst, parts <= 6)
   spheres3d(sbst$coords, col = "green", radius = 1.5)
   scan()
}

fls <- list.files("/Volumes/Face3D/Glasgow-controls/Liberty/Liberty-Controls-dmp/Liberty-Controls-dmp-final/Data-Liberty/",
                  full.names = TRUE)
load(fls[1])
plot(face)
cat(paste0(rep("-", length(fls)), collapse = ""), "\n")
for (i in 1:length(fls)) {
   cat(".")
   load(fls[i])
   face$lmks <- NULL
   face <- initiallandmarks.face3d(face)
   plot(face, new = FALSE)
   spheres3d(face$lmks, col = "green", radius = 4)
   snapshot3d(paste("temp/temp-", i, ".png", sep = ""))
}
cat("\n")

load(fls[44])
face$lmks <- NULL
face <- initiallandmarks.face3d(face)
plot(face, new = FALSE)
spheres3d(face$lmks, radius = 3)

# Load the template
load("~/owncloud/Shared/visualisation-shape/examples/data/template_male_symmetric.dmp")
load("~/owncloud/Shared/visualisation-shape/examples/data/template_female_symmetric.dmp")
template <- template.male
template <- index.face3d(template)
plot(template, col = "shape index", new = FALSE)
plot(template, col = 1+template$kappa1, new = FALSE)
plot(template, type = "mesh", new = FALSE)
spheres3d(template$lmks, col = "yellow")
spheres3d(template$curves, col = "red")

ind  <- match(rownames(face$lmks), rownames(template$lmks))
# Warp only landmarks
mat1 <- rbind(template$lmks[ind, ], template$lmks[-ind, ])
wp   <- warp.face3d(mat1, face$lmks, subset = length(ind) + 1:nrow(template$lmks[-ind, ]),
                    general = TRUE, project = FALSE)
spheres3d(wp, col = "green", radius = 3)

# Warp the template as well
mat1 <- rbind(template$lmks[ind, ], template$lmks[-ind, ], template$curves, template$coords)
wp   <- warp.face3d(mat1, face$lmks, subset = (length(ind) + 1):nrow(mat1), general = TRUE, project = FALSE)
nnew <- nrow(template$lmks) - nrow(face$lmks)
lmks.new <- wp[1:nnew, ]
model    <- list(curves = wp[nnew + 1:nrow(template$curves), ],
                 coords = wp[nnew + nrow(template$curves) + 1:nrow(template$coords), ],
                 triples = template$triples)
class(model) <- "face3d"
model <- normals.face3d(model)
plot(face, new = FALSE, type = "mesh")
plot(model, add = TRUE, type = "mesh")
spheres3d(face$lmks, radius = 3)
spheres3d(lmks.new, col = "blue", radius = 3)
spheres3d(model$curves, col = "yellow", radius = 2)

template   <- normals.face3d(template)
wp.ind     <- apply(rdist(wp, template$coords), 1, which.min)
wp.normals <- template$normals[ind, ]
for (i in 1:nrow(wp)) {
   prj     <- c(sweep(face$coords, 2, wp[i, ]) %*% wp.normals[i, ])
   dst     <- apply((face$coords - outer(prj, wp.normals[i, ]))^2, 1, sum)
   wp[i, ] <- face$coords[which.min(dst)]
}
spheres3d(wp, col = "green", radius = 3)

plot(template, new = FALSE)
spheres3d(template$coords)
spheres3d(template$meshes, col = "green")
spheres3d(template$full.mesh, col = "green", radius = 1.5)
spheres3d(template$meshes, col = "red", radius = 2)

