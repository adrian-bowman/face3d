#     Test initiallandmarks.face3d

# setwd("/Volumes/Face3D/Glasgow-controls/Liberty/Liberty-Controls-dmp/Liberty-Controls-dmp-final/Data-Liberty/")
setwd("/Volumes/Face3D/Smile Train/Brazil/captures-R/controls/")

library(Face3D)
library(fields)
library(rgl)
library(geometry)

# fls <- list.files(full.names = TRUE)
fls <- list.files(pattern = "resolution-0.002")

# Remove problematic images
"2019-03-11-10-30-23-resolution-0.002.Rdata"
"2019-03-12-20-16-43-resolution-0.002.Rdata"
"2019-03-12-20-16-43-resolution-0.002.Rdata"
"2019-03-13-09-04-00-resolution-0.002.Rdata"
"2019-03-13-09-04-56-resolution-0.002.Rdata"
"2019-03-16-09-47-00-resolution-0.002.Rdata"

# Shape index issues
"2019-03-14-06-12-47-resolution-0.002.Rdata"
"2019-03-14-06-22-52-resolution-0.002.Rdata"
"2019-03-15-10-04-09-resolution-0.002.Rdata"
"2019-03-15-13-57-32-resolution-0.002.Rdata"
"2019-03-16-09-40-35-resolution-0.002.Rdata"
"2019-03-16-09-46-17-resolution-0.002.Rdata"

# See test-ball.R for initial detection of the face

for (i in 49:length(fls)) {
   cat(i, fls[i], "\n")
   load(fls[i])
   plot(face, new = FALSE)
   
   face <- index.face3d(face, overwrite = TRUE)
   plot(face, new = FALSE, col = "shape index")
   gc   <- face$kappa1 * face$kappa2
   ind  <- (face$kappa1 < 0) | (face$kappa2 < 0)
   gc[ind] <- 0
   plot(face, new = FALSE, col = gc, key = TRUE)

   # face$lmks <- NULL
   # face <- initiallandmarks.face3d(face, "pn", monitor = 2, graphics = "add")
   # # face <- rotate.face3d(face, c("pn", "se", "enL", "enR"))
   # plot(face, new = FALSE)
   # spheres3d(face$lmks, radius = 2, col = "red")
   scan()
}

# Look at integral of gc over parts with shape index > 0.75
sbst  <- subset(face, face$shape.index < -0.75)
sbst  <- subset(face, face$shape.index >  0.75)
gc    <- sbst$kappa1 * sbst$kappa2
plot(sbst, col = gc, breaks.type = "linear", key = TRUE, new = FALSE)
parts <- connected.face3d(sbst)
wt    <- rep(0, length(unique(parts)))
for (i in unique(parts)) wt[i] <- sum(gc[parts == i])
sbst1 <- subset(sbst, parts == which.max(wt))
gc1   <- sbst1$kappa1 * sbst1$kappa2
plot(sbst1, display = "spheres", col = gc1, range.colour = range(gc), add = TRUE)



i <- 5
load(fls[i])
plot(face, new = FALSE)
face <- index.face3d(face, overwrite = TRUE)
plot(face, new = FALSE, col = "shape index")
gc   <- face$kappa1 * face$kappa2
ind  <- (face$kappa1 > 0) | (face$kappa2 > 0)
gc[ind] <- 0
plot(face, new = FALSE, col = gc, breaks.type = "quantile", key = TRUE)
plot(face, new = FALSE, col = gc, breaks.type = "linear", key = TRUE)

edg <- edges.face3d(face)
lapply(edg, function (x) rgl::lines3d(face$coords[x, ], lwd = 2, col = "blue"))
spheres3d(face$coords[unlist(edg), ], col = "blue")

gc <- face$kappa1 * face$kappa2
plot(face, new = FALSE, colour = "shape index")
plot(face, new = FALSE, colour = gc)
plot(face, new = FALSE, colour = "grey")


i     <- 1
load(fls[i])
dst   <- c(rdist(t(face$lmks["pn", ]), face$coords))
face1 <- subset(face, dst < 70)
plot(face1, new = FALSE)
face1 <- index.face3d(face1, distance = 20, overwrite = TRUE)
plot(face1, colour = "shape index", new = FALSE)
gc    <- face1$kappa1 * face1$kappa2
gc[(face1$kappa1 < 0) | (face1$kappa2 < 0)] <- 0
plot(face1, colour = gc, new = FALSE)


# Look at positive and negative Gaussian curvature separately
gc    <- face$kappa1 * face$kappa2

face1 <- subset(face, (face$kappa1 < 0) & (face$kappa2 < 0) & (gc > median(gc[gc > 0])))
gc1   <- face1$kappa1 * face1$kappa2
plot(face1, new = FALSE, col = gc1, breaks.type = "linear", key = TRUE)
tr    <- (abs(gc1))^0.25 * sign(gc1)
plot(face1, new = FALSE, col = tr, breaks.type = "linear", key = TRUE)

face2 <- subset(face, (face$kappa1 > 0) & (face$kappa2 > 0) & (gc > quantile(gc[gc > 0], 0.75)))
gc2   <- face2$kappa1 * face2$kappa2
plot(face2, new = FALSE, col = gc2, breaks.type = "linear", key = TRUE)
tr <- (abs(gc2))^0.25 * sign(gc2)
plot(face2, new = FALSE, col = tr, breaks.type = "linear", key = TRUE)







# Orient and compute shape index for convenience
# fls <- list.files("~/Desktop/Data-Liberty/", full.names = TRUE)
# for (i in 1:length(fls)) {
#    load(fls[i])
#    ed     <- edist.face3d(face$coords, face$lmks["pn", ])
#    ind    <- (ed < 70)
#    # face   <- orient.face3d(face)
#    face   <- index.face3d(face, distance = 5, subset = ind, overwrite = TRUE)
#    # plot(face, new = FALSE)
#    # plot(face, col = "shape index", add = TRUE)
#    save(face, file = fls[i])
#    cat(i, "")
# }
# cat("\n")

# Find and store initial landmarks
# for (i in 1:length(fls)) {
#    load(fls[i])
#    lmk.nms        <- c("pn", "se", "enL", "enR")
#    lmks           <- face$lmks
#    face$lmks      <- NULL
#    face           <- initiallandmarks.face3d(face, lmk.nms)
#    face$lmks.orig <- lmks
#    save(face, file = fls[i])
#    cat(i, "")
# }

# Check the initial landmarks are ok and look at the curvature information
for (i in (1:length(fls))) {
   load(fls[i])
   face   <- subset(face, edist.face3d(face$coords, face$lmks["pn", ]) < 100)
   thresh <- quantile(face$kappa2[face$kappa2 < 0], 0.35, na.rm = TRUE)
   sbst1  <- subset(face,  face$kappa2 < thresh)
   sbst2  <- subset(sbst1, connected.face3d(sbst1) == 1)
   ind1   <- as.numeric(rownames(sbst1$coords))
   ind2   <- as.numeric(rownames(sbst2$coords))
   ind    <- ind1[ind2]
   # sbst2  <- index.face3d(sbst2, distance = 5, overwrite = TRUE)
   # face1  <- subset(face, -ind)
   # plot(face1, new = FALSE)
   plot(sbst2, new = FALSE, col = sbst2$kappa2)
   spheres3d(face$lmks[c("se", "pn"), ], col = "red")
   ppv1 <- planepath.face3d(sbst2, face$lmks["se", ], face$lmks["pn", ],
                            rotation.range = pi/8, bridge.gaps = TRUE)
   # ppv2 <- planepath.face3d(sbst2, face$lmks["pn", ],
   #                          rotation = "maximise", rotation.range = pi/8,
   #                          graphics = FALSE,
   #                          direction = face$lmks["pn", ] - face$lmks["se", ])
   spheres3d(ppv1$path, col = "red")
   # spheres3d(ppv2$path, col = "red")
   # nrm  <- apply(ppv1$shape$normals[ppv1$pts1, ], 2, mean)
   nrm <- ppv1$normal
   drn  <- c(crossproduct(face$lmks["se", ] - face$lmks["pn", ], nrm))
   pph1 <- planepath.face3d(sbst2, face$lmks["pn", ], si.target = 1,
                           direction = drn, normal = nrm, rotation.range = pi/8)
   drn  <- c(crossproduct(face$lmks["pn", ] - face$lmks["se", ], nrm))
   pph2 <- planepath.face3d(sbst2, face$lmks["pn", ], si.target = 1,
                           direction = drn, normal = nrm, rotation.range = pi/8)
   spheres3d(pph1$path, col = "red")
   spheres3d(pph2$path, col = "red")
   cat(i, "")
   snapshot3d(paste("temp/image-", i, ".png", sep = ""))
}

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

