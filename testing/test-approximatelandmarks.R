#     Approximate location of landmarks on the face

setwd("~/OneDrive - University of Glasgow/research/face3d")

install.packages("face3d", repos = NULL, type = "source")
library(face3d)
library(rgl)
library(fields)
library(geometry)
library(shapes)
library(MASS)

# Test on controls

fls <- list.files("~/Desktop/Glasgow-controls", recursive = TRUE, full.names = TRUE)
fls <- fls[-grep("details", fls)]
fls <- fls[-grep("lmks", fls)]
fls <- fls[-182]      # main area of positive curvature is the clothing
fls <- fls[-61]       # eye ridge instead of pn

landmark.names <- c("pn", "enL", "enR", "se")
# landmark.names <- "none"
# landmark.names <- "pn"

# Find pn, se (and enL/R)

dst <- matrix(NA, nrow = length(fls), ncol = 4, dimnames = list(NULL, landmark.names))
for (i in 1:length(fls)) {
   cat(i, "")
   load(fls[i])
   class(face) <- "face3d"
   face$landmarks <- NULL
   s.spacing  <- 10
   s.distance <- 45
   trim       <- 30
   face <- approximatelandmarks.face3d(face, landmark.names,
                                       sample.spacing = s.spacing, trim = trim,
                                       distance = 10, sample.distance = s.distance,
                                       monitor = 3)
   plot(face)
   clr <- rep("blue", nrow(face$landmarks))
   clr[grep("R", rownames(face$landmarks))] <- "green"
   clr[grep("L", rownames(face$landmarks))] <- "orange"
   spheres3d(face$landmarks, radius = 3, col = clr)

   if ("lmks" %in% names(face))
      dst[i, ] <- apply(face$landmarks - face$lmks[landmark.names, ], 1,
                              function(x) sqrt(sum(x^2)))
   save(face, file = fls[i])
   
   # snapshot3d(paste("~/Desktop/temp/temp_", i, ".png", sep = ""))
   # snapshot3d(paste("~/Desktop/temp/temp_crv", i, ".png", sep = ""))
}
save(dst, file = "~/Desktop/discrepancies.Rda")

# load("lmks-liberty.rda")
# load(fls[81])
# rownames(lmks.liberty) <- rownames(face$lmks)
# save(lmks.liberty, file = "lmks-liberty.rda")

# Find acL/R

load("lmks-liberty.rda")
lmks.liberty <- lmks.liberty[-match(c("tL", "tR", "oiL", "oiR"), rownames(lmks.liberty)), , ]
landmark.names <- c("pn", "se", "acL", "acR")
lmks   <- lmks.liberty[landmark.names, , ]

gpa     <- gpa.face3d(lmks, scale = FALSE)
mn.popn <- gpa$mean
n.lmks  <- nrow(lmks)
tan     <- apply(sweep(gpa$aligned, 1:2, gpa$mean), 3, c)
mnt     <- apply(tan, 1, mean)
covt    <- cov(t(tan))

i <- 106

# Look at the curvature characteristics at acL and acR
for (i in 107:length(fls)) {
   cat(i, "")
   load(fls[i])
   if ("coords" %in% names(face)) {
      face$vertices  <- face$coords
      face$coords    <- NULL
   }
   if ("triples" %in% names(face)) {
      face$triangles <- matrix(face$triples, ncol = 3, byrow = TRUE)
      face$triples   <- NULL
   }
   rnms <- rownames(face$landmarks)
   if (all(c("acL", "acR") %in% rnms))
      face$landmarks <- face$landmarks[-match(c("acL", "acR"), rnms), ]
   face <- findac(face, monitor = 1)
   face$landmarks
   plot(face)
   spheres3d(face$landmarks, col = "yellow", radius = 2)
   save(face, file = fls[i])
}

# Inspect the results
for (i in 1:length(fls)) {
   cat(i, "")
   load(fls[i])
   plot(face, col = "kappa2")
   spheres3d(face$landmarks, col = "yellow", radius = 2)
   scan()
}


# Construct a likelihood function for se

for (i in 62:length(fls)) {
   cat(i, "")
   load(fls[i])
   lmk.initial <- face$landmarks["se", ]
   patch <- lmklik(face, "se", "pn", monitor = 0)
   face$landmarks["se", ] <- patch$landmark
   plot(face)
   plot(patch$subset, display = "spheres", col = patch$subset$pattern, add = TRUE)
   spheres3d(lmk.initial, radius = 2)
   spheres3d(face$lmks["se", ], col = "blue", radius = 2)
   spheres3d(face$landmarks["se", ], col = "red", radius = 2)
   spheres3d(face$landmarks, col = "yellow", radius = 1.5)
   snapshot3d(paste("~/Desktop/temp/temp_", i, ".png", sep = ""))
}

# Check the position of se in manual and approximate landmarks
for (i in 114:length(fls)) {
   cat(i, "")
   load(fls[i])
   plot(face)
   dst   <- c(rdist(t(face$landmarks["se", ]), face$vertices))
   face  <- index.face3d(face, subset = dst < 20, distance = 10,
                         directions = TRUE, overwrite = TRUE)
   sbst0 <- subset(face, dst < 20)
   sbst  <- subset(face, dst < 10)
   
   # plot(sbst, col = "shape index")
   # plot(sbst, col = sbst$kappa1)
   # plot(sbst, col = sbst$kappa2)
   # plot(sbst, col = sbst$kappa1 * sbst$kappa2)
   # plot(sbst, display = "direction 1", add = TRUE)
   # plot(sbst, display = "direction 2", add = TRUE)
   
   nvert <- nrow(sbst$vertices)
   crv1  <- numeric(0)
   crv2  <- numeric(0)
   pmax  <- matrix(nrow = 0, ncol = 3)
   # edge  <- edges.face3d(sbst)
   # dst   <- rdist(sbst$vertices, sbst$vertices[edge[[1]], ])
   # dst   <- apply(dst, 1, min)
   sbst0$directions <- NULL
   for (j in 1:nvert) {
      sbst1 <- sbst0
      ppath <- planepath.face3d(sbst1, sbst$vertices[j, ], direction = sbst$directions[ , 1, j],
                                rotation = 0, bothways = TRUE)
      gcrv  <- gcurvature.face3d(ppath$path, 3)
      jdst  <- c(rdist(t(sbst$vertices[j, ]), gcrv$resampled.curve))
      jarc  <- which.min(jdst)
      gcrvj <- gcrv$gcurvature[jarc]
      # plot(gcrv$arclength, gcrv$gcurvature)
      pmax  <- rbind(pmax, gcrv$pos.max)
      crv1  <- c(crv1, gcrvj / gcrv$gcurvature[gcrv$ind.max])
      ppath <- planepath.face3d(sbst1, sbst$vertices[j, ], direction = sbst$directions[ , 2, j],
                                rotation = 0, bothways = TRUE)
      gcrv  <- gcurvature.face3d(ppath$path, 3)
      jdst  <- c(rdist(t(sbst$vertices[j, ]), gcrv$resampled.curve))
      jarc  <- which.min(jdst)
      gcrvj <- gcrv$gcurvature[jarc]
      # plot(gcrv$arclength, gcrv$gcurvature)
      pmax  <- rbind(pmax, gcrv$pos.max)
      # crv2  <- c(crv2, gcrvj / gcrv$gcurvature[gcrv$ind.max])
      crv2  <- c(crv2, gcrvj)
   }
   crv <- crv1 * crv2
   dst <- c(rdist(t(face$landmarks["se", ]), sbst0$vertices))
   plot(subset(sbst0, dst > 10))
   plot(sbst, col = crv2, add = TRUE)
   plot(sbst, display = "principal 1", add = TRUE)
   spheres3d(face$landmarks["se", ], col = "red")
   mode <- mode.face3d(sbst, crv, 5)$mode
   # mode <- sbst$vertices[which.max(crv), ]
   spheres3d(mode, col = "blue")
   
   ppath <- planepath.face3d(face, face$landmarks["pn", ], face$landmarks["se", ], si.target = 1)
   plot(face)
   spheres3d(ppath$path)
   spath <- smoothpath.face3d(face, face$landmarks["pn", ], face$landmarks["se", ],
                              penalty = 0.02, ppath = ppath, si.target = 1)
   plot(face)
   spheres3d(face$landmarks[c("se", "pn"), ], col = "red", radius = 2)
   spheres3d(spath)
    
   plot(sbst, col = "kappa2")
   plot(sbst, display = "principal 2", add = TRUE)
   
   snapshot3d(paste("~/Desktop/temp/temp_", i, ".png", sep = ""))
}




for (i in 1:length(fls)) {
   load(fls[i])
   dst   <- c(rdist(t(face$landmarks["se", ]), face$coords))
   face1 <- index.face3d(face, distance = 5, subset = dst < 10, overwrite = TRUE)
   sbst <- subset(face1, dst < 10)
   sbst$gc <- sbst$kappa1 * sbst$kappa2
   plot(sbst, col = sbst$gc)
   spheres3d(face$landmarks["se", ], col = "yellow")
   mode <- mode.face3d(sbst, -sbst$gc, threshold = 8)$mode
   spheres3d(mode, col = "red")
   scan()
}




# Optimise over location and orientation - no likelihood yet
dst  <- c(rdist(t(face$landmarks["pn", ]), face$coords))
face <- subset(face, dst < 100)



face <- index.face3d(face, overwrite = TRUE)
plot(face, col = "shape index")
plot(face, col = face$kappa1)
plot(face, col = face$kappa2)
plot(face, col = face$kappa1 * face$kappa2)

fit <- function()

# Prior information

plot(face)
spheres3d(face$landmarks, col = "yellow", radius = 3)
n.lmks <- nrow(lmks.liberty)
for (i in 1:n.lmks) {
   ind <- c(i, i + n.lmks, i + 2 * n.lmks)
   plot3d(ellipse3d(covt[ind, ind], centre = gpa$mean[i, ]), 
          col = "lightblue", alpha = 0.5, add = TRUE)
}

# Procrustes registration on the four landmarks
opa <- opa.face3d(gpa$mean[landmark.names, ], face$landmarks[landmark.names, ],
                  gpa$mean, scale = FALSE, return.parameters = TRUE)
for (i in 1:n.lmks) pop3d()
for (i in 1:n.lmks) {
   ind <- c(i, i + n.lmks, i + 2 * n.lmks)
   plot3d(ellipse3d(opa$rotation %*% covt[ind, ind] %*% t(opa$rotation), centre = opa$opa[i, ]),
          col = "lightblue", alpha = 0.5, add = TRUE)
}

# Conditional distribution of the remaining landmark.names
ind.cond <- match(landmark.names, rownames(lmks.liberty))
ind.cond <- c(ind.cond, ind.cond + n.lmks, ind.cond + 2 * n.lmks)
mnt.cond <- mnt[-ind.cond] + covt[-ind.cond, ind.cond] %*%
   ginv(covt[ind.cond, ind.cond]) %*%
   (c(face$landmarks[landmark.names, ] - gpa$mean[ind.cond]))
spheres3d(matrix(mnt.cond + gpa$mean[-ind.cond], ncol = 3), col = "red")


# Other landmarks

# sn

dst      <- rdist(t(face$lmks["pn", ]), face$coords)
drn      <- face$landmarks["pn", ] - face$landmarks["se", ]
nrm      <- face$normals[which.min(dst), ]
drn      <- crossproduct(drn, nrm)
prp      <- sweep(sbst$coords, 2, face$landmarks["pn", ]) %*% drn
# Adjust the distance from the prior information: more than to sn but less than to lips?
sbst     <- subset(face, dst < 40 & abs(prp) < 5)
crv      <- sbst$kappa1 * sbst$kappa2
mode     <- mode.face3d(sbst, -crv, 5)$mode
spheres3d(mode, col = "yellow")

sbst     <- subset(sbst, sbst$crv < 0)
plot(sbst, col = -sbst$crv, new = TRUE)
sbst <- subset(sbst, sbst$crv < median(sbst$crv))
plot(sbst, col = -sbst$crv)
plot(sbst, col = pmax(0, -sbst$kappa1 * sbst$kappa2))
plot(sbst, col = "shape index")
ppath <- planepath.face3d(sbst, face$landmarks["pn", ],
                          direction = face$landmarks["pn", ] - face$landmarks["se", ],
                          rotation = 0)$path
spheres3d(ppath)


for (i in 1:length(fls)) {
   cat(i, "")
   load(fls[i])
   
   sbst  <- subset(face, face$shape.index > 0.2)
   parts <- connected.face3d(sbst)
   sbst  <- subset(sbst, parts == 1)
   chull <- convhulln(sbst$coords)
   cpts  <- unique(c(chull))
   
   plot(sbst, col = sbst$kappa2)
   # spheres3d(sbst$coords[cpts, ], radius = 1, col = "white")
   
   # sbst1 <- subset(sbst, sbst$kappa2 < quantile(sbst$kappa2, 0.10))
   # plot(sbst1, col = sbst1$kappa2)
   # 
   # edges <- edges.face3d(sbst)
   # lapply(edges, function(x) lines3d(sbst$coords[x, ], lwd = 2, col = "blue"))
   # 
   # dst  <- rdist(sbst$coords[cpts, ], sbst$coords[cpts, ])
   # ind  <- which(apply(dst, 1, min) > 10)
   # 
   # crv <- -(sbst$kappa1[cpts] + sbst$kappa2[cpts])
   # dst <- rdist(sbst$coords[cpts, ], sbst$coords)
   # val <- apply(dst, 1, function(x) length(which(x < 5)))
   # dst <- rdist(sbst$coords[cpts, ], sbst$coords[cpts, ])
   # val <- apply(dst, 1, function(x) length(which(x < 5))) / val
   # val <- apply(dst, 1, function(x) mean(crv[x < 5])) * val
   # pn  <- sbst$coords[cpts, ][which.max(val), ]
   # spheres3d(pn, col = "yellow", radius = 3)
   
   # plot(sbst, col = sbst$kappa2)
   
   snapshot3d(paste("~/Desktop/temp_pos/temp_", i, ".png", sep = ""))
}


# What other approximate landmarks can we find?

for (i in 1:length(fls)) {
   
   cat(i, "")
   load(fls[i])
   plot(face)
   spheres3d(face$landmarks, col = "green", radius = 2)

   pn   <- face$landmarks["pn", ]
   dst  <- c(rdist(face$coords, t(pn)))
   ind  <- which(dst < 50)
   sbst <- subset(face, ind)
   plot(sbst)
   spheres3d(pn, col = "green", radius = 2)
   sbst <- index.face3d(sbst, overwrite = TRUE)

   nrm   <- face$normals[which.min(dst), ]
   drn   <- c(crossproduct(face$landmarks["se", ] - pn, nrm))
   pathL <- planepath.face3d(sbst, pn, direction =  drn, rotation.range = pi / 4, si.target = 1)
   pathR <- planepath.face3d(sbst, pn, direction = -drn, rotation.range = pi / 4, si.target = 1)
   acL   <- gcurvature.face3d(pathL$path, 10)$pos.max
   acR   <- gcurvature.face3d(pathR$path, 10)$pos.max
   
   drn   <- pn -face$landmarks["se", ]
   pathM <- planepath.face3d(sbst, pn, direction =  drn, rotation.range = pi / 4, si.target = 1)
   sn    <- gcurvature.face3d(pathM$path, 10)$pos.max
   
   plot(sbst, col = sbst$kappa2)
   spheres3d(pn, col = "green", radius = 2)
   spheres3d(pathL$path, col = "yellow")
   spheres3d(pathR$path, col = "yellow")
   spheres3d(pathM$path,  col = "yellow")
   spheres3d(acL, col = "green", radius = 2)
   spheres3d(acR, col = "green", radius = 2)
   spheres3d(sn,  col = "green", radius = 2)
   
   snapshot3d(paste("~/Desktop/temp_nose/temp_", i, ".png", sep = ""))
}


face     <- template_male
selected <- sample(1:nrow(face$coords), 1000)
face     <- index.face3d(face, subset = selected, distance = 50, overwrite = TRUE)
wp       <- warp.face3d(face$coords[selected, ], face$kappa1[selected], face$coords)
plot(face, col = wp)

# face.1  <- index.face3d(face, distance = 10,  overwrite = TRUE)
# face.5  <- index.face3d(face, distance = 50,  overwrite = TRUE)
# face.10 <- index.face3d(face, distance = 100, overwrite = TRUE)
# face.15 <- index.face3d(face, distance = 150, overwrite = TRUE)
# 
# plot(face.1,  col = "shape index")
# plot(face.5,  col = "shape index", new = TRUE)
# plot(face.10, col = "shape index", new = TRUE)
# plot(face.15, col = "shape index", new = TRUE)
