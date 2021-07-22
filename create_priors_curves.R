#     Create priors for curves

# install.packages("~/research/face3d/face3d", repos = NULL, type = "source")
# library(face3d)
# library(rgl)
# library(fields)
# library(geometry)
# library(shapes)
# library(MASS)
# library(colorspace)

# Liberty's controls
fls <- list.files("../Glasgow-controls/Liberty", full.names = TRUE)

loadface <- function(x) {
   load(x)
   if ("coords" %in% names(face)) {
      face$vertices  <- face$coords
      face$triangles <- matrix(face$triples, ncol = 3, byrow = TRUE)
      face$coords    <- NULL
      face$triples   <- NULL
   }
   face$landmarks <- NULL
   class(face)    <- "face3d"
   invisible(face)
}

# Create an array of curves
curves <- template_male$curves
curves <- curves[curves[ , 3] > 0, ]
nms    <- rownames(curves)
nms    <- sub("upper-lip", "upper lip", nms)
nms    <- nms[-grep("cheek-eye", nms)]
nms    <- nms[-grep("philtrum lip", nms)]
pn.id  <- match("mid-line nasal profile 20", rownames(curves))
curves <- array(0, dim = c(length(nms), 3, length(fls)),
                dimnames = list(nms))
for (i in 1:length(fls)) {
   face <- loadface(fls[i])
   crvs <- face$curves
   curves[ , , i] <- crvs[nms, ]
}

gpa <- gpa.face3d(curves, scale = FALSE)
curves <- gpa$aligned
clear3d()
apply(curves, 3, function(x) spheres3d(x))

pca <- pca.face3d(curves)
names(pca)
head(pca$percent)

pc.high <- pca$mean + 2 * pca$sd[1] * matrix(pca$evecs[ , 1], ncol = 3)
pc.low  <- pca$mean - 2 * pca$sd[1] * matrix(pca$evecs[ , 1], ncol = 3)
clear3d()
spheres3d(pca$mean)
spheres3d(pc.high, col = "red")
spheres3d(pc.low,  col = "blue")

# Find the indices for differencing.
# This currently assumes equally-spaced points with consistent rowname ordering in each curve.
rnms <- rownames(pca$mean)
cnms <- rnms
for (i in 0:9) cnms <- gsub(as.character(i), "", cnms)
cnms <- unique(cnms)
dff.fn <- function(x) {
   nms <- grep(x, rnms)
   lng <- length(nms)
   cbind(nms[c(2:lng, lng)], nms[c(1, 1:(lng - 1))])
}
diffmat <- matrix(ncol = 2, nrow = 0)
for (nm in cnms) diffmat <- rbind(diffmat, dff.fn(nm))
cbind(rnms[diffmat[ , 1]], rnms[diffmat[ , 2]])
