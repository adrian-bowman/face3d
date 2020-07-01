edist.face3d  <- function(x, x0)
   sqrt(rowSums(sweep(x, 2, x0)^2))

arclength.face3d <- function(curve)
   cumsum(c(0, apply(diff(curve)^2, 1, function(x) sqrt(sum(x)))))

projectline.face3d <- function(shape, x1, x2, normal) {
   unit  <- (x2 - x1) / sqrt(sum((x2 - x1)^2))
   proj  <- c(sweep(shape$coords, 2, x1) %*% unit)
   dst   <- apply(outer(proj, unit) - sweep(shape$coords, 2, x1), 1,
                          function(x) sqrt(sum(x^2)))
   if (!missing(normal)) {
      unity <- crossproduct(unit, normal)
      projy <- c(sweep(shape$coords, 2, x1) %*% unity)
      dst   <- dst * sign(projy)
   }
   invisible(list(x = proj, y = dst))
}

is.face3d <- function(object, report = FALSE) {
   ind <- (class(object) == "face3d")
   if (!ind & report) cat("the object does not have class face3d.\n")
   ind <- ind && is.list(object)
   if (!ind & report) cat("the object is not a list.\n")
   ind <- ind && ("coords" %in% names(object))
   if (!ind & report) cat("the object does not have a coords component.\n")
   ind <- ind && ("triples" %in% names(object))
   if (!ind & report) cat("the object does not have a triples component.\n")
   ind <- ind && is.matrix(object$coords)
   if (!ind & report) cat("the coords component is not a matrix.\n")
   ind <- ind && (dim(object$coords)[2] == 3)
   if (!ind & report) cat("the coords component does not have three columns.\n")
   ind
}

area.face3d <- function(shape) {
   n              <- length(shape$triples) / 3
   triangles      <- array(t(shape$coords[shape$triples,]), c(3, 3, n))
   triangle.areas <- apply(triangles, 3, function(x) {
      1/2 * sqrt(sum((crossproduct(x[,2] - x[,1], x[,3] - x[,1], scale = FALSE))^2))})
   if (is.null(rownames(shape$coords))) stop("'shape' must have rownames.")
   triangles  <- matrix(shape$triples, ncol = 3, byrow = TRUE)
   a.wts      <- tapply(rep(triangle.areas / 3, 3), c(triangles), sum)
   wts        <- numeric(nrow(shape$coords))
   ind.a      <- as.numeric(names(a.wts))
   wts[ind.a] <- a.wts
   isol       <- summary(shape, checks = TRUE, print = FALSE)$isolated
   if (!is.null(isol)) {
      l.wts     <- numeric(0)
      crvs      <- c("mid-line columella", "upper face right", "upper face left")
      for (crv in crvs) {
         ind        <- grep(crv, rownames(shape$coords))
         if (grepl("upper face", crv)) ind <- isol[isol %in% ind]
         nind       <- length(ind)
         del        <- apply(diff(shape$coords[ind, ]), 1, function(x) sqrt(sum(x^2)))
         ind1       <- rep(1:nind, each = 2)[-c(1, 2 * nind)]
         new        <- tapply(rep(del / 2, each = 2), ind1, sum)
         names(new) <- as.character(ind)
         l.wts      <- c(l.wts, new)
      }
      wts[as.numeric(names(l.wts))] <- l.wts * sqrt(mean(a.wts))
   }
   invisible(list(area = sum(triangle.areas), triangles = triangle.areas, points = wts))
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Package `face3d', version 0.1-1: type help(face3d) for summary information")
}


if(getRversion() >= "2.15.1")
   utils::globalVariables(c("use.colours", "i", "imageplot", "imagematrix",
               "phi", "theta", "speed", "extremes.showing", "spin.animate",
               "pc", "closest.showing", "control.mean.showing",
               "data.showing", "mean.showing", "get.lam", "sparseMatrix",
               "tsearch", "delaunayn", "principal.curve", "colour", "type",
               "x1.arclength", "drns.showing", "path.showing", "template",
               "rdist", "result", "pcurve", "insp.case"))
