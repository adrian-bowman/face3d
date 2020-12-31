"edist.face3d"  <- function(x, x0) sqrt(rowSums(sweep(x, 2, x0)^2))

if(getRversion() >= "2.15.1")
   utils::globalVariables(c("use.colours", "i", "imageplot", "imagematrix",
               "phi", "theta", "speed", "extremes.showing", "spin.animate",
               "pc", "closest.showing", "control.mean.showing",
               "data.showing", "mean.showing", "get.lam", "sparseMatrix",
               "tsearch", "delaunayn", "principal.curve", "colour", "type",
               "x1.arclength", "drns.showing", "path.showing", "template",
               "rdist", "result", "pcurve", "insp.case"))
