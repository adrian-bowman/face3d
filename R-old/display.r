"display.face3d" <- function (shape, type = "surface", colour = "texture", colour.range,
                              new = TRUE, add = FALSE, subset,
                              windowRect = c(300, 300, 800, 800), ...) {
                              	
    if (!require(rgl)) stop("the rgl package is required.")
                              	
    if (!(is.list(shape) && ("coords" %in% names(shape)) && ("triples" %in% names(shape))))
        stop("this is not a face3d object.")
    if ("lit" %in% names(list(...)))
        stop("the lit argument is set by default.")
    if (!missing(subset)) {
        if (length(subset) != nrow(shape$coords))
            stop("subset does not match the shape coords.")
        shape <- subset.face3d(shape, subset)
    }
    if (add)
        new <- FALSE
    if (new) {
        open3d(windowRect = windowRect)
        view3d()
    }
    
    missing.colour.range <- missing(colour.range)
    
    if (!is.na(colour[1]) && colour[1] == "texture") {
        if ("colour" %in% names(shape))
            colour <- shape$colour
        else colour <- "grey"
    }
    
    if (all(colour == "shape index")) {
        if ("shape.index" %in% names(shape)) {
            # clr <- rgb(c(rep(0, 3), 0.5, rep(1, 5)), c(rep(1, 
                # 7), 0.5, 0), c(0, 0.5, rep(1, 3), 0.5, rep(0, 
                # 3)))
            # clr.brks <- c(-1, seq(-7, 7, by = 2)/8, 1)
            # colour <- as.character(cut(shape$shape.index[, colname], 
                # breaks = clr.brks, labels = clr))
            si      <- shape$shape.index
            clr.r   <- c(rep(0, 3), 0.5, rep(1, 5))
            clr.g   <- c(rep(1, 7), 0.5, 0)
            clr.b   <- c(0, 0.5, rep(1, 3), 0.5, rep(0, 3))
            clr.pts <- seq(-1, 1, by = 0.25)
            colour  <- rep("grey", length(shape$shape.index))
            ind     <- !is.na(si)
            colour[ind] <- rgb(approx(clr.pts, clr.r, xout = si[ind])$y,
                               approx(clr.pts, clr.g, xout = si[ind])$y,
                               approx(clr.pts, clr.b, xout = si[ind])$y)
        }
        else {
            colour <- "grey"
            cat("Information on shape index is not present - reverting to grey.\n")
        }
    }
    
    if ((length(colour) == 1) && (substr(colour, 1, 6) == "normal")) {
       if (!("normals" %in% names(shape))) {
          colour <- "grey"
          cat("Normals are not present - reverting to grey.\n")
       }
       else {
          clr <- switch(colour,
                       "normal-x" = abs(shape$normals[ , 1]),
                       "normal-y" = shape$normals[ , 2],
                       "normal-z" = shape$normals[ , 3])
          if (colour == "normal-x")
          	 brks <- seq(0, 1, length = 21)
          else
       	     brks <- seq(-1, 1, length = 21)
       	  colour <- topo.colors(20)[cut(clr, brks, labels = FALSE)]
       }
    }
    
    if (is.numeric(colour)) {
        if (missing.colour.range) colour.range <- range(colour, na.rm = TRUE)
        if (missing.colour.range & (colour.range[1] >= -1) & (colour.range[2] <= 1)) {
           clr.r       <- c(rep(0, 3), 0.5, rep(1, 5))
           clr.g       <- c(rep(1, 7), 0.5, 0)
           clr.b       <- c(0, 0.5, rep(1, 3), 0.5, rep(0, 3))
           clr.pts     <- seq(-1, 1, by = 0.25)
           ind         <- !is.na(colour)
           colour[ind] <- rgb(approx(clr.pts, clr.r, xout = colour[ind])$y,
                              approx(clr.pts, clr.g, xout = colour[ind])$y,
                              approx(clr.pts, clr.b, xout = colour[ind])$y)
        }
        else {
           clr         <- rep(0, nrow(shape$coords))
           ind         <- !is.na(colour)
           clr[ind]    <- cut(colour[ind], seq(colour.range[1], colour.range[2], length = 20),
                                  labels = FALSE, include.lowest = TRUE)
           colour      <- rep("grey", length(colour))
           colour[ind] <- topo.colors(20)[clr[ind]]
        }
        if (length(colour) == length(shape$triples) / 3)
           colour <- rep(colour, each = 3)
        else if (length(colour) != nrow(shape$coords)) {
           print(length(colour))
           print(nrow(shape$coords))
           colour <- "grey"
           cat("The length of the colour vector is inappropriate - reverting to grey.\n")
        }
     }
     
   single.colour <- (length(colour) == 1)
   if (single.colour) colour <- rep(colour, nrow(shape$coords))
   if ((type == "points") & (length(colour) != nrow(shape$coords)))
   	     stop("mismatch between type and the length of colour.")
   if (type %in% c("mesh", "surface")) {
   	  if (length(colour) == nrow(shape$coords)) colour <- colour[shape$triples]
   	  else if (length(colour) != length(shape$triples))
   	     stop("mismatch between type and the length of colour.")
   }
      
   if (!new & !add) ids <- rgl.ids()$id
   if (type == "points")
      points3d(shape$coords, col = colour, ...)
   else if (type == "mesh")
      triangles3d(shape$coords[shape$triples, ], col = colour, 
            front = "line", back = "line", lit = FALSE, ...)
   else if (type == "surface")
      triangles3d(shape$coords[shape$triples, ], col = colour, lit = single.colour, ...)
   if (!new & !add) pop3d(id = ids)
   # for (i in 1:length(ids)) pop3d(id = ids)
        
    invisible()

}

plot.face3d <- display.face3d
