plot.face3d <- function (shape, display = c("surface", "lines"), plot.isolated = TRUE,
                              colour = "texture", palette, breaks, scaling,
                              range.colour, key = FALSE, clabel = NULL, cex.axis = 1,
                              new = FALSE, add = FALSE, subset, vector.scale = 1,
                              windowRect = c(300, 300, 800, 800), theta, phi, zoom,
                              ...) {
                              	
   if (!requireNamespace("rgl", quietly = TRUE)) stop("the rgl package is required.")
                        
   breaks.missing  <- missing(breaks)
   palette.missing <- missing(palette)
   if (!breaks.missing && !palette.missing && length(breaks) != length(palette) + 1)
      stop("length(breaks) is not length(palette) + 1.")
   if (missing(breaks)) {
      breaks  <- NA
      nbreaks <- if (!palette.missing && length(palette) != 1) length(palette) + 1 else 50
   }
   else if (length(breaks) == 1) {
      nbreaks <- breaks
      breaks  <- NA
   }
   else
      nbreaks <- length(breaks)
   palette.missing      <- missing(palette)
   col.palette          <- if (palette.missing) NA else palette
   scaling.missing  <- missing(scaling)
   range.colour.missing <- missing(range.colour)

   if (!(is.list(shape) && ("vertices" %in% names(shape)) && ("triangles" %in% names(shape))))
		stop("this is not a face3d object.")

   pars <- list(...)
   if ("lit" %in% names(pars))
      stop("the lit argument is set by default.")

   if (!missing(subset)) {
      # if (!is.logical(subset) & !)
      # if (is.logical(subset) & (length(subset) != nrow(shape$vertices)))
      #    stop("the length of subset does not match the shape vertices.")
      # if (is.integer(subset) & all(subset >= 1) & all(subset <= nrow(shape$vertices)))
         shape <- subset.face3d(shape, subset)
      # else
      #    stop("subset contains inappropriate values.")
   }

   if (length(rgl::rgl.dev.list()) == 0) new <- TRUE
   else if (add) new <- FALSE
   if (!new & !add) ids <- rgl::rgl.ids()$id
   if (new) rgl::open3d(windowRect = windowRect)
   if (!add) {
      if (missing(theta)) theta <- 0
      if (missing(phi))   phi   <- 0
      if (missing(zoom))  zoom  <- 0.7
      rgl::view3d(theta, phi, zoom  = zoom)
   }
   
   if (all(colour == "shape.index")) colour <- "shape index"
   if (length(colour) == 1 && is.character(colour) &&
       !(substr(colour, 1, 7) %in% c("texture", "normal-", "shape i"))) {
      if (colour %in% names(shape)) colour <- shape[[colour]]
   }
  
   if ("normal" %in% display) {
      if (!("normals" %in% names(shape))) shape <- normals.face3d(shape)
      clr  <- if (length(colour) == 1 && is.character(colour) &&
                  !(substr(colour, 1, 7) %in% c("texture", "normal-"))) colour else "black"
      crds <- cbind(shape$vertices, shape$vertices + vector.scale * shape$normals)
      crds <- matrix(c(t(crds)), ncol = 3, byrow = TRUE)
      rgl::segments3d(crds, col = clr)
   }
  
   if ("principal" %in% substr(display, 1, 9)) {
      if (!("directions" %in% names(shape)))
         stop("principal directions not present.  These can be computed by index.face3d with the argument 'directions = TRUE'.")
      clr  <- if (length(colour) == 1 && is.character(colour) &&
                  !(substr(colour, 1, 7) %in% c("texture", "normal-"))) colour else "black"
      display1 <- display[match("principal", substr(display, 1, 9))]
      drn  <- substr(display1, nchar(display1), nchar(display1))
      if (!(drn %in% c("1", "2"))) stop("the final character of display should be '1' or '2'.")
      drn  <- as.numeric(drn)
      crds <- cbind(shape$vertices, shape$vertices + vector.scale * t(shape$directions[ , drn, ]))
      crds <- matrix(c(t(crds)), ncol = 3, byrow = TRUE)
      rgl::segments3d(crds, col = clr)
   }

    if (!is.na(colour[1]) && colour[1] == "texture") {
       if ("colour" %in% names(shape))
          colour <- shape$colour
        else
          colour <- "grey"
    }

    if (length(colour) == 1 && !is.na(colour) && colour == "shape index") {
        if ("shape.index" %in% names(shape)) {
            colour      <- sicolour.face3d(shape$shape.index)
            breaks      <- seq(-1, 1, by = 0.05)
            scaling     <- "linear"
            brks.r      <- c(rep(0, 3), 0.5, rep(1, 5))
            brks.g      <- c(rep(1, 7), 0.5, 0)
            brks.b      <- c(0, 0.5, rep(1, 3), 0.5, rep(0, 3))
            brks        <- seq(-1, 1, by = 0.25)
            col.palette <- rgb(approx(brks, brks.r, xout = breaks[-1] - 0.025)$y,
                               approx(brks, brks.g, xout = breaks[-1] - 0.025)$y,
                               approx(brks, brks.b, xout = breaks[-1] - 0.025)$y)
            # col.palette <- rgb(clr.r, clr.g, clr.b)
            # breaks      <- c(-1, seq(-7, 7, by = 2)/8, 1)
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
          col.palette <- topo.colors(20)
       	 colour <- col.palette[cut(clr, brks, labels = FALSE)]
       }
    }

    if (is.numeric(colour) & length(colour) > 1) {
       if (range.colour.missing)
          range.colour <- range(colour, na.rm = TRUE)
       if (diff(range.colour) < sqrt(.Machine$double.eps)) {
          warning("the colour values are all identical.")
          range.colour <- mean(colour) + c(-1, 1)
       }
       bxp <- boxplot(colour, plot = FALSE)
       lmt <- c(bxp$stats[c(1, 5), 1])
       ex1 <- length(which(colour < lmt[1])) > 0
       ex2 <- length(which(colour > lmt[2])) > 0
       if (scaling.missing) {
          scaling <- if (range.colour.missing & (ex1 | ex2)) "trimmed" else "linear"
       }
       if (all(is.na(breaks))) {
          if (prod(range.colour) < 0 & all(is.na(col.palette))) {
             range.colour <- c(-1, 1) * max(abs(range.colour))
             col.palette  <- "diverge"
          }
          if (scaling == "quantile") {
             if (!all(is.na(col.palette)) && substr(col.palette[1], 1, 7) == "diverge") {
                q <- quantile(abs(colour), seq(0, 1, length = (nbreaks %/% 2) + 1), na.rm = TRUE)
                q <- q[-1]
                breaks <- unique(c(rev(-q), 0, q))
             }
             else
                breaks <- unique(quantile(colour, seq(0, 1, length = nbreaks), na.rm = TRUE))
          }
          else if (scaling == "trimmed") {
             if (!all(is.na(col.palette)) && substr(col.palette[1], 1, 7) == "diverge") {
                almt   <- max(abs(lmt))
                breaks <- c(min(colour),
                            seq(-almt, almt, length = nbreaks - 2),
                            max(colour))
             }
             else
               breaks <- c(min(colour),
                           seq(lmt[1], lmt[2], length = nbreaks - 2),
                           max(colour))
          }
          else
             breaks <- seq(range.colour[1], range.colour[2], length = nbreaks)
       }
       nbreaks         <- length(breaks)
       breaks[1]       <- breaks[1] - 0.001 * diff(range(breaks))
       breaks[nbreaks] <- breaks[nbreaks] + 0.001 * diff(range(breaks))

       if (all(is.na(col.palette)))
          col.palette <- topo.colors(nbreaks - 1)
       else if (col.palette[1] == "diverge") {
          if (requireNamespace("colorspace", quietly = TRUE)) 
             col.palette <- colorspace::diverge_hcl(nbreaks - 1) 
          else
             # col.palette  <- topo.colors(nbreaks - 1)
             col.palette  <- diverge.topo(nbreaks - 1)
       }
       else if (col.palette[1] == "diverge_topo")
          col.palette  <- diverge.topo(nbreaks - 1)
 
       clr         <- rep(0, nrow(shape$vertices))
       ind         <- !is.na(colour)
       clr[ind]    <- cut(colour[ind], breaks, labels = FALSE, include.lowest = TRUE)
       colour      <- rep("white", length(colour))
       colour[ind] <- col.palette[clr[ind]]
       if (length(colour) == length(c(shape$triangles)) / 3)
          colour <- rep(colour, each = 3)
       else if (length(colour) != nrow(shape$vertices)) {
          colour <- "grey"
          cat("The length of the colour vector is inappropriate - reverting to grey.\n")
       }
    }

   if (key) colourkey.face3d(col.palette, breaks, scaling,
                             clabel = clabel, cex.axis = cex.axis)

   single.colour <- (length(table(colour)) == 1)
   if (single.colour & (length(colour) == 1)) colour <- rep(colour, nrow(shape$vertices))
   if (("points" %in% display) & (length(colour) != nrow(shape$vertices)))
   	  stop("mismatch between display and the length of colour.")
   if (any(c("mesh", "surface") %in% display)) {
   	  if (length(colour) == nrow(shape$vertices)) colour <- colour[c(t(shape$triangles))]
   	  else if (length(colour) != length(c(shape$triangles)))
   	     stop("mismatch between display and the length of colour.")
   }
   
   if ("points" %in% display)
      rgl::points3d(shape$vertices, col = colour, ...)
   if ("spheres" %in% display)
      rgl::spheres3d(shape$vertices, col = colour, ...)
   if ("mesh" %in% display)
      rgl::triangles3d(shape$vertices[c(t(shape$triangles)), ], col = colour,
            front = "line", back = "line", lit = FALSE, ...)
   if ("surface" %in% display)
      rgl::triangles3d(shape$vertices[c(t(shape$triangles)), ], col = colour, lit = single.colour, ...)
   if (("lines" %in% display) & ("lines" %in% names(shape)))
      rgl::segments3d(shape$vertices[c(t(shape$lines)), ], ...)
   # Plot isolated points
   if (!("points" %in% display) & any(display %in% c("mesh", "surface")) & plot.isolated) {
      ind <- which(!(1:nrow(shape$vertices) %in% c(shape$triangles)))
      if (("lines" %in% display) & ("lines" %in% names(shape)))
         ind <- ind[!(ind %in% c(shape$lines))]
      rgl::spheres3d(shape$vertices[ind, ], ...)
   }

   if (!new & !add) rgl::pop3d(id = ids)

    invisible(list(palette = col.palette, breaks = breaks))

}

sicolour.face3d <- function(x) {
   outside  <- (!is.na(x) && abs(x) > 1)
   if (any(outside)) stop("some values are outside the range -1 to 1.")
   brks.r   <- c(rep(0, 3), 0.5, rep(1, 5))
   brks.g   <- c(rep(1, 7), 0.5, 0)
   brks.b   <- c(0, 0.5, rep(1, 3), 0.5, rep(0, 3))
   brks     <- seq(-1, 1, by = 0.25)
   clr      <- rep("white", length(x))
   ind      <- which(!is.na(x) & !outside)
   clr[ind] <- rgb(approx(brks, brks.r, xout = x[ind])$y,
                   approx(brks, brks.g, xout = x[ind])$y,
                   approx(brks, brks.b, xout = x[ind])$y)
   clr
}

sicolor.face3d <- function(x)
   sicolour.face3d(x)

diverge.topo <- function(n, alpha = 1) {
   n   <- n %/% 2
   neg <- hsv(h = seq.int(from = 43/60, to = 31/60, length.out = n), alpha = alpha)
   j   <- n %/% 3
   k   <- n %/% 3
   i   <- n - j - k
   pos <-  c(hsv(h = seq.int(from = 23/60, to = 11/60, length.out = i), alpha = alpha), 
             hsv(h = seq.int(from = 10/60, to = 6/60, length.out = (j + k)), alpha = alpha,
                 s = seq.int(from = 1, to = 0.3, length.out = (j + k)), v = 1))
   c(neg, pos)
}

colourkey.face3d <- function(cols, breaks, scaling, cex.axis = 1, par.mar = c(2, 1, 2, 3) + 0.1,
                             margin = FALSE, clabel = NULL, rgl = FALSE)  {

   par.old <- par()[c("mgp", "mar", "tcl")]
   ngrid   <- length(cols)
   xvec    <- rep(0, ngrid)
   if (length(breaks) == 2)
      breaks <- seq(breaks[1], breaks[2], length = ngrid + 1)
   else if (length(breaks) != ngrid + 1)
      stop("inappropriate length of breaks in colourkey.")
   if (scaling == "linear") {
      zlim      <- range(breaks)
      yaxs      <- "r"
   }
   else if (scaling == "trimmed") {
      dbrks     <- breaks[3] - breaks[2]
      zlim      <- c(breaks[2] - dbrks, breaks[length(breaks) - 1] + dbrks)
      yaxs      <- "r"
   }
   else {
      zlim      <- c(0, ngrid)
      brks.orig <- breaks
      breaks    <- 0:ngrid
      yaxs      <- "i"
   }

   if (rgl) {
      xlim <- diff(zlim) / 10
      crds <- cbind(cbind(0,    breaks[1:ngrid], 0), cbind(0,    breaks[-1],      0),
                    cbind(xlim, breaks[-1],      xlim), cbind(xlim, breaks[1:ngrid], xlim))
      crds <- matrix(c(t(crds)), ncol = 3, byrow = TRUE)
      clr  <- rep(cols, each = 4)
      view3d(135, 0)
      rgl.material(lit = FALSE)
      quads3d(crds, col = clr)
      # axes3d(col = "black")
      axis3d('y', pos = c(0, NA, 0), col = "black")
      mtext3d("Something", 'y++', pos = c(0.5, 0, zlim[2]), col = "black")
      par3d(mouseMode = rep("none", 4))
   }
   else {
      par(mar = par.mar, mgp = c(1.5, 0.2 + 0.5 * (cex.axis - 1), 0), tcl = -0.2)
      xrange <- if (margin) c(-1, 1) else c(0, 1)
      plot(xrange, zlim, type = "n", axes = FALSE, xaxs = "i", yaxs = yaxs, xlab = " ", ylab = " ")
      if (scaling == "linear")
         axis(4, cex.axis = cex.axis)
      else if (scaling == "trimmed") {
         ticks <- pretty(breaks[2:(length(breaks) - 1)])
         axis(4, at = ticks, labels = TRUE, cex.axis = cex.axis)
      }
      else {
         ticks <- pretty(0:ngrid)
         lbls  <- as.character(signif(brks.orig, 4))[match(ticks, 0:ngrid)]
         lbls[lbls == "Inf"] <- NA
         axis(4, at = ticks, labels = lbls, cex.axis = cex.axis)
      }
      nbrks <- length(breaks)
      breaks[c(1, nbrks)] <- par()$usr[3:4]
      rect(xvec, breaks[-nbrks], xvec + 1, breaks[-1], col = cols, border = NA)
      lines(c(0, 0, 1, 1, 0), breaks[c(1, nbrks, nbrks, 1, 1)])
      if (scaling == "trimmed") {
         # abline(h =      breaks[1:2],       col = "white")
         # abline(h =      rev(breaks)[1:2],  col = "white")
         segments(xrange[1], breaks[1], xrange[2], breaks[2], col = "white", lwd = 2)
         segments(xrange[1], breaks[2], xrange[2], breaks[1], col = "white", lwd = 2)
         segments(xrange[1], rev(breaks)[1], xrange[2], rev(breaks)[2], col = "white", lwd = 2)
         segments(xrange[2], rev(breaks)[1], xrange[1], rev(breaks)[2], col = "white", lwd = 2)
      }
      if (!is.null(clabel)) mtext(clabel, 4, 1.5, cex = cex.axis)
   }
   par(par.old)
   invisible()  
}
