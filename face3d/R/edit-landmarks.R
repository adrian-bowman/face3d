editlandmarks.face3d <- function(shape, lmk.names = rownames(shape$landmarks), panel = TRUE,
                                    lmk.name = "none", directions = FALSE) {
   
   if (!requireNamespace("rpanel", quitely = TRUE)) stop("the rpanel package is not available.")
   
   radius         <- 0.5
   delta          <- 0.2
   angle          <- 0
   if (is.null(lmk.names)) lmk.names <- paste("lmk", 1:3, sep = "")
   id.lmks        <- integer(length = length(lmk.names))
   names(id.lmks) <- lmk.names
   env            <- new.env()
   
   rotate <- function(panel) {
   	  if (panel$lmk.name == "none" | !panel$path.showing) return(panel)
   	  shape <- if (panel$zoom) panel$sbst else panel$shape
   	  # if (panel$drns.showing) {
   	  	 # if (!("directions" %in% names(shape))) shape <- curvatures(shape, directions = TRUE)
   	     # drns             <- array(NA, dim = c(3, 3, dim(shape$directions)[3]))
         # drns[ , 1:2, ]   <- shape$directions
         # drns[ ,   3, ]   <- t(shape$normals)
         # shape$directions <- drns
   	  # }
      path  <- planepath.face3d(shape, panel$shape$landmarks[panel$lmk.name, ],
                                direction = c(sin(panel$angle), cos(panel$angle), 0),
                                rotation = 0, monitor = FALSE, boundary = c(Inf, Inf),
                                directions = panel$drns.showing)
      panel$arclength    <- path$arclength
      panel$path         <- path$path
      panel$x1.arclength <- path$x1.arclength
      panel$shape        <- path$shape
      id.line.old        <- panel$id.line
      panel$id.line      <- rgl::lines3d(panel$path, col = "green", lwd = 2)
      if (!is.na(id.line.old) && (id.line.old %in% rgl::rgl.ids()$id)) rgl::pop3d(id = id.line.old)
      if (panel$drns.showing) {
         panel$directions <- path$directions[ , 1:3, ]
         panel$kappa1     <- path$kappa1
         panel$kappa2     <- path$kappa2
         panel            <- show.drns(panel)
      }
      panel
   }
   
   along <- function(panel) {
   	  if (panel$lmk.name == "none" | !panel$path.showing) return(panel)
      panel$x1.arclength      <- min(max(panel$arclength), panel$x1.arclength)
      lmk                     <- c(approx(panel$arclength, panel$path[ , 1], panel$x1.arclength)$y,
                                   approx(panel$arclength, panel$path[ , 2], panel$x1.arclength)$y,
                                   approx(panel$arclength, panel$path[ , 3], panel$x1.arclength)$y)
      ind                     <- panel$lmk.name
      panel$shape$landmarks[ind, ] <- lmk
      old                     <- panel$id.lmks[ind]
      panel$id.lmks[ind]      <- rgl::spheres3d(t(lmk), radius = panel$radius, col = "red")
      rgl::pop3d(id = old)
      if (panel$drns.showing) panel <- show.drns(panel)
      panel
   }
   
   show.drns <- function(panel) {
      drn1          <- c(approx(panel$arclength, panel$directions[1, 1, ], panel$x1.arclength)$y,
                         approx(panel$arclength, panel$directions[2, 1, ], panel$x1.arclength)$y,
                         approx(panel$arclength, panel$directions[3, 1, ], panel$x1.arclength)$y)
      drn2          <- c(approx(panel$arclength, panel$directions[1, 2, ], panel$x1.arclength)$y,
                         approx(panel$arclength, panel$directions[2, 2, ], panel$x1.arclength)$y,
                         approx(panel$arclength, panel$directions[3, 2, ], panel$x1.arclength)$y)
      drn3          <- c(approx(panel$arclength, panel$directions[1, 3, ], panel$x1.arclength)$y,
                         approx(panel$arclength, panel$directions[2, 3, ], panel$x1.arclength)$y,
                         approx(panel$arclength, panel$directions[3, 3, ], panel$x1.arclength)$y)
      kappa1        <-   approx(panel$arclength, panel$kappa1,             panel$x1.arclength)$y
      kappa2        <-   approx(panel$arclength, panel$kappa2,             panel$x1.arclength)$y
      id.drn1.old   <- panel$id.drn1
      id.drn2.old   <- panel$id.drn2
      id.drn3.old   <- panel$id.drn3
      x             <- panel$shape$landmarks[panel$lmk.name, ]
      xgrid         <- seq(-panel$shape$si.distance, panel$shape$si.distance, length = 20)
      zgrid         <- sweep(xgrid %o% drn1 + kappa1 * xgrid^2 %o% drn3, 2, x, "+")
      panel$id.drn1 <- rgl::lines3d(zgrid, col = "blue")
      zgrid         <- sweep(xgrid %o% drn2 + kappa2 * xgrid^2 %o% drn3, 2, x, "+")
      panel$id.drn2 <- rgl::lines3d(zgrid, col = "red")
      # panel$id.drn1 <- lines3d(rbind(x, x + 10 * drn1), col = "blue")
      # panel$id.drn2 <- lines3d(rbind(x, x + 10 * drn2), col = "red")
      panel$id.drn3 <- rgl::lines3d(rbind(x, x + 10 * drn3), col = "yellow")
      if (!is.na(id.drn1.old) && (id.drn1.old %in% rgl::rgl.ids()$id)) rgl::pop3d(id = id.drn1.old)
      if (!is.na(id.drn2.old) && (id.drn2.old %in% rgl::rgl.ids()$id)) rgl::pop3d(id = id.drn2.old)
      if (!is.na(id.drn3.old) && (id.drn3.old %in% rgl::rgl.ids()$id)) rgl::pop3d(id = id.drn3.old)
      panel
   }
   
   new.lmk <- function(panel) {
      if (panel$lmk.name.old != "none") {
         rgl::pop3d(id = panel$id.lmks[panel$lmk.name.old])
         panel$id.lmks[panel$lmk.name.old] <- rgl::spheres3d(t(panel$shape$landmarks[panel$lmk.name.old, ]),
                                                        radius = panel$radius, col = "green")
      }
      if (panel$lmk.name != "none") {
      	 if (!(panel$lmk.name %in% rownames(panel$shape$landmarks))) {
      	 	origin <- apply(summary.face3d(panel$shape, print = FALSE)$ranges, 2, mean)
      	   umat   <- solve(rgl::par3d("userMatrix")[1:3, 1:3])
            crds   <- sweep(panel$shape$vertices, 2, origin) %*% umat
            ind    <- order(edist.face3d(crds[ , 1:2], rep(0, 2)))[1:10]
            ind1   <- which.max(crds[ind, 3])
            panel$shape$landmarks <- rbind(panel$shape$landmarks, panel$shape$vertices[ind[ind1], ])
            rownames(panel$shape$landmarks)[nrow(panel$shape$landmarks)] <- panel$lmk.name
      	 }
         else
            rgl::pop3d(id = panel$id.lmks[panel$lmk.name])
         panel$id.lmks[panel$lmk.name] <- rgl::spheres3d(t(panel$shape$landmarks[panel$lmk.name, ]),
                                                    radius = panel$radius, col = "red")
         panel <- rotate(panel)
      }
      else {
         rgl::pop3d(id = panel$id.line)
         panel$id.line <- NA
      }
      panel$lmk.name.old               <- panel$lmk.name
      panel
   }
   
   save.shape <- function(panel) {
      assign("shape", panel$shape, pos = panel$env)
      panel
   }
   
   display.shape <- function(panel) {
   	  shape <- if (panel$zoom) panel$sbst else panel$shape
   	  if (panel$colour == "shape index") shape <- curvatures(shape)
      plot(shape, colour = panel$colour, display = panel$type, new = panel$first)
      panel$first <- FALSE
      if (("landmarks" %in% names(panel$shape)) & !panel$zoom)  {
         for (i in 1:nrow(panel$shape$landmarks))
            if (all(!is.na(panel$shape$landmarks[i, ])))
               panel$id.lmks[i] <- rgl::spheres3d(t(panel$shape$landmarks[i, ]), col = "green", radius = panel$radius)
         if (!(panel$lmk.name == "none")) rgl::pop3d(id = panel$id.lmks[panel$lmk.name])
      }
      if (!(panel$lmk.name == "none")) {
         panel$id.lmks[panel$lmk.name] <- rgl::spheres3d(t(panel$shape$landmarks[panel$lmk.name, ]),
                                                    radius = panel$radius, col = "red")
         panel <- rotate(panel)
      }
      panel
   }
   
   zoom <- function(panel) {
   	  if (panel$lmk.name == "none") {
   	     rpanel::rp.messagebox("zooming requires a landmark to be specified.", title = "Face3D message")
   	     return(panel)
   	  }
   	  panel$zoom <- !panel$zoom
   	  if (panel$zoom) {
   	     lmk <- panel$shape$landmarks[panel$lmk.name, ]
         panel$sbst <- subset(panel$shape, edist.face3d(shape$vertices, lmk) < 20)
      }
      panel$id.line <- NA
      panel <- display.shape(panel)
      panel
   }
   
   perpend <- function(panel) {
   	  panel$perpendicular <- !panel$perpendicular
   	  if (panel$perpendicular)
   	     panel$angle <- panel$angle + pi / 2
   	  else
   	     panel$angle <- panel$angle - pi / 2
      panel <- rotate(panel)
      panel
   }
   
   show.path <- function(panel) {
      if (!panel$path.showing) {
         if (!is.na(panel$id.line)) rgl::pop3d(id = panel$id.line)
         panel$id.line <- NA
      }
      else
         panel <- rotate(panel)
      panel
   }
   
   if (panel) {
      panel <- rpanel::rp.control(shape = shape, angle = angle, id.lmks = id.lmks, id.line = NA, first = TRUE,
                          id.drn1 = NA, id.drn2 = NA, id.drn3 = NA,
                          delta = delta, env = env, lmk.name = lmk.name, lmk.name.old = lmk.name,
                          path = NA, radius = radius, zoom = FALSE, colour = "texture", type = "surface",
                          path.showing = FALSE, drns.showing = FALSE, perpendicular = FALSE)
      rpanel::rp.do(panel, display.shape)
      rpanel::rp.combo(panel, lmk.name, "Landmark name", c("none", lmk.names), action = new.lmk)
      rpanel::rp.button(panel, zoom, "Zoom in/out")
      rpanel::rp.checkbox(panel, path.showing, show.path, "Show path")
      rpanel::rp.checkbox(panel, drns.showing, show.path, "Show directions")
      rpanel::rp.doublebutton(panel, angle,        pi/48, "direction", action = rotate)
      rpanel::rp.doublebutton(panel, x1.arclength, delta, "distance",  action = along)
      rpanel::rp.button(panel, perpend, "Perpendicular direction")
      rpanel::rp.radiogroup(panel, type,   c("points", "mesh", "surface"), action = display.shape)
      rpanel::rp.radiogroup(panel, colour, c("texture", "shape index"), action = display.shape)
      rpanel::rp.button(panel, save.shape, "Save shape")
      rpanel::rp.block(panel)
      if (exists("shape", env, inherits = FALSE))
         shape <- get("shape", env, inherits = FALSE)
      else
         cat("The shape has not been saved and the original shape has been returned.")
   }
   else {
      panel <- list(shape = shape, angle = angle, id.lmks = id.lmks, id.line = NA, first = TRUE,
                    id.drn1 = NA, id.drn2 = NA, id.drn3 = NA,
                    delta = delta, env = env, lmk.name = lmk.name, lmk.name.old = lmk.name,
                    path = NA, radius = radius, zoom = FALSE, colour = "texture", type = "surface",
                    path.showing = TRUE, drns.showing = directions, perpendicular = FALSE)
      panel <- display.shape(panel)
      rotate(panel)
   }

   invisible(shape) 
}
