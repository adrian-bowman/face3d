#   Display of 3d shape pca using rgl.

plotpca.face3d <- function(gpa, npc = 6, connect = list(),
        group = NA, control = NA, title = '3D PCA display', rgl.size = 2) {

rotate <- function(panel) {
      with(panel, {
         if (phi < -90) phi <- -90
         if (phi >  90) phi <-  90
         rgl.viewpoint(theta = theta, phi = phi, fov = fov)
         })
      panel
      }

pc.animate <- function(panel) {
   with(panel,{
      shapes3d.pca.plot(panel)
      cmult.vec <- c(seq(-0.2, -3, by = -0.2), seq(-2.8, 3, by = 0.2), 
                     seq(2.8, 0, by = -0.2))
      pc        <- as.numeric(pc)
      if (!mean.showing)
         shapes3d.plot.case(gpa$mshape, connect, col = "black", size = rgl.size)
      for (cmult in cmult.vec) {
         temp <- inv.tangcoord(v(cmult, pc, gpa), gpa$mshape)
         for (i in 1:length(connect)) pop3d()
         shapes3d.plot.case(temp, connect, col = "black", size = rgl.size)
         if (spin.animate) {
            rgl.viewpoint(theta = theta, phi = phi, fov = fov)
            theta <- theta + 1
            }
         hold(speed)
         }
      for (i in 1:length(connect)) pop3d()
      if (mean.showing)
         shapes3d.plot.case(gpa$mshape, connect, col = "black", size = rgl.size)
      })
   if (panel$spin.animate) panel$theta <- panel$theta + 60
   panel
   }
   
hold <- function(speed) {
   if (speed < 100) {
      tim <- proc.time()[3]
      while (tim + 0.1 * (1 - speed / 100) > proc.time()[3]) proc.time()[3]
      }
   }

shapes3d.pca.plot <- function(panel) {
   with(panel, {
      clear3d()
      shapes3d.plot.case(pc.limits, col = 'white', alpha = 0)
      if (mean.showing)
         shapes3d.plot.case(gpa$mshape, connect, col = "black", size = rgl.size)
      if (control.mean.showing)
         shapes3d.plot.case(control.mean, connect, col = "black", size = rgl.size)
      if (data.showing) {
         shapes3d.plot.case(gpa$rotated[, , insp.case], connect = connect,
            col = group.code[insp.case] + 3, size = rgl.size)
         # text.coords <- c(mean(range(pc.limits[,1])), 
         #                  # max(pc.limits[,2]) + diff(range(pc.limits[,2]))/4,
         #                  max(pc.limits[,2]),
         #                  mean(range(pc.limits[,3])))
         # rgl.texts(text.coords[1], text.coords[2], text.coords[3],
         #    paste('Observation', insp.case), 
         #    col = group.code[insp.case] + 3, justify = 'center')
         if ((closest.showing) & (insp.case %in% case.ind))  {
            closest.ind <- which(insp.case == case.ind)
            if (p.ctrl[closest.ind] < 1)
               shapes3d.plot.case(closest.ctrl[, , closest.ind], 
                  connect = connect, col = 3, size = rgl.size)
            }
         }
      if (extremes.showing) {
         pc <- as.numeric(pc)
         shapes3d.plot.case(lower[, , pc], connect, col = "blue", size = rgl.size)
         shapes3d.plot.case(upper[, , pc], connect, col = "red",  size = rgl.size)
         }
      })
   # if (!panel$scores.showing) panel <- shapes3d.scores.plot(panel)
   panel
   }
   
shapes3d.plot.case <- function(shape, connect = list(), ...) {
   if (length(connect) == 0) connect$points <- 1:nrow(shape)
   cpts   <- connect$points
   clines <- connect$lines
   if (length(cpts) > 0)
      rgl.points(shape[cpts,1], shape[cpts, 2], shape[cpts, 3], ...)
   if (length(clines) > 0)
      rgl.lines(shape[clines, 1], shape[clines, 2], shape[clines, 3], ...)
   }
   
shapes3d.tour <- function(panel) {
      with(panel,{
      shapes3d.pca.plot(panel)
      npos   <- 5
      c.vals <- matrix(rnorm(npos * 6, sd = 1.5), ncol = npos)
      c.vals <- cbind(rep(0, 6), c.vals, rep(0, 6))
      shapes3d.plot.case(gpa$mshape, connect, col = 'darkgreen', 
         size = rgl.size)
      for (i in 1:(npos + 1)) {
      for (j in 1:20) {
         c.temp   <- (j / 20) * c.vals[, i + 1] + ((20 - j) / 20) * c.vals[,i ]
         temp.tan <- as.vector(gpa$pcar[, 1:6] %*% diag(gpa$pcasd[1:6]) %*% 
                          c.temp)
         temp.tan <- apply(gpa$tan, 1, mean) + temp.tan
         temp     <- inv.tangcoord(temp.tan, gpa$mshape)
         for (l in 1:length(connect)) pop3d()
         shapes3d.plot.case(temp, connect, col = 2, size = rgl.size)
         if (spin.animate) {
            rgl.viewpoint(theta = theta, phi = phi, fov = fov)
            theta <- theta + 1
            }
         hold(speed)
         }}
      for (l in 1:length(connect)) rgl.pop()
      })
   panel
   }
  
shapes3d.scores.plot <- function(panel) {
   with(panel, {
      if (!scores.showing) dev.new()       # x11(width = 3.7, height = 4)
      if (all(group.code == 1))
         plot(gpa$scores[, 1] * gpa$pcasd[1], gpa$scores[, 2] * gpa$pcasd[2],
             xlab = "Score 1", ylab = "Score 2")
      else {
         plot(gpa$scores[, 1] * gpa$pcasd[1], gpa$scores[, 2] * gpa$pcasd[2],
             xlab = "Score 1", ylab = "Score 2", 
             pch = group.code, col = group.code + 3)
         }
      })
   panel$scores.showing <- TRUE
   panel
   } 

shapes3d.show.data <- function(panel) {
   # if (!panel$scores.showing) x11(width = 3.7, height = 4)
   panel <- shapes3d.scores.plot(panel)
   with(panel, {
      if (data.showing) 
         text(gpa$scores[insp.case, 1] * gpa$pcasd[1], 
              gpa$scores[insp.case, 2] * gpa$pcasd[2], insp.case)
      shapes3d.pca.plot(panel)
      })
   panel
   }

shapes3d.inspect <- function(panel) {
   if (!panel$scores.showing) dev.new()      # x11(width = 3.7, height = 4)
   panel         <- shapes3d.scores.plot(panel)
   with(panel, {
      text(gpa$scores[insp.case, 1] * gpa$pcasd[1], 
           gpa$scores[insp.case, 2] * gpa$pcasd[2], insp.case)

      shapes3d.pca.plot(panel)
      })
   panel
   }

shapes3d.scores.identify <- function(panel) {
   if (!panel$scores.showing) {
      dev.new()      # x11(width = 3.7, height = 4)
      panel <- shapes3d.scores.plot(panel)
      }
   ind <- identify(panel$gpa$scores[, 1] * panel$gpa$pcasd[1], 
                   panel$gpa$scores[, 2] * panel$gpa$pcasd[2], n = 1)
   if (length(ind) > 0) {
      panel$insp.case <- ind[1]
      shapes3d.pca.plot(panel)
      }
   panel
   } 

closest.control <- function(gpa, npc, group, ctrl = 'Control') {

     ctrl.ind  <- which(group == ctrl)
     n.ctrls   <- length(ctrl.ind)
     case.ind  <- which(group != ctrl)
     n.cases   <- length(case.ind)

     x     <- t(gpa$tan - apply(gpa$tan, 1, mean)) %*% gpa$pcar[, 1:npc]
     x.c   <- apply(x[ctrl.ind, ], 2, mean)
     Sigma <- cov(x[ctrl.ind, ])

     p.y                       <- rep(NA, n.cases)
     closest.control.scores    <- array(NA, dim = c(n.cases, npc))
     closest.control.tan.coord <- array(NA, dim = dim(gpa$tan[, case.ind]))
     closest.ctrls             <- array(NA, dim = dim(gpa$rotated[, , case.ind]))

     for (i in 1:n.cases) {
         yi   <- x[case.ind[i], ]
         p.yi <- sqrt(qchisq(0.95, npc) / mahalanobis(yi, x.c, Sigma))
         if (p.yi < 1) {
            closest.control.scores[i,]    <- p.yi * yi + (1 - p.yi) * x.c
            closest.control.tan.coord[,i] <- apply(gpa$tan, 1, mean)
            for (j in 1:npc) closest.control.tan.coord[,i] <-
               closest.control.tan.coord[,i] + 
               closest.control.scores[i,j] * gpa$pcar[,j]
            p.y[i] <- p.yi
            closest.ctrls[, , i] <-
               inv.tangcoord(closest.control.tan.coord[, i], gpa$mshape)
            }
         }
     control.mean.tan <- apply(gpa$tan, 1, mean)
     for (j in 1:npc) control.mean.tan <- control.mean.tan +
                              x.c[j] * gpa$pcar[, j]
     control.mean     <- inv.tangcoord(control.mean.tan, gpa$mshape)
     invisible(list(closest.ctrl = closest.ctrls, p.ctrl = p.y,
                    case.ind = case.ind, control.mean = control.mean))
     }

v <- function(c, j, gpapca) apply(gpapca$tan, 1, mean) + 
                                 c * gpapca$pcasd[j] * gpapca$pcar[, j]

inv.tangcoord <- function(A, mu, approximate = TRUE) {
   # A is a 'configuration' in tangent coordinates km*1
   # mu is the mean around which the original tangent coordinates
   # were calculated.
   # Usually, A will correponds to the effect of one of the PC of
   # variation Eqn 5.25 pp 96, i.e. A = v(c, j).
   # The definition of inv.TC corresponds to Eqn 4.36 on pp37 with
   # identity rotation matrix (i.e. no rotation)
   #### inv.TC <- sqrt(1 - A %*% A) * c(mu) + A
   
   m      <- dim(mu)[2]
   if (approximate) inv.tc <- matrix(A, ncol = m) + mu
   else {
      inv.tc <- sqrt(centroid.size(mu)^2 - A %*% A) * c(mu) + A
      inv.tc <- matrix(inv.tc, ncol = m)
      }
   invisible(inv.tc)
   }


   if (gpa$m != 3) stop('three-dimensional data were expected.')
   npc    <- max(npc, 2)
   pcmult <- 3
   
   pc.limits <- rbind(apply(gpa$rotated, 2, min), apply(gpa$rotated, 2, max))
   for (i in 1:npc) {
      pc.limits <- rbind(pc.limits, inv.tangcoord(v( pcmult, i, gpa), gpa$mshape))
      pc.limits <- rbind(pc.limits, inv.tangcoord(v(-pcmult, i, gpa), gpa$mshape))
      }
   pc.limits <- rbind(apply(pc.limits, 2, max), apply(pc.limits, 2, min))
   
   lower <- array(0, dim = c(gpa$k, 3, npc))
   upper <- lower
   for (i in 1:npc) {
      lower[, , i] <- inv.tangcoord(v(-3, i, gpa), gpa$mshape)
      upper[, , i] <- inv.tangcoord(v( 3, i, gpa), gpa$mshape)
      }
      
   if (!all(is.na(group))) {
      group      <- as.factor(group)
      group.code <- match(group, levels(group))
      if (!is.na(control)) {
         closest    <- closest.control(gpa, npc, group, control)
         ctrl.code  <- which(levels(group) == control)
         if (ctrl.code != 1) {
            ind.ctrl  <- which(group.code == ctrl.code)
            ind.black <- which(group.code == 1)
            group.code[ind.ctrl]  <- 1
            group.code[ind.black] <- ctrl.code
            }
         }
      else
         closest    <- list(closest.ctrl = NA, p.ctrl = NA, case.ind = NA, control.mean = NA)
      }
   else {
      group.code <- rep(1, dim(gpa$rotated)[3])
      closest    <- list(closest.ctrl = NA, p.ctrl = NA, case.ind = NA, control.mean = NA)
      }
   
   cpts <- which(!(1:gpa$k %in% unlist(connect)))
   if (length(connect) > 0) {
      ind <- numeric(0)
      for (i in 1:length(connect)) {
         ci  <- connect[[i]]
         ni  <- length(ci)
         ind <- c(ind, rep(ci, rep(2, ni))[-c(1, 2 * ni)])
         }
      connect <- list(lines = ind)
      }
   if (length(cpts) > 0) connect$points <- cpts

   pht <- if ((!all(is.na(group))) & (!is.na(control))) 365 else 330
   pca.display <- rp.control(title, 
                     theta = -30, phi = 30, pc.limits = pc.limits, fov = 1,
                     prev.pc = "Mean", prev.insp = 0, gpa = gpa, npc = npc,
                     connect = connect, scores.showing = FALSE,
                     lower = lower, upper = upper, pcmult = pcmult,
                     control = control, closest.ctrl = closest$closest.ctrl,
                     p.ctrl = closest$p.ctrl, case.ind = closest$case.ind,
                     group.code = group.code, control.mean.showing = FALSE,
                     closest.showing = FALSE, rgl.size = rgl.size,
                     control.mean = closest$control.mean)
   rp.checkbox(pca.display, mean.showing, shapes3d.pca.plot,
                     initval = TRUE, title = "Show mean")
   rp.checkbox(pca.display, data.showing, shapes3d.show.data,
                     title = "Show data")
   rp.doublebutton(pca.display, insp.case, 1, title = "",
                     range = c(1, dim(gpa$rotated)[3]), initval = 1,
                     action = shapes3d.inspect)
   if ((!all(is.na(group))) & (!is.na(control))) {
      rp.checkbox(pca.display, control.mean.showing, 
                     shapes3d.pca.plot,
                     title = "Control mean")
      rp.checkbox(pca.display, closest.showing, shapes3d.pca.plot,
                     title = "Closest control")
      }
   rp.button(pca.display, action = shapes3d.scores.plot, "Scores plot")
   rp.button(pca.display, action = shapes3d.scores.identify, "Identify points")
   # print(paste(1:npc, ' (', round(gpa$percent[1:npc]), '%)', sep = ''))
   rp.radiogroup(pca.display, pc, 1:npc,
         labels = paste(1:npc, ' (', round(gpa$percent[1:npc]), '%)', sep = ''),
         action = shapes3d.pca.plot, title = "Principal components")
   rp.button(pca.display, action = pc.animate, "Animate")
   rp.checkbox(pca.display, spin.animate, title = 'Spin animation')
   rp.checkbox(pca.display, extremes.showing, shapes3d.pca.plot, title = 'Show extremes')
   rp.button(pca.display, action = shapes3d.tour, "Tour")
   rp.slider(pca.display, speed, 0, 100, initval = 100, title = "Animation speed")
   rp.doublebutton(pca.display, theta, -2, "Theta", action = rotate)
   rp.doublebutton(pca.display, phi, -2,  "Phi", action = rotate)

   rgl.open()
   rgl.bg(color = c('white', 'black'))
   rgl.viewpoint(0, 0, fov = 1)
   rp.do(pca.display, shapes3d.pca.plot)
      
   invisible(pca.display)

   }
