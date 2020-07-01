rp.shapeindex <- function(panel = TRUE, s = -1) {

   if (!require(rpanel)) stop("the rpanel package is not available.")
   if (!require(rgl))    stop("the rgl package is not available.")

   surface.draw <- function(panel) {
      with(panel, {
         tn <- tan(s * pi / 2)
         if (s <= 0) {
            k1 <- 1
            k2 <- -k1 * (1 + tn) / (1 - tn)
         }
         else {
            k2 <- -1
            k1 <- -k2 * (1 - tn) / (1 + tn)
         }
         k  <- c(k1, k2)
         s  <- atan((k[2] + k[1]) / (k[2] - k[1])) * 2 / pi
         z  <- 0.5 * (k[1] * x^2 + k[2] * y^2)
         save <- par3d(skipRedraw=TRUE)
         on.exit(par3d(save))
         if (use.colours) {
            i <- findInterval(s, clr.brks, all.inside = TRUE)
            colr <- clr[i]
         }
         else
            colr <- "green"
         persp3d(x, y, z, col = colr, aspect = "iso",
                   box = FALSE, axes = FALSE,
                   xlab = "", ylab = "", zlab = "", alpha = 1)
         persp3d(x, y, z, col = "black", aspect = "iso",
                   box = FALSE, axes = FALSE,
                   front = "lines", back = "lines", lwd = 2,
                   xlab = "", ylab = "", zlab = "", add = TRUE)
         # text3d(0, 1, 0.5 - 0.4 * (s > 0), "s = 1")
         # desc <- c("spherical cup", "trough", "rut", "saddle rut",
         #           "saddle", "saddle ridge", "ridge", "dome", 
         #           "spherical cap")
         # text3d(0, 0, max(0 - i/5, -1), desc[i])
      })
      panel
   }

   ngrid <- 20
   r     <- matrix(seq(0, 1, length = ngrid), ngrid, ngrid, byrow = TRUE)
   theta <- matrix(seq(0, 2 * pi, length = ngrid), ngrid, ngrid)
   x     <- r * cos(theta)
   y     <- r * sin(theta)
   clr   <- rgb(c(rep(0, 3), 0.5, rep(1, 5)), c(rep(1, 7), 0.5, 0),
                c(0, 0.5, rep(1, 3), 0.5, rep(0, 3)))
   clr.brks <- c(-1, seq(-7, 7, by = 2) / 8, 1)

   open3d()
   
   if (panel) {
      pnl <- rp.control(x = x, y = y, use.colours = TRUE,
                        clr = clr, clr.brks = clr.brks)
      # rp.slider(panel, k, rep(-1, 2), rep(1, 2), surface.draw, showvalue = TRUE)
      rp.slider(pnl, s, -1, 1, surface.draw, showvalue = TRUE,
                labels = "Shape index")
      rp.checkbox(pnl, use.colours, surface.draw, "use colours")
      rp.do(pnl, surface.draw)
   }
   else {
      pnl <- list(x = x, y = y, use.colours = TRUE, clr = clr, clr.brks = clr.brks,
                  s = s)
      surface.draw(pnl)
   }

   invisible()
}
