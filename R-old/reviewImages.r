reviewImages.face3d <- function(path = ".", pattern = ".jpg", recursive = TRUE, hscale = 1, vscale = 0.8) {
	
   if (!require(jpeg))   stop("the jpeg package is required.")
   if (!require(rpanel)) stop("the rpanel package is required.")
	
   files <- list.files(path = path, pattern = pattern, recursive = recursive, full.names = TRUE)
   fls   <- list.files(path = path, pattern = pattern, recursive = recursive)
   if (length(files) == 0) stop("no files found")
   
   plot.image <- function(panel) {
   	  with(panel, {
         picture <- readJPEG(files[i])
         par(mar = c(1, 1, 2, 1))
         plot(c(1, ncol(picture)), c(1, nrow(picture)), type = "n", xlab = "", ylab = "",
               axes = FALSE)
         rasterImage(picture, 1, 1, ncol(picture),nrow(picture))
         title(fls[i])
      })
      panel
   }
   
   replot.image <- function(panel) {
      rp.tkrreplot(panel, imageplot)
      panel
   }
   
   panel <- rp.control(files = files, i = 1)
   rp.tkrplot(panel, imageplot, plot.image, pos = "right", hscale = hscale, vscale = vscale)
   rp.doublebutton(panel, i, 1, "image", replot.image)

   invisible(files)
}
