# Create a movie of an animation from one shape to another

animate.face3d <- function(shapes, ngrid = 4, image.folder, ...) {
   
   if (!is.list(shapes)) stop("shapes is not a list.")
   if (!all(sapply(shapes, is.face3d)))
      stop("some elements of shapes are not face3d objects.")
   
   snap <- !missing(image.folder)
   if (snap && !(image.folder %in% list.files()))
      stop(paste("the image folder (", image.folder, ") does not exist.", sep = ""))

   i <- 1
   plot(shapes[[1]], ...)
   if (snap) snapshot3d(paste(image.folder ,"/temp-", i, ".png", sep = ""))

   for (k in 1:(length(shapes) - 1)) {
      shape12 <- shapes[[k]]
      for (j in 1:ngrid) {
         beta <- j / ngrid
         shape12$coords <- beta * shapes[[k + 1]]$coords + (1 - beta) * shapes[[k]]$coords
         sv <- par3d(skipRedraw = TRUE)
         plot(shape12, ...)
         par3d(sv)
         i <- i + 1
         if (snap) snapshot3d(paste(image.folder ,"/temp-", i, ".png", sep = ""))
      }
   }
   
   invisible()
}

#    system(paste("ffmpeg -start_number 1 -i ", storage.folder, "/temp-", 
#                 "%d.png -vframes ", i,
#                 " -c:v libx264 -pix_fmt yuv420p -vf 'scale=trunc(iw/2)*2:trunc(ih/2)*2' -y ",
#                 output, sep = ""))


# system(paste("ffmpeg -start_number 1 -i temp/temp-", 
#              "%d.png -vframes ", 4 * ngrid,
#              " -c:v libx264 -pix_fmt yuv420p -vf 'scale=trunc(iw/2)*2:trunc(ih/2)*2' -y ",
#              "animation.mp4", sep = ""))
