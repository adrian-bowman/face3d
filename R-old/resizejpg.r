resizejpg.face3d <- function(filename, recursive = TRUE, xsize = 1024, ysize = 1024,
                                        file.extension = ".jpg", file.addition = "-small",
                                        overwrite = FALSE, verbose = 1) {

#     Function to create .jpg files of reduced size

   path <- NULL
   if (missing(filename))
      path <- "."
   else if (length(filename) > 1)
      files <- filename
   else if (substr(filename, nchar(filename) + 1 - nchar(file.extension), nchar(filename)) == 
                          file.extension)
      files <- filename
   else
      path  <- filename

   if (!is.null(path)) {
      files <- list.files(path, recursive = recursive, full.names = TRUE)
      ind   <- substr(files, nchar(files) + 1 - nchar(file.extension), nchar(files)) == 
                          file.extension
      files <- files[ind]
   }
   if (length(files) == 0) {
     cat("No files found.\n")
     return(invisible())
   }

   for (i in 1:length(files)) {
      end.name <- substr(files[i],
   	                     nchar(files[i]) + 1 - nchar(file.extension) - nchar(file.addition),
   	                     nchar(files[i]))
      smallname <- paste(substr(files[i], 1, nchar(files[i]) - nchar(file.extension)),
                         file.addition, file.extension, sep = "")
   	  flag <- TRUE
   	  if (end.name == paste(file.addition, file.extension, sep = "")) {
   	     if (verbose > 1)
   	        cat("\nThis file has the ending", end.name, "already:", files[i], "\n")
   	     flag <- FALSE
   	  }
   	  else if (smallname %in% files) {
   	  	 print("here")
   	     if (verbose > 1)
   	        cat("\nThe file", smallname, "already exists.\n")
   	     flag <- FALSE
   	  }
   	  if (overwrite) flag <- TRUE
   	  if (flag) {    
         system(paste("convert -geometry ", as.character(xsize), "x", as.character(ysize),
                            " ", files[i], " ", smallname, sep = ""))
   	     if (verbose > 0) cat(files[i], "\n")
   	  }
   }

   invisible()
}
