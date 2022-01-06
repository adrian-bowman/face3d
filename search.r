#     Search for text in all files

target    <- "tsearch"
diry      <- "."
diry      <- "~/research/face3d/face3d"
recursive <- FALSE
recursive <- TRUE

files <- list.files(diry, full.names = TRUE, recursive = recursive)
ind   <- (!grepl(".R", files) & !grep(".r", files))
if(length(ind) > 0) files <- files[!ind]
ind   <-        grep(".Rcheck", files)
ind   <- c(ind, grep(".Rproj",  files))
ind   <- c(ind, grep(".Rda",    files))
ind   <- c(ind, grep(".rda",    files))
ind   <- c(ind, grep("archive", files))
ind   <- c(ind, grep("R-old",   files))
if (length(ind) > 0) files <- files[-ind]

for (ifl in files) {
   # cat(ifl, "\n")
   file <- readLines(ifl)
   grp <- grep(target, file)
   if (length(grp) > 0) {
      cat(ifl, "\n")
      for (jfl in grp) cat(jfl, file[jfl], "\n")
   }
}
