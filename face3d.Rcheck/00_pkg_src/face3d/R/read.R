"read.face3d" <- function (filename, colour = TRUE, jpgfile, convert.gray = FALSE,
                           quality = 100, monitor = TRUE)  {

    nc <- nchar(filename)

    #-----------------------------------------------------------------
    #                            obj files
    #-----------------------------------------------------------------
  
    if (substr(filename, nc - 3, nc) == ".obj") {
        coords  <- readLines(filename, warn = FALSE)
        coordsv <- coords[substr(coords, 1, 2) == "v "]
        coordsv <- strsplit(coordsv, " ")
        coordsv <- lapply(coordsv, function(x) x[x != ""])
        ncrdv   <- sapply(coordsv, length)
        if (diff(range(ncrdv)) != 0) stop("the rows for vertices have different lengths.")
        ncrdv   <- ncrdv[1]
        coordsv <- t(matrix(unlist(coordsv), nrow = length(coordsv[[1]])))
        vcol    <- NULL
        if (ncrdv == 7) {
           vcol <- cbind(as.numeric(coordsv[ , 5]), as.numeric(coordsv[ , 6]), as.numeric(coordsv[ , 7]))
           vcol <- if (all(vcol >= 0 & vcol <= 1)) vcol <- rgb(vcol[ , 1], vcol[ , 2], vcol[ , 3]) else NA
        }
        coordsv <- cbind(as.numeric(coordsv[ , 2]), as.numeric(coordsv[ , 3]), as.numeric(coordsv[ , 4]))
        coordsf <- coords[substr(coords, 1, 2) == "f "]
        if (length(coordsf) == 0)
          stop("there are no faces present.\n")
        coordsf <- strsplit(coordsf, " ")
        if (length(grep("/", unlist(coordsf))) == 0) {
            coordsf <- t(matrix(unlist(coordsf), nrow = 4))
            coordsf <- cbind(as.numeric(coordsf[ , 2]), as.numeric(coordsf[ , 3]), as.numeric(coordsf[ , 4]))
            coordsf <- c(t(coordsf))
            result  <- list(coords = coordsv, triples = coordsf)
            #return(invisible(result))
        }
        else {
            n.edges <- sapply(coordsf, length) - 1
            if (any(n.edges > 5)) {
                print(table(n.edges))
                stop("some faces have more than 5 edges.")
            }
        if (length(table(n.edges)) == 1) {
            coordsf  <- coords[substr(coords, 1, 2) == "f "]
            coordsf  <- strsplit(coordsf, " ")
            coordsf  <- unlist(coordsf)
            coordsf  <- matrix(coordsf, ncol = 4, byrow = TRUE)
            coordsf  <- c(t(coordsf[, 2:4]))
            coordsf  <- strsplit(coordsf, "/")
            coordsfv <- as.numeric(sapply(coordsf, function(x) x[1]))
            coordsfv <- as.numeric(coordsfv)
            if (length(coordsf[[1]]) > 1)
               coordsft <- as.numeric(sapply(coordsf, function(x) x[2]))
            result <- list(coords = coordsv, triples = coordsfv)
        }
        if (length(table(n.edges)) == 2) {
            coordsf <- coords[substr(coords, 1, 2) == "f "]
            coordsf.three <- coordsf[n.edges == 3]
            coordsf.three <- strsplit(coordsf.three, " ")
            coordsf.three <- unlist(coordsf.three)
            coordsf.three <- t(matrix(unlist(coordsf.three), nrow = 4))
            coordsf.three <- c(t(coordsf.three[, 2:4]))
            coordsf.three <- strsplit(coordsf.three, "/")
            coordsfv.three <- c(matrix(unlist(coordsf.three), nrow = 9)[c(1, 4, 7), ])
            coordsfv.three <- as.numeric(coordsfv.three)
            coordsft.three <- c(matrix(unlist(coordsf.three), nrow = 9)[c(2, 5, 8), ])
            coordsft.three <- as.numeric(coordsft.three)
            coordsf.four <- coordsf[n.edges == 4]
            coordsf.four <- strsplit(coordsf.four, " ")
            coordsf.four <- unlist(coordsf.four)
            coordsf.four <- t(matrix(unlist(coordsf.four), nrow = 5))
            coordsf.four <- c(t(coordsf.four[, 2:5]))
            coordsf.four <- strsplit(coordsf.four, "/")
            coordsfv.four <- c(matrix(unlist(coordsf.four), nrow = 12)[c(1, 4, 7, 10), ])
            coordsfv.four <- as.numeric(coordsfv.four)
            coordsft.four <- c(matrix(unlist(coordsf.four), nrow = 12)[c(2, 5, 8, 11), ])
            coordsft.four <- as.numeric(coordsft.four)
            coordsfv.four <- matrix(coordsfv.four, ncol = 4, byrow = TRUE)
            coordsfv.four <- c(c(t(coordsfv.four[, 1:3])), c(t(coordsfv.four[ , c(1, 3, 4)])))
            coordsfv <- c(coordsfv.three, coordsfv.four)
            coordsfv <- as.numeric(coordsfv)
            coordsft.four <- matrix(coordsft.four, ncol = 4, byrow = TRUE)
            coordsft.four <- c(c(t(coordsft.four[, 1:3])), c(t(coordsft.four[ , c(1, 3, 4)])))
            coordsft <- c(coordsft.three, coordsft.four)
            coordsft <- as.numeric(coordsft)
            result <- list(coords = coordsv, triples = coordsfv)
        }
        if (length(table(n.edges)) == 3) {
            coordsf <- coords[substr(coords, 1, 2) == "f "]
            coordsf.three <- coordsf[n.edges == 3]
            coordsf.three <- strsplit(coordsf.three, " ")
            coordsf.three <- unlist(coordsf.three)
            coordsf.three <- t(matrix(unlist(coordsf.three),
                nrow = 4))
            coordsf.three <- c(t(coordsf.three[, 2:4]))
            coordsf.three <- strsplit(coordsf.three, "/")
            coordsfv.three <- c(matrix(unlist(coordsf.three),
                nrow = 9)[c(1, 4, 7), ])
            coordsfv.three <- as.numeric(coordsfv.three)
            coordsft.three <- c(matrix(unlist(coordsf.three),
                nrow = 9)[c(2, 5, 8), ])
            coordsft.three <- as.numeric(coordsft.three)
            coordsf.four <- coordsf[n.edges == 4]
            coordsf.four <- strsplit(coordsf.four, " ")
            coordsf.four <- unlist(coordsf.four)
            coordsf.four <- t(matrix(unlist(coordsf.four), nrow = 5))
            coordsf.four <- c(t(coordsf.four[, 2:5]))
            coordsf.four <- strsplit(coordsf.four, "/")
            coordsfv.four <- c(matrix(unlist(coordsf.four), nrow = 12)[c(1,
                4, 7, 10), ])
            coordsfv.four <- as.numeric(coordsfv.four)
            coordsft.four <- c(matrix(unlist(coordsf.four), nrow = 12)[c(2,
                5, 8, 11), ])
            coordsft.four <- as.numeric(coordsft.four)
            coordsfv.four <- matrix(coordsfv.four, ncol = 4,
                byrow = TRUE)
            coordsfv.four <- c(c(t(coordsfv.four[, 1:3])), c(t(coordsfv.four[,
                c(1, 3, 4)])))
            coordsf.five <- coordsf[n.edges == 5]
            coordsf.five <- strsplit(coordsf.five, " ")
            coordsf.five <- unlist(coordsf.five)
            coordsf.five <- t(matrix(unlist(coordsf.five), nrow = 6))
            coordsf.five <- c(t(coordsf.five[, 2:6]))
            coordsf.five <- strsplit(coordsf.five, "/")
            coordsfv.five <- c(matrix(unlist(coordsf.five), nrow = 15)[c(1,
                4, 7, 10, 13), ])
            coordsfv.five <- as.numeric(coordsfv.five)
            coordsft.five <- c(matrix(unlist(coordsf.five), nrow = 15)[c(2,
                5, 8, 11, 14), ])
            coordsft.five <- as.numeric(coordsft.five)
            coordsfv.five <- matrix(coordsfv.five, ncol = 5,
                byrow = TRUE)
            coordsfv.five <- c(c(t(coordsfv.five[, 1:3])), c(t(coordsfv.five[,
                c(1, 3, 4)])),c(t(coordsfv.five[,c(1, 4, 5)])))
            coordsfv <- c(coordsfv.three, coordsfv.four, coordsfv.five)
            coordsfv <- as.numeric(coordsfv)

            coordsft.five <- matrix(coordsft.five, ncol = 5,
                byrow = TRUE)
            coordsft.five <- c(c(t(coordsft.five[, 1:3])), c(t(coordsft.five[,
                c(1, 3, 4)])),c(t(coordsft.five[,c(1, 4, 5)])))
            coordsft <- c(coordsft.three, coordsft.four, coordsft.five)
            coordsft <- as.numeric(coordsft)
            result <- list(coords = coordsv, triples = coordsfv)
        }
        }
        colour.present <- any(substr(coords, 1, 3) == "vt ")
        if (colour && requireNamespace("jpeg", quietly = TRUE) & colour.present) {
            coordst <- coords[substr(coords, 1, 3) == "vt "]
            coordst <- strsplit(coordst, " ")
            coordst <- t(matrix(unlist(coordst), nrow = 3))
            coordst <- cbind(as.numeric(coordst[, 2]), as.numeric(coordst[, 3]))
            coordst1 <- rep(NA, nrow(coordsv))
            coordsfv <- as.numeric(coordsfv)
            coordst1[coordsfv] <- coordsft
            coordst1 <- as.numeric(coordst1)
            coordst2 <- coordst[coordst1, ]

            if (missing(jpgfile))
               jpgfile <- paste(substr(filename, 1, nchar(filename) - 4), ".jpg", sep = "")
            image.mat <- try(jpeg::readJPEG(jpgfile), silent = TRUE)
            if (class(image.mat) != "try-error") {
            	m <- dim(image.mat)[1]
            	n <- dim(image.mat)[2]
            	IC1 <- 1 + round((1 - coordst2[, 2]) * (m - 1))
            	IC2 <- 1 + round((coordst2[, 1]) * (n - 1))
            	r <- image.mat[cbind(IC1, IC2, 1)]
            	g <- image.mat[cbind(IC1, IC2, 2)]
            	b <- image.mat[cbind(IC1, IC2, 3)]
            	ind <- which(is.na(r + g + b))
            	r[ind] <- 0.5
            	g[ind] <- 0.5
            	b[ind] <- 0.5
            	clr <- rgb(r, g, b)
               result$colour <- clr
            }
        }
        if (colour && !colour.present && !is.null(vcol) && length(vcol) > 1) result$colour <- vcol
        class(result) <- "face3d"
    }
    
    #-----------------------------------------------------------------
    #                           ply files
    #-----------------------------------------------------------------
    
    else if (substr(filename, nc - 3, nc) == ".ply") {
       
       # Open file, identify format and read header
       fl <- file(filename, "rb")
       # open(fl)
       content <- readLines(fl, 2, warn = FALSE)
       if (content[1] != "ply") stop("this is not a ply file.")
       fmt <- content[2]
       fmt <- strsplit(fmt, " ")[[1]]
       fmt <- fmt[fmt != ""]
       fmt <- fmt[2]
       if (fmt != "ascii") fmt <- strsplit(fmt, "_")[[1]][2]
       n <- 2
       while (content[n] != "end_header") {
          n       <- n + 1
          content <- c(content, readLines(fl, 1, warn = FALSE))
       }

       # Identify vertex and face elements and their properties
       elements <- character(0)
       els      <- grep("element", content)
       for (i in els) {
          el <- strsplit(content[i], " ")[[1]]
          el <- el[el != ""]
          elements <- c(elements, el[2])
       }
       ind      <- match(c("vertex", "face"), elements)
       if (any(is.na(ind))) stop("vertex and face elements are not both present.")
       props    <- list()
       type     <- list()
       n        <- rep(NA, 2)
       names(n) <- elements
       for (el in elements) {
          ind   <- grep(paste("element", el), content)
          # if (length(ind) == 0) stop("an", el, "element is not present.")
          m     <- strsplit(content[ind], " ")[[1]]
          m     <- m[m != ""]
          n[el] <- as.integer(m[3])
          prps  <- character(0)
          i     <- ind + 1
          while (substr(content[i], 1, 8) == "property") {
             prps <- c(prps, content[i])
             i    <- i + 1
          }
          prps  <- strsplit(prps, " ")
          prps  <- lapply(prps, function(x) x[x != ""])
          typ   <- unlist(lapply(prps, function(x) x[2]))
          prps  <- unlist(lapply(prps, function(x) x[3]))
          props[[el]] <- prps
          type[[el]]  <- typ
       }
       result  <- list()
       
       # Read the content
       if (fmt == "ascii") {
          if (monitor)
             cat("Reading", n["vertex"], "vertex co-ordinates and",
                            n["face"],   "triangle indices", "...")
          rl1     <- readLines(fl, n[1], warn = FALSE)
          rl2     <- readLines(fl, n[2], warn = FALSE)
          coords  <- if (names(n)[1] == "vertex") rl1 else rl2
          triples <- if (names(n)[1] == "face")   rl1 else rl2
          coords  <- strsplit(coords, " ")
          coords  <- lapply(coords, function(x) x[x != ""])
          ncrds   <- sapply(coords, length)
          if (diff(range(ncrds)) != 0) stop("the rows for vertices have different lengths.")
          ncrds   <- ncrds[1]
          coords  <- t(matrix(as.numeric(unlist(coords)), nrow = ncrds))
          colnames(coords) <- props$vertex
          if (all(c("red", "green", "blue") %in% props$vertex))
             colour <- coords[, c("red", "green", "blue")]
          if (!all(c("x", "y", "z") %in% props$vertex))
             stop("x, y, z not all present.")
          result$coords <- coords[, c("x", "y", "z")]
          triples <- strsplit(triples, " ")
          triples <- lapply(triples, function(x) x[x != ""])
          nvals   <- sapply(triples, function(x) x[1])
          if (!(all(nvals == "3"))) stop("not all faces are triangles.")
          triples <- t(matrix(as.numeric(unlist(triples)), nrow = 4))
          result$triples <- c(t(triples[ , 2:4])) + 1
          if (monitor) cat(" completed.\n")
       }
       else {
          coords <- matrix(nrow = n["vertex"], ncol = 3,
                           dimnames = list(NULL, c("x", "y", "z")))
          colour <- matrix(nrow = n["vertex"], ncol = 3,
                           dimnames = list(NULL, c("red", "green", "blue")))
          triples <- integer(3 * n["face"])
          for (el in elements) {
             if (monitor) {
                txt <- if (el == "vertex") "vertex co-ordinates" else "triangle indices"
                cat("Reading", n[el], txt, "...")
             }
             # thcount <- 0
             for (i in 1:n[el]) {
               if (el == "vertex") {
                  for (j in 1:length(props[[el]])) {
                     typ  <- type[[el]][j]
                     sgnd <- if (substr(typ, 1, 1) == "u") FALSE else TRUE
                     siz  <- if (typ == "uchar") 1 else NA_integer_
                     if (typ == "uchar") typ <- "integer"
                     val <- readBin(fl, typ, 1, signed = sgnd, size = siz, endian = fmt)
                     vbl <- props[[el]][j]
                     if (vbl %in% c("x", "y", "z")) coords[i, vbl] <- val
                     if (vbl %in% c("red", "green", "blue")) colour[i, vbl] <- val
                  }
               }
               if (el == "face") {
                  m <- readBin(fl, "integer", 1, size = 1, signed = FALSE, endian = fmt)
                  if (m != 3) stop("not all faces are triangles.")
                  triples[(i - 1) *3 + 1:3] <-
                     readBin(fl, "integer", 3, endian = fmt)
               }
               # if (i %/% 1000 > thcount) {
               #    thcount <- thcount + 1
               #    cat(thcount, "")
               # }
             }
             if (monitor) cat(" completed.\n")
          }
          result$coords  <- coords
          result$triples <- triples + 1
       }
       
       if (all(c("red", "green", "blue") %in% props$vertex)) {
             mx <- if (max(colour) > 1) 255 else 1
             result$colour <- rgb(colour[ , "red"], colour[ , "green"],
                                  colour[ , "blue"], maxColorValue = mx)
       }
       class(result)  <- "face3d"

       close(fl)
    }
    
    #-----------------------------------------------------------------
    #                             pts files
    #-----------------------------------------------------------------

    else if (substr(filename, nc - 3, nc) == ".pts") {
        coords <- readLines(filename)
        coordss <- coords[substr(coords, 1, 1) == "S"]
        coordss <- strsplit(coordss, " ")
        coordssnew <- list()
        k1 <- length(coordss)
        for (j in 1:k1) coordssnew[[j]] <- coordss[[j]][coordss[[j]] != 
            ""][-1]
        lmks <- matrix(as.numeric(unlist(coordssnew)), k1, 3, 
            byrow = TRUE)
        result <- list(lmks = lmks)
        curves.present <- any(substr(coords, 1, 1) == "C")
        if (curves.present) {
            coordss <- coords[substr(coords, 1, 1) == "C"]
            coordss <- strsplit(coordss, " ")
            coordssnew <- list()
            k2 <- length(coordss)
            for (j in 1:k2) coordssnew[[j]] <- coordss[[j]][coordss[[j]] != 
                ""][-1]
            curves <- matrix(as.numeric(unlist(coordssnew)), 
                k2, 3, byrow = TRUE)
            result$curves <- curves
        }
    }
    
    #-----------------------------------------------------------------
    #                            dilm files
    #-----------------------------------------------------------------
 
    else if (substr(filename, nc - 4, nc) == ".dilm") {
        coords <- readLines(filename, warn = FALSE)
        k <- length(coords) - 1
        coords <- coords[4:k]
        coordss <- strsplit(coords, " ")
        coordssnew <- list()
        k1 <- length(4:k)
        lmks <- matrix(0, k1, 3)
        for (j in 1:k1) coordssnew[[j]] <- coordss[[j]][coordss[[j]] != 
            ""][-1]
        for (j in 1:k1) {
            n <- nchar(coordssnew[[j]][1])
            lmks[j, 2] <- as.numeric(strsplit(coordssnew[[j]][1], 
                "\"")[[1]][2])
            lmks[j, 3] <- as.numeric(strsplit(coordssnew[[j]][3], 
                "\"")[[1]][2])
            lmks[j, 1] <- as.numeric(strsplit(coordssnew[[j]][9], 
                "\"")[[1]][2])
        }
        result <- list(lmks = lmks)
    }
    

    #-----------------------------------------------------------------
    #                             tps files
    #-----------------------------------------------------------------
    
    else if (substr(filename, nc - 3, nc) == ".tps") {
        coords <- readLines(filename, warn = FALSE)
        coords <- coords[coords != ""]
        k1 <- length(coords)
        k <- as.numeric(substr(coords[1], 4, 5))
        ID.lm <- which(substr(coords, 1, 3) == "LM=")
        ID.id <- which(substr(coords, 1, 3) == "ID=")
        ID.image <- which(substr(coords, 1, 6) == "IMAGE=")
        index <- matrix(coords[ID.id],length(ID.id),1)
        index <- as.numeric(apply(index,1,id.string))
        images <- matrix(coords[ID.image],length(ID.image),1)
        images <- apply(images,1,image.string)
        ID.out <- c(ID.lm,ID.id,ID.image)
        k2 <- k1 - length(ID.out)
        coords <- matrix(coords[-ID.out],k2,1)
        coords <- matrix(as.numeric(apply(coords,1,split.string)),2,k2)
        n <- k2/k
        coords <- array(coords, c(2, k, n))
        coords <- aperm(coords, c(2, 1, 3))
        result <- list()
        result$coords <- coords
        result$images <- images
        result$index <- index
        cat("Number of landmarks:", k, "\n")
        cat("Number of dimensions:", dim(coords)[2], "\n")
        cat("Number of images:", n, "\n")
    }

    #-----------------------------------------------------------------
    #                        landmarkAscii files
    #-----------------------------------------------------------------
    
    else if (substr(filename, nc - 13, nc) == ".landmarkAscii") {
            coords <- readLines(filename, warn = FALSE)
            coords <- coords[coords != ""]
            k <- length(coords)
            k1 <- length(10:k)
            coords <- matrix(coords[10:k],k1,1)
            coords <- matrix(as.numeric(apply(coords,1,split.string)),k1,3, byrow=TRUE)
            cat("Number of landmarks:", k1, "\n")
            result <- list()
            result$coords <- coords
    }

    #-----------------------------------------------------------------
    #                              jpg files
    #-----------------------------------------------------------------

    else if (substr(filename, nc - 3, nc) == ".jpg") {
            if (!requireNamespace("jpeg", quietly = TRUE)) stop("the jpeg package is not available.")
            jpgfile <- jpeg::readJPEG(filename)
            jpgfile <- jpeg::imagematrix(jpgfile, type="rgb")
            if (convert.gray == FALSE) jpgfile <- jpeg::imagematrix(jpgfile,type="rgb")
            if (convert.gray == TRUE)  jpgfile <- jpeg::imagematrix(jpgfile,type="grey")
            result <- list()
            result$jpgfile <- jpgfile
    }

    # else if (substr(filename, nc - 3, nc) == ".tif") {
            # if (!require(rtiff)) stop("the rtiff package is not available.")
            # tiffile <- readTiff(filename, reduce = 1 - quality/100)
            # if (convert.gray){
            # colgrey <- (tiffile@red + tiffile@green + tiffile@blue)/3
            # tiffile <- newPixmapRGB(red = colgrey, green = colgrey, blue = colgrey)
            # }
            # result <- list()
            # result$tiffile <- tiffile
     # }
    
    else stop("file extension not recognised.")
        
    invisible(result)
}

"sub.string"   <- function(x,nch) substr(x, nch - 3, nch)
"image.string" <- function(x) {
                    nch <- nchar(x)
                    jpgfile <- substr(x, 7, nch)
                    jpgfile
}
"id.string"    <- function(x) {
                    nch <- nchar(x)
                    jpgfile <- substr(x, 4, nch)
                    jpgfile
}
"split.string" <- function(x) strsplit(x, " ")[[1]]
