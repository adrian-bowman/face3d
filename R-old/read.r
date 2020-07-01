"read.face3d" <- function (filename, colour = TRUE, jpgfile.addition = "-small", convert.gray = FALSE,
                           quality = 100)  {

    nc <- nchar(filename)
    if (substr(filename, nc - 3, nc) == ".obj") {
        coords <- readLines(filename,warn = FALSE)
        coordsv <- coords[substr(coords, 1, 2) == "v "]
        coordsv <- strsplit(coordsv, " ")
        ########################################################################
        if(length(coordsv[[1]]) == 5) coordsv <- lapply(coordsv,function(x) x[x != ""])
        ########################################################################
        coordsv <- t(matrix(unlist(coordsv), nrow = 4))
        coordsv <- cbind(as.numeric(coordsv[ , 2]), as.numeric(coordsv[ , 3]), as.numeric(coordsv[ , 4]))
        coordsf <- coords[substr(coords, 1, 2) == "f "]
        coordsf <- strsplit(coordsf, " ")
        if (length(grep("/", unlist(coordsf))) == 0) {
            coordsf <- t(matrix(unlist(coordsf), nrow = 4))
            coordsf <- cbind(as.numeric(coordsf[ , 2]), as.numeric(coordsf[ , 3]), as.numeric(coordsf[ , 4]))
            coordsf <- c(t(coordsf))
            result <- list(coords = coordsv, triples = coordsf)
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
        colour.present <- any(substr(coords, 1, 3) == "vt ")
        if (colour && require(jpeg) & colour.present) {
            coordst <- coords[substr(coords, 1, 3) == "vt "]
            coordst <- strsplit(coordst, " ")
            coordst <- t(matrix(unlist(coordst), nrow = 3))
            coordst <- cbind(as.numeric(coordst[, 2]), as.numeric(coordst[, 3]))
            coordst1 <- rep(NA, nrow(coordsv))
            coordsfv <- as.numeric(coordsfv)
            coordst1[coordsfv] <- coordsft
            coordst1 <- as.numeric(coordst1)
            coordst2 <- coordst[coordst1, ]
            jpgfile <- paste(substr(filename, 1, nchar(filename) - 4), jpgfile.addition, ".jpg", sep = "")
            image.mat <- readJPEG(jpgfile)
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
            class(result) <- "face3d"
        } 
    }
    }   
    
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
   else if (substr(filename, nc - 3, nc) == ".jpg") {
            if (!require(jpeg)) stop("the jpeg package is not available.")
            jpgfile <- readJPEG(filename)
            jpgfile <- imagematrix(jpgfile,type="rgb")
            if (convert.gray == FALSE) jpgfile <- imagematrix(jpgfile,type="rgb")
            if (convert.gray == TRUE)  jpgfile <- imagematrix(jpgfile,type="grey")
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

"sub.string" <- function(x,nch) substr(x, nch - 3, nch)
"image.string" <- function(x) {
                    nch <- nchar(x)
                    jpgfile <- substr(x, 7, nch)
                    jpgfile
}
"id.string" <- function(x) {
                    nch <- nchar(x)
                    jpgfile <- substr(x, 4, nch)
                    jpgfile
}
"split.string" <- function(x) strsplit(x, " ")[[1]]

