"write.face3d" <- function (filename, face, jpgfile, tiffile, quality = 100) 
{
    nc <- nchar(filename)
    if (substr(filename, nc - 3, nc) == ".ply") {
              plyfile <- paste(substr(filename, 1, nchar(filename) - 4), ".ply", sep = "")
              tripl <- readLines(plyfile)
              tripl <- tripl[substr(tripl, 1, 2) == "3 "]
              tripl <- strsplit(tripl, " ")
              tripl <- t(matrix(unlist(tripl), nrow=4))
              tripl <- cbind(as.numeric(tripl[,1]),as.numeric(tripl[,2]),as.numeric(tripl[,3]),as.numeric(tripl[,4]))
              ## COL-to-RGB tranfromation  
              CLR <- round(t(col2rgb(face$colour, alpha = FALSE)),6)
              ## pulling all together
              ALL <- cbind(face$coords,CLR,rep(255,dim(CLR)[1]))
              ALL1 <- apply(ALL,1,function(x) paste(x, collapse = " "))
              tripl1 <- apply(tripl,1,function(x) paste(x, collapse = " "))
              ## saving final results
              filesavename <- paste(substr(plyfile, 1, nchar(filename) - 4),"-colour",".ply",sep="")
              xx <- file(filesavename,"w")
              cat("ply",
                 "format ascii 1.0", file=xx,sep="\n") 
                 cat("element vertex",dim(ALL)[1],file=xx,sep=" ")
                 cat("","property float x","property float y","property float z",
                 #"property float nx","property float ny","property float nz",
                 "property uchar red","property uchar green","property uchar blue","property uchar alpha",
                 file=xx,sep="\n")
                 cat("element face",dim(tripl)[1],file=xx,sep=" ")
                 cat("","property list uchar int vertex_indices","end_header",file=xx,sep="\n")
                 cat(ALL1,file=xx,sep="\n")
                 cat(tripl1,file=xx,sep="\n")
              close(xx)     
    }                
    else if (substr(filename, nc - 3, nc) == ".dmp") {
            dmpfile <- paste(substr(filename, 1, nchar(filename) - 4), ".dmp", sep = "")
            if ( class(face) == "face3d") save(face, file = dmpfile)
            if ( class(face) == "array") save(face, file = dmpfile)
    }
    else if (substr(filename, nc - 3, nc) == ".jpg") {
             if (quality < 100)  jpgfilename <- paste(substr(filename, 1, nchar(filename) - 4),"-small",".jpg",sep="")
             if (quality == 100) jpgfilename <- paste(substr(filename, 1, nchar(filename) - 4),".jpg",sep="")
             jpeg(jpgfilename,quality = 100 - quality, width = dim(jpgfile)[2], height = dim(jpgfile)[1], units = "px")
             par(mar=c(0,0,0,0),mai = c(0,0,0,0))
             plot(jpgfile)
             dev.off()
    }
    else if (substr(filename, nc - 3, nc) == ".tif") {
             dms <- dim(tiffile@red)
             if (quality < 100)  tiffilename <- paste(substr(filename, 1, nchar(filename) - 4),"-",quality,".tif",sep="")
             if (quality == 100)  tiffilename <- paste(substr(filename, 1, nchar(filename) - 4),".tif",sep="")
             tiff(tiffilename, width = dms[2], height = dms[1], units = "px")
             par(mar=c(0,0,0,0),mai = c(0,0,0,0))
             plot(tiffile)
             dev.off()
    }
    else stop("file extension not recognised.")
}
