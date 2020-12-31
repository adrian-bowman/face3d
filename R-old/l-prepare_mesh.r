prepare_mesh.face3d<- function(face){
#Removing Duplicated Points

a           <- grep("mid-face right", rownames(face$meshes))
remove      <- a[c(115:120)]
face$meshes <- face$meshes[-remove, ]

aa          <- grep("mid-face right", rownames(face$meshes))
remove.2    <- seq(0,length(aa), by=6)	 
remove.2    <- remove.2[-1]
aa          <- aa[remove.2]
face$meshes <- face$meshes[-aa, ]

b           <- grep("mid-face left", rownames(face$meshes))
remove      <- b[c(115:120)]
face$meshes <- face$meshes[-remove, ]
bb          <- grep("mid-face left", rownames(face$meshes))
remove.2    <- seq(0,length(bb), by=6)	 
remove.2    <- remove.2[-1]
bb          <- bb[remove.2]
face$meshes <- face$meshes[-bb, ]

c           <- grep("lower face right", rownames(face$meshes))
remove      <- seq(0,length(c), by=5)	 
remove      <- remove[-1]
c           <- c[remove]
face$meshes <- face$meshes[-c, ]

cc           <- grep("lower face right", rownames(face$meshes))
cc           <- cc[109:112]
face$meshes <- face$meshes[-cc, ]


d           <- grep("lower face left", rownames(face$meshes))
remove      <- seq(0,length(d), by=5)	 
remove      <- remove[-1]
d           <- d[remove]
face$meshes <- face$meshes[-d, ]

dd           <- grep("lower face left", rownames(face$meshes))
dd           <- dd[109:112]
face$meshes <- face$meshes[-dd, ]


i           <- grep("upper mid face right", rownames(face$meshes))
remove      <- i[c(96:100)]
face$meshes <- face$meshes[-remove, ]

k           <- grep("upper mid face left", rownames(face$meshes))
remove      <- k[c(96:100)]
face$meshes <- face$meshes[-remove, ]



n           <- grep("upper face right 3", rownames(face$meshes))
remove      <- n[c(19:21)]
face$meshes <- face$meshes[-remove, ]

n           <- grep("upper face left 3", rownames(face$meshes))
remove      <- n[c(19:21)]
face$meshes <- face$meshes[-remove, ]

zz          <- grep("upper face right 3", rownames(face$meshes))
remove      <- seq(0,length(zz), by=3)
remove      <- remove[-1]
zz          <- zz[remove]
face$meshes <- face$meshes[-zz,]

zz          <- grep("upper face left 3", rownames(face$meshes))
remove      <- seq(0,length(zz), by=3)
remove      <- remove[-1]
zz          <- zz[remove]
face$meshes <- face$meshes[-zz,]

l           <- grep("nose right", rownames(face$meshes))
remove      <- l[c(31:35)]
face$meshes <- face$meshes[-remove, ]

m           <- grep("nose left", rownames(face$meshes))
remove      <- m[c(31:35)]
face$meshes <- face$meshes[-remove, ]




e           <- grep("lower lip right", rownames(face$meshes))
remove      <- seq(0,length(e), by=3)	 
remove      <- remove[-1]
e           <- e[remove]
face$meshes <- face$meshes[-e, ]

g          <- grep("lower lip right", rownames(face$meshes))
remove      <- g[c(13:14)]
face$meshes <- face$meshes[-remove, ]


f           <- grep("lower lip left", rownames(face$meshes))
remove      <- seq(0,length(f), by=3)	 
remove      <- remove[-1]
f           <- f[remove]
face$meshes <- face$meshes[-f, ]

g          <- grep("lower lip left", rownames(face$meshes))
remove      <- g[c(13:14)]
face$meshes <- face$meshes[-remove, ]



g          <- grep("upper lip right", rownames(face$meshes))
remove      <- g[c(19:21)]
face$meshes <- face$meshes[-remove, ]

g          <- grep("upper lip left", rownames(face$meshes))
remove      <- g[c(19:21)]
face$meshes <- face$meshes[-remove, ]








j           <- grep("philtrum right", rownames(face$meshes))
remove      <- j[c(49:54)]
face$meshes <- face$meshes[-remove, ]

j           <- grep("philtrum left", rownames(face$meshes))
remove      <- j[c(49:54)]
face$meshes <- face$meshes[-remove, ]


ii           <- grep("philtrum right", rownames(face$meshes))
remove      <- seq(1,length(ii), by=6) 
ii           <- ii[remove]
face$meshes <- face$meshes[-ii, ]

jj           <- grep("philtrum left", rownames(face$meshes))
remove      <- seq(1,length(jj), by=6) 
jj           <- jj[remove]
face$meshes <- face$meshes[-jj, ]

jjj           <- grep("philtrum right", rownames(face$meshes))
remove      <- jjj[c(6:10)]
face$meshes <- face$meshes[-remove, ]

ll           <- grep("philtrum left", rownames(face$meshes))
remove      <- ll[c(6:10)]
face$meshes <- face$meshes[-remove, ]



M   <- grep("philtrum right", rownames(face$meshes))
path <- face$meshes[M, ]
ind <- (substr(rownames(face$meshes), 1, nchar("philtrum right")) == "philtrum right")
        face$meshes    <- face$meshes[!ind, ]
        rownames(path) <- paste("philtrum right", rep(1:7, each=5), rep(2:6, 7))
        face$meshes    <- rbind(face$meshes, path) 

M   <- grep("philtrum left", rownames(face$meshes))
path <- face$meshes[M, ]
ind <- (substr(rownames(face$meshes), 1, nchar("philtrum left")) == "philtrum left")
        face$meshes    <- face$meshes[!ind, ]
        rownames(path) <- paste("philtrum left", rep(1:7, each=5), rep(2:6, 7))
        face$meshes    <- rbind(face$meshes, path) 



zz          <- grep("upper face right 1", rownames(face$meshes))
remove      <- seq(0,length(zz), by=3)
remove      <- remove[-1]
zz          <- zz[remove]
face$meshes <- face$meshes[-zz,]

zz          <- grep("upper face left 1", rownames(face$meshes))
remove      <- seq(0,length(zz), by=3)
remove      <- remove[-1]
zz          <- zz[remove]
face$meshes <- face$meshes[-zz,]



# Remove upper face right/left 2 and upper eye socket from curves, but add in brow ridge to mesh 

right        <- grep("upper face right 2", rownames(face$meshes))
left         <- grep("upper face left 2", rownames(face$meshes))
remove.right <- right[c(2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27,29,30)]
remove.left  <- left[c(2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27,29,30)]
face$meshes  <- face$meshes[-c(remove.right, remove.left), ]

#add mid-lines
X           <- face$curves[grep("mid-line", rownames(face$curves)),]
X           <- X[-c(52:55,70:87),]
face$meshes <- rbind(face$meshes,X)

M <- grep("mid-line nasal profile", rownames(face$meshes))
path <- resample.face3d(face$meshes[M,], n=5)
ind            <- (substr(rownames(face$meshes), 1, nchar("mid-line nasal profile")) == "mid-line nasal profile")
        face$meshes    <- face$meshes[!ind, ]
        rownames(path) <- paste("mid-line nasal profile ", 1:nrow(path), sep="")
        face$meshes    <- rbind(face$meshes, path) 

M <- rev(grep("mid-line chin", rownames(face$meshes)))
path <- resample.face3d(face$meshes[M,], n=4)
ind            <- (substr(rownames(face$meshes), 1, nchar("mid-line chin")) == "mid-line chin")
        face$meshes    <- face$meshes[!ind, ]
        rownames(path) <- paste("mid-line chin ", 1:nrow(path), sep="")
        face$meshes    <- rbind(face$meshes, path) 

M <- grep("mid-line nasal root", rownames(face$meshes))
path <- resample.face3d(face$meshes[M,], n=3)
ind     <- (substr(rownames(face$meshes), 1, nchar("mid-line nasal root")) == "mid-line nasal root")
        face$meshes    <- face$meshes[!ind, ]
        rownames(path) <- paste("mid-line nasal root ", 1:nrow(path), sep="")
        face$meshes    <- rbind(face$meshes, path) 

M <- grep("mid-line philtral", rownames(face$meshes))
path <- resample.face3d(face$meshes[M,], n=6)
ind     <- (substr(rownames(face$meshes), 1, nchar("mid-line philtral")) == "mid-line philtral")
        face$meshes    <- face$meshes[!ind, ]
        rownames(path) <- paste("mid-line philtral ", 1:nrow(path), sep="")
        face$meshes    <- rbind(face$meshes, path) 

M <- grep("mid-line bottom lip", rownames(face$meshes))
path <- resample.face3d(face$meshes[M,], n=3)
ind     <- (substr(rownames(face$meshes), 1, nchar("mid-line bottom lip")) == "mid-line bottom lip")
        face$meshes    <- face$meshes[!ind, ]
        rownames(path) <- paste("mid-line bottom lip ", 1:nrow(path), sep="")
        face$meshes    <- rbind(face$meshes, path) 

M <- grep("mid-line upper lip", rownames(face$meshes))
path <- resample.face3d(face$meshes[M,], n=3)
ind     <- (substr(rownames(face$meshes), 1, nchar("mid-line upper lip")) == "mid-line upper lip")
        face$meshes    <- face$meshes[!ind, ]
        rownames(path) <- paste("mid-line upper lip ", 1:nrow(path), sep="")
        face$meshes    <- rbind(face$meshes, path) 


#add lip curves
ULR <- face$curves[grep("upper lip right", rownames(face$curves))[1:2], ]
ULL <- face$curves[grep("upper lip left", rownames(face$curves))[1:2], ]
LLR <- face$curves[grep("lower lip right", rownames(face$curves))[1:2], ]
LLL <- face$curves[grep("lower lip left", rownames(face$curves))[1:2], ]
MLL <- face$curves[grep("mid-line lip left", rownames(face$curves))[1:2], ]
MLR <- face$curves[grep("mid-line lip right", rownames(face$curves))[1:2], ]
rownames(ULR) <- paste("ULR", 1:nrow(ULR))
rownames(ULL) <- paste("ULL", 1:nrow(ULL))
rownames(LLR) <- paste("LLR", 1:nrow(LLR))
rownames(LLL) <- paste("LLL", 1:nrow(LLL))
rownames(MLR) <- paste("MLR", 1:nrow(MLR))
rownames(MLL) <- paste("MLL", 1:nrow(MLL))

face$meshes    <- rbind(face$meshes, ULR, ULL, LLR, LLL, MLR, MLL) 


#organizing so all in the same order    

meshes <- c(grep("mid-face right", rownames(face$meshes)),
            grep("mid-face left", rownames(face$meshes)),
            grep("upper mid face right", rownames(face$meshes)),
            grep("upper mid face left", rownames(face$meshes)),
            grep("upper lip left", rownames(face$meshes)),
            grep("lower lip left", rownames(face$meshes)),
            grep("philtrum right", rownames(face$meshes)),
            grep("philtrum left", rownames(face$meshes)),
            grep("upper lip right", rownames(face$meshes)),
            grep("lower face right", rownames(face$meshes)),
            grep("nose left", rownames(face$meshes)),
	        grep("nose right", rownames(face$meshes)),
            grep("lower lip right", rownames(face$meshes)),
	        grep("upper face right 1", rownames(face$meshes)),
	        grep("upper face right 2", rownames(face$meshes)),
	        grep("upper face right 3", rownames(face$meshes)),
	        grep("upper face left 1", rownames(face$meshes)),
	        grep("upper face left 2", rownames(face$meshes)),
	        grep("upper face left 3", rownames(face$meshes)),
	        grep("lower face left", rownames(face$meshes)),
	        grep("LLR", rownames(face$meshes)),
	        grep("LLL", rownames(face$meshes)),
	        grep("ULL", rownames(face$meshes)),
	        grep("ULR", rownames(face$meshes)),
	        grep("MLR", rownames(face$meshes)),
	        grep("MLL", rownames(face$meshes)),
	        grep("mid-line", rownames(face$meshes)))
	      
	          
curves   <-  c(grep("lower eye socket right", rownames(face$curves)),
              grep("lower eye socket left", rownames(face$curves)),
              grep("upper eye socket right", rownames(face$curves)),
              grep("upper eye socket left", rownames(face$curves)),
              grep("nasal root right", rownames(face$curves)),
              grep("nasal root left", rownames(face$curves)),
              grep("nasal boundary right", rownames(face$curves)),
              grep("nasal boundary left", rownames(face$curves)),
              grep("nasal bridge left", rownames(face$curves)),
              grep("nasal bridge right", rownames(face$curves)),
              grep("nasolabial right", rownames(face$curves)),
              grep("nasolabial left", rownames(face$curves)),
              grep("nasal base left", rownames(face$curves)),
              grep("nasal base right", rownames(face$curves)),                
              grep("mid-line nasal profile", rownames(face$curves)),
              grep("mid-line nasal root", rownames(face$curves)),
              grep("mid-line columella", rownames(face$curves)),                 
              grep("mid-line philtral", rownames(face$curves)),
              grep("mid-line upper lip", rownames(face$curves)),
              grep("mid-line bottom lip", rownames(face$curves)),
              grep("mid-line mentolabial", rownames(face$curves)),
              grep("mid-line chin", rownames(face$curves)),
              grep("lower lip left", rownames(face$curves)),
              grep("lower lip right", rownames(face$curves)),
              grep("upper lip left", rownames(face$curves)),
              grep("upper lip right", rownames(face$curves)),
              grep("philtrum ridge right", rownames(face$curves)),
              grep("philtrum ridge left", rownames(face$curves)),
              grep("brow ridge left", rownames(face$curves)),
              grep("brow ridge right", rownames(face$curves)),
              grep("mid-line lip left", rownames(face$curves)),
              grep("mid-line lip right", rownames(face$curves)),
              grep("mandible left", rownames(face$curves)),
              grep("mandible right", rownames(face$curves)),
              grep("philtrum lip left", rownames(face$curves)),
              grep("philtrum lip right", rownames(face$curves)),
              grep("cheek-eye left", rownames(face$curves)),
              grep("cheek-eye right", rownames(face$curves)),
              grep("cheek-nose left", rownames(face$curves)),
              grep("cheek-nose right", rownames(face$curves)),
              grep("cheek-lip left", rownames(face$curves)),
              grep("cheek-lip right", rownames(face$curves)))
face$meshes <- face$meshes[meshes, ]
face$curves <- face$curves[curves, ]

invisible(face)
}
    
