prepare_mesh.face3d<- function(face){
#Removing Duplicated Points

a           <- grep("mid-face right", rownames(face$mesh))
remove      <- a[c(115:120)]
face$mesh <- face$mesh[-remove, ]

aa          <- grep("mid-face right", rownames(face$mesh))
remove.2    <- seq(0,length(aa), by=6)	 
remove.2    <- remove.2[-1]
aa          <- aa[remove.2]
face$mesh <- face$mesh[-aa, ]

b           <- grep("mid-face left", rownames(face$mesh))
remove      <- b[c(115:120)]
face$mesh <- face$mesh[-remove, ]
bb          <- grep("mid-face left", rownames(face$mesh))
remove.2    <- seq(0,length(bb), by=6)	 
remove.2    <- remove.2[-1]
bb          <- bb[remove.2]
face$mesh <- face$mesh[-bb, ]

c           <- grep("lower face right", rownames(face$mesh))
remove      <- seq(0,length(c), by=5)	 
remove      <- remove[-1]
c           <- c[remove]
face$mesh <- face$mesh[-c, ]

cc           <- grep("lower face right", rownames(face$mesh))
cc           <- cc[109:112]
face$mesh <- face$mesh[-cc, ]


d           <- grep("lower face left", rownames(face$mesh))
remove      <- seq(0,length(d), by=5)	 
remove      <- remove[-1]
d           <- d[remove]
face$mesh <- face$mesh[-d, ]

dd           <- grep("lower face left", rownames(face$mesh))
dd           <- dd[109:112]
face$mesh <- face$mesh[-dd, ]


i           <- grep("upper mid face right", rownames(face$mesh))
remove      <- i[c(96:100)]
face$mesh <- face$mesh[-remove, ]

k           <- grep("upper mid face left", rownames(face$mesh))
remove      <- k[c(96:100)]
face$mesh <- face$mesh[-remove, ]



n           <- grep("upper face right 3", rownames(face$mesh))
remove      <- n[c(19:21)]
face$mesh <- face$mesh[-remove, ]

n           <- grep("upper face left 3", rownames(face$mesh))
remove      <- n[c(19:21)]
face$mesh <- face$mesh[-remove, ]

zz          <- grep("upper face right 3", rownames(face$mesh))
remove      <- seq(0,length(zz), by=3)
remove      <- remove[-1]
zz          <- zz[remove]
face$mesh <- face$mesh[-zz,]

zz          <- grep("upper face left 3", rownames(face$mesh))
remove      <- seq(0,length(zz), by=3)
remove      <- remove[-1]
zz          <- zz[remove]
face$mesh <- face$mesh[-zz,]

l           <- grep("nose right", rownames(face$mesh))
remove      <- l[c(31:35)]
face$mesh <- face$mesh[-remove, ]

m           <- grep("nose left", rownames(face$mesh))
remove      <- m[c(31:35)]
face$mesh <- face$mesh[-remove, ]




e           <- grep("lower lip right", rownames(face$mesh))
remove      <- seq(0,length(e), by=3)	 
remove      <- remove[-1]
e           <- e[remove]
face$mesh <- face$mesh[-e, ]

g          <- grep("lower lip right", rownames(face$mesh))
remove      <- g[c(13:14)]
face$mesh <- face$mesh[-remove, ]


f           <- grep("lower lip left", rownames(face$mesh))
remove      <- seq(0,length(f), by=3)	 
remove      <- remove[-1]
f           <- f[remove]
face$mesh <- face$mesh[-f, ]

g          <- grep("lower lip left", rownames(face$mesh))
remove      <- g[c(13:14)]
face$mesh <- face$mesh[-remove, ]



g          <- grep("upper lip right", rownames(face$mesh))
remove      <- g[c(19:21)]
face$mesh <- face$mesh[-remove, ]

g          <- grep("upper lip left", rownames(face$mesh))
remove      <- g[c(19:21)]
face$mesh <- face$mesh[-remove, ]








j           <- grep("philtrum right", rownames(face$mesh))
remove      <- j[c(49:54)]
face$mesh <- face$mesh[-remove, ]

j           <- grep("philtrum left", rownames(face$mesh))
remove      <- j[c(49:54)]
face$mesh <- face$mesh[-remove, ]


ii           <- grep("philtrum right", rownames(face$mesh))
remove      <- seq(1,length(ii), by=6) 
ii           <- ii[remove]
face$mesh <- face$mesh[-ii, ]

jj           <- grep("philtrum left", rownames(face$mesh))
remove      <- seq(1,length(jj), by=6) 
jj           <- jj[remove]
face$mesh <- face$mesh[-jj, ]

jjj           <- grep("philtrum right", rownames(face$mesh))
remove      <- jjj[c(6:10)]
face$mesh <- face$mesh[-remove, ]

ll           <- grep("philtrum left", rownames(face$mesh))
remove      <- ll[c(6:10)]
face$mesh <- face$mesh[-remove, ]



M   <- grep("philtrum right", rownames(face$mesh))
path <- face$mesh[M, ]
ind <- (substr(rownames(face$mesh), 1, nchar("philtrum right")) == "philtrum right")
        face$mesh    <- face$mesh[!ind, ]
        rownames(path) <- paste("philtrum right", rep(1:7, each=5), rep(2:6, 7))
        face$mesh    <- rbind(face$mesh, path) 

M   <- grep("philtrum left", rownames(face$mesh))
path <- face$mesh[M, ]
ind <- (substr(rownames(face$mesh), 1, nchar("philtrum left")) == "philtrum left")
        face$mesh    <- face$mesh[!ind, ]
        rownames(path) <- paste("philtrum left", rep(1:7, each=5), rep(2:6, 7))
        face$mesh    <- rbind(face$mesh, path) 



zz          <- grep("upper face right 1", rownames(face$mesh))
remove      <- seq(0,length(zz), by=3)
remove      <- remove[-1]
zz          <- zz[remove]
face$mesh <- face$mesh[-zz,]

zz          <- grep("upper face left 1", rownames(face$mesh))
remove      <- seq(0,length(zz), by=3)
remove      <- remove[-1]
zz          <- zz[remove]
face$mesh <- face$mesh[-zz,]



# Remove upper face right/left 2 and upper eye socket from curves, but add in brow ridge to mesh 
# 
# right        <- grep("upper face right 2", rownames(face$mesh))
# left         <- grep("upper face left 2", rownames(face$mesh))
# remove.right <- right[c(2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27,29,30)]
# remove.left  <- left[c(2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27,29,30)]
# face$mesh  <- face$mesh[-c(remove.right, remove.left), ]

#add mid-lines
X           <- face$curves[grep("mid-line", rownames(face$curves)),]
X           <- X[-c(52:55,70:87),]
face$mesh <- rbind(face$mesh,X)

M <- grep("mid-line nasal profile", rownames(face$mesh))
path <- resample.face3d(face$mesh[M,], n=5)
ind            <- (substr(rownames(face$mesh), 1, nchar("mid-line nasal profile")) == "mid-line nasal profile")
        face$mesh    <- face$mesh[!ind, ]
        rownames(path) <- paste("mid-line nasal profile ", 1:nrow(path), sep="")
        face$mesh    <- rbind(face$mesh, path) 

M <- rev(grep("mid-line chin", rownames(face$mesh)))
path <- resample.face3d(face$mesh[M,], n=4)
ind            <- (substr(rownames(face$mesh), 1, nchar("mid-line chin")) == "mid-line chin")
        face$mesh    <- face$mesh[!ind, ]
        rownames(path) <- paste("mid-line chin ", 1:nrow(path), sep="")
        face$mesh    <- rbind(face$mesh, path) 

M <- grep("mid-line nasal root", rownames(face$mesh))
path <- resample.face3d(face$mesh[M,], n=3)
ind     <- (substr(rownames(face$mesh), 1, nchar("mid-line nasal root")) == "mid-line nasal root")
        face$mesh    <- face$mesh[!ind, ]
        rownames(path) <- paste("mid-line nasal root ", 1:nrow(path), sep="")
        face$mesh    <- rbind(face$mesh, path) 

M <- grep("mid-line philtral", rownames(face$mesh))
path <- resample.face3d(face$mesh[M,], n=6)
ind     <- (substr(rownames(face$mesh), 1, nchar("mid-line philtral")) == "mid-line philtral")
        face$mesh    <- face$mesh[!ind, ]
        rownames(path) <- paste("mid-line philtral ", 1:nrow(path), sep="")
        face$mesh    <- rbind(face$mesh, path) 

M <- grep("mid-line bottom lip", rownames(face$mesh))
path <- resample.face3d(face$mesh[M,], n=3)
ind     <- (substr(rownames(face$mesh), 1, nchar("mid-line bottom lip")) == "mid-line bottom lip")
        face$mesh    <- face$mesh[!ind, ]
        rownames(path) <- paste("mid-line bottom lip ", 1:nrow(path), sep="")
        face$mesh    <- rbind(face$mesh, path) 

M <- grep("mid-line upper-lip", rownames(face$mesh))
path <- resample.face3d(face$mesh[M,], n=3)
ind     <- (substr(rownames(face$mesh), 1, nchar("mid-line upper-lip")) == "mid-line upper-lip")
        face$mesh    <- face$mesh[!ind, ]
        rownames(path) <- paste("mid-line upper-lip ", 1:nrow(path), sep="")
        face$mesh    <- rbind(face$mesh, path) 


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

face$mesh    <- rbind(face$mesh, ULR, ULL, LLR, LLL, MLR, MLL) 


#organizing so all in the same order    

mesh <- c(grep("mid-face right", rownames(face$mesh)),
            grep("mid-face left", rownames(face$mesh)),
            grep("upper mid face right", rownames(face$mesh)),
            grep("upper mid face left", rownames(face$mesh)),
            grep("upper lip left", rownames(face$mesh)),
            grep("lower lip left", rownames(face$mesh)),
            grep("philtrum right", rownames(face$mesh)),
            grep("philtrum left", rownames(face$mesh)),
            grep("upper lip right", rownames(face$mesh)),
            grep("lower face right", rownames(face$mesh)),
            grep("nose left", rownames(face$mesh)),
	        grep("nose right", rownames(face$mesh)),
            grep("lower lip right", rownames(face$mesh)),
	        grep("upper face right 1", rownames(face$mesh)),
	        grep("upper face right 2", rownames(face$mesh)),
	        grep("upper face right 3", rownames(face$mesh)),
	        grep("upper face left 1", rownames(face$mesh)),
	        grep("upper face left 2", rownames(face$mesh)),
	        grep("upper face left 3", rownames(face$mesh)),
	        grep("lower face left", rownames(face$mesh)),
	        grep("LLR", rownames(face$mesh)),
	        grep("LLL", rownames(face$mesh)),
	        grep("ULL", rownames(face$mesh)),
	        grep("ULR", rownames(face$mesh)),
	        grep("MLR", rownames(face$mesh)),
	        grep("MLL", rownames(face$mesh)),
	        grep("mid-line", rownames(face$mesh)))
	      
	          
curves   <-  c(grep("lower eye socket right", rownames(face$curves)),
              grep("lower eye socket left", rownames(face$curves)),
              #grep("upper eye socket right", rownames(face$curves)),
              #grep("upper eye socket left", rownames(face$curves)),
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
              grep("mid-line upper-lip", rownames(face$curves)),
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
face$mesh <- face$mesh[mesh, ]
face$curves <- face$curves[curves, ]


invisible(face)
}
    
