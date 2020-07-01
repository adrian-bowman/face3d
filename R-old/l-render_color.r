 render_color.face3d<- function(face, reference, target, landmark.names, nbrk = c(-30, seq(-10,-1,1), seq(1,10,1),30)){
  
  
  
  

################################################################
#rendering original object
rownames(reference)     <- landmark.names
rownames(target)        <- landmark.names

UFR.R     <- reference[grep("upper face right 2", rownames(reference)), ]
UFL.R     <- reference[grep("upper face left 2", rownames(reference)), ]
MLCOLUM.R <- reference[grep("mid-line columella", rownames(reference)), ]
UFR.T     <- target[grep("upper face right 2", rownames(target)), ]
UFL.T     <- target[grep("upper face left 2", rownames(target)), ]
MLCOLUM.T <- target[grep("mid-line columella", rownames(target)), ]

landmark.names <- landmark.names[-grep("upper face right 2", landmark.names)]
landmark.names <- landmark.names[-grep("upper face left 2", landmark.names)]
landmark.names <- landmark.names[-grep("mid-line columella ", landmark.names)]
reference <- reference[-grep("upper face right 2", rownames(reference)), ]
reference <- reference[-grep("upper face left 2", rownames(reference)), ]
reference <- reference[-grep("mid-line columella", rownames(reference)), ]
target    <- target[-grep("upper face right 2", rownames(target)), ]
target    <- target[-grep("upper face left 2", rownames(target)), ]
target    <- target[-grep("mid-line columella", rownames(target)), ]

rownames(reference)     <- landmark.names
rownames(target)        <- landmark.names
rend.reference          <- render.face3d(reference)
rend.target             <- render.face3d(target)

# creating painted face by normals +/-

normals               <- normals.face3d(rend.reference)$normal
difference            <- rend.reference$coords - rend.target$coords 
euc.dist              <- sqrt(apply(difference^2, 1, sum))
ind                   <- apply(difference * normals, 1, sum)
pts                   <- euc.dist * sign(ind)
#nbrk                  <- c(-30, seq(-10,-1,1), seq(1,10,1),30)
cc                    <- topo.colors(length(nbrk))
colours               <- cc[cut(pts, nbrk)]
if (any(is.na(colours))==TRUE){
corners               <- which(is.na(colours)==TRUE)
# print(length(corners))
colours[corners]      <- colours[corners + 1 ]
}
rend.reference$colour <- colours
rend.target$colour    <- colours



#calling up the brow and columella


#Find normals of brows in original face in order to paint them properly
if(!("normals" %in% names(face))) face <- normals.face3d(face)
face.brow.R              <- face$meshes[grep("upper face right 2", rownames(face$meshes)),]
brow.R.ids               <- rep(NA)
for(k in 1:length(face.brow.R[,1])){
brow.R.ids[k]            <- closest.face3d(face.brow.R[k,], face)$id[1]
}
brow.R.normals.reference <- face$normals[brow.R.ids, ]


face.brow.L              <- face$meshes[grep("upper face left 2", rownames(face$meshes)),]
brow.L.ids               <- rep(NA)
for(k in 1:length(face.brow.L[,1])){
brow.L.ids[k]            <- closest.face3d(face.brow.L[k,], face)$id[1]
}
brow.L.normals.reference <- face$normals[brow.L.ids, ]


reference.brow.R       <- UFR.R
reference.brow.L       <- UFL.R
columella.post.ref     <- MLCOLUM.R
target.brow.R          <- UFR.T
target.brow.L          <- UFL.T
columella.post.tar     <- MLCOLUM.T



# adding on brow ridge
# left side
normals <- brow.L.normals.reference 
difference <- reference.brow.L -  target.brow.L
euc.dist        <- sqrt(apply(difference^2, 1, sum))
ind             <- apply(difference * normals, 1, sum)
pts             <- euc.dist * sign(ind)
#nbrk            <- c(seq(-30,-1,1), seq(1,30,1))
cc   <- topo.colors(length(nbrk))
colours.l <- cc[cut(pts, nbrk)]


#right side
normals <- brow.R.normals.reference 
difference <- reference.brow.R - target.brow.R
euc.dist        <- sqrt(apply(difference^2, 1, sum))
ind             <- apply(difference * normals, 1, sum)
pts             <- euc.dist * sign(ind)
#nbrk            <- c(seq(-30,-1,1), seq(1,30,1))
cc   <- topo.colors(length(nbrk))
colours.r <- cc[cut(pts, nbrk)]

   
brow.ref.coords <- rbind(reference.brow.R, reference.brow.L) 
brow.ref.colour <- c(colours.r, colours.l)
brow.tar.coords <- rbind(target.brow.R, target.brow.L)  
  
rend.reference$brow <- brow.ref.coords
rend.reference$brow.colour <- brow.ref.colour
rend.target$brow    <- brow.tar.coords



rend.reference$extra <- columella.post.ref
rend.target$extra    <- columella.post.tar

# display.face3d(rend.reference)
# lines3d(rend.reference$brow[1:10, ], col = rend.reference$brow.colour[1:10])
# lines3d(rend.reference$brow[11:20, ], col = rend.reference$brow.colour[11:20])
# lines3d(rend.reference$extra)

      invisible(list(rend.reference = rend.reference, rend.target = rend.target))
}
