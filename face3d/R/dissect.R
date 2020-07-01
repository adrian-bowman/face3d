dissect.face3d <- function(shape,  curve, endedge = F){
  
  if (!is.face3d(shape)) stop("this is not a face3d shape.")
  ##################check duplicates#############################################
  if (any(duplicated(curve))) stop("there are duplicates on the path")
  ##################check if the path is on the shape############################
  distshape <- apply(curve, 1, function(x)closest.face3d(x, shape)$distance)
  if(!all(distshape<0.001)) stop("the curve is not on the shape")
  
  shape <- index.face3d(shape)
  ##################check in-between edges###########################################
  triangles <-  matrix(shape$triples,  ncol = 3,  byrow = TRUE)
  inbetween <- FALSE
  for (i in 1:(nrow(curve)-1)) {
    if (check.edges(curve[i, ], curve[i+1,], triangles, shape)){
      inbetween <- TRUE
      # stop(paste("point",i,"and point",i+1,"have edges in between.",sep = " "))
      break
    }
  }
  if(inbetween){
    curve2 <- rbind(curve, curve[1,])
    curve <- matrix(nrow=0,ncol=3)
    for(i in 1:(nrow(curve2)-1)){
      pp <- planepath.face3d(shape, curve2[i,], curve2[i+1,])$path
      curve <- rbind(curve,pp[-nrow(pp),])
    }
  }
  ##################find edges which are crossed by the path#########################
  edge.mat <- matrix(nrow=0, ncol = 2)
  cids     <- vector()
  for(i in 1:nrow(curve)){
    ids <- pospts(curve[i, ], shape)$id
    if (length(ids)==2) {
      edge.mat <- rbind(edge.mat, ids)
      cids     <- c(cids, i)
    }
  }
  
  # plot(shape, colour="grey", display = "mesh")
  # spheres3d(curve, radius = 0.05)
  # spheres3d(shape$coords[edge.mat[1, ], ], radius = 0.05, col=2:3)
  # spheres3d(curve[cids, ], radius = 0.06, col=2)
  
  #####################account,  for each edge,  that how many path points lying on################################
  edge.ts <- rep(0, nrow(edge.mat))
  for (i in 1:nrow(edge.mat)) {
    for (j in (1:nrow(edge.mat))[-i]) {
      if (length(which(edge.mat[i, ]%in%edge.mat[j, ]))==2) {
        edge.ts[i] <- j
      }
    }
  }
  
  # table(edge.ts)
  # spheres3d(curve[cids[c(221, 223)], ], radius = 0.06, col=2:3)
  
  ###############################specify inner and outer vertices to the path###############################################
  norm.mat <- matrix(ncol = 3, nrow=nrow(edge.mat))
  dir.mat  <- matrix(ncol = 3,  nrow = nrow(edge.mat))
  cor.vec  <- vector()
  # inn.ind  <- vector()
  # out.ind  <- vector()
  inn.ind  <- rep(0, nrow(curve))
  out.ind  <- rep(0, nrow(curve))
  ecurve   <- rbind(curve, curve[1, ])
  
  for (i in 1:nrow(edge.mat)) {
    if (!edge.ts[i]==0) {
      
      names(inn.ind)[cids[i]] <- "Multi"
      names(out.ind)[cids[i]] <- "Multi"
      
      dir.vec  <- ecurve[cids[i]+1, ]-ecurve[cids[i], ]
      nor.vec  <- (shape$normals[edge.mat[i, 1], ]+shape$normals[edge.mat[i, 2], ])/2
      axi.vec  <- crossproduct(nor.vec, dir.vec)
      proj.pts <- c(curve[cids[i], ]%*%t(axi.vec), curve[cids[edge.ts[i]], ]%*%t(axi.vec))
      proj.vtc <- c(shape$coords[edge.mat[i,  1], ]%*%t(axi.vec), shape$coords[edge.mat[i,  2], ]%*%t(axi.vec))
      ide      <- which.min(abs(proj.vtc-proj.pts[1]))
      if (which.min(abs(proj.pts-proj.vtc[ide]))==1) {
        if (proj.vtc[ide] < proj.pts[1]) {
          out.ind[cids[i]] <- edge.mat[i, ide]
        } else if (!proj.vtc[ide] < proj.pts[1]) {
          inn.ind[cids[i]] <- edge.mat[i, ide]
        }
      } else if (which.min(abs(proj.pts-proj.vtc[-ide]))==1) {
        if(proj.vtc[-ide] < proj.pts[1]){
          out.ind[cids[i]] <- edge.mat[i, -ide]
        } else if(!proj.vtc[-ide] < proj.pts[1]){
          inn.ind[cids[i]] <- edge.mat[i, -ide]
        }
      }
    } else {
      cond <- TRUE

      if(cids[i]==1){
        scids <- unlist(c(apply(curve[c(1, 2, nrow(curve)), ], 1, function(x)pospts(x, shape)$id)))
        if(length(unique(scids))==3){
          cond    <- FALSE
          cor.vec <- c(cor.vec ,i)
        }
      } else if(cids[i]==nrow(curve)){
        scids <- unlist(c(apply(curve[c(1, cids[i], cids[i]-1), ], 1, function(x)pospts(x, shape)$id)))
        if(length(unique(scids))==3){
          cond    <- FALSE
          cor.vec <- c(cor.vec ,i)
        }
      } else {
        scids <- unlist(c(apply(curve[c(cids[i]+1, cids[i], cids[i]-1), ], 1, function(x)pospts(x, shape)$id)))
        if(length(unique(scids))==3){
          cond    <- FALSE
          cor.vec <- c(cor.vec ,i)
        }
      }
      if(cond){
        norm.mat[i, ] <- (shape$normals[edge.mat[i, 1], ]+shape$normals[edge.mat[i, 2], ])/2
        dir.mat[i, ]  <- ecurve[cids[i]+1, ]-ecurve[cids[i], ]
        axi.vec       <- crossproduct(norm.mat[i, ], dir.mat[i, ])
        pg1           <- shape$coords[edge.mat[i, 1], ]%*%t(axi.vec)
        pg2           <- shape$coords[edge.mat[i, 2], ]%*%t(axi.vec)
        
        inn.ind[cids[i]]        <- edge.mat[i, which.max(c(pg1, pg2))]
        out.ind[cids[i]]        <- edge.mat[i, which.min(c(pg1, pg2))]
        names(inn.ind)[cids[i]] <- "Edge"
        names(out.ind)[cids[i]] <- "Edge"
      } else {
        names(inn.ind)[cids[i]] <- "Corner"
        names(out.ind)[cids[i]] <- "Corner"
      }
    }
  }
  names(inn.ind)[which(is.na(names(inn.ind)))] <- "VI"
  names(out.ind)[which(is.na(names(out.ind)))] <- "VI"
  ##################################delete the triangles and cut the surface#####################################
  
  
  idtri <- vector()
  for(i in 1:nrow(edge.mat)){
    if(!names(inn.ind)[cids[i]]=="Corner"){
      for(j in 1:nrow(triangles)){
        if(length(which(edge.mat[i, ]%in%triangles[j, ]))%in%2:3){  
          idtri <- c(idtri, j)
          # print(j)
        }
      }
    }
  }
  s1 <- shape
  if(endedge){
    if(pospts(curve[1,],shape)$pos=="vertex"){
      idtri <- c(idtri ,c(which(triangles[,1]==pospts(curve[1,],shape)$id)
                         ,which(triangles[,2]==pospts(curve[1,],shape)$id),which(triangles[,3]==pospts(curve[1,],shape)$id)))
    }
    if(pospts(curve[nrow(curve),],shape)$pos=="vertex"){
      idtri <- c(idtri ,c(which(triangles[,1]==pospts(nrow(curve),shape)$id)
                         ,which(triangles[,2]==pospts(nrow(curve),shape)$id),which(triangles[,3]==pospts(nrow(curve),shape)$id)))
    }
  }
  s1$triples <- c(t(triangles[-idtri, ]))
  path.s     <- connected.face3d(s1)
  
  ##################################inner surface###################################################################
  
  subshape <- subset.face3d(s1, (!path.s==(path.s[out.ind[which(names(out.ind)=="Edge")[1]]]))&(!1:(nrow(shape$coords))%in%out.ind[which(names(out.ind)=="Edge")]), remove.singles = F)
  
  ##################################New triangulation#################################
  orow       <- rownames(subshape$coords)
  lrow       <- nrow(subshape$coords)
  inn.ind3   <- sapply(inn.ind, function(x)if(!x==0){return(which(orow==x))}else{return(x)})
  newcoords  <- rbind(subshape$coords, curve)
  iid3       <- c(inn.ind3 ,inn.ind3[1])
  newtriples <- vector()
  for(i in 1:nrow(curve)){
    if(i==nrow(curve) & endedge) break
    if(names(iid3)[i]=="Corner"){
      if(i>1){
        if(pospts(curve[i,],subshape)$pos=="edge"){
          idcor <- pospts(curve[i,],subshape)$id
          if(newtriples[length(newtriples)-1]%in%idcor){
            newtriples <- c(newtriples, i+lrow, iid3[i+1], idcor[which(!idcor==newtriples[length(newtriples)-1])],i+1+lrow ,iid3[i+1],i+lrow)
            # if(!(length(newtriples)%%3)==0) print(i)
          } else {
            idmul <- sapply(pospts(curve[i+1,],shape)$id ,function(x)which(orow==x))
            if(length(unlist(idmul))==length(idmul)){
              newtriples <- c(newtriples, i+lrow, idcor[which(!idcor%in%idmul)], iid3[i-1], i+1+lrow ,iid3[i+1],i+lrow)
            } else {
              newtriples <- c(newtriples, i+1+lrow, newtriples[length(newtriples)-1], i+lrow)
              # if(!(length(newtriples)%%3)==0) print(i)
            }
          }
        } else {
          newtriples <- c(newtriples, i+1+lrow, newtriples[length(newtriples)-1], i+lrow)
          # if(!(length(newtriples)%%3)==0) print(i)
        } 
      } else {
        if(pospts(curve[i,],subshape)$pos=="edge"){
          idcor <- pospts(curve[i,],subshape)$id
          nonids <- iid3[which(iid3>0)]
          if(nonids[length(nonids)]%in%idcor){
            newtriples <- c(newtriples, i+lrow, iid3[i+1], idcor[which(!idcor==nonids[length(nonids)])],i+1+lrow ,iid3[i+1],i+lrow)
            # if(!(length(newtriples)%%3)==0) print(i)
          } else {
            idmul <- sapply(pospts(curve[i+1,],shape)$id ,function(x)which(orow==x))
            newtriples <- c(newtriples, i+lrow, idcor[which(!idcor%in%idmul)], nonids[length(nonids)], i+1+lrow ,iid3[i+1],i+lrow)
            # if(!(length(newtriples)%%3)==0) print(i)
          }
        } else {
          newtriples <- c(newtriples, i+1+lrow, nonids[length(nonids)], i+lrow)
          # if(!(length(newtriples)%%3)==0) print(i)
        } 
      }
    } else if(i==1 & all(iid3[i:(i+1)]==0)){
      nonids <- iid3[which(iid3>0)]
      newtriples <- c(newtriples, i+1+lrow, nonids[length(nonids)], i+lrow)
      # if(!(length(newtriples)%%3)==0) print(i)
    } else if(iid3[i]==0){
      if(all(pospts(curve[i,],shape)$id%in%orow)){
        if(pospts(curve[i,],shape)$pos=="inside" & names(iid3[i+1])=="Multi"){
          if(i>1){
            if(names(iid3[i-1])=="Multi"){
              idcor      <- sapply(pospts(curve[i,],shape)$id ,function(x)which(orow==x))
              idmul      <- sapply(pospts(curve[i+1,],shape)$id ,function(x)which(orow==x))
              newtriples <- c(newtriples, i+lrow, newtriples[length(newtriples)-1], idcor[which(!idcor%in%idmul)], i+lrow, idcor[which(!idcor%in%idmul)], iid3[i+1], i+lrow, iid3[i+1], i+1+lrow)
              # if(!(length(newtriples)%%3)==0) print(i)
            }
          } else {
            if(names(iid3[length(inn.ind3)])=="Multi"){
              idcor      <- sapply(pospts(curve[i,],shape)$id ,function(x)which(orow==x))
              idmul      <- sapply(pospts(curve[i+1,],shape)$id ,function(x)which(orow==x))
              newtriples <- c(newtriples, i+lrow, iid3[length(inn.ind3)], idcor[which(!idcor%in%idmul)], i+lrow, idcor[which(!idcor%in%idmul)], iid3[i+1], i+lrow, iid3[i+1], i+1+lrow)
              # if(!(length(newtriples)%%3)==0) print(i)
            }
          }
        }
      } else if(iid3[i+1]==0){
        newtriples <- c(newtriples, i+1+lrow, newtriples[length(newtriples)-1], i+lrow)
        # if(!(length(newtriples)%%3)==0) print(i)
      } else {
        if(i==1){
          if(!iid3[length(inn.ind3)]==0){
            newtriples <- c(newtriples, i+1+lrow, iid3[i+1], i+lrow, i+lrow, iid3[i+1], iid3[length(inn.ind3)])
            if(!(length(newtriples)%%3)==0) print(i)
          } else if(endedge){
            newtriples <- c(newtriples, i+1+lrow, iid3[i+1], i+lrow)
          }
        } else if(!iid3[i-1]==0){
          newtriples <- c(newtriples, i+1+lrow, iid3[i+1], i+lrow, i+lrow, iid3[i+1], iid3[i-1])
          # if(!(length(newtriples)%%3)==0) print(i)
        } else {
          newtriples <- c(newtriples, i+1+lrow, iid3[i+1], i+lrow)
          # if(!(length(newtriples)%%3)==0) print(i)
        }
      }
    } else {
      if(iid3[i+1]==0){
        newtriples <- c(newtriples, i+1+lrow, iid3[i], i+lrow)
        # if(!(length(newtriples)%%3)==0) print(i)
      } else {
        if(iid3[i+1]==iid3[i]){
          newtriples <- c(newtriples, i+1+lrow, iid3[i+1], i+lrow)
          # if(!(length(newtriples)%%3)==0) print(i)
        } else {
          newtriples <- c(newtriples, i+1+lrow, iid3[i+1], i+lrow, i+lrow, iid3[i+1], iid3[i])
          # if(!(length(newtriples)%%3)==0) print(i)
        }
      } 
    }
    # if(length(which(newtriples==0))>0){
    #   print(i)
    # }
  }
  newtriples[which(newtriples==(nrow(newcoords)+1))] <- nrow(subshape$coords)+1
  newlist                                            <- list(coords=newcoords, triples=c(subshape$triples, newtriples))
  shapeL                                             <- as.face3d(newlist)
  
  
  
  ##################################outer surface###################################################################  
  
  # subshape2 <- subset.face3d(s1, (!path.s==path.s[inn.ind[1]])&(!1:(nrow(shape$coords))%in%inn.ind), remove.singles = F)
  subshape2 <- subset.face3d(s1, (!path.s==(path.s[inn.ind[which(names(inn.ind)=="Edge")[1]]]))&(!1:(nrow(shape$coords))%in%inn.ind[which(names(inn.ind)=="Edge")]), remove.singles = F)
  # subshape2 <- subset.face3d(s1, (path.s==path.s[out.ind[1]])|(1:(nrow(shape$coords))%in%out.ind), remove.singles = F)
  
  ##################################New triangulation#################################
  orow        <- rownames(subshape2$coords)
  lrow        <- nrow(subshape2$coords)
  out.ind3    <- sapply(out.ind, function(x)if(!x==0){return(which(orow==x))}else{return(x)})
  newcoords   <- rbind(subshape2$coords, curve)
  oid3        <- c(out.ind3 ,out.ind3[1])
  newtriples2 <- vector()
  for(i in 1:nrow(curve)){
    if(i==nrow(curve) & endedge) break
    if(names(oid3)[i]=="Corner"){
      if(i>1){
        if(pospts(curve[i,],subshape2)$pos=="edge"){
          idcor <- pospts(curve[i,],subshape2)$id
          if(newtriples2[length(newtriples2)-1]%in%idcor){
            newtriples2 <- c(newtriples2, i+lrow, idcor[which(!idcor==newtriples2[length(newtriples2)-1])], oid3[i+1],i+lrow ,oid3[i+1],i+1+lrow)
            # if(length(newtriples2)%%3)
          } else {
            idmul <- sapply(pospts(curve[i+1,],shape)$id ,function(x)which(orow==x))
            newtriples2 <- c(newtriples2, i+lrow, oid3[i-1], idcor[which(!idcor%in%idmul)], i+lrow ,oid3[i+1],i+1+lrow)
          }
        } else {
          newtriples2 <- c(newtriples2, i+lrow, newtriples2[length(newtriples2)-1], i+1+lrow)
        } 
      } else {
        if(pospts(curve[i,],subshape2)$pos=="edge"){
          idcor <- pospts(curve[i,],subshape2)$id
          nonids <- oid3[which(oid3>0)]
          if(nonids[length(nonids)]%in%idcor){
            newtriples2 <- c(newtriples2, i+lrow, idcor[which(!idcor==nonids[length(nonids)])], oid3[i+1],i+lrow ,oid3[i+1],i+1+lrow)
          } else {
            idmul <- sapply(pospts(curve[i+1,],shape)$id ,function(x)which(orow==x))
            newtriples2 <- c(newtriples2, i+lrow, nonids[length(nonids)], idcor[which(!idcor%in%idmul)], i+lrow ,oid3[i+1],i+1+lrow)
          }
        } else {
          newtriples2 <- c(newtriples2, i+lrow, nonids[length(nonids)], i+1+lrow)
        } 
      }
    } else if(i==1 & all(oid3[i:(i+1)]==0)){
      nonids      <- oid3[which(oid3>0)]
      newtriples2 <- c(newtriples2, i+lrow, nonids[length(nonids)], i+1+lrow)
    } else if(oid3[i]==0){
      if(all(pospts(curve[i,],shape)$id%in%orow)){
        if(pospts(curve[i,],shape)$pos=="inside" & names(oid3)[i+1]=="Multi"){
          if(i>1){
            if(names(oid3[i-1])=="Multi"){
              idcor       <- sapply(pospts(curve[i,],shape)$id ,function(x)which(orow==x))
              idmul       <- sapply(pospts(curve[i+1,],shape)$id ,function(x)which(orow==x))
              newtriples2 <- c(newtriples2, i+lrow, newtriples2[length(newtriples2)-1], idcor[which(!idcor%in%idmul)], i+lrow, idcor[which(!idcor%in%idmul)], oid3[i+1], i+lrow, oid3[i+1], i+1+lrow)
            }
          } else {
            if(names(oid3[length(out.ind3)])=="Multi"){
              idcor       <- sapply(pospts(curve[i,],shape)$id ,function(x)which(orow==x))
              idmul       <- sapply(pospts(curve[i+1,],shape)$id ,function(x)which(orow==x))
              newtriples2 <- c(newtriples2, i+lrow, oid3[length(out.ind3)], idcor[which(!idcor%in%idmul)], i+lrow, idcor[which(!idcor%in%idmul)], oid3[i+1], i+lrow, oid3[i+1], i+1+lrow)
            }
          }
        }
      } else if(oid3[i+1]==0){
        newtriples2 <- c(newtriples2, i+lrow, newtriples2[length(newtriples2)-1], i+1+lrow)
      } else {
        if(i==1){
          if(!oid3[length(out.ind3)]==0){
            newtriples2 <- c(newtriples2, i+lrow, oid3[i+1], i+1+lrow, oid3[length(out.ind3)], oid3[i+1], i+lrow)
          } else if(endedge){
            newtriples2 <- c(newtriples2, i+lrow, oid3[i+1], i+1+lrow)
          }
        } else if(!oid3[i-1]==0){
          newtriples2 <- c(newtriples2, i+lrow, oid3[i+1], i+1+lrow, oid3[i-1], oid3[i+1], i+lrow)
        } else {
          newtriples2 <- c(newtriples2, i+lrow, oid3[i+1], i+1+lrow)
        }
      }
    } else {
      if(oid3[i+1]==0){
        newtriples2 <- c(newtriples2, i+lrow, oid3[i], i+1+lrow)
      } else {
        if(oid3[i+1]==oid3[i]){
          newtriples2 <- c(newtriples2, i+lrow, oid3[i+1], i+1+lrow)
        } else {newtriples2 <- c(newtriples2, i+lrow, oid3[i+1], i+1+lrow, oid3[i], oid3[i+1], i+lrow)}
      } 
    }
    # if(length(which(newtriples2==0))>0){
    #   print(i)
    # }
  }
  
  newtriples2[which(newtriples2==(nrow(newcoords)+1))] <- nrow(subshape2$coords)+1
  newlist2                                             <- list(coords=newcoords, triples=c(subshape2$triples, newtriples2))
  shapeR                                               <- as.face3d(newlist2)
  
  return(list(shapeL=shapeL , shapeR=shapeR))
}

shape.curve <- function(fish, b1, b2, lk1, lk2, dista=15){
  boundary <- c(b1, b2)
  
  lmk1 <- lk1
  lmk2 <- lk2
  
  rng   <- sqrt(sum((lmk2 - lmk1)^2))
  bndry <- boundary * rng
  unit  <- (lmk2 - lmk1) / rng
  prjn  <- c(sweep(fish$coords,  2,  lmk1) %*% unit)
  near  <- function(x,  pts,  distance) sqrt(.rowSums((sweep(pts,  2,  x))^2,  nrow(pts),  3)) < distance
  ind1  <- (prjn > - bndry[1]) & (prjn < rng + bndry[1])
  prjn2 <- outer(prjn,  unit)
  ind2  <- apply((sweep(fish$coords,  2,  lmk1) - prjn2)^2,  1,  function(x) sqrt(sum(x)) < bndry[2])
  shape <- subset.face3d(fish,  ind1 & ind2,  remove.singles = TRUE)
  
  shape <- index.face3d(shape,  distance = dista,  overwrite = TRUE)
  path  <- connected.face3d(shape)
  shape <- subset.face3d(shape, path==1)
}

pospts <- function(x, shp){
  ind <- closest.face3d(x, shp)$id
  if(length(ind)==1){
    return(list(id=ind, pos="vertex"))
  } else if(length(ind)==2){
    return(list(id=ind, pos="edge"))
  } else {
    vec_1 <- shp$coords[ind[1], ]-x
    vec_2 <- shp$coords[ind[2], ]-x
    vec_3 <- shp$coords[ind[3], ]-x
    vec_1 <- vec_1/Eucdist(vec_1, c(0, 0, 0))
    vec_2 <- vec_2/Eucdist(vec_2, c(0, 0, 0))
    vec_3 <- vec_3/Eucdist(vec_3, c(0, 0, 0))
    if(vec_1%*%vec_2 < -0.99999){
      return(list(id=ind[1:2], pos="edge"))
    } else if (vec_1%*%vec_3 < -0.99999){
      return(list(id=ind[c(1, 3)], pos="edge"))
    } else if (vec_2%*%vec_3 < -0.99999) {
      return(list(id=ind[2:3], pos="edge"))
    } else {return(list(id=closest.face3d(x, shp)$id, pos="inside"))}
  }
} 

Eucdist    <- function(x1, x2){sqrt(sum((x1-x2)^2))}

intri.fun <- function(x ,mat ,edge=F){
  cn <- FALSE
  if(edge){
    for(i in 1:nrow(mat)){
      if(length(which(x%in%mat[i,]))==2){
        cn <- TRUE
        break
      }
    }
  } else {
    for(i in 1:nrow(mat)){
      if(all(x%in%mat[i,])){
        cn <- TRUE
        break
      }
    }
  }
  return(cn)
}

check.edges <- function(pt1, pt2, triangles, shape){
  conn <- TRUE
  ##################check in-between edges###########################################
    p1 <- pospts(pt1 ,shape)
    p2 <- pospts(pt2 ,shape)
    if(all(c(p1$pos ,p2$pos)%in%"edge")){
      if(all(unique(c(p1$id ,p2$id))%in%intersect(p1$id ,p2$id))){
        conn <- FALSE
      } else if (intri.fun(unique(c(p1$id ,p2$id)),triangles)){
        conn <- FALSE
      }
    } else if(any(c(p1$pos ,p2$pos)%in%"edge") & any(c(p1$pos ,p2$pos)%in%"vertex")){
      ide <- which(c(p1$pos ,p2$pos)%in%"edge")
      idv <- which(c(p1$pos ,p2$pos)%in%"vertex")
      if(list(p1 ,p2)[[idv]]$id%in%list(p1 ,p2)[[ide]]$id){
        conn <- FALSE
      } else if(intri.fun(unique(c(p1$id ,p2$id)),triangles)){
        conn <- FALSE
      }
    } else if((any(c(p1$pos ,p2$pos)%in%"edge") & any(c(p1$pos ,p2$pos)%in%"inside")) | (any(c(p1$pos ,p2$pos)%in%"vertex") & any(c(p1$pos ,p2$pos)%in%"inside"))){
      ide <- which.min(c(length(p1$id),length(p2$id)))
      idi <- which.max(c(length(p1$id),length(p2$id)))
      if(all(list(p1 ,p2)[[ide]]$id%in%list(p1 ,p2)[[idi]]$id)){
        conn <- FALSE
      }
    } else if(all(c(p1$pos ,p2$pos)%in%"vertex")){
      if (intri.fun(unique(c(p1$id ,p2$id)),triangles ,edge = T)){
        conn <- FALSE
      }
    } else if(all(c(p1$pos ,p2$pos)%in%"inside")){
      if (all(p1$id%in%p2$id) & all(p2$id%in%p1$id)){
        conn <- FALSE
      }
    } 
    return(conn)
}

antic.fun <- function(pt1 ,pt2 ,pt3 ,shape){
  #requires normal vectors on the shape
  dir <- pt2-pt1
  ppt <- pospts(pt1 ,shape)
  if(ppt$pos=="vertex"){
    nor <- shape$normals[ppt$id,]
  } else {
    nor <- colMeans(shape$normals[ppt$id,])
  }
  axi <- crossproduct(nor ,dir)
  pg2 <- pt2%*%t(axi)
  pg3 <- pt3%*%t(axi)
  if(pg2<pg3){
    return(c(1 ,2 ,3))
  } else {return(c(1 ,3 ,2))}
}
