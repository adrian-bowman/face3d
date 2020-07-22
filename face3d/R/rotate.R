translate.face3d <- function(shape, shift, exclude = NULL) {

    if (!any(c("face3d", "matrix") %in% class(shape)))
        stop("shape must be a matrix or a face3d object.")
    if (missing(shift))
        stop("shift must be supplied.")
    if (!is.null(exclude) && !is.character(exclude))
        stop("exclude must be a character vector")
    
    translate.fn     <- function(x) sweep(x, 2, shift, "+")

    if (is.face3d(shape)) {
        ind        <- sapply(shape, function(x) is.matrix(x) && ncol(x) == 3)
        ind        <- which(ind)
        exclude    <- c("triangles", exclude)
        ind        <- ind[-match(exclude, names(ind))]
        shape[ind] <- lapply(shape[ind], translate.fn) 
    }
    else
        shape <- translate.fn(shape)

    return(shape)
}

rotate.face3d <- function(shape, angle, raxis, centre = rep(0, 3), center,
                              landmarks = c("sn", "n", "exR", "exL"), rotation = "coronal",
                              exclude = NULL) {
        
    if (!any(c("face3d", "matrix") %in% class(shape)))
        stop("shape must be a matrix or a face3d object.")
    if (!missing(center)) centre <- center
    if (!is.null(exclude) && !is.character(exclude))
        stop("exclude must be a character vector")

    # Nominated angle
    if (!missing(angle)) {
        if (missing(raxis))
            stop("raxis must be specified when an angle is specified.")
        else if (is.face3d(shape)) {
            ind     <- sapply(shape, function(x) is.matrix(x) && ncol(x) == 3)
            ind     <- which(ind)
            exclude <- c("triangles", exclude)
            ind     <- ind[-match(exclude, names(ind))]
            fn      <- function(x) {
                x <- sweep(x, 2, centre)
                x <- rgl::rotate3d(x, angle, raxis[1] , raxis[2], raxis[3])
                x <- sweep(x, 2, centre, "+")
            }
            shape[ind] <- lapply(shape[ind], fn) 
        }
        else {
            shape <- sweep(shape, 2, centre)
            shape <- rgl::rotate3d(shape, angle, raxis[1], raxis[2], raxis[3])
            shape <- sweep(shape, 2, centre, "+")
        }
        return(shape)
    }
    
    # Rotation to landmarks
    id.lndms      <- landmarks
    lndms.two.mid <- id.lndms[1:2]
    lndms.two.lat <- id.lndms[3:4]
    shape1 <- shape
    if ("face3d" %in% class(shape1)) shape <- shape$landmarks

    if (rotation == "coronal" | rotation == "sagittal sinister" | rotation == "sagittal dexter" | rotation == "transversal"){
       angles <- numeric(3)
       ROT.shape <- shape
       ROTlndms <- ROT.shape[lndms.two.mid,]
       ROTlndms.vec <- apply(ROTlndms,2,diff)
       ANGLE.xz.mid <- as.numeric(atan2(ROTlndms.vec[2],ROTlndms.vec[1]))
       if (ANGLE.xz.mid <= pi/2) ANGLE.xz.mid <- -1*(pi/2 - ANGLE.xz.mid) #clockwise
       if (ANGLE.xz.mid > pi/2)  ANGLE.xz.mid <- ANGLE.xz.mid - pi/2 #clockwise
       angles[1] <- ANGLE.xz.mid
       ROT.shape <- rgl::rotate3d(ROT.shape,ANGLE.xz.mid,0,0,1)
       ROTlndms <- ROT.shape[lndms.two.mid,]
       ROTlndms.vec <- apply(ROTlndms,2,diff)
       ANGLE.xz.mid <- as.numeric(atan2(ROTlndms.vec[2],ROTlndms.vec[3]))
       if (ANGLE.xz.mid <=  pi/2) ANGLE.xz.mid <- (pi/2 - ANGLE.xz.mid) #anticlockwise
       if (ANGLE.xz.mid > pi/2)  ANGLE.xz.mid <- -1*(ANGLE.xz.mid - pi/2) #anticlockwise
       angles[2] <- ANGLE.xz.mid
       ROT.shape <- rgl::rotate3d(ROT.shape,ANGLE.xz.mid,1,0,0)
       ROTlndms <- ROT.shape[lndms.two.lat,]
       ROTlndms.vec <- apply(ROTlndms,2,diff)
       ANGLE.xz.lat <- -1*as.numeric(atan2(ROTlndms.vec[3],ROTlndms.vec[1]))
       angles[3] <- ANGLE.xz.lat
       ROT.shape <- rgl::rotate3d(ROT.shape,ANGLE.xz.lat,0,1,0)
       names(angles) <- c("coronal z","coronal x","coronal y")
       ROT <- list(ROT.shape = ROT.shape, angles = angles)
       if ("face3d" %in% class(shape1)) {
          shape1$vertices <- rgl::rotate3d(shape1$vertices,angles[1],0,0,1)
          shape1$vertices <- rgl::rotate3d(shape1$vertices,angles[2],1,0,0)
          shape1$vertices <- rgl::rotate3d(shape1$vertices,angles[3],0,1,0)
          shape1$landmarks <- ROT.shape
          if ("curves" %in% names(shape1)) {
             shape1$curves <- rgl::rotate3d(shape1$curves,angles[1],0,0,1)
             shape1$curves <- rgl::rotate3d(shape1$curves,angles[2],1,0,0)
             shape1$curves <- rgl::rotate3d(shape1$curves,angles[3],0,1,0)
          }
          if ("mesh" %in% names(shape1)) {
             shape1$mesh <- rgl::rotate3d(shape1$mesh,angles[1],0,0,1)
             shape1$mesh <- rgl::rotate3d(shape1$mesh,angles[2],1,0,0)
             shape1$mesh <- rgl::rotate3d(shape1$mesh,angles[3],0,1,0)
             }
       }

    if (rotation == "sagittal sinister"){
       ROT.shape <- rgl::rotate3d(ROT.shape,pi/2,0,1,0)
       ROT <- list(ROT.shape = ROT.shape, angles = angles)
         if ("face3d" %in% class(shape1)) {
               shape1$vertices <- rgl::rotate3d(shape1$vertices,pi/2,0,1,0)
               shape1$landmarks <- ROT.shape
               if ("curves" %in% names(shape1)){
               shape1$curves <- rgl::rotate3d(shape1$curves,pi/2,0,1,0)
               }
               if ("mesh" %in% names(shape1)){
               shape1$mesh <- rgl::rotate3d(shape1$mesh,pi/2,0,1,0)
               }
               }
       }

    if (rotation == "sagittal dexter") {
       ROT.shape <- rgl::rotate3d(ROT.shape,-pi/2,0,1,0)
       ROT <- list(ROT.shape = ROT.shape, angles = angles)
         if ("face3d" %in% class(shape1)) {
               shape1$vertices <- rgl::rotate3d(shape1$vertices,-pi/2,0,1,0)
               shape1$landmarks <- ROT.shape
               if ("curves" %in% names(shape1)){
               shape1$curves <- rgl::rotate3d(shape1$curves,-pi/2,0,1,0)
               }
               if ("mesh" %in% names(shape1)){
               shape1$mesh <- rgl::rotate3d(shape1$mesh,-pi/2,0,1,0)
               }
               }
    }

    if (rotation == "transversal") {
       ROT.shape <- rgl::rotate3d(ROT.shape,pi/2,1,0,0)
       ROT <- list(ROT.shape = ROT.shape, angles = angles)
        if ("face3d" %in% class(shape1)) {
               shape1$vertices <- rgl::rotate3d(shape1$vertices,-pi/2,1,0,0)
               shape1$landmarks <- ROT.shape
               if ("curves" %in% names(shape1)){
               shape1$curves <- rgl::rotate3d(shape1$curves,-pi/2,1,0,0)
               }
               if ("mesh" %in% names(shape1)){
               shape1$mesh <- rgl::rotate3d(shape1$mesh,-pi/2,1,0,0)
               }
         }
    }
    
    if ("face3d" %in% class(shape1))
       return(shape1)
    else if ("matrix" %in% class(shape1)) {
       ROT <- list(ROT.shape = ROT.shape, angles = angles)
       return(ROT)
       }
    }
}
