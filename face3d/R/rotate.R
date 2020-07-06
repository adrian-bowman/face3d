"rotate.face3d" <- function(shape, id.lndms = c("sn","n","exR","exL"), rotation = "coronal") {
    # if (missing(rotation))  {
    #    cat("Type of rotation not specified.\n Those currently supported are: \n coronal (frontal) \n sagittal sinister (lateral sinister) \n sagittal sinister (lateral sinister) \n transversal (vertical/basal)\n")
    #    return()
    #    }
    lndms.two.mid <- id.lndms[1:2]
    lndms.two.lat <- id.lndms[3:4]
    shape1 <- shape
    if ("face3d" %in% class(shape1)) {
                  shape <- shape$landmarks
                  }
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
    if (rotation == "sagittal dexter"){
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

    if (rotation == "transversal"){
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
        #class(shape1) <- "face3d"
       if ("face3d" %in% class(shape1)) {
           #class(shape1) <- "face3d"
        return(shape1)
       }
    if ("matrix" %in% class(shape1)) {
       ROT <- list(ROT.shape = ROT.shape, angles = angles)
       return(ROT)
       }
    }
}
