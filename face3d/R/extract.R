"extract.face3d" <- function(model) {

   if (!requireNamespace("prodlim", quietly = TRUE)) stop("the prodlim package is required.")

   mesh <- if (class(model) == "face3d") face$meshes else model
           mesh.nms         <- c( "mid-line columella 1",        "philtrum left 1 6" ,               "philtrum right 1 6",     
                                             "mid-line columella 8" ,       "mid-line nasal profile 1",     "mid-line nasal root 1",    
                                             "upper mid face left 12 5",  "upper mid face right 12 5", "nose left 1 5",           
                                             "nose right 1 5",                  "mid-line philtral 6",             "upper lip left 6 1",       
                                             "upper lip right 6 1",            "LLL 1",                               "LLR 1",                    
                                             "mid-line bottom lip 1" ,       "mid-line bottom lip 3",       "mid-line chin 4" ,         
                                             "mid-line chin 1")
           lmks.names        <- c( "pn",   "acL" , "acR",  "sn",   "se",   "n",    "exL",  "exR",  "enL",  "enR", 
                                               "ls",   "cphL", "cphR", "chL",  "chR",  "st",   "li",   "sl",   "gn" )  
           lmks                     <- mesh[mesh.nms, ]
           rownames(lmks)  <- lmks.names
  
   if (class(model) == "face3d"){
               curves                                     <- data.frame(face$curves)
               mesh                                       <- data.frame(face$meshes)
               ind                                           <- row.match(curves, mesh)
               cfc                                           <- which(is.na(ind) == "TRUE")
               curves.from.curves                  <- curves[cfc, ]
               cfm                                          <- which(!is.na(ind) == "TRUE")
               curves.from.mesh                   <- mesh[ind[cfm], ]   
               curves.names                          <- rownames(curves[cfm, ])      
               rownames(curves.from.mesh) <- curves.names                           
               curves.extracted                      <- rbind(curves.from.curves, curves.from.mesh)
         }
     else
        curves.extracted <- NA
 
      invisible(list(curves= curves.extracted, lmks = lmks))
}
