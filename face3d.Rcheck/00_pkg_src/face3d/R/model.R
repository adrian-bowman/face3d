model.face3d <- function(face, type.face = "standard", sex = "female", graphics = FALSE){
  #creating all curves
 # if (type.face == "standard"){ #create for standard face
  
    face <- facecurves.face3d(face, "nasal root left")
    face <- facecurves.face3d(face, "nasal root right")
    face <- facecurves.face3d(face, "nasal boundary left")
    face <- facecurves.face3d(face, "nasal boundary right")
    face <- facecurves.face3d(face, "nasal bridge left")
    face <- facecurves.face3d(face, "nasal bridge right")
    face <- facecurves.face3d(face, "mid-line philtral")
    face <- facecurves.face3d(face, "nasal base left")
    face <- facecurves.face3d(face, "nasal base right")
    face <- facecurves.face3d(face, "mid-line lip")
    face <- facecurves.face3d(face, "lower lip right")
    face <- facecurves.face3d(face, "lower lip left")
    face <- facecurves.face3d(face, "upper lip right")
    face <- facecurves.face3d(face, "upper lip left")
    face <- facecurves.face3d(face, "philtrum ridge right")
    face <- facecurves.face3d(face, "philtrum ridge left")
    face <- facecurves.face3d(face, "philtrum lip right")
    face <- facecurves.face3d(face, "philtrum lip left")
    face <- facecurves.face3d(face, "mid-line nasal profile")
    face <- facecurves.face3d(face, "mid-line nasal root")
    face <- facecurves.face3d(face, "mid-line columella")
    face <- facecurves.face3d(face, "mid-line upper-lip")
    face <- facecurves.face3d(face, "mid-line bottom lip")
    face <- facecurves.face3d(face, "mid-line mentolabial")
    face <- facecurves.face3d(face, "mid-line chin")
    face <- facecurves.face3d(face, "nasolabial right")
    face <- facecurves.face3d(face, "nasolabial left")
    face <- facecurves.face3d(face, "cheek-lip left")
    face <- facecurves.face3d(face, "cheek-lip right")
    face <- facecurves.face3d(face, "cheek-eye left")
    face <- facecurves.face3d(face, "cheek-eye right")
    face <- facecurves.face3d(face, "cheek-nose left")
    face <- facecurves.face3d(face, "cheek-nose right")
   # face <- facecurves.face3d(face, "mandible right")
   # face <- facecurves.face3d(face, "mandible left")
   # face <- facecurves.face3d(face, "Mandible")
    #face <- facecurves.face3d(face, "upper eye socket right")
    #face <- facecurves.face3d(face, "upper eye socket left")
    face <- facecurves.face3d(face, "lower eye socket right")
    face <- facecurves.face3d(face, "lower eye socket left")
    face <- facecurves.face3d(face, "brow ridge left")
    face <- facecurves.face3d(face, "brow ridge right")
  # 
  # } else if (type.face == "open.mouth") { #create for open mouth face
  #   
  #   face <- facecurves.face3d(face, "nasal root left")
  #   face <- facecurves.face3d(face, "nasal root right")
  #   face <- facecurves.face3d(face, "nasal boundary left")
  #   face <- facecurves.face3d(face, "nasal boundary right")
  #   face <- facecurves.face3d(face, "nasal bridge left")
  #   face <- facecurves.face3d(face, "nasal bridge right")
  #   face <- facecurves.face3d(face, "mid-line philtral")
  #   face <- facecurves.face3d(face, "nasal base left")
  #   face <- facecurves.face3d(face, "nasal base right")
  #   face <- facecurves.face3d(face, "OM mid-line lip")
  #   face <- facecurves.face3d(face, "OM upper lip right")
  #   face <- facecurves.face3d(face, "OM upper lip left")
  #   face <- facecurves.face3d(face, "OM lower lip right")
  #   face <- facecurves.face3d(face, "OM lower lip left")
  #   
  #   X    <- face$curves[grep("OM mid-line lip", rownames(face$curves)),]
  #   ind  <- (substr(rownames(face$curves),1,nchar("OM mid-line lip"))=="OM mid-line lip")
  #   face$curves <- face$curves[!ind, ]
  #   rownames(X) <- paste("mid-line lip", 1:nrow(X))
  #   face$curves <- rbind(face$curves,X)
  #   
  #   X    <- face$curves[grep("OM upper lip right", rownames(face$curves)),]
  #   ind  <- (substr(rownames(face$curves),1,nchar("OM upper lip right"))=="OM upper lip right")
  #   face$curves <- face$curves[!ind, ]
  #   rownames(X) <- paste("upper lip right", 1:nrow(X))
  #   face$curves <- rbind(face$curves,X)
  #   
  #   X    <- face$curves[grep("OM upper lip left", rownames(face$curves)),]
  #   ind  <- (substr(rownames(face$curves),1,nchar("OM upper lip left"))=="OM upper lip left")
  #   face$curves <- face$curves[!ind, ]
  #   rownames(X) <- paste("upper lip left", 1:nrow(X))
  #   face$curves <- rbind(face$curves,X)
  #   
  #   X    <- face$curves[grep("OM lower lip right", rownames(face$curves)),]
  #   ind  <- (substr(rownames(face$curves),1,nchar("OM lower lip right"))=="OM lower lip right")
  #   face$curves <- face$curves[!ind, ]
  #   rownames(X) <- paste("lower lip right", 1:nrow(X))
  #   face$curves <- rbind(face$curves,X)
  #   
  #   X    <- face$curves[grep("OM lower lip left", rownames(face$curves)),]
  #   ind  <- (substr(rownames(face$curves),1,nchar("OM lower lip left"))=="OM lower lip left")
  #   face$curves <- face$curves[!ind, ]
  #   rownames(X) <- paste("lower lip left", 1:nrow(X))
  #   face$curves <- rbind(face$curves,X)
  #   
  #   face <- facecurves.face3d(face, "philtrum ridge right")
  #   face <- facecurves.face3d(face, "philtrum ridge left")
  #   face <- facecurves.face3d(face, "philtrum lip right")
  #   face <- facecurves.face3d(face, "philtrum lip left")
  #   face <- facecurves.face3d(face, "mid-line nasal profile")
  #   face <- facecurves.face3d(face, "mid-line nasal root")
  #   face <- facecurves.face3d(face, "mid-line columella")
  #   face <- facecurves.face3d(face, "mid-line upper-lip")
  #   face <- facecurves.face3d(face, "mid-line bottom lip")
  #   face <- facecurves.face3d(face, "mid-line mentolabial")
  #   face <- facecurves.face3d(face, "mid-line chin")
  #   face <- facecurves.face3d(face, "nasolabial right")
  #   face <- facecurves.face3d(face, "nasolabial left")
  #   face <- facecurves.face3d(face, "cheek-lip left")
  #   face <- facecurves.face3d(face, "cheek-lip right")
  #   face <- facecurves.face3d(face, "cheek-eye left")
  #   face <- facecurves.face3d(face, "cheek-eye right")
  #   face <- facecurves.face3d(face, "cheek-nose left")
  #   face <- facecurves.face3d(face, "cheek-nose right")
  #   face <- facecurves.face3d(face, "mandible right")
  #   face <- facecurves.face3d(face, "mandible left")
  #   face <- facecurves.face3d(face, "Mandible")
  #   #face <- facecurves.face3d(face, "upper eye socket right")
  #   #face <- facecurves.face3d(face, "upper eye socket left")
  #   face <- facecurves.face3d(face, "lower eye socket right")
  #   face <- facecurves.face3d(face, "lower eye socket left")
  #   face <- facecurves.face3d(face, "brow ridge left")
  #   face <- facecurves.face3d(face, "brow ridge right")
    
#   }else if (type.face == "child"){  #for child face
#     
#     
#   }
#   
#   
# #cutting face down at sides etc. 
#   curves.full      <- face$curves
#   face$curves.full <- curves.full
#   face             <- cutface.face3d(face)
# #sliding and resampling curves 
#   #BUT instead of sliding and resampling, we are just resampling the curves- Done Here. with a new argument
#   if (sex == "female"){ curves.new <- slidingcurves.face3d(template.female, face, type = "resample")
#   }else{                curves.new <- slidingcurves.face3d(template.male, face, type = "resample")}
#   face$curves         <- curves.new
#   
  
# #create the mesh 
#   face <- mesh.face3d(face, "mid-face right")
#   face <- mesh.face3d(face, "mid-face left")
#   face <- mesh.face3d(face, "upper mid face right")
#   face <- mesh.face3d(face, "upper mid face left")
#   face <- mesh.face3d(face, "upper lip right")
#   face <- mesh.face3d(face, "upper lip left")
#   face <- mesh.face3d(face, "lower lip right")
#   face <- mesh.face3d(face, "lower lip left")
#   face <- mesh.face3d(face, "philtrum right")
#   face <- mesh.face3d(face, "philtrum left")
#   face <- mesh.face3d(face, "lower face right")
#   face <- mesh.face3d(face, "lower face left")
#   face <- mesh.face3d(face, "nose right")
#   face <- mesh.face3d(face, "nose left")
#   face <- mesh.face3d(face, "upper face right 1")
#   face <- mesh.face3d(face, "upper face left 1")
#   #face <- mesh.face3d(face, "upper face right 2")
#   #face <- mesh.face3d(face, "upper face left 2")
#   face <- mesh.face3d(face, "upper face right 3")
#   face <- mesh.face3d(face, "upper face left 3")
#   face$mesh.full  <- face$mesh
#   
# #sliding and resampling mesh 
#   #BUT instead of sliding and resampling, we are just resampling the mesh- Done Here. with a new argument
#  # here need to fix 1156 resample versus all else#############ß
#   if (sex == "female"){ mesh.new <- slidingmesh.face3d(template.female, face, type = "resample")
#    }else{               mesh.new <- slidingmesh.face3d(template.male, face, type = "resample")}
# 
#   face$mesh       <- mesh.new
# 
#   #prepare the mesh
#   
#    face             <- prepare_mesh.face3d(face)
   invisible(face)
  
}
  

