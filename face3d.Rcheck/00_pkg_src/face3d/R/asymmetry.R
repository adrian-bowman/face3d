#     Asymmetry

asymmetry.face3d <- function(X, sides = c("left", "right"), reference.size,
                             model.mesh = FALSE, region, match.subset,
                             controls, labels, labels.size = 1) {
  
   if (is.null(rownames(X))) stop("X has no rownames.")
   missing.region   <- missing(region)
   missing.controls <- missing(controls)
   if (!missing.controls && is.list(controls) &&
               (all(c("asymmetry", "reference.size") %in% names(controls)))) {
      reference.size <- controls$reference.size
      if (missing.region & ("region" %in% names(controls))) {
         region         <- controls$region
         missing.region <- FALSE
      }
      controls <- controls$asymmetry
   }
   if (missing(match.subset)) match.subset <- 1:nrow(X)
   rnms <- rownames(X)
   if (missing.region) {
      region <- list("global"     = 1:nrow(X),
                     "lower face" = c(grep("lower face", rnms), grep("chin", rnms)),
                     "lips"       = grep("lip", rnms),
                     "philtrum"   = grep("philtr", rnms),
                     "nose"       = c(grep("nose", rnms), grep("nasal", rnms), grep("columella", rnms)),
                     "mid-face"   = grep("mid-face", rnms),
                     "brow"       = grep("upper face", rnms),
                     "upper face" = grep("upper mid face", rnms))
   }
   else {
      if (length(region) == 1) {
         region      <- list(region)
         names(list) <- deparse(substitute(region))
      }
      else if (!is.list(region))
         stop("'region' should be a list.")
   }
   if (!("global" %in% names(region))) region$global <- 1:nrow(X)

   if (length(dim(X)) == 2)
      X <- array(c(X), dim = c(dim(X), 1), dimnames = list(rownames(X), NULL, NULL))
   n <- dim(X)[3]
   missing.labels <- missing(labels)
   if (missing.labels) labels <- 1:(dim(X)[3])
   
   asymmetry <- matrix(nrow = n, ncol = length(names(region)))
   size      <- matrix(nrow = n, ncol = length(names(region)))
   for (i in 1:n) {

      rownames(X) <- rnms
      if (model.mesh) areas <- area.face3d(as.face3d(X[ , , i], model.mesh = TRUE))$points
      rownames(X) <- sub("LLL", "LL left",  rownames(X))
      rownames(X) <- sub("LLR", "LL right", rownames(X))
      rownames(X) <- sub("ULL", "UL left",  rownames(X))
      rownames(X) <- sub("ULR", "UL right", rownames(X))
      rownames(X) <- sub("MLL", "ML left",  rownames(X))
      rownames(X) <- sub("MLR", "ML right", rownames(X))
      
      for (j in 1:length(region)) {
         Xreg   <- X[region[[j]], , i]
         nms    <- rownames(Xreg)
         ind.l  <- grep(sides[1], nms)
         ind.r  <- grep(sides[2], nms)
         nms.m  <- sub(sides[1], "", nms)
         nms.m  <- sub(sides[2], "", nms.m)
         ind.rl <- match(nms.m[ind.r], nms.m[ind.l])
         ind.lr <- match(nms.m[ind.l], nms.m[ind.r])

         flag.r <- which(is.na(ind.rl))
         flag.l <- which(is.na(ind.lr))
         if (length(flag.l) + length(flag.r) > 0) {
            cat(paste("Some", sides[1], "and", sides[2], "points do not match.\n"))
            if (length(flag.r) > 0) print(nms[ind.r][flag.r])
            if (length(flag.l) > 0) print(nms[ind.l][flag.l])
         }
         
         XR               <- Xreg
         XR[ind.l, ]      <- Xreg[ind.r[ind.lr], ]
         XR[ind.r, ]      <- Xreg[ind.l[ind.rl], ]
         XR[ , 1]         <- -XR[ , 1]

         mtch.sbst <- if (names(region)[j] == "global") match.subset else 1:nrow(Xreg)

         if (model.mesh) {
            areas.reg        <- areas[region[[j]]]
            areas.R          <- areas.reg
            areas.R[ind.l]   <- areas.reg[ind.r[ind.lr]]
            areas.R[ind.r]   <- areas.reg[ind.l[ind.rl]]
            areas.reg        <- (areas.reg + areas.R) / 2
            XR               <- opa.face3d(XR, Xreg, scale = FALSE, match.subset = mtch.sbst,
                                           weights = areas.reg, model.mesh = TRUE)
         }
         else
            XR <- opa.face3d(XR, Xreg, scale = FALSE, match.subset = mtch.sbst)

         XR              <- (Xreg + XR) / 2

         if (model.mesh) {
            size[i, j]      <- sum(areas.reg)
            asymmetry[i, j] <- sqrt(sum(rowSums((Xreg - XR)^2) * areas.reg) / size[i, j])
         }
         else {
            centre          <- apply(Xreg, 2, mean)
            size[i, j]      <- sqrt(sum(sweep(Xreg, 2, centre)^2) / nrow(Xreg))
            asymmetry[i, j] <- sqrt(sum((Xreg - XR)^2) / nrow(Xreg)) / size[i, j]
         }
         if ((names(region)[j] == "global") & (n == 1)) symmetric.shape <- XR
      }
   }
   colnames(asymmetry) <- names(region)
   colnames(size)      <- names(region)
   if (missing(reference.size)) reference.size <- apply(size, 2, mean)
   if (!model.mesh) asymmetry <- sweep(asymmetry, 2, reference.size, "*")

   if (!missing.controls) {
      ngrid <- 500
      x     <- seq(1/ngrid, 1.1 * max(c(controls)), length = ngrid)
      dens  <- apply(controls, 2, function(z) density(z, n = ngrid, from = min(x), to = max(x))$y)
      dat   <- data.frame(x = rep(x, ncol(controls)), dens = c(dens),
                          gp = factor(rep(colnames(controls), each = length(x)), levels = colnames(controls)))
      reg   <- factor(rep(colnames(asymmetry), each = nrow(asymmetry)), levels = colnames(asymmetry))
      dfrm  <- data.frame(asy = c(asymmetry), reg = reg,
                          clr = factor(1 + rep(1:n, length(colnames(asymmetry)))),
                          label = rep(labels, n))
      ht    <- 0.6
      wd    <- 1.2 * ht / 2
      plt   <- ggplot(dat, aes(x = x, y = gp)) + theme_classic() +
                  geom_tile(aes(fill = dens), height = ht) +
                  scale_fill_gradient(low = "white", high = "black") +
                  geom_segment(aes(x = asy, xend = asy,
                                   y = as.numeric(reg) + wd, yend = as.numeric(reg) - wd,
                                   col = clr), data = dfrm, size = 1) +
                  theme(axis.title.y = element_blank(), legend.position  = "none",
                        axis.line.y = element_blank(), axis.ticks.y = element_blank()) +
                  xlab("asymmetry (mm.)")
      if (!is.null(labels) & !any(is.na(labels)))
         plt <- plt + geom_label(aes(x = asy, y = as.numeric(reg) + wd, label = label, col = clr),
                                       data = dfrm, size = labels.size)
      print(plt)
   }
   
   result    <- list(asymmetry = asymmetry, region = region, size = size, reference.size = reference.size)
   if (n == 1) {
      asymmetry <- c(asymmetry)
      size      <- c(size)
      result$symmetric <- symmetric.shape
   }
   if (!missing(controls)) result$plot <- plt
   
   invisible(result)
}
