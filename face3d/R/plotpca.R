plotpca.face3d <- function(x, group, high, display, npc, pc, p,
                            magnification = 2, template, template.match, range.colour,
                            n.animation = 12, image.folder, ...) {
   
   # Check the arguments and set appropriate defaults
   horizontal      <- TRUE
   npc.missing     <- missing(npc)
   p.missing       <- missing(p)
   group.missing   <- missing(group)
   display.missing <- missing(display)
   if(display.missing) display <- "scores"
   if ("animation" %in% display) {
      save.animation <- !missing(image.folder)
      if (!(image.folder %in% list.files()))
         stop(paste("the image folder (", image.folder, ") does not exist.", sep = ""))
   }
   if (all(c("pca", "group", "npc", "p", "weights") %in% names(x))) {
      pca                  <- x$pca
      if (npc.missing) npc <- x$npc
      if (p.missing)   p   <- x$p
      npc.missing          <- FALSE
      p.missing            <- FALSE
      if (!missing(group) && !all(group == x$group))
         stop("the group variable in x does not match the group argument.")
      group                <- x$group
      group.missing        <- FALSE
   } 
   else
      pca <- x
   if (!("weights" %in% names(pca))) stop("the pca information does not have a 'weights' component.")

   # if ("ci" %in% display) {
   #    if (group.missing) stop("ci has been requested but group has not been set.")
   #    if (length(unique(group)) != 2)
   #       stop("ci has been requested but there are more than two groups specified in group.")
   # }
   if (group.missing) group <- rep(1, nrow(pca$scores))
   if (npc.missing)  npc  <- if ("shapes" %in% display) 1
                             else min(which(cumsum(pca$percent) > 80))
   if (any(c("shapes", "x", "y", "z", "normal", "s-normal", "animation") %in% display)) {
      template.missing       <- missing(template)
      template.match.missing <- missing(template.match)
      if (template.missing)
         stop("a template must be specified when display is set to `shapes' or `animation'.")
      if (template.match.missing) {
         if      (nrow(pca$mean) == nrow(template$mesh))   template.match <- template$mesh
         else if (nrow(pca$mean) == nrow(template$curves)) template.match <- template$curves
         else if (nrow(pca$mean) == nrow(template$lmks))   template.match <- template$lmks
         else stop("the dimensionality of the shapes does not match the mesh, curves or lmks components of the template.")
      }
   }
   range.colour.missing <- missing(range.colour)
   
   result <- list()
   
   # Ensure that the PCs tend to have higher values for the nominated `high' group
   if (!group.missing) {
      if (missing(high)) high <- levels(factor(as.character(group)))[2]
      for (ipc in 1:ncol(pca$scores)) {
         mns                <- c(tapply(pca$scores[ , ipc], group, mean))
         ind                <- match(high, names(mns))
         sgn                <- if (order(mns)[ind] < order(-mns)[ind]) -1 else 1
         pca$scores[ , ipc] <- sgn * pca$scores[ , ipc]
         pca$evecs[  , ipc] <- sgn * pca$evecs[  , ipc]
      }
   }

   if ("scores" %in% display) {
      clr    <- rainbow_hcl(length(unique(group)))
      pclev  <- paste("PC ", 1:npc, " ", round(pca$pct[1:npc], 1), "%", sep = "")
      if (!p.missing) pclev <- paste(pclev, paste("p=", round(p[-1], 3), sep = ""))
      # if (horizontal) pclev <- rev(pclev)
      dfrm   <- data.frame(scores = c(pca$scores[ , 1:npc]),
                           group  = as.factor(rep(group, npc)),
                           pc     = factor(rep(pclev, each = nrow(pca$scores)), levels = pclev))
      gg.scores <- ggplot(dfrm, aes(group, scores, fill = group)) +
                   geom_boxplot() +
                   scale_fill_manual(values = clr) +
                   theme(legend.position = "right")
      if (horizontal)
         gg.scores <- gg.scores + facet_grid(. ~ pc, labeller = label_wrap_gen(width = 6)) +
                      # theme(,axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1))
                      theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                            axis.ticks.x = element_blank())
       else
         gg.scores <- gg.scores + facet_grid(pc ~ ., labeller = label_wrap_gen(width = 6)) +
                      coord_flip() + theme(axis.title.y = element_blank())
      if (group.missing) {
         if (horizontal)
            gg.scores <- gg.scores + theme(axis.title.x = element_blank(),
                                           axis.text.x  = element_blank(),
                                           axis.ticks.x = element_blank())
         else
            gg.scores <- gg.scores + theme(axis.title.y = element_blank(),
                                           axis.text.y  = element_blank(),
                                           axis.ticks.y = element_blank())
         gg.scores <- gg.scores + theme(legend.position = "none")
      }
      print(gg.scores)
      result <- gg.scores
   }

   # if ("ci" %in% display) {
   #    ci     <- t(apply(pca$scores[ , 1:npc], 2, function(x) -t.test(x ~ group, conf.level = 1 - 0.05/npc)$conf.int))
   #    ci.clr <- clr[1 + as.numeric(ci[ , 1] * ci[ , 2] > 0)]
   #    dfrm1  <- data.frame(x1 = ci[ , 1], x2 = ci[ , 2],  pc = 1:npc)
   #    gg.ci  <- ggplot(dfrm1, aes(range(x1, x2), range(1:npc))) +
   #                      geom_segment(aes(x = x1, y = pc, xend = x2, yend = pc), size = 1.5, col = ci.clr) +
   #                      geom_hline(yintercept = 1:(npc - 1) + 0.5, col = "white", size = 1.5) +
   #                      geom_vline(xintercept = 0, linetype = "dashed", col = clr[1], size = 1.5) +
   #                      scale_y_continuous(breaks = 1:npc) + 
   #                      ylab("PC") + xlab("CI")
   #    # gg.ci <- gg.ci + coord_flip()
   #    print(gg.ci)
   #    result$ci <- gg.ci
   # }
   
   if (any(c("shapes", "x", "y", "z", "normal", "s-normal") %in% display)) {
      if (missing(pc)) stop("the indices of the pcs to be displayed must be specified in the variable pc.")
      adj <- matrix(rep(0, nrow(pca$evecs)), ncol = 3)
      for (i in 1:length(pc))
         adj <- adj + matrix(magnification * pca$sd[pc[i]] * pca$evecs[ , pc[i]], ncol = 3)
      adj    <- sweep(adj, 1, sqrt(pca$weights), "/") / sqrt(length(pc))
      shape1 <- warp.face3d(template.match, pca$mean - adj, template)
      shape2 <- warp.face3d(template.match, pca$mean + adj, template)
      if ("shapes" %in% display) {
         clr <- rainbow_hcl(2)
         plot(shape1, col = clr[1], new = FALSE)
         plot(shape2, col = clr[2], add = TRUE)
      }
      if (any(c("x", "y", "z", "normal", "s-normal") %in% display)) {
         dst <- distance.face3d(shape1, shape2)
         dst <- if (display == "s-normal") dst$xyz * sign(dst$normal) else dst[[display]]
         if (range.colour.missing) range.colour <- range(dst, na.rm = TRUE)
         plot(shape1, col = dst, new = FALSE, range.colour = range.colour)
      }
   }

   if ("animation" %in% display) {
      if (missing(pc)) stop("the indices of the pcs to be displayed nust be specified in the variable pc.")
      cat("Preparing ...")
      wp <- list()
      mg <- c(0, magnification, -magnification)
      for (i in 1:length(mg)) {
         adj <- matrix(rep(0, nrow(pca$evecs)), ncol = 3)
         for (j in 1:length(pc))
            adj <- adj + matrix(mg[i] * pca$sd[pc[j]] * pca$evecs[ , pc[j]], ncol = 3)
         adj     <- sweep(adj, 1, sqrt(pca$weights), "/") / sqrt(length(pc))
         wp[[i]] <- warp.face3d(template.match, pca$mean + adj, template)
         cat(".")
      }
      cat(" complete\n")
      animate.face3d(list(wp[[1]], wp[[2]], wp[[1]], wp[[3]], wp[[1]]), n.animation,
                     image.folder, ...)
      # sq <- c(1, 2, 1, 3, 1)
      # plot(wp[[1]], col = "grey", new = FALSE)
      # for (i in 1:(length(sq) - 1)) {
      #    for (j in 1:n.animation) {
      #       shp        <- wp[[sq[i]]]
      #       wt         <- j / n.animation
      #       shp$coords <- (1 - wt) * wp[[sq[i]]]$coords + wt * wp[[sq[i + 1]]]$coords
      #       p3d <- rgl::par3d(skipRedraw = TRUE)
      #       rgl::pop3d()
      #       plot(shp, col = "grey", add = TRUE)
      #       rgl::par3d(p3d)
      #       if (save.animation) snapshot3d(paste(image.folder ,"/temp-", i, ".png", sep = ""))
      #    }
      # }
   }
   
   if (length(result) == 0) result <- NULL 
   invisible(result)
}
