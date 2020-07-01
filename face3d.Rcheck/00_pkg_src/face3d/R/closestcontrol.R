closestcontrol.face3d <- function(cases, gpa, pca, npc, pc.proportion = 0.8,
                                  display = c("unusual", "histogram"), labels) {
   
   if (all(c("all", "unusual") %in% display))
      display <- display[-which(display == "all")]
   if (!any(c("all", "unusual") %in% display))
      display <- c(display, "unusual")
   if (all(c("histogram", "density strip") %in% display))
      display <- display[-which(display == "density strip")]
   if (!any(c("histogram", "density strip") %in% display))
      display <- c(display, "histogram")
   if (is.matrix(cases))
      cases <- array(c(cases), dim = c(dim(cases), 1))
   if (missing(labels)) labels <- 1:(dim(cases)[3])
   
   n.cases   <- dim(cases)[3]
   n.ctrls   <- dim(gpa$rotated)[3]
   k         <- dim(gpa$rotated)[1]
   if (missing(pca)) pca <- pca.face3d(gpa)
   weights   <- if ("weights" %in% names(gpa)) gpa$weights else rep(1, k)
   
   if (missing(npc))
      npc <- min(which(cumsum(pca$percent) > 100 * pc.proportion))
   else if (npc == 1) {
      cat("npc should be at least 2 - resetting.\n")
      npc <- 2
   }
   
   # Register the controls against the control mean (to be even-handed)
   ctrls <- gpa$rotated
   # for (i in 1:n.ctrls)
   #    ctrls[ , , i] <- opa.face3d(gpa$rotated[ , , i], gpa$mean,
   #                                     weights = gpa$weights, model.mesh = TRUE)
   X            <- t(apply(ctrls, 3, c))
   X            <- sweep(X, 2, c(gpa$mean))
   X            <- sweep(X, 2, rep(sqrt(gpa$weights), 3), "*")
   scores.ctrls <- X %*% pca$evecs[ , 1:npc]

   # Register the cases against the control mean
   cases.proc <- cases
   for (i in 1:n.cases)
      cases.proc[ , , i] <- opa.face3d(cases[ , , i], gpa$mean,
                                       weights = gpa$weights, model.mesh = TRUE)
   X            <- t(apply(cases.proc, 3, c))
   X            <- sweep(X, 2, c(gpa$mean))
   X            <- sweep(X, 2, rep(sqrt(gpa$weights), 3), "*")
   scores.cases <- X %*% pca$evecs[ , 1:npc]

   # Find the closest controls
   mn        <- apply(scores.ctrls[ , 1:npc], 2, mean)
   covmat    <- diag(pca$sd[1:npc]^2)
   mhd       <- mahalanobis(scores.cases, mn, covmat)
   atyp      <- pchisq(mhd, npc)
   wt.norm   <- pmin(sqrt(qchisq(0.95, npc) / mhd), 1)
   wt.mat    <- if (n.cases == 1) matrix(wt.norm)     else diag(wt.norm)
   wt.mat1   <- if (n.cases == 1) matrix(1 - wt.norm) else diag(1 - wt.norm)
   cc.scores <- wt.mat  %*% scores.cases +
                wt.mat1 %*% matrix(rep(rep(0, npc), n.cases), byrow = TRUE, ncol = npc)
   cc.tan    <- cc.scores %*% t(pca$evecs[ , 1:npc])
   cc.tan    <- sweep(cc.tan, 2, rep(sqrt(gpa$weights), 3), "/")
   cc        <- array(0, dim = c(nrow(gpa$mean), 3, n.cases))
   for (i in 1:n.cases)
      cc[ , , i] <- gpa$mean + matrix(cc.tan[i, ], ncol = 3)

   # Plot the Mahalanobis distances in the PC space
   mhd.ctrl <- mahalanobis(scores.ctrls[ , 1:npc], mn, covmat)
   mhd.case <- mahalanobis(scores.cases, mn, covmat)

   # plot(mhd.ctrl, mhd.case)
   # abline(0, 1)

   if (any(c("histogram", "density strip") %in% display)) {
      ind      <- if ("unusual" %in% display) which(wt.norm < 1) else 1:length(mhd.case)
      nind     <- length(ind)
      clr      <- factor(as.numeric(wt.norm[ind] >= 1), levels = c("0", "1"))
      ttl      <- paste("Mahalanobis distance (npc = ", as.character(npc), ")", sep = "")
   }
   
   if ("histogram" %in% display) {
      hst      <- hist(mhd.ctrl, plot = FALSE)
      xmax     <- max(hst$breaks, mhd.case, npc + 3 * sqrt(2 * npc))
      xgrid    <- seq(0, xmax, length = 100)
      ygrid    <- dchisq(xgrid, npc) * length(mhd.ctrl) * diff(hst$breaks[1:2])
      ymax     <- max(hst$counts, ygrid)
      df.ctrl  <- data.frame(mhd.ctrl)
      df.grid  <- data.frame(xgrid, ygrid)
      line.col <- if (requireNamespace("scales", quietly = TRUE)) scales::hue_pal()(2)[2] else "blue"
      plt.pc <- ggplot(df.ctrl) +
         geom_histogram(aes(mhd.ctrl), breaks = hst$breaks, fill = "grey") +
         xlim(0, xmax) + ylim(0, ymax) + xlab(ttl) +
         geom_line(aes(x = xgrid, y = ygrid), data = df.grid, col = line.col, size = 1) +
         theme(legend.position = "none")
      if (nind > 0) {
         ord     <- order(mhd.case[ind])
         tht     <- ymax - (1:nind * ymax) / (nind + 1)
         df.case <- data.frame(x = mhd.case[ind][ord], y = rep(0, nind), yend = tht, label = labels[ind[ord]],
                               clr = clr[ord])
         plt.pc  <- plt.pc +
            geom_segment(aes(x = x, y = y,  xend = x, yend = yend, col = clr), data = df.case, size = 1)
         if (!is.null(labels) & !any(is.na(labels)))
            plt.pc <- plt.pc + geom_label(aes(x = x, y = yend, label = label, col = clr), data = df.case)
      }
   }
   else if ("density strip" %in% display) {
      ngrid <- 500
      x     <- seq(1/ngrid, 1.1 * max(mhd.ctrl), length = ngrid)
      dens  <- density(mhd.ctrl, n = ngrid, from = min(x), to = max(x))$y
      dat   <- data.frame(x, dens, gp = factor(rep("PC", ngrid)))
      ind   <- if ("unusual" %in% display) which(wt.norm < 1) else 1:length(mhd.case)
      nind  <- length(ind)
      ht    <- 0.6
      wd    <- 1.2 * ht / 2
      plt.pc <- ggplot(dat, aes(x = x, y = gp)) + theme_classic() +
         geom_tile(aes(fill = dens), height = ht) +
         scale_fill_gradient(low = "white", high = "black") +
         theme(axis.title.y = element_blank(), legend.position  = "none",
               axis.line.y = element_blank(), axis.text.y = element_blank(),
               axis.ticks.y = element_blank()) +
         xlab(ttl)
      if (nind > 0) {
         df.case <- data.frame(x = mhd.case[ind], y = rep(0, nind), label = labels[ind], clr = clr)
         plt.pc <- plt.pc +
            geom_segment(aes(x = x, xend = x, y = 1 + wd, yend = 1 - wd, col = clr),
                         data = df.case, size = 1)
         if (!is.null(labels) & !any(is.na(labels)))
            plt.pc <- plt.pc + geom_label(aes(x = x, y = 1 + wd, label = label, col = clr), data = df.case)
      }
   }
   
   # Find the case and control pca approximations
   cases.tan <- scores.cases %*% t(pca$evecs[ , 1:npc])
   cases.tan <- sweep(cases.tan, 2, rep(sqrt(gpa$weights), 3), "/")
   cases.pca <- array(0, dim = c(k, 3, n.cases))
   for (i in 1:n.cases)
      cases.pca[ , , i] <- gpa$mean + matrix(cases.tan[i, ], ncol = 3)
   ctrls.tan <- scores.ctrls[ , 1:npc] %*% t(pca$evecs[ , 1:npc])
   ctrls.tan <- sweep(ctrls.tan, 2, rep(sqrt(gpa$weights), 3), "/")
   ctrls.pca <- array(0, dim = c(k, 3, n.ctrls))
   for (i in 1:n.ctrls)
      ctrls.pca[ , , i] <- gpa$mean + matrix(ctrls.tan[i, ], ncol = 3)

   ctrls.res  <- ctrls - ctrls.pca
   cases.res  <- cases.proc - cases.pca
   ctrls.dist <- apply(ctrls.res, c(1, 3), function(x) sqrt(sum(x^2)))
   cases.dist <- apply(cases.res, c(1, 3), function(x) sqrt(sum(x^2)))
   ctrls.sd   <- apply(ctrls.dist, 1, sd)
   cases.sd   <- apply(cases.dist, 1, sd)
   ctrls.dist <- sweep(ctrls.dist, 1, ctrls.sd, "/")
   cases.dist <- sweep(cases.dist, 1, ctrls.sd, "/")
   ctrls.dev  <- apply(ctrls.dist, 2, mean)
   cases.dev  <- apply(cases.dist, 2, mean)

   threshold <- quantile(ctrls.dev, 0.95)
   scl       <- pmin(threshold / cases.dev[ind], 1)
   cc.pca    <- cc
   for (i in 1:n.cases)
      cc[ , , i] <- cc[ , , i] + scl[i] * cases.res[ , , i]
   
   if (any(c("histogram", "density strip") %in% display)) {
      ind       <- if ("unusual" %in% display) which(cases.dev > threshold) else 1:n.cases
      nind      <- length(ind)
      clr       <- factor(as.numeric(cases.dev[ind] <= threshold))
      ttl       <- paste("Residual scale (npc = ", as.character(npc), ")", sep = "")
   }
   
   # Plot the residual scales
   if ("histogram" %in% display) {
      hst         <- hist(ctrls.dev, plot = FALSE)
      xmax        <- max(hst$breaks, cases.dev)
      ymax        <- max(hst$counts)
      tht         <- ymax - (1:nind) * ymax / (nind + 1)
      df.ctrl     <- data.frame(ctrls.dev)
      plt.res     <- ggplot(df.ctrl) +
         geom_histogram(aes(ctrls.dev), breaks = hst$breaks, fill = "grey") +
         xlim(min(hst$breaks - diff(range(hst$breaks)) / 20, cases.dev[ind]), xmax) +
         xlab(ttl)
      if (nind > 0) {
         cdev        <- cases.dev[ind]
         names(cdev) <- as.character(ind)
         ord      <- order(cases.dev[ind])
         df.case  <- data.frame(x = cases.dev[ind][ord], y = rep(0, nind), yend = tht,
                                label = labels[ind[ord]], clr = clr[ord])
         plt.res  <- plt.res +
            geom_segment(aes(x = x, y = y, xend = x, yend = yend, col = clr), data = df.case, size = 1) +
            theme(legend.position = "none")
         if (!is.null(labels) & !any(is.na(labels)))
            plt.res <- plt.res + geom_label(aes(x = x, y = yend, label = label, col = clr), data = df.case)
            
      }
   }
   else if ("density strip" %in% display) {
      x     <- seq(1/ngrid, 1.1 * max(ctrls.dev), length = ngrid)
      dens  <- density(ctrls.dev, n = ngrid, from = min(x), to = max(x))$y
      dat   <- data.frame(x, dens, gp = factor(rep("Residuals", ngrid)))
      ht    <- 0.6
      wd    <- 1.2 * ht / 2
      plt.res <- ggplot(dat, aes(x = x, y = gp)) + theme_classic() +
         geom_tile(aes(fill = dens), height = ht) +
         scale_fill_gradient(low = "white", high = "black") +
         theme(axis.title.y = element_blank(), legend.position  = "none",
               axis.line.y = element_blank(), axis.ticks.y = element_blank(),
               axis.text.y = element_blank()) +
         xlab(ttl)
      if (nind > 0) {
         cdev        <- cases.dev[ind]
         names(cdev) <- as.character(ind)
         df.case     <- data.frame(x = cases.dev[ind], y = rep(0, nind), label = labels[ind], clr = clr)
         plt.res <- plt.res +
            geom_segment(aes(x = x, xend = x, y = 1 + wd, yend = 1 - wd, col = clr),
                         data = df.case, size = 1)
         if (!is.null(labels) & !any(is.na(labels)))
            plt.res <- plt.res + geom_label(aes(x = x, y = 1 + wd, label = label, col = clr), data = df.case)
      }
   }

   if (!("none") %in% display) {
      if (requireNamespace("gridExtra", quietly = TRUE))
         grid.arrange(plt.pc, plt.res, nrow = 2)
      else {
         print(plt.pc)
         print(plt.res)
      }
   }

   result <- list(closest.control = cc, case.pca = cases.pca, case.procrustes = cases.proc,
                  atypicality = atyp, cases.scl = scl, cases.residuals = cases.res,
                  cc.pca = cc.pca, npc = npc)
   if (any(c("histogram", "density strip") %in% display)) {
      result$plot.pc = plt.pc
      result$plot.residual = plt.res
   }

   invisible(result)
}
