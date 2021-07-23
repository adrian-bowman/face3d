#     Permutation test for the comparison of PCs across two groups

compare.pca <- function(x, group, group.pca = FALSE, weights, npc, nsim, monitor = TRUE) {

   if (is.list(x)) {
      nms  <- c("pca", "group", "group.pca", "weights", "p", "npc", "nsim",
                "t2", "t2.sim")
      mtch <- match(nms, names(x))
      if (any(is.na(mtch)))
         stop(paste("x is missing the component(s):",
                    paste(nms[is.na(mtch)], collapse = ", ")))
      pca       <- x$pca
      group     <- x$group
      group.pca <- x$group.pca
      p         <- x$p
      npc       <- x$npc
      weights   <- x$weights
      t2        <- x$t2
      t2.sim    <- x$t2.sim
      nsim      <- x$nsim
   }
   else {
      rotated <- x
      if (missing(group)) stop("groups must be specified.")
      if (missing(weights)) weights <- rep(1, nrow(rotated))
   }
   
   lvls <- levels(factor(as.character(group)))
   if (missing(nsim)) {
      nsim <- if (group.pca) 200 else 500
   }
   
   hotstat <- function(x, group, group.pca) {
      x1  <- x[group == lvls[1], ]
      x2  <- x[group == lvls[2], ]
      n1  <- nrow(x1)
      n2  <- nrow(x2)
      m1  <- apply(x1, 2, mean)
      m2  <- apply(x2, 2, mean)
      cv  <- ((n1 - 1) * cov(x1) + (n2 - 1) * cov(x2)) / (n1 + n2 - 2)
      eig <- eigen(cv)
      d   <- abs(t(eig$vectors) %*% (m1 - m2) / sqrt(eig$values))
      # d   <- d^2 * n1 * n2 / (n1 + n2)
      # ans <- mean(d)
      # if (!group.pca)
      #    d <- apply(x, 2, function(z) t.test(z ~ group)$statistic^2)
      d   <- d * sqrt(n1 * n2 / (n1 + n2))
      ans <- sqrt(mean(d^2))
      if (!group.pca)
         d <- apply(x, 2, function(z) abs(t.test(z ~ group)$statistic))
      ans <- c(ans, d)
      names(ans) <- c("global", paste("pc", 1:npc))
      ans
   }
   
   if (!is.list(x)) {
      pca  <- if (group.pca) pca.face3d(rotated, group, weights) else pca.face3d(rotated, weights = weights)
      if (missing(npc)) npc <- min(which(pca$cumpct > 80))

      t2 <- hotstat(pca$scores[ , 1:npc], group, group.pca)

      # Generate the null hypthesis distribution by random permutations
      t2.sim <- matrix(0, nrow = nsim, ncol = 1 + npc)
      for (isim in 1:nsim) {
         gp <- sample(group)
         pca.sim <- if (group.pca) pca.face3d(rotated, gp, weights = weights) else pca
         t2.sim[isim, ] <- hotstat(pca.sim$scores[ , 1:npc], gp, group.pca)
         if (group.pca) cat(isim, "")
      }
      if (group.pca) cat("\n")
   }
   
   # Print results
   p.fn <- function(x, tstat) length(which(x > tstat)) / length(x)
   p    <- rep(0, 1 + npc)
   for (i in 1:(1 + npc)) p[i] <- p.fn(t2.sim[ , i], t2[i])
   names(p) <- c("global", as.character(1:npc))
   cat("Permutation p-values:\n")
   print(p)
   
   # Plot 
   if (monitor) {
      pclev <- c(paste("global  p=", p[1], sep = ""),
                 paste("PC ", 1:npc, " ", round(pca$pct[1:npc], 1), "%", " p=", round(p[-1], 3), sep = ""))
      pc    <- factor(rep(pclev, each = nsim), levels = pclev)
      dfrm  <- data.frame(t2.sim = c(t2.sim), gp = factor(rep(1, nsim)), pc = pc)
      plt   <- ggplot(dfrm, aes(gp, t2.sim)) +
                      ylab("|t| statistic") +
                      theme(legend.position = "none", axis.title.x = element_blank(),
                            axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
                      geom_boxplot() +
                      facet_grid(. ~ pc, labeller = label_wrap_gen(width = 6)) +
                      geom_point(aes(x, y, col = clr, shape = g, size = 1.5),
                                 data = data.frame(x = rep(1, length(t2)), y = t2,
                                                   clr = p <= 0.05 / c(1, rep(npc, npc)),
                                                   g = (names(t2) == "global"), pc = pclev)) +
                      scale_colour_manual(values = c("blue", "red"))
      print(plt)
   }
   
   # Plot density strips
   # del  <- diff(range(t2.sim)) / 10
   # evp  <- seq(min(t2.sim) - del, max(t2.sim) + del, length = 100)
   # sm.f <- function(x) sm.density(x, eval.points = evp, display = "none")$estimate
   # dens <- apply(t2.sim, 2, sm.f)
   # gp   <- factor(rep(c("global", paste("pc", 1:npc)), each = nsim),
   #                levels = c(paste("pc", npc:1), "global"))
   # dfrm <- data.frame(dens = c(dens), gp = gp, x = rep(evp, 1 + npc))
   # plt  <- ggplot(dfrm, aes(x = x, y = gp)) +
   #                theme_bw() + ylab("") + xlab("test statistic") +
   #                theme(legend.position = "none", axis.title.y = element_blank()) +
   #                geom_tile(aes(fill = dens)) +
   #                scale_fill_continuous(low="white", high="black") +
   #                geom_point(aes(x, y, col = "red"), data = data.frame(x = t2, y = names(t2)))
   #                # facet_grid(gp ~ ., scales = 'free_y', space = 'free')
   # print(plt)
   
   # Returned objects
   result <- list(pca = pca, group = group, group.pca = group.pca, p = p, npc = npc, weights = weights,
                  t2 = t2, t2.sim = t2.sim, nsim = nsim)
   if (is.list(x) & monitor) result <- plt
   invisible(result)
}
