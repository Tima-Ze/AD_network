library(tidyverse)
library(Rcpp)

sourceCpp("bin/wTO_correlation_snakemake.cpp")



wTO_edit <- function(A, sign = c("abs", "sign")) {
  if (sign %in% c("abs", "absolute")) {
    A <- abs(A)
  }
  A_TF <- as.data.frame(subset(A, select = row.names(A)))
  C <- as.matrix(A) %*% t(A)
  W <- C + A_TF
  KI <- rowSums(abs(A), na.rm = T)

  KII <- KI
  K <- sapply(KII, pmin, y = KI)
  WTO <- round(W / (K + 1 - abs(A_TF)), 3)
  return(WTO)
}


CorrelationOverlap_edit <- function(Data, Overlap, method) {
  if (method == "p") {
    COR <- cp_cor(t(Data))
    # COR <- cp_cor(t(Data))
  } else if (method == "s") {
    COR <- spearman_cor(t(Data))
  }
  rownames(COR) <- rownames(Data)
  colnames(COR) <- rownames(Data)
  # diag(COR) <- 0
  # COR <- t(COR)
  COR[!is.finite(COR)] <- 0
  COR[is.na(COR)] <- 0
  A <- subset(COR, row.names(COR) %in% Overlap)
  return(A)
}

wTO.in.line_edit <- function(d) {
  upperTriangle <- upper.tri(d, diag = F)
  d.upperTriangle <- d
  d.upperTriangle[!upperTriangle] <- NA
  d_melted <- data.table::as.data.table(stats::na.omit(reshape2::melt(as.matrix(d.upperTriangle),
    value.name = "correlationCoef"
  )))
  names(d_melted) <- c("Node.1", "Node.2", "wTO")
  return(d_melted)
}




wTOFast_edit <- function(Data, Overlap = row.names(Data), method = "p", sign = "sign",
                         delta = 0.2, n = 10, method_resampling = "Bootstrap", lag = NULL,
                         ID = NULL) {
  Overlap <- unique(as.character(Overlap))
  `%ni%` <- Negate(`%in%`)
  if (is.numeric(n) == F) {
    stop("n must be numeric.")
  }
  if (n <= 0) {
    stop("n must be greater than 0.")
  }
  if (is.data.frame(Data) == F) {
    stop("Data must be a data.frame.")
  }
  if (method %ni% c("s", "spearman", "p", "pearson")) {
    stop("Method must be: \"s\", \"spearman\", \"p\" or \"pearson\".")
  }
  if (method_resampling %ni% c("Bootstrap", "BlockBootstrap")) {
    stop("Method must be: \"Bootstrap\" or \"BlockBootstrap\".")
  }
  if (method_resampling %in% "BlockBootstrap") {
    if (is.null(lag) & is.null(ID)) {
      stop("If you want to use the \"BlockBootstrap\" please give a lag or the indivuals ID.")
    }
    if (!is.null(lag) & !is.null(ID)) {
      stop("If you want to use the \"BlockBootstrap\" please give a lag OR the indivuals ID.")
    }
  }
  DIM_Overlap <- nrow(subset(Data, row.names(Data) %in% Overlap))
  if (DIM_Overlap == 0) {
    stop("There is no overlapping nodes. Please check your input \"Overlap\"")
  }
  if (!is.null(DIM_Overlap)) {
    message(paste(
      "There are", DIM_Overlap, "overlapping nodes,",
      dim(Data)[1], "total nodes and", dim(Data)[2], "individuals."
    ))
  }
  message("This function might take a long time to run. Don't turn off the computer.")
  wtomelt0 <- CorrelationOverlap_edit(
    Data = Data, Overlap = Overlap,
    method = method
  ) %>% wTO_edit(., sign)
  `%>%` <- magrittr::`%>%`
  . <- NULL
  method_resampling <- "Bootstrap"
  for (i in 1:n) {
    message(" ", i, " ", appendLF = FALSE)
    if (method_resampling == "BlockBootstrap") {
      if (!is.null(lag)) {
        nsampl <- ifelse(ncol(Data) %% lag == 0, ncol(Data) %/% lag,
          ncol(Data) %/% lag + 1
        )
        Y <- sample(1:nsampl, size = nsampl, replace = T)
        Vect <- Y * lag
        j <- lag - 1
        while (j > 0) {
          Vect <- cbind(Vect, Y * lag - j)
          j <- j - 1
        }
        SAMPLES <- c(Vect)
        SAMPLES[SAMPLES > ncol(Data)] <- NA
        SAMPLE <- stats::na.exclude(SAMPLES)
        Data_boot <- Data[, SAMPLE]
      }
      if (!is.null(ID)) {
        ID %<>% as.factor
        bootID <- sample(levels(ID), replace = TRUE)
        Data_boot <- subset(Data, select = ID %in% bootID[1])
        for (k in 2:length(bootID)) {
          Data_boot <- cbind(Data_boot, subset(Data,
            select = ID %in% bootID[k]
          ))
        }
      }
      #      res = wTO::CorrelationOverlap(Data = Data_boot,
      #                                    Overlap = Overlap, method = method) %>% wTO::wTO(.,
      #                                                                                     sign)
      res <- CorrelationOverlap_edit(
        Data = Data_boot,
        Overlap = Overlap, method = method
      ) %>% wTO_edit(
        .,
        sign
      )
    } else if (method_resampling != "BlockBootstrap") {
      res <- CorrelationOverlap_edit(Data = Data[, sample(1:ncol(Data), replace = TRUE)], Overlap = Overlap, method = method) %>%
        wTO_edit(., sign)
    }
    U <- (res < wtomelt0 - delta) + (res > wtomelt0 + delta)
    if (i == 1) {
      out <- U
    }
    if (i != 1) {
      out <- out + U
    }
    rm(res)
    rm(U)
  }
  wtomelt0 <- wTO.in.line_edit(wtomelt0)
  cor <- wTO.in.line_edit(out)
  adj.pval <- p.adjust(cor$wTO / n, method = "BH")
  pval <- data.table::data.table(wtomelt0,
    pval = cor$wTO / n,
    pval.adj = adj.pval
  )
  message("Done!")
  return(pval)
}


wTO.Complete_edit <- function(k = 1, n = 100, Data, Overlap = row.names(Data),
                              method = "p", method_resampling = "Bootstrap",
                              pvalmethod = "BH", savecor = F,
                              expected.diff = 0.20, lag = NULL, ID = NULL,
                              normalize = F, plot = T, plotName = "plot") {
  N <- k
  Overlap <- unique(as.character(Overlap))
  `%ni%` <- Negate(`%in%`)
  ##### Messages
  if (is.numeric(k) == F) {
    stop("k must be numeric.")
  }
  if (k <= 0) {
    stop("k must be greater than 0.")
  }
  if (is.numeric(n) == F) {
    stop("n must be numeric.")
  }
  if (n <= 0) {
    stop("n must be greater than 0.")
  }
  if (is.data.frame(Data) == F) {
    stop("Data must be a data.frame.")
  }

  if (method %ni% c("s", "spearman", "p", "pearson")) {
    stop('Method must be: "s", "spearman", "p" or "pearson".')
  }

  if (method_resampling %ni% c("Bootstrap", "Reshuffle", "BlockBootstrap")) {
    stop('Method must be: "Bootstrap", "BlockBootstrap" or "Reshuffle".')
  }
  if (method_resampling %in% "BlockBootstrap") {
    if (is.null(lag) & is.null(ID)) {
      stop('If you want to use the "BlockBootstrap" please give a lag or the indivuals ID.')
    }
    if (!is.null(lag) & !is.null(ID)) {
      stop('If you want to use the "BlockBootstrap" please give a lag OR the indivuals ID.')
    }
  }
  if (pvalmethod %ni% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")) {
    stop("pvalmethod must be:  'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr' or 'none'")
  }

  if (normalize %ni% c(T, F)) {
    stop("normalize must be:  TRUE or FALSE.")
  }

  if (normalize == T) {
    Data.n <- as.data.frame(som::normalize(Data))
    row.names(Data.n) <- row.names(Data)
    Data <- Data.n
  }

  DIM_Overlap <- nrow(subset(Data, row.names(Data) %in% Overlap))
  if (DIM_Overlap == 0) {
    stop('There is no overlapping nodes. Please check your input "Overlap"')
  }
  if (!is.null(DIM_Overlap)) {
    message(paste(
      "There are", DIM_Overlap, "overlapping nodes,", dim(Data)[1],
      "total nodes and", dim(Data)[2], "individuals."
    ))
  }


  message("This function might take a long time to run. Don't turn off the computer.")
  PAR <- par()

  ## For the original data
  # real_Genes = Data
  Saving <- CorrelationOverlap_edit(Data = Data, Overlap = Overlap, method = method)
  WTO_abs <- wTO_edit(A = Saving, sign = "abs")
  WTO_sign <- wTO_edit(A = Saving, sign = "sign")
  Cor_real <- wTO.in.line_edit(WTO_sign)
  Cor_real_abs <- wTO.in.line_edit(WTO_abs)
  names(Cor_real) <- names(Cor_real_abs) <- c("Node.1", "Node.2", "wTO_0")
  idcol <- c("Node.1", "Node.2")
  rm("WTO_abs")
  rm("WTO_sign")
  data.table::setkeyv(Cor_real, c("Node.1", "Node.2"))
  data.table::setkeyv(Cor_real_abs, c("Node.1", "Node.2"))

  Orig <- cbind(Rep = 0, Cor_real[Cor_real_abs])
  names(Orig) <- c("Rep", "Node.1", "Node.2", "wTO_sign", "wTO_abs")
  reps_rest <- n
  ### If only one node
  if (k == 1) {
    a <- 0
    while (reps_rest > 0) {
      message(" ", reps_rest, " ", appendLF = FALSE)
      # message(a)
      K <- 1:min(N, reps_rest)

      # K = 1:n
      OUTPUT <- lapply(K, wTO.aux.each_edit,
        Data = Data,
        Overlap = Overlap, method = method, ID, lag = lag, method_resampling = method_resampling
      )

      ALL <- data.table::rbindlist(OUTPUT, idcol = idcol)
      names(ALL) <- names(Orig) <- c("Rep", "Node.1", "Node.2", "wTO_sign", "wTO_abs")

      ALL_DT_sig <- data.table::dcast(ALL, Node.1 + Node.2 ~ Rep, value.var = "wTO_sign")
      ALL_DT_abs <- data.table::dcast(ALL, Node.1 + Node.2 ~ Rep, value.var = "wTO_abs")

      if (a == 0) {
        Ps1 <- rowSums(ALL_DT_sig[, -c(1:2)] < Orig$wTO_sign - expected.diff)
        Ps2 <- rowSums(ALL_DT_sig[, -c(1:2)] > Orig$wTO_sign + expected.diff)
        Ps <- Ps1 + Ps2

        Pa1 <- rowSums(ALL_DT_abs[, -c(1:2)] < Orig$wTO_abs - expected.diff)
        Pa2 <- rowSums(ALL_DT_abs[, -c(1:2)] > Orig$wTO_abs + expected.diff)
        Pa <- Pa1 + Pa2

        TAB_SIGN <- as.data.frame(table(unlist(round(ALL_DT_sig[, -c(1:2)], 2))))
        TAB_ABS <- as.data.frame(table(unlist(round(ALL_DT_abs[, -c(1:2)], 2))))
      }
      if (a > 0) {
        Ps1 <- rowSums(ALL_DT_sig[, -c(1:2)] < Orig$wTO_sign - expected.diff)
        Ps2 <- rowSums(ALL_DT_sig[, -c(1:2)] > Orig$wTO_sign + expected.diff)
        Ps <- Ps + Ps1 + Ps2

        Pa1 <- rowSums(ALL_DT_abs[, -c(1:2)] < Orig$wTO_abs - expected.diff)
        Pa2 <- rowSums(ALL_DT_abs[, -c(1:2)] > Orig$wTO_abs + expected.diff)
        Pa <- Pa + Pa1 + Pa2

        TAB_SIGN_aux <- as.data.frame(table(unlist(round(ALL_DT_sig[, -c(1:2)], 2))))
        TAB_ABS_aux <- as.data.frame(table(unlist(round(ALL_DT_abs[, -c(1:2)], 2))))

        TAB_SIGN <- plyr::join(TAB_SIGN, TAB_SIGN_aux, by = "Var1")
        TAB_SIGN <- data.frame(
          Var1 = TAB_SIGN$Var1,
          Sum = rowSums(TAB_SIGN[, -1])
        )
        TAB_ABS <- plyr::join(TAB_ABS, TAB_ABS_aux, by = "Var1")
        TAB_ABS <- data.frame(
          Var1 = TAB_ABS$Var1,
          Sum = rowSums(TAB_ABS[, -1])
        )
      }
      rm("ALL_DT_sig", "ALL_DT_abs", "ALL", "OUTPUT")

      reps_rest <- (reps_rest - N)
      a <- a + 1
    }
  } else if (k > 1) {
    WTO <- new.env()
    assign("Data", Data, envir = WTO)
    assign("Overlap", Overlap, envir = WTO)
    assign("method", method, envir = WTO)
    assign("CorrelationOverlap_edit", CorrelationOverlap_edit, envir = WTO)
    assign("wTO_edit", wTO_edit, envir = WTO)
    assign("wTO.in.line_edit", wTO_edit, envir = WTO)
    assign("wTO.aux.each_edit", wTO.aux.each_edit, envir = WTO)
    assign("method_resampling", method_resampling, envir = WTO)
    assign("sample_ind", sample_ind, envir = WTO)
    assign("lag", lag, envir = WTO)
    assign("ID", ID, envir = WTO)
    cl <- parallel::makeCluster(k)
    parallel::clusterExport(cl, "Data", envir = WTO)
    parallel::clusterExport(cl, "wTO.in.line_edit", envir = WTO)
    parallel::clusterExport(cl, "lag", envir = WTO)
    parallel::clusterExport(cl, "Overlap", envir = WTO)
    parallel::clusterExport(cl, "method", envir = WTO)
    parallel::clusterExport(cl, "CorrelationOverlap_edit", envir = WTO)
    parallel::clusterExport(cl, "wTO_edit", envir = WTO)
    parallel::clusterExport(cl, "wTO.aux.each_edit", envir = WTO)
    parallel::clusterExport(cl, "method_resampling", envir = WTO)
    parallel::clusterExport(cl, "sample_ind", envir = WTO)
    # message("cluster")
    # K = 1:n

    a <- 0
    while (reps_rest > 0) {
      # message(a)
      K <- 1:min(N, reps_rest)

      OUTPUT <- parallel::clusterApply(cl, K, wTO.aux.each_edit,
        Data = Data,
        Overlap = Overlap, ID, lag = lag, method = method, method_resampling = method_resampling
      )
      ALL <- data.table::rbindlist(OUTPUT, idcol = idcol)
      names(ALL) <- names(Orig) <- c("Rep", "Node.1", "Node.2", "wTO_sign", "wTO_abs")

      ALL_DT_sig <- data.table::dcast(ALL, Node.1 + Node.2 ~ Rep, value.var = "wTO_sign")
      ALL_DT_abs <- data.table::dcast(ALL, Node.1 + Node.2 ~ Rep, value.var = "wTO_abs")


      if (a == 0) {
        Ps <- rowSums(ALL_DT_sig[, -c(1:2)] < Orig$wTO_sign - expected.diff) +
          rowSums(ALL_DT_sig[, -c(1:2)] > Orig$wTO_sign + expected.diff)

        Pa <- rowSums(ALL_DT_abs[, -c(1:2)] < Orig$wTO_abs - expected.diff) +
          rowSums(ALL_DT_abs[, -c(1:2)] > Orig$wTO_abs + expected.diff)


        TAB_SIGN <- as.data.frame(table(unlist(round(ALL_DT_sig[, -c(1:2)], 2))))
        TAB_ABS <- as.data.frame(table(unlist(round(ALL_DT_abs[, -c(1:2)], 2))))
      }
      if (a > 0) {
        Ps <- Ps + rowSums(ALL_DT_sig[, -c(1:2)] < Orig$wTO_sign - expected.diff) +
          rowSums(ALL_DT_sig[, -c(1:2)] > Orig$wTO_sign + expected.diff)
        Pa <- Pa + rowSums(ALL_DT_abs[, -c(1:2)] < Orig$wTO_abs - expected.diff) +
          rowSums(ALL_DT_abs[, -c(1:2)] > Orig$wTO_abs + expected.diff)

        # message(Pa)
        # message(Ps)

        TAB_SIGN_aux <- as.data.frame(table(unlist(round(ALL_DT_sig[, -c(1:2)], 2))))
        TAB_ABS_aux <- as.data.frame(table(unlist(round(ALL_DT_abs[, -c(1:2)], 2))))

        TAB_SIGN <- plyr::join(TAB_SIGN, TAB_SIGN_aux, by = "Var1")
        TAB_SIGN <- data.frame(
          Var1 = TAB_SIGN$Var1,
          Sum = rowSums(TAB_SIGN[, -1])
        )
        TAB_ABS <- plyr::join(TAB_ABS, TAB_ABS_aux, by = "Var1")
        TAB_ABS <- data.frame(
          Var1 = TAB_ABS$Var1,
          Sum = rowSums(TAB_ABS[, -1])
        )
      }
      rm("ALL_DT_sig", "ALL_DT_abs", "ALL", "OUTPUT")

      reps_rest <- (reps_rest - N)
      a <- a + 1
    }

    parallel::stopCluster(cl)
  }

  message("Simulations are done.")

  message("Computing p-values")
  Orig$pval_sig <- Ps / n
  Orig$pval_abs <- Pa / n


  if (method_resampling == "Reshuffle") {
    Orig$pval_sig <- 1 - Orig$pval_sig
    Orig$pval_abs <- 1 - Orig$pval_abs
  }

  Orig$Padj_sig <- (stats::p.adjust(Orig$pval_sig, method = pvalmethod))
  Orig$Padj_abs <- (stats::p.adjust(Orig$pval_abs, method = pvalmethod))


  ## Running the correlation
  if (savecor == T) {
    Total_Correlation <- as.data.frame(stats::cor(t(Data), method = method))
    Total_Correlation <- wTO.in.line_edit(Total_Correlation)
    names(Total_Correlation) <- c("Node.1", "Node.2", "Cor")
  }
  if (savecor == F) {
    Total_Correlation <- NULL
  }

  TAB_SIGN_aux <- as.data.frame(table(round(Orig$wTO_sign, 2)))
  TAB_ABS_aux <- as.data.frame(table(round(Orig$wTO_abs, 2)))
  TAB_SIGN <- plyr::join(TAB_SIGN, TAB_SIGN_aux, by = "Var1")
  TAB_ABS <- plyr::join(TAB_ABS, TAB_ABS_aux, by = "Var1")

  message("Computing cutoffs")
  if (plot == TRUE) {
    png(paste0("wTO.complete_", plotName, "_.png"), width = 1200, height = 900)
    graphics::par(mar = c(5.1, 4.1, 4.1, 8.1), xpd = TRUE, mfrow = c(3, 2))
  }

  Cutoffs <- Cut.off(TAB_SIGN, "wTO - Resampling", plot = plot)
  Cutoffs_abs <- Cut.off(TAB_ABS, "|wTO| - Resampling", plot = plot)


  Orig <- Orig[, -"Rep"]


  Orig$wTO_abs <- as.numeric(Orig$wTO_abs)
  Orig$wTO_sign <- as.numeric(Orig$wTO_sign)
  Orig$pval_abs <- as.numeric(Orig$pval_abs)
  Orig$pval_sig <- as.numeric(Orig$pval_sig)
  Orig$Padj_abs <- as.numeric(Orig$Padj_abs)
  Orig$Padj_sig <- as.numeric(Orig$Padj_sig)

  Quantiles <- rbind(
    Cutoffs$Empirical.Quantile,
    Cutoffs$Quantile,
    Cutoffs_abs$Empirical.Quantile,
    Cutoffs_abs$Quantile
  )
  row.names(Quantiles) <- c(
    "Empirical.Quantile",
    "Quantile",
    "Empirical.Quantile.abs",
    "Quantile.abs"
  )

  tQ <- as.data.frame(t(Quantiles))
  output <- list(
    wTO = Orig,
    Correlation = Total_Correlation,
    Quantiles = Quantiles
  )
  col <- ifelse(Orig$pval_sig < 0.05 & Orig$pval_abs < 0.05, "red",
    ifelse(Orig$pval_sig < 0.05, "orange",
      ifelse(Orig$pval_abs < 0.05, "yellow", "black")
    )
  )
  if (plot == T) {
    # png(paste0("wTO.complete_Plot1.png"),width = 1200,height = 900)

    # par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE, mfrow = c(3,1))
    graphics::plot(Orig$wTO_sign, Orig$wTO_abs,
      axes = F,
      xlab = "|wTO|", ylab = "wTO",
      main = "|wTO| vs wTO", pch = ".", xlim = c(-1, 1), ylim = c(0, 1),
      col.main = "steelblue2", col.lab = "steelblue2", col = col
    )
    graphics::axis(1,
      las = 1, cex.axis = 0.8, col = "steelblue",
      col.ticks = "steelblue3", col.axis = "steelblue"
    )
    graphics::axis(2, las = 1, cex.axis = 0.8, col = "steelblue", col.ticks = "steelblue3", col.axis = "steelblue")


    graphics::legend(c(0.9, 0), c(
      "p-value < 0.05", "wTO sign & |wTO|",
      "wTO sign", "|wTO|"
    ),
    col = c("transparent", "red", "orange", "yellow"), pch = 16, bty = "n",
    inset = c(-0.8, 0), cex = 0.8
    )


    graphics::par(xpd = FALSE)
    graphics::abline(h = 0, lty = 2, col = "gray50")
    graphics::abline(v = 0, lty = 2, col = "gray50")
    # dev.off()

    # png(paste0("wTO.complete_Plot2.png"),width = 1200,height = 900)
    graphics::plot(Orig$wTO_sign, Orig$pval_sig,
      axes = F,
      xlab = "wTO", ylab = "p-value", ylim = c(0, 1), xlim = c(-1, 1), col.main = "steelblue2", col.lab = "steelblue2",
      main = "wTO vs p-value",
      pch = 16
    )
    graphics::axis(1,
      las = 1, cex.axis = 0.8, col = "steelblue",
      col.ticks = "steelblue3", col.axis = "steelblue"
    )
    graphics::axis(2, las = 1, cex.axis = 0.8, col = "steelblue", col.ticks = "steelblue3", col.axis = "steelblue")

    graphics::par(xpd = FALSE)
    graphics::abline(v = tQ$Empirical.Quantile, lty = 2, col = c("red", "orange", "yellow", "yellow", "orange", "red"))
    graphics::par(xpd = T)
    graphics::legend(c(0.9, 0), c("Empirical Quantiles", "0.1%", "2.5%", "10%", "90%", "97.5%", "99.9%"),
      col = c("white", "red", "orange", "yellow", "yellow", "orange", "red"), lwd = 2, bty = "n",
      inset = c(-0.8, 0), cex = 0.8
    )
    # dev.off()

    # png(paste0("wTO.complete_Plot3.png"),width = 1200,height = 900)
    graphics::par(xpd = FALSE)
    graphics::plot(Orig$wTO_abs, Orig$pval_abs,
      axes = F,
      xlab = "|wTO|", ylab = "p-value", ylim = c(0, 1), xlim = c(0, 1),
      main = "|wTO| vs p-value",
      pch = 16, col.main = "steelblue2", col.lab = "steelblue2"
    )
    graphics::axis(1,
      las = 1, cex.axis = 0.8, col = "steelblue",
      col.ticks = "steelblue3", col.axis = "steelblue"
    )
    graphics::axis(2, las = 1, cex.axis = 0.8, col = "steelblue", col.ticks = "steelblue3", col.axis = "steelblue")


    graphics::abline(v = tQ$Empirical.Quantile.abs, lty = 2, col = c("red", "orange", "yellow", "yellow", "orange", "red"))
    graphics::par(xpd = T)
    graphics::legend(c(0.9, 0), c("Empirical Quantiles", "0.1%", "2.5%", "10%", "90%", "97.5%", "99.9%"),
      col = c("white", "red", "orange", "yellow", "yellow", "orange", "red"), lwd = 2, bty = "n",
      inset = c(-0.8, 0), cex = 0.8
    )
    # dev.off()
    dev.off()
  }


  class(output) <- append("wTO", class(output))
  message("Done!")
  return(output)
}


wTO.aux.each_edit <- function(n, Data, Overlap, method, method_resampling, lag, ID) {
  if (method_resampling == "Bootstrap") {
    real_Genes <- sample(Data, replace = T)
  }
  if (method_resampling == "BlockBootstrap") {
    if (!is.null(lag)) {
      nsampl <- ifelse(ncol(Data) %% lag == 0, ncol(Data) %/% lag, ncol(Data) %/% lag + 1)
      Y <- sample(1:nsampl, size = nsampl, replace = T)
      Vect <- Y * lag
      i <- lag - 1
      while (i > 0) {
        Vect <- cbind(Vect, Y * lag - i)
        i <- i - 1
      }

      SAMPLES <- c(Vect)
      SAMPLES[SAMPLES > ncol(Data)] <- NA
      SAMPLE <- stats::na.exclude(SAMPLES)
      real_Genes <- Data[, SAMPLE]
      row.names(real_Genes) <- row.names(Data)
    }

    if (!is.null(ID)) {
      ID %<>% as.factor
      bootID <- sample(levels(ID), replace = TRUE)


      Data_boot <- subset(Data, select = ID %in% bootID[1])

      for (k in 2:length(bootID)) {
        real_Genes <- cbind(
          Data_boot,
          subset(Data, select = ID %in% bootID[k])
        )
      }
    }
  } else if (method_resampling == "Reshuffle") {
    real_Genes <- as.data.frame(lapply(1:ncol(Data), FUN = sample_ind, dfExpression = Data))
    names(real_Genes) <- names(Data)
    row.names(real_Genes) <- row.names(Data)
  }


  Saving <- CorrelationOverlap_edit(Data = real_Genes, Overlap = Overlap, method = method)
  WTO_abs <- wTO_edit(A = Saving, sign = "abs")
  WTO_sig <- wTO_edit(A = Saving, sign = "sign")

  # message(".", appendLF = F)
  Cor_star <- wTO.in.line_edit(WTO_sig)
  Cor_star_abs <- wTO.in.line_edit(WTO_abs)
  names(Cor_star) <- c("Node.1", "Node.2", "wTo_sign")
  names(Cor_star_abs) <- c("Node.1", "Node.2", "wTo_abs")
  data.table::setkeyv(Cor_star, c("Node.1", "Node.2"))
  data.table::setkeyv(Cor_star_abs, c("Node.1", "Node.2"))

  return(Cor_star[Cor_star_abs])
}


Cut.off <- function(wTO_value, type, plot) {
  `%ni%` <- Negate(`%in%`)
  wTO_value <- plyr::arrange(wTO_value, wTO_value$Var1)
  wTO_value[is.na(wTO_value)] <- 0
  wTO_value$relstar <- wTO_value[, 2] / sum(wTO_value[, 2], na.rm = T)
  wTO_value$relreal <- wTO_value[, 3] / sum(wTO_value[, 3], na.rm = T)


  quantile.from.freq <- function(vals, freq, quant) {
    ord <- order(vals)
    cs <- cumsum(freq[ord])

    if (length(which(cs < quant)) > 0) {
      return(vals[max(which(cs < quant)) + 1])
    }
    if (length(which(cs < quant)) == 0) {
      return(min(vals))
    }
  }


  wTO_value$Var1 <- as.numeric(as.matrix(wTO_value$Var1))
  quantile_star <- data.frame(
    quantile.from.freq(wTO_value$Var1, wTO_value$relstar, 0.001),
    quantile.from.freq(wTO_value$Var1, wTO_value$relstar, 0.025),
    quantile.from.freq(wTO_value$Var1, wTO_value$relstar, 0.1),
    quantile.from.freq(wTO_value$Var1, wTO_value$relstar, 0.9),
    quantile.from.freq(wTO_value$Var1, wTO_value$relstar, 0.975),
    quantile.from.freq(wTO_value$Var1, wTO_value$relstar, 0.999)
  )

  quantile_real <- data.frame(
    quantile.from.freq(wTO_value$Var1, wTO_value$relreal, 0.001),
    quantile.from.freq(wTO_value$Var1, wTO_value$relreal, 0.025),
    quantile.from.freq(wTO_value$Var1, wTO_value$relreal, 0.1),
    quantile.from.freq(wTO_value$Var1, wTO_value$relreal, 0.9),
    quantile.from.freq(wTO_value$Var1, wTO_value$relreal, 0.975),
    quantile.from.freq(wTO_value$Var1, wTO_value$relreal, 0.999)
  )


  names(quantile_real) <- names(quantile_star) <- c("0.1%", "2.5%", "10%", "90%", "97.5%", "99.9%")
  if (plot == TRUE) {
    PLOT <- function(wTO_value) {
      graphics::par(xpd = FALSE)
      graphics::plot(wTO_value$relstar ~ wTO_value$Var1,
        type = "l",
        xlim = c(floor(min(wTO_value$Var1)), 1),
        main = type,
        ylim = c(0, 1), axes = F,
        # ylim = c(0, max(wTO_value$relstar)), axes = F,
        xlab = "wTO", ylab = "Density", col.main = "steelblue2", col.lab = "steelblue2"
      )
      graphics::lines(wTO_value$Var1, wTO_value$relreal, type = "l", col = "violet")
      graphics::abline(h = 0, col = "gray", lty = 4)


      graphics::abline(v = c(quantile_real), col = c("red", "orange", "yellow", "yellow", "orange", "red"), lty = 2)

      graphics::axis(1,
        las = 1, cex.axis = 0.8, col = "steelblue",
        col.ticks = "steelblue3", col.axis = "steelblue"
      )
      graphics::axis(2, las = 1, cex.axis = 0.8, col = "steelblue", col.ticks = "steelblue3", col.axis = "steelblue")
      graphics::par(xpd = T)
      graphics::legend(c(0.9, max(wTO_value$relstar)), c(
        "wTO - Data set",
        "wTO - Reshuffle",
        "99.9%",
        "95%",
        "80%"
      ),
      inset = c(-0.8, 0), lwd = 2,
      lty = 1, col = c(
        "violet",
        "black",
        "yellow", "orange", "red"
      ), bty = "n", cex = 0.8
      )
    }
    res <- try(PLOT(wTO_value))
    if (!methods::is(res, "try-error")) {
      res
    }
  }
  return(list(Empirical.Quantile = quantile_star, Quantile = quantile_real))
}


sample_ind <- function(x, dfExpression) {
  z <- base::sample(dfExpression[, x], replace = F)
  return(z)
}
