
create_joined <- function(dnds_path_eden = "",
                          dnds_path_alt = "",
                                 gap_path = "",
                                 only_legend = FALSE,
                                 window_size = 10,
                                 gap_threshold = 0.6,
                                 gapcolor = TRUE,
                                 points = TRUE) {
  # process alt
  require(zoo)
  gap_data <- try(read.table(gap_path, header = F))
  dnds_data <- try(read.table(dnds_path_alt, header = T, sep=','))
  dN <- dnds_data$dN
  dS <- dnds_data$dS
  dNdS <- rep(0, length(dN))
  dNdS[which(dS > 0)] <- dN[which(dS > 0)] / dS[which(dS > 0)]
  # process the dnds data with sliding window
  pvalue_window <- sliding_window(dN,
                                  dS,
                                  gap_data =  gap_data,
                                  w_size = window_size,
                                  g_threshold = gap_threshold)
  
  # generate a matrix with start and end positions for plotting of significant windows
  epitopes_alt <-
    cluster_matrix(pvalues = pvalue_window,
                   w_size = window_size,
                   significance_level = 0.05)
  epitopes_alt$method <- "SLAC"
  
  # create data frame with all informations
  df_alt <- data.frame(
    position = dnds_data$Codon,
    selection = dNdS,
    selection_smooth = smooth.spline(dnds_data$Codon, dNdS, spar =
                                       0.35)$y,
    gap = gap_data$V1,
    gap_smooth = smooth.spline(dnds_data$Codon, gap_data$V1, spar =
                                 0.35)$y,
    stringsAsFactors = FALSE
  )
  df_alt$method <- "SLAC"
  # process eden
  
  gap_data <- try(read.table(gap_path, header = F))
  dnds_data <- try(read.table(dnds_path_eden, header = T))
  dN <- dnds_data$Nd_tip + dnds_data$Nd
  dS <- dnds_data$Sd_tip + dnds_data$Sd
  dNdS <- rep(0, length(dN))
  dNdS[which(dS > 0)] <- dN[which(dS > 0)] / dS[which(dS > 0)]
  # process the dnds data with sliding window
  pvalue_window <- sliding_window(dN,
                                  dS,
                                  gap_data =  gap_data,
                                  w_size = window_size,
                                  g_threshold = gap_threshold)
  
  # generate a matrix with start and end positions for plotting of significant windows
  epitopes_eden <-
    cluster_matrix(pvalues = pvalue_window,
                   w_size = window_size,
                   significance_level = 0.05)
  epitopes_eden$method <- "EDEN"
  
  # create data frame with all informations
  df_eden <- data.frame(
    position = dnds_data$pos,
    selection = dNdS,
    selection_smooth = smooth.spline(dnds_data$pos, dNdS, spar =
                                       0.35)$y,
    gap = gap_data$V1,
    gap_smooth = smooth.spline(dnds_data$pos, gap_data$V1, spar =
                                 0.35)$y,
    stringsAsFactors = FALSE
  )
  df_eden$method <- "EDEN"
 
  epitopes <- rbind(epitopes_alt, epitopes_eden)
  
  df_b <- rbind(df_eden, df_alt)

    fig_a <- ggplot(df_b, aes(position, selection, color=method))
   # fig_a <- fig_a + geom_point(size = 1.7, alpha = 3 / 4)
    fig_a <- fig_a + geom_line(aes(y = selection_smooth))
    fig_a <- fig_a + theme_classic() + xlab("position in alignment") + ylab("dN/dS ratio")
    fig_a <- fig_a + ggtitle(paste(dnds_path_eden, sep = ""))
    
    # highlight significant aereas
    fig_a <-
      fig_a + geom_rect(data = epitopes_eden,
        aes(NULL, NULL, xmin = start, xmax = end),
        ymin = -Inf,
        ymax = Inf,
        fill = "red"
      )
    fig_a <-
      fig_a + geom_rect(data = epitopes,
                        aes(NULL, NULL, xmin = start, xmax = end),
                        ymin = -Inf,
                        ymax = Inf,
                        fill = "blue"
      )

    fig_a
}








create_msa_plot_slac <- function(dnds_path = "",
                            gap_path = "",
                            only_legend = FALSE,
                            window_size = 10,
                            gap_threshold = 0.6,
                            gapcolor = TRUE,
                            points = TRUE) {
  require(zoo)
  gap_data <- try(read.table(gap_path, header = F))
  dnds_data <- try(read.table(dnds_path, header = T, sep=','))
  dN <- dnds_data$dN
  dS <- dnds_data$dS
  dNdS <- rep(0, length(dN))
  dNdS[which(dS > 0)] <- dN[which(dS > 0)] / dS[which(dS > 0)]
  # process the dnds data with sliding window
  pvalue_window <- sliding_window(dN,
                                  dS,
                                  gap_data =  gap_data,
                                  w_size = window_size,
                                  g_threshold = gap_threshold)
  
  # generate a matrix with start and end positions for plotting of significant windows
  epitopes <-
    cluster_matrix(pvalues = pvalue_window,
                   w_size = window_size,
                   significance_level = 0.05)
  
  # create data frame with all informations
  df <- data.frame(
    position = dnds_data$Codon,
    selection = dNdS,
    selection_smooth = smooth.spline(dnds_data$Codon, dNdS, spar =
                                       0.35)$y,
    gap = gap_data$V1,
    gap_smooth = smooth.spline(dnds_data$Codon, gap_data$V1, spar =
                                 0.35)$y,
    stringsAsFactors = FALSE
  )
  
  #### genereate plot ####
  # create sequence plot figure with dnds and gap informations
  fig_a <- ggplot(df, aes(position, selection))
  fig_a <-
    fig_a + theme_classic() + xlab("position in alignment") + ylab("dN/dS ratio")
  fig_a <- fig_a + ggtitle(paste(dnds_path, sep = ""))
  # highlight significant aereas
  fig_a <-
    fig_a + geom_rect(
      data = epitopes,
      aes(NULL, NULL, xmin = start, xmax = end),
      ymin = -Inf,
      ymax = Inf,
      fill = "red"
    )
  fig_a <-
    fig_a + geom_ribbon(aes(x = position, ymax = selection_smooth, ymin = 1), fill =
                          "grey80")
  if (points) {
    if (gapcolor) {
      fig_a <- fig_a + geom_point(aes(colour = gap), size = 1.7, alpha = 3 / 4)
    } else {
      fig_a <- fig_a + geom_point(size = 1.7, alpha = 3 / 4)
    }
  }
  
  fig_a <-
    fig_a + scale_colour_gradient(low = "#4db898", high = "#ea5420")
  
  
  fig_a <- fig_a + geom_line(aes(y = selection_smooth))
  #fig_a <- fig_a + labs(title=paste(sample_name,"; g_threshold=", gap_threshold,"; window_size=", window_size,sep=""))
  # fig_a <- fig_a + geom_hline(aes(yintercept=1))
  
  return(fig_a)
}



create_msa_plot <- function(dnds_path = "",
                            gap_path = "",
                            only_legend = FALSE,
                            window_size = 10,
                            gap_threshold = 0.6,
                            gapcolor = TRUE,
                            points = TRUE) {
  require(zoo)
  gap_data <- try(read.table(gap_path, header = F))
  dnds_data <- try(read.table(dnds_path, header = T))
  dN <- dnds_data$Nd_tip + dnds_data$Nd
  dS <- dnds_data$Sd_tip + dnds_data$Sd
  dNdS <- rep(0, length(dN))
  dNdS[which(dS > 0)] <- dN[which(dS > 0)] / dS[which(dS > 0)]
  # process the dnds data with sliding window
  pvalue_window <- sliding_window(dN,
                                  dS,
                                  gap_data =  gap_data,
                                  w_size = window_size,
                                  g_threshold = gap_threshold)
  
  # generate a matrix with start and end positions for plotting of significant windows
  epitopes <-
    cluster_matrix(pvalues = pvalue_window,
                   w_size = window_size,
                   significance_level = 0.05)
  
  # create data frame with all informations
  df <- data.frame(
    position = dnds_data$pos,
    selection = dNdS,
    selection_smooth = smooth.spline(dnds_data$pos, dNdS, spar =
                                       0.35)$y,
    gap = gap_data$V1,
    gap_smooth = smooth.spline(dnds_data$pos, gap_data$V1, spar =
                                 0.35)$y,
    stringsAsFactors = FALSE
  )
  
  #### genereate plot ####
  # create sequence plot figure with dnds and gap informations
  fig_a <- ggplot(df, aes(position, selection))
  fig_a <-
    fig_a + theme_classic() + xlab("position in alignment") + ylab("dN/dS ratio")
  fig_a <- fig_a + ggtitle(paste(dnds_path, sep = ""))
  # highlight significant aereas
  fig_a <-
    fig_a + geom_rect(
      data = epitopes,
      aes(NULL, NULL, xmin = start, xmax = end),
      ymin = -Inf,
      ymax = Inf,
      fill = "red"
    )
  fig_a <-
    fig_a + geom_ribbon(aes(x = position, ymax = selection_smooth, ymin = 1), fill =
                          "grey80")
  if (points) {
    if (gapcolor) {
      fig_a <- fig_a + geom_point(aes(colour = gap), size = 1.7, alpha = 3 / 4)
    } else {
      fig_a <- fig_a + geom_point(size = 1.7, alpha = 3 / 4)
    }
  }
  
  fig_a <-
    fig_a + scale_colour_gradient(low = "#4db898", high = "#ea5420")
  
  
  fig_a <- fig_a + geom_line(aes(y = selection_smooth))
  #fig_a <- fig_a + labs(title=paste(sample_name,"; g_threshold=", gap_threshold,"; window_size=", window_size,sep=""))
  # fig_a <- fig_a + geom_hline(aes(yintercept=1))
  
  return(fig_a)
}


# This function creates a p-value matrix indicating significant positions based
sliding_window <-
  function(dN,
           dS,
           gap_data,
           w_size = 10,
           g_threshold = 0.2) {
    require(zoo)
    dNdS_window <- rep(0, length(dN))
    if (w_size > 0) {
      dN_window <- rollsum(dN, w_size, fill = list(NA, NULL, NA))
      dS_window <- rollsum(dS, w_size, fill = list(NA, NULL, NA))
      dNdS_window[which(dS_window > 0)] <-
        dN_window[which(dS_window > 0)] /
        dS_window[which(dS_window > 0)]
      dNdS_gap <- rollmean(gap_data$V1, w_size, fill = list(NA, NULL, NA))
      dNdS_gap[which(is.na(dNdS_gap))] <- 1
    }
    
    # iterate over window and find significant windows with FDR correction
    dN_all <- sum(dN, na.rm = T)
    dS_all <- sum(dS, na.rm = T)
    window <- rep(1, length(dN))
    for (i in 1:length(dNdS_window)) {
      if (dNdS_gap[i] < g_threshold) {
        if (!is.na(dN_window[i]) && !is.na(dS_window[i])) {
          test_matrix <-  matrix(c(dN_window[i], dN_all,
                                   dS_window[i], dS_all), nrow = 2)
          pval <-
            fisher.test(test_matrix, alternative = c("greater"))[["p.value"]]
          window[i] <- p.adjust(pval, method = "fdr", n = length(window))
        } else {
          window[i] <- 1
        }
      } else {
        # g_threshold
        window[i] <- 1
      }
    }
    return(window)
  }


# This function creates start and end coordinates of clusters froma list of pvalues
cluster_matrix <- function(pvalues,
                           significance_level = 0.05,
                           w_size = 10) {
  iteration <- 1
  i <- 1
  clusters <- NULL
  
  while (!is.na (which(pvalues[i:length(pvalues)] <= significance_level)[1])) {
    # check if there are some significant elements
    if (i == 1) {
      i <-
        which(pvalues[i:length(pvalues)] <= significance_level)[1]  # first sig elem
      clusters$start[iteration] <- i
      i <-
        which(pvalues[i:length(pvalues)] > significance_level)[1] + clusters$start[iteration] - 1 # last sig. elem.
      clusters$end[iteration] <- i
      iteration <- iteration + 1
    } else {
      i <-
        which(pvalues[i:length(pvalues)] <= significance_level)[1] + clusters$end[iteration -
                                                                                    1]
      clusters$start[iteration] <- i
      i <-
        which(pvalues[i:length(pvalues)] > significance_level)[1] + clusters$start[iteration] - 1
      clusters$end[iteration] <- i
      iteration <- iteration + 1
    }
  }
  if (is.null(clusters)) {
    # add some zero values to prevent plotting issues
    clusters$start[1] <- 0
    clusters$end[1] <- 0
    clusters <- as.data.frame(clusters)
  } else {
    clusters <- as.data.frame(clusters)
    clusters$start <- clusters$start - w_size / 2
    clusters$end <- clusters$end + w_size / 2
  }
  return(clusters)
}