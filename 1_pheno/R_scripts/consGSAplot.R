# Description:
#     Makes a plot showing the most significant change genesets, as identified
#     from consensus geneset analysis using Piano. For the highest ranking
#     genesets, a barchart is plotted showing the number of up- and down
#     regulated genes. A p-value cutoff is used to signify what genes have
#     significantly changed expression.
#
# Input:      resList   result list from consensus GSA, ideally from consGSA.R
#             rankScore cutoff for non-directional consensus score, whether
#                       genesets should be included in barchart.
#             Pcutoff   cutoff for adjusted P-values of differentially expressed
#                       genes. Recommended is 0.05.
#             distinct  'distinct' if distinct-up and -down should be used for rank
#                       cutoffs. 'distinct' will give less results, but more clearly
#                       up or down regulated. 'mixed' will give more results, but
#                       this will also include GO-terms with mixed regulation.
#             savePlot  true if plot should be saved as file, false if not
#             title     short title used in plot and filename
# Output:     plot      the plot
#
# 2016-11-28 EJK: Separate running of GSA from plotting. Use consGSA.R first.
# 2016-08-22 EJK: Made some adjustments to automatically resize PDF depending
#                 on number of significant GO terms.
# 2016-05-16 Eduard Kerkhoven (eduardk@chalmers.se)

addline_format <- function(x,...){
  
}

consGSAplot <-
  function(resList,
           rankScore,
           Pcutoff,
           distinct,
           savePlot,
           title) {
    library(ggplot2)
    library(piano)
    library(parallel)
    library(snowfall)
    library(tidyr)
    library(scales)
    library(stringr)
    
    # Extract non-directional, distinct up and down genesets.
    non <- consensusScores(resList, class = "non", plot = F)
    df <- data.frame(Name=rownames(non$rankMat[non$rankMat[,1] < rankScore + 1,])) # Select any non-directional with rank =< rankScore.
    
    Pval <- resList[[1]]$geneLevelStats
    FC <- resList[[1]]$directions
    GS <- resList[[1]]$gsc

    
if (distinct=='distinct') {
  up <-
    consensusScores(resList,
                    class = "distinct",
                    direction = "up",
                    plot = F)
  dn <-
    consensusScores(resList,
                    class = "distinct",
                    direction = "down",
                    plot = F)
  non <-
    non$rankMat[non$rankMat[, "ConsScore"] < rankScore + 1 , , drop = F]
  dn <-
    dn$rankMat[dn$rankMat[, "ConsScore"] < rankScore + 1, , drop = F]
  up <-
    up$rankMat[up$rankMat[, "ConsScore"] < rankScore + 1, , drop = F]
  dn <- dn[rownames(dn) %in% rownames(non), , drop = F]
  up <- up[rownames(up) %in% rownames(non), , drop = F]

  df <- data.frame(Name = c(rownames(up), rownames(dn)))
}

    sumTable<-data.frame(Name=names(resList[[1]]$gsc),up=resList[[1]]$nGenesUp,dn=resList[[1]]$nGenesDn,tot=resList[[1]]$nGenesTot)
    df <- merge(df, sumTable)
    
    colnames(df) <- c("geneset", "up", "dn", "tot")
    
    df$highup <- 0
    df$highdown <- 0
    df$lowup <- 0
    df$lowdown <- 0
    
    tmp <- GS[names(GS) %in% df$geneset]
    tmp <- tmp[order(names(tmp))]
    
    for (gset in 1:length(tmp)) {
      df[gset,]$highup <- sum(tmp[[gset]] %in% names(Pval[Pval <
                                                            Pcutoff,]) &
                                tmp[[gset]] %in% names(FC[FC > 0,]))  # Highly significant up
      df[gset,]$highdown <- sum(tmp[[gset]] %in% names(Pval[Pval <
                                                              Pcutoff,]) &
                                  tmp[[gset]] %in% names(FC[FC < 0,]))  # Highly significant down
      df[gset,]$lowup <- sum(tmp[[gset]] %in% names(Pval[Pval >
                                                           Pcutoff,]) &
                               tmp[[gset]] %in% names(FC[FC > 0,]))  # Low significant up
      df[gset,]$lowdown <- sum(tmp[[gset]] %in% names(Pval[Pval >
                                                             Pcutoff,]) &
                                 tmp[[gset]] %in% names(FC[FC < 0,]))  # Low significant down
    }
    
    df$geneset <- factor(df$geneset, levels = df[order(df$up / df$dn),
                                                 1])  # Order from mostly up to mostly down.
    df$geneset<-gsub('(.{50})(.+)','\\1...',df$geneset) # I don't want no damn ellipses // Erik B
    # df$geneset <- str_wrap(df$geneset, width = 50)
    
    # Add ALL genes
    df3 <- data.frame(geneset = 'All')
    df3$geneset <- 'All'
    df3$up <- sum(FC > 0)
    df3$dn <- sum(FC < 0)
    df3$tot <- length(FC)
    df3$highup <- sum(Pval < Pcutoff & FC > 0) 
    df3$highdown <- sum(Pval < Pcutoff & FC < 0) 
    df3$lowup <- sum(Pval > Pcutoff & FC > 0) 
    df3$lowdown <- sum(Pval > Pcutoff & FC < 0)
    df<-rbind(df,df3)
    
    df2 <- gather(df[, c("geneset", "highup", "highdown", "lowup",
                         "lowdown")], "directAndSignif", "genes", 2:5)
    
    df2$geneset <-
      factor(df2$geneset, levels = df2[order(df$up / df$dn),
                                       1])  # Order from mostly up to mostly down.
    df2$directAndSignif <-
      factor(df2$directAndSignif,
             levels = c("highup",
                        "lowup", "lowdown", "highdown")) # Pointless, order is not maintained by ggplot if stat="identity" is used

    
    
        
    # Very inelegant way to order significance and direction of changed genes, instead of order mentioned above
    df2$order <- 0
    df2$order[df2$directAndSignif == "highup"] <- 1
    df2$order[df2$directAndSignif == "lowup"] <- 2
    df2$order[df2$directAndSignif == "lowdown"] <- 3
    df2$order[df2$directAndSignif == "highdown"] <- 4
    df2 <- df2[order(df2$order),]

    grht <-
      3 + (dim(df)[1] * 0.5) # Determine height of output graph in cm, 3 cm and additional 0.5 per GO term
    
    # df[nrow(df)+1, ] <- NA
    # df2[nrow(df2)+1, ] <- NA
    # levels(df2$geneset) <- c(levels(df2$geneset), "")
    # df2$geneset[is.na(df2$geneset)] <- ""
    # levels(df$geneset) <- c(levels(df$geneset), "")
    # df$geneset[is.na(df$geneset)] <- ""
    
    
    plot <-
      
      
      ggplot(df2, aes(x = geneset, y = genes, fill = directAndSignif)) +
      geom_bar(position = "fill", stat = "identity") + coord_flip() +
      scale_y_continuous(labels = percent_format(), breaks = c(0,
                                                               0.25, 0.5, 0.75, 1)) +
      expand_limits(x = length(levels(df2$geneset)) + 1, y = 1.1) +
      scale_fill_manual(
        values = c("#963836", "#b57372", "#abbdd2", "#4682b4"),
        name = "Direction",
        labels = c(
          paste0("Up (p<", Pcutoff, ")"),
          paste0("Up (p>", Pcutoff, ")"),
          paste0("Down (p>", Pcutoff, ")"),
          paste0("Down (p<", Pcutoff, ")")
        )
      ) + labs(x = "GO terms", y = "Genes") +
      theme_set(theme_bw()) + theme(
        text = element_text(size = 8),
        legend.position = c(0,1.01),
        legend.direction = "horizontal",
        legend.justification = "left",
        #legend.background = element_rect(fill = alpha("white", 0.1)),
        #legend.key.size = unit(7, "points"),
        axis.text = element_text(colour = "black"),
        panel.grid = element_blank(),
        line = element_line(linewidth = 0.25),
        plot.margin = unit(c(10,10,10, 10), "points"),
        legend.key=element_blank(),
        panel.border=element_blank(),
        axis.ticks.y=element_blank()
      ) + geom_text(
        data = df,
        aes(x = geneset,
            y = 1.06, label = tot),
        inherit.aes = FALSE,
        size = 3
      ) + annotate(
        "text", x = length(df$geneset) + 0.8, y = 1.06, label = "Total genes",
        size = 3
      ) +
      ggtitle(title)
    
    
    if (savePlot == T) {
      plot + ggsave(
        file = paste0("consGSAplot_", title, ".pdf"),
        device = "pdf",
        height = grht,
        units = "cm"
      )
    }
    return(plot)
  }