# Description:
#     Plots the DESeq2 results table in the form of a volcano plot, with adjustable
#     log-fold and padj cutoffs. Top significant genes are labeled (hidden if they
#     overlap).
#     
# Input:  input             Tibble with contrasted results table from DESeq2
#         gene_description  Dataframe with gene names and ID
#         title             String with the plot title, default is name of input
#         pcut              Numeric with the padj cutoff
#         foldcut           Numeric with the log2fold change cutoff
#         labelcount        Numeric with how many of the top significant genes are
#                           labeled.
# Output: p                 ggplot2 object with the custom volcano plot.
# 
# 2023-03-24 Erik Bjurström (erikbju@chalmers.se)

volcano_custom <- function(input, gene_description, title = deparse(substitute(input)), pcut = 0.05, foldcut = 0.58, labelcount = 10){
  volc2 <- input %>%
    data.frame()
  
  volc2<- volc2[abs(volc2$log2FoldChange) < 10,]
  if (!("Gene_name" %in% colnames(volc2))) {
    volc2 <- nicefy_gene_results(volc2, gene_description)
  }
  volc2 <- volc2 %>% 
    filter(!is.na(padj)) %>% 
    mutate(genelabel = ifelse(Gene_name == "" | is.na(Gene_name),
                              Gene_stable_ID, 
                              paste("italic('", Gene_name,"')", sep = ""))) %>% 
    mutate(top = NA) %>% 
    arrange(padj)
  
  
  
  
  volc2$diffexpressed <- "NO"
  volc2$diffexpressed[volc2$log2FoldChange >= foldcut & volc2$padj < pcut] <- "UP"
  volc2$diffexpressed[volc2$log2FoldChange <= -foldcut & volc2$padj < pcut] <- "DOWN"
  
  noes <- sum(volc2$diffexpressed[1:labelcount] == "NO")
  siglabels <- labelcount

  
  
  while (siglabels - noes < labelcount) {
    siglabels = labelcount + noes
    noes <- sum(volc2$diffexpressed[1:siglabels] == "NO")
  }
  if (labelcount == 0) {
    volc2$top[1:siglabels] <- c("")
  } else {
    volc2$top[1:siglabels] <- volc2$genelabel[1:siglabels]
  }
  
  
  
  volc2 <- within(volc2, {
    f <- diffexpressed == "NO"
    top[f] <- NA
  })
  
  colorlist <- c("#CCCCCC", "#CCCCCC","#CCCCCC")
  if("UP" %in% volc2$diffexpressed) {
    colorlist[3] <- "#FF3333"
  }
  if("DOWN" %in% volc2$diffexpressed) {
    colorlist[1] <- "#0000CC"
  }

  
  maxfoldchange <- abs(volc2[which.max(abs(volc2$log2FoldChange)),]$log2FoldChange)
  
  p <- ggplot(volc2 %>% arrange(match(diffexpressed, c("NO", "UP", "DOWN"))), aes(x=log2FoldChange, y=-log10(padj), label = top, col=diffexpressed, te)) +
    geom_point(alpha = 0.7) + 
    theme_minimal() +
    theme(legend.title = element_text(size = 20),legend.text = element_text(size = 20), axis.text = element_text(size=20), axis.title = element_text(size = 20)) +
    labs(color='Differentially \nexpressed') +
    geom_text_repel(parse = T, show.legend =F, max.overlaps = 10, size = 7) +
    ggtitle(title) +
    xlim(-maxfoldchange, maxfoldchange) +
    scale_color_manual(values=colorlist) +
    geom_vline(xintercept=c(-foldcut, foldcut), linetype = "dashed") +
    geom_hline(yintercept=-log10(pcut), linetype = "dashed")
  
  return(p)
}