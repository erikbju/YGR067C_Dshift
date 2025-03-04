# Saves a venn diagram

save_ven <- function(g, e, w, t){
  
  
  ven <- venndetail(list("Wildtype vs. tda∆1 \n on glucose" = g$Gene_stable_ID, "Wildtype vs. tda∆1 \n on ethanol" = e$Gene_stable_ID
                         , "Glucose vs. Ethanol \n on Wildtype" = w$Gene_stable_ID, "Glucose vs. Ethanol \n for tda∆1" = t$Gene_stable_ID))
  png(filename = "verbose_venn.png", res = 300, height = 2000, width = 2000)
  print(plot(ven, cat.cex = 1))
  dev.off()
  print("done")
  return(ven)
}