# Description:
#   Attaches 'Gene name' to 'Gene stable ID' if it exists in bioMart's database. 
#   Also attaches description, taken from SGD.
#
# Input:  input       Tibble with contrast results table form DESeq2
#         genelist    Dataframe with gene names and ID
# Output: out_table   Tibble with contrast results table with gene names + descr.
#
# 2023-03-24 Erik Bjurström (erikbju@chalmers.se)

nicefy_gene_results <- function(input, genelist) {
  temp_table <- input
  genes_in_input <- subset(genelist, (genelist$Gene_stable_ID %in% input$Gene_stable_ID))
  
  out_table <- merge(temp_table, genes_in_input, by.x = c("Gene_stable_ID"), by.y = c("Gene_stable_ID"), all = TRUE) %>% 
    relocate(Gene_name, .after = Gene_stable_ID)
  
  return(out_table)
}