# Description:
#     Generates a dataframe with genes included in given datasets but
#     excludes genes which are part of the sets specified in the exclude arg.
#     
# Input:  myset             Data frame with DE genes
#         venn              Venn-set of the datasets. Keeps track of overlaps
#         normalized_counts Data frame with normalized counts
#         meta              Data frame with meta data.
#         mysubset          Optional arg, string specifiying set from venn
#         
# Output: out               Data frame with DE genes that are found in the 
#                           included data sets but not in the excluded.
# 
# 2023-03-24 Erik Bjurström (erikbju@chalmers.se)

heatmap_custom <- function(myset, venn, normalized_counts, meta, mysubset = "") {
  if (mysubset == ""){
    tmp <- myset
  } else {
    tmp <- subset(myset, myset$Gene_stable_ID %in% getSet(venn, subset = mysubset)$Detail) 
  }
  filtered_tmp <- tmp %>% subset(startsWith(Gene_stable_ID, "Y") & baseMean > 100) %>%   mutate(genelabel = ifelse(Gene_name == "" | Gene_name == "NA" 
                                                                                                                   | is.na(Gene_name),Gene_stable_ID, Gene_name))
  
  mat <- normalized_counts[filtered_tmp$Gene_stable_ID,] %>% apply(1,scale) %>% t()
  colnames(mat) <- rownames(meta)
  
  h <- pheatmap(mat, cluster_rows = T, cluster_cols = T, labels_col =  colnames(mat), labels_row = filtered_tmp$genelabel, scale = "row", fontsize = 15)
  return(h)
}