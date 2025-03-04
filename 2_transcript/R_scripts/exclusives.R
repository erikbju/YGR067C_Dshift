# Description:
#     Generates a dataframe with genes included in given datasets but
#     excludes genes which are part of the sets specified in the exclude arg.
#     
# Input:  ven               Venn-set of the datasets. Keeps track of overlaps
#         sig               Data frame with significant DE genes.
#         comb_list         String list with the combinations of set names
#         include           Which sets to include. String or list.
#         exclude           Which sets to exclude. String or list.
#                           labeled.
# Output: out               Data frame with DE genes that are found in the 
#                           included data sets but not in the excluded.
# 
# 2023-03-24 Erik Bjurström (erikbju@chalmers.se)

exclusives <- function(ven, sig, comb_list, include = "", exclude =""){
  tmp <- comb_list[grep(paste(c(include),collapse="|"), comb_list)]
  
  if (exclude == ""){
    tmp <-  c(tmp, "Shared")
  } else {
    tmp <- tmp[-grep(paste(c(exclude),collapse="|"), tmp)]
  }
  out <- sig %>% subset(Gene_stable_ID %in% getSet(ven, subset = tmp)$Detail)
  return(out)
}