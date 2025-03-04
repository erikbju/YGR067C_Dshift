# Description:
#     Performs consensus GSA using five common (and relatively fast) methods that are
#     included in Piano: "mean", "median", "sum", "stouffer" and "tailStrength".
#     Runs on maximum number of CPU nodes. Gives resList as output, which can be
#     used to make a plot using consGSAplot.
#
# Input:      Pval      named vector of gene-level adjusted P-values
#             FC        named vector of gene-level log2 fold-changes
#             GS        gene-set, as loaded usign Piano
# Output:     resList   results list
#
# 2016-11-28 EJK: Separate funtion to run consensus GSA.
# 2016-05-16 Eduard Kerkhoven (eduardk@chalmers.se)



consGSA <-
  function(Pval,
           FC,
           GS, gsSizeDn, gsSizeUp) {
    require(piano)
    require(parallel)
    require(snowfall)
    require(tidyr)
    #  if (!exists("rankScore")) # Attempt to set default settings, not sure how to do this...
    
    # Find out number of processor cores, run GSA on n-1 cores.
    cores <- as.numeric(detectCores())
    cat("Running GSA 1/8")
    gsaRes1 <- runGSA(
      Pval,
      FC,
      geneSetStat = "mean",
      gsc = GS,
      nPerm = round(1000 / cores) * cores,
      gsSizeLim = c(gsSizeDn, gsSizeUp),
      ncpus = cores
    )
    cat("Running GSA 2/8")
    gsaRes2 <- runGSA(
      Pval,
      FC,
      geneSetStat = "median",
      gsc = GS,
      nPerm = round(1000 / cores) * cores,
      gsSizeLim = c(gsSizeDn, gsSizeUp),
      ncpus = cores
    )
    cat("Running GSA 3/8")
    gsaRes3 <- runGSA(
      Pval,
      FC,
      geneSetStat = "sum",
      gsc = GS,
      nPerm = round(1000 / cores) * cores,
      gsSizeLim = c(gsSizeDn, gsSizeUp),
      ncpus = cores
    )
    cat("Running GSA 4/8")
    gsaRes4 <- runGSA(
      Pval,
      FC,
      geneSetStat = "stouffer",
      gsc = GS,
      nPerm = round(1000 / cores) * cores,
      gsSizeLim = c(gsSizeDn, gsSizeUp),
      ncpus = cores
    )
    cat("Running GSA 5/8")
    gsaRes5 <- runGSA(
      Pval,
      FC,
      geneSetStat = "tailStrength",
      gsc = GS,
      nPerm = round(1000 / cores) * cores,
      gsSizeLim = c(gsSizeDn, gsSizeUp),
      ncpus = cores
    )
    cat("Running GSA 6/8")
    gsaRes6 <- runGSA(
      Pval,
      FC,
      geneSetStat = "wilcoxon",
      gsc = GS,
      nPerm = round(1000 / cores) * cores,
      gsSizeLim = c(gsSizeDn, gsSizeUp),
      ncpus = cores
    )
    cat("Running GSA 7/8")
    gsaRes7 <- runGSA(
      Pval,
      FC,
      geneSetStat = "fisher",
      gsc = GS,
      nPerm = round(1000 / cores) * cores,
      gsSizeLim = c(gsSizeDn, gsSizeUp),
      ncpus = cores
    )
    cat("Running GSA 8/8")
    gsaRes8 <- runGSA(
      FC,
      geneSetStat = "maxmean",
      gsc = GS,
      nPerm = round(1000 / cores) * cores,
      gsSizeLim = c(gsSizeDn, gsSizeUp),
      ncpus = cores
    )
    # No maxmean and fisher: don't support direction. No wilcoxon, too slow.
    cat("Reorganizing data and prepare for plotting")
    # Combine results in list
    resList <- list(gsaRes1, gsaRes2, gsaRes3, gsaRes4, gsaRes5, gsaRes6, gsaRes7,
                    gsaRes8)
    resList <-
      setNames(resList,
               c("mean", "median", "sum", "stouffer",
                 "tailStrength", "wilcoxon", "fisher", "maxmean"))
    
return(resList)
  }