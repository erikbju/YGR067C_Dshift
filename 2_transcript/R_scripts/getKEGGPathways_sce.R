## Copyright 2010 Laurent Jacob, Pierre Neuvial and Sandrine Dudoit.

## This file is part of DEGraph.

## DEGraph is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## DEGraph is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with DEGraph.  If not, see <http://www.gnu.org/licenses/>.

#########################################################################/**
## @RdocFunction getKEGGPathways
##
## @title "Builds a graph for each of the KEGG pathways"
##
## \description{
##  @get "title".
## }
##
## @synopsis
##
## \arguments{
##   \item{path}{A @character value, the local _full_ path of KGML data.}
##   \item{rootPath}{A @character value, the local _root_ path of KGML data.}
##   \item{organism}{A @character value specifying the organism whose
##      pathways should be considered. Defaults to "hsa" (Homo Sapiens).}
##   \item{metaTag}{A @character value, specifying the type of pathways to
##     be considered ("metabolic" or "non-metabolic"). Defaults to "non-metabolic".}
##   \item{pattern}{An optional @character value specifying a file name pattern to
##     look for.}
##   \item{verbose}{If @TRUE, extra information is output.}
## }
##
## \value{
##   A @list containing a \code{\link[=graph-class]{graph}} object for each KEGG pathway with at least one edge.
## }
##
## \details{If 'path' is supplied, KGML files in this directory are loaded.
##   Otherwise, KGML files are assumed to be in
##   <rootPath>/<metaTag>/"organisms"/<organism>, which mirrors the
##  structure of the KEGG KGML file repository.
## }
##
## @author
##
## \seealso{
##   @see "parseKGML"
##   @see "KEGGpathway2Graph"
## }
##
## @examples "../incl/getKEGGPathways.Rex"
##
##*/###########################################################################

getKEGGPathways_sce <- function(path=NULL, rootPath="networkData/ftp.genome.jp/pub/kegg/xml/kgml", organism="sce", metaTag=c("non-metabolic", "metabolic"), pattern=NULL, verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  if (is.null(path)) {
    rootPath <- Arguments$getReadablePath(rootPath)
    metaTag <- match.arg(metaTag)
    organism <- match.arg(organism)
    path <- file.path(rootPath, metaTag, "organisms", organism)
  }
  path <- Arguments$getReadablePath(path)
  # Argument 'pattern':
  pattern <- Arguments$getCharacter(pattern)
  
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    cat <- R.utils::cat
    pushState(verbose)
    on.exit(popState(verbose))
  } 
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Setup
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  if (length(pattern)) {
    patt <- sprintf("(%s)", pattern)
  } else {
    patt <- paste("^", organism, "([0-9]+).xml$", sep="")
  }
  
  filenames <- list.files(path, pattern=patt)
  if (!length(filenames)) {
    throw("No pathway found with pattern :'", patt, "' in directory ", path)
  }
  print(patt)
  pIds <- gsub(patt, "\\1", filenames)  ## pathway IDs
  verbose && cat(verbose, "Pathway IDs:")
  verbose && str(verbose, pIds)
  
  pathnames <- file.path(path, filenames)
  names(pathnames) <- pIds
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## get all KEGG pathways
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  pb <- ProgressBar(stepLength=100/length(pathnames))
  reset(pb)
  
  grList <- lapply(pathnames, FUN=function(pathname) {
    pw <- parseKGML(pathname)
    pwInfo <- getPathwayInfo(pw)
    gr <- KEGGpathway2Graph(pw, genesOnly=TRUE, expandGenes=TRUE)
    #attr(gr, "info") <- pwInfo
    gr@graphData$info <- pwInfo
    #attr(gr, "label") <- getTitle(pwInfo)
    gr@graphData$label <- KEGGgraph::getTitle(pwInfo)
    increase(pb)
    gr
  })
  
  verbose && cat(verbose, "KEGG pathways retrieved:")
  verbose && str(verbose, names(grList))
  
  ## remove those with 0 edges
  ne <- sapply(grList, numEdges)
  idxs <- which(ne==0)
  if (length(idxs)) {
    verbose && cat(verbose, "Discarding ", length(idxs), " pathways with 0 edges")
    grList <- grList[-idxs]
  }
  verbose && str(verbose, names(grList))
  
  return(grList)
}

############################################################################
# HISTORY:
## 2010-10-08
## o Now validating argument 'verbose'.
# 2010-09-23
# o Added a 'path' argument to load KGML files sitting in non-standard
##  directories.
# 2010-09-20
# o Added a 'pattern' argument.
# 2010-09-14
# o Clean-ups.
# o Now returning only a list of 'graph' elements (labels are passed as
#   attributes.
# o pathway labels are now inferred from the XML file, not from the
#   (unreliable) "KEGGPATHID2NAME" environment.
############################################################################
