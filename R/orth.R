## ----------- orth.R ------------ ##
#                                   #
#      orthG                        #
#      orthP                        #
#      getseqGS                     #
#      subsetGS                     #
#      speciesGS                    #
#      orthology                    #
#                                   #
## ------------------------------- ##


## --------------------------------------- ##
##                   orthG                 ##
## --------------------------------------- ##
#' Infer GS OrthoGroups Within a Set of Species
#' @description Infers GS orthogroups using tree reconciliation
#' @usage orthG(set = "all")
#' @param set set of species of interest provided as a character vector either with the binomial or short code of the species (see data(sdf)).
#' @details When set = "all", all the species in the database will be included.
#' @return  A list with two elements. The first one is the adjacency matrix (1 for orthologous, 0 for paralogous). The second element is an orthogroup graph.
#' @examples orthG(set = c("Pp", "Psy", "Psm", "Ap"))
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph as_data_frame
#' @importFrom utils data
#' @export

orthG <- function(set = "all"){

  A_selected <- A_selected
  A <- A_selected
  A[is.na(A)] <- 0
  a <- t(A)
  A <- A + a
  sdf <- sdf
  if (set[1] == "all"){
    subA <- A
  } else {
    gsprot <- c()
    for (i in 1:length(set)){
      t <- set[i]
      if (nchar(t) > 5){
        t <- sdf$short[which(sdf$species == t)][1]
      }
      gsprot <- c(gsprot, sdf$Sec.Name_[which(sdf$short == t)])
    }
    subA <- A[which(rownames(A) %in% gsprot), which(colnames(A) %in% gsprot)]

  }

  g <- graph_from_adjacency_matrix(subA, mode = "undirected")
  plot(g)

  return(list(subA, g))
}

## --------------------------------------- ##
##                   orthP                 ##
## --------------------------------------- ##
#' Search Orthologous of a Given Protein
#' @description Searchs orthologous of a given protein within a set of selected species
#' @usage orthP(phylo_id, set = "all")
#' @param phylo_id phylo_id of the query protein
#' @param set set of species of interest provided as a character vector, either with the binomial or short code of the species (see details).
#' @details When set = "all", the search will be carry out against all the species in the database.
#' @return A list with thee elements: 1. subtree of the relevant proteins; 2. vector color; 3. phylo_ids of the orthologous found.
#' @examples orthP(phylo_id = "Pp_GS1a", set = c("Pp", "Psy", "Psm", "Ap"))
#' @importFrom ape read.tree
#' @importFrom ape getMRCA
#' @importFrom TreeTools Subtree
#' @importFrom TreeTools Preorder
#' @importFrom utils data
#' @export

orthP <- function(phylo_id, set = "all"){

  if (set[1] == "all"){
    sdf <- sdf
    set <- unique(sdf$short)
  } else {
    ## Make sure to include the query into the set:
    query <- strsplit(phylo_id, split = "_")[[1]][1]
    set <- unique(c(query, set))
  }

  selected_tr <- selected_tr
  tr <- Preorder(selected_tr)

  A <- orthG(set = set)[[1]]
  A <- A[ , colnames(A) == phylo_id]

  setProt <- names(A)[which(A == 1)]
  if (length(setProt) == 0){
    return(paste("No orthologs of ", phylo_id, " have been detected", sep = ""))
  }
  setProt <- c(setProt, phylo_id)

  lca <- getMRCA(tr, setProt)
  str <- Subtree(tr, lca)
  col <- rep("black", length(str$tip.label))
  col[which(str$tip.label %in% setProt)] <- "blue"
  col[which(str$tip.label == phylo_id)] <- "red"

  output <- list(str, col, setProt)
  attr(output, "phylo_id") <- phylo_id
  plot(str, tip.color = col)

  return(output)

}


## --------------------------------------- ##
##                 getseqGS                ##
## --------------------------------------- ##
#' Get the GS Sequence
#' @description Provides the requested GS sequence
#' @usage getseqGS(phylo_id, molecule = "Prot")
#' @param phylo_id the unique sequence identifier
#' @param molecule either "Prot" or "CDS"
#' @details The identifier should be one of the 'phylo_id' from data(agf).
#' @return  The requested sequence as a character string.
#' @examples getseqGS("Pp_GS1b_2")
#' @importFrom utils data
#' @export

getseqGS <- function(phylo_id, molecule = "Prot"){

  agf <- agf
  if (! phylo_id %in% agf$phylo_id){
    return("Sorry, the id you have provided has not been found in our database")
  } else {
    if (molecule == "CDS"){
      output <- agf$dna[which(agf$phylo_id == phylo_id)]
    } else {
      output <- agf$prot[which(agf$phylo_id == phylo_id)]
    }
    if (!is.na(output)){
      attr(output, "phylo_id") <- phylo_id
      return(output)
    } else {
      return("Sorry, no sequence could be retrieved!")
    }
  }
}

## --------------------------------------- ##
##                 subsetGS                ##
## --------------------------------------- ##
#' GS Proteins Report
#' @description Assembles a report regarding the GS proteins found in the indicated subset of species
#' @usage subsetGS(sp)
#' @param sp set of species of interest (either binomial or short code name)
#' @details This function returns the protein and DNA sequences of the different isoforms found in each species, along with other relevant data.
#' @return  A dataframe with the information for the requested species.
#' @examples subsetGS(c("Pinus pinaster", "Ath"))
#' @importFrom utils data
#' @export

subsetGS <- function(sp){

  agf <- agf
  output <- agf[which(agf$species %in% sp | agf$short %in% sp), ]
  absent <- c()
  for (t in sp){
    if (t %in% agf$short | t %in% agf$species) {
      # do nothing
    } else {
      absent <- c(absent, t)
    }
  }

  attr(output, "absent") <- absent
  if (length(absent) == 1){
    warning(paste("The following species has not been found in our database:  ", absent))
  } else if (length(absent) > 1){
    species_not_found = paste(absent, collapse = " , ")
    warning(paste("The following species have not been found in our database:  ", species_not_found))
  }
  return(output)
}

## --------------------------------------- ##
##                 speciesGS               ##
## --------------------------------------- ##
#' Map Species Names
#' @description Map binomial species name to short code species name and vice versa
#' @usage speciesGS(sp)
#' @param sp set of species of interest (either binomial or short code name)
#' @details The species set should be given as a character vector (see example)
#' @return  A datafrane containing the information for the requested species.
#' @examples speciesGS(c("Pinus pinaster", "Ath"))
#' @export

speciesGS <- function(sp){
  agf <- agf
  usp <- unique(agf$species)
  names(usp) <- unique(agf$short)
  l <- length(sp)
  output <- data.frame(input = sp, output = NA)
  for (i in 1:length(sp)){
    if (sp[i] %in% usp){
      output$output[i] <- names(usp[which(usp == sp[i])])
    } else if (sp[i] %in% names(usp)){
      output$output[i] <- usp[which(names(usp) == sp[i])]
    } else {
      output$output[i] <- "species not found in out database"
    }
  }
  return(output)
}

## --------------------------------------- ##
##               orthology                 ##
## --------------------------------------- ##
#' Infer OrthoGroups Using Tree Reconciliation
#' @description Infer orthogroups using species and gene trees reconciliation
#' @usage orthology(trees, invoke, d = 2, t = 10, l = 1, plot = TRUE, saverec = FALSE)
#' @param trees path to a single file containing first the species tree, followed by a  single gene/protein tree (see details).
#' @param invoke character string representing the way in which the executable of RANGER-DTL (see details) is invoked.
#' @param d cost assigned to gene duplication.
#' @param t cost assigned to gene transfer.
#' @param l cost assigned to gene loss.
#' @param plot when TRUE, the orthology network graph is plotted.
#' @param saverec path to the directory where to save the reconcilation file. If not provided the file is not saved (default)
#' @details The executable of RANGER-DTL (https://compbio.engr.uconn.edu/software/RANGER-DTL) should be installed. All input trees must be expressed using the Newick format terminated by a semicolon, and they must be fully binary (fully resolved) and rooted. Species names in the species tree must be unique. E.g, E.g., (((speciesA_gene1, speciesC_gene1), speciesB_geneX), speciesC_gene2); and (((speciesA, speciesC), speciesB), speciesC); are both valid gene tree inputs and, in fact, represent the same gene tree. This gene tree contains one copy of the gene from speciesA and speciesB, and two copies from speciesC.
#' @return  A list with four elements. The first one is a 'phylo' object where the nodelabels indicate the event: D, duplication or T transfer. If no label is shown is because the event correspond to speciation. The second element is a dataframe (the first column is the label of the internal nodes in the gene tree; the second column is the label of the internal nodes in the species tree, and the third and fourth columns label each internal node according to the inferred event). The third element of the list is an adjacency matrix: 1 when two proteins are orthologous, 0 if they are paralogous. The last element of the list is an orthogroup graph.
#' @examples orthology(trees = system.file("extdata", "input.trees", package = "orthGS"))
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph as_data_frame
#' @importFrom utils data
#' @export

orthology <- function(trees, invoke = "Ranger-DTL.mac", d = 2, t = 10, l = 1, plot = TRUE, saverec = FALSE){

  if (!is.logical(saverec)){
    cmd <- paste(invoke, " -i ", trees, " -D ", d, " -T ", t, " -L ", l, " -o ", saverec, sep = "")
    system(cmd)
    m <- mapTrees(saverec)
  } else {
    recfile <- paste(tempdir(), 'tempFile', sep = "")
    on.exit(unlink(recfile))
    cmd <- paste(invoke, " -i ", trees, " -D ", d, " -T ", t, " -L ", l, " -o ", recfile, sep = "")
    system(cmd)
    m <- mapTrees(recfile)
  }


  A <- m[[3]]
  A[is.na(A)] <- 0
  a <- t(A)
  A <- A + a
  g <- graph_from_adjacency_matrix(A, mode = "undirected")
  if (plot){
    plot(g)
  }

  output <- list(m[[1]], m[[2]], m[[3]], g)
  names(output) <- c("phylo", "map", "adjacency", "graph")
  return(output)
}
