## ---- parseRec.R--------------- ##
#                                  #
#      coltips                     #
#      mapTrees                    #
#                                  #
## ------------------------------ ##

## --------------------------------------- ##
##                coltips                  ##
## --------------------------------------- ##
#' Colouring Tree Tips
#' @description Make a color vector for colouring tree tips
#' @usage coltips(phy)
#' @param phy tree as a phylo object
#' @details Each tip is given a color according to the nature of the isoform: green (GS2), blue (GS1a), brown (GS1b Gym), salmon (GS1b Ang), purple (other).
#' @return a color vector as long as the number of tips
#' @examples coltips(ape::read.tree(text = "((Bdi, Sly), (Pp, Ap));"))
#' @importFrom phangorn as.phyDat
#' @importFrom utils data
#' @export

coltips <- function(phy){
  # Make sure that load("./AngGymFern.Rda") has been executed
  agf <- agf
  tr <- phy
  col <- numeric(length(tr$tip.label))
  for (i in 1:length(tr$tip.label)){
    if (grepl("GS2", tr$tip.label[i])){
      col[i] <- "darkgreen"
    } else if (grepl("GS1a", tr$tip.label[i])){
      col[i] <- "blue"
    } else if (grepl("GS1b", tr$tip.label[i]) &
               strsplit(tr$tip.label[i], split = "_")[[1]][1] %in% agf$short[which(agf$taxon == "Acrogymnospermae")]){
      col[i] <- "brown"
    } else if (grepl("GS1b", tr$tip.label[i]) &
               strsplit(tr$tip.label[i], split = "_")[[1]][1] %in% agf$short[which(agf$taxon == "Angiospermae")]){
      col[i] <- "salmon"
    } else {
      col[i] <- "purple"
    }
  }
  return(col)
}

## --------------------------------------- ##
##                mapTrees                 ##
## --------------------------------------- ##
#' Map Gene Tree into Species Tree
#' @description Maps a gene/protein tree into a species tree
#' @usage mapTrees(path2rec)
#' @param path2rec path to the file containing the reconciliation output.
#' @details Mapping gene tree into species tree allow to infer the sequence of events (Duplication, Speciation, Transfer).
#' @return A list with three elements. The first one is a 'phylo' object where the nodelabels indicate the event: D, duplication or T transfer. If no label is shown is because the event correspond to speciation. The second element is a dataframe (the first column is the label of the internal nodes in the gene tree; the second column is the label of the internal nodes in the species tree, and the third and fourth columns label each internal node according to the inferred event). The third element of the list is an adjacency matrix: 1 when two proteins are orthologous, 0 if they are paralogous.
#' @examples mapTrees(fs::path_package("extdata", "representatives", package = "orthGS"))
#' @importFrom ape read.tree
#' @importFrom castor get_pairwise_mrcas
#' @export

mapTrees <- function(path2rec){

  ## --- Parsing reconciliation file
  con <- file(path2rec, "r")
  nodes <- c()
  while(TRUE){
    line = readLines(con, n = 1)
    if (length(line) == 0){
      break
    }
    nodes <- c(nodes, line)
  }
  close(con)
  lca <- nodes[grepl("LCA\\[", nodes)]

  ## --- Gene family tree
  gtr <- read.tree(text = nodes[grep("Gene Tree: ", nodes) + 1])

  ## --- Maping nodes
  m <- unlist(lapply(lca, function(x) strsplit(x, split = " ")[[1]][1]))
  z <- as.numeric(substr(m, 2, nchar(m)))
  z <- z[order(z)]

  ndf <- data.frame(n = (gtr$Nnode+2):(2*gtr$Nnode+1),
                    m = paste("m", z, " ", sep = ""),
                    event = NA,
                    label = NA)

  for (i in 1:nrow(ndf)){
    at <- grep(ndf$m[i], lca)
    if (grepl("Speciation", lca[at])){
      ndf$event[i] <- "Speciation"
      ndf$label[i] <- ""
    } else if (grepl("Duplication", lca[at])){
      ndf$event[i] <- "Duplication"
      ndf$label[i] <- "D"
    } else if (grepl("Transfer", lca[at])){
      ndf$event[i] <- "Transfer"
      ndf$event[i] <- "T"
    }
  }

  ## --- Adjacency matrix
  n <- length(gtr$tip.label)
  A <- matrix(rep(NA,n*n), ncol = n)
  rownames(A) <- colnames(A) <- gtr$tip.label

  for (i in 1:(n-1)){
    t <- rownames(A)[i]
    t_ <- strsplit(t, "_")[[1]][1]

    for (j in (i+1):n){
      vs <- colnames(A)[j]
      vs_ <- strsplit(vs, "_")[[1]][1]
      tnode <- get_pairwise_mrcas(gtr, t, vs)
      tevent <- ndf$event[which(ndf$n == tnode)]

      if (tevent == "Duplication"){
        A[i,j] <- 0
      } else if (tevent == "Speciation"){
        A[i,j] <- 1
      } else if (tevent == "Transfer"){ # Se adscribirá como ortólogo (ortólogo por transferencia horizontal)
        A[i,j] <- 1
      }
    }
  }
  gtr$node.label <- ndf$label

  return(list(gtr, ndf, A))
}
