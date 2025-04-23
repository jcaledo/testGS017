
## ----------- ancillary.R ------------ ##
#                                        #
#      msa                               #
#      mltree                            #
#      gapless_msa                       #
#      madRoot                           #
#                                        #
## ------------------------------------ ##


## ---------------------------------------------------------------- ##
#     msa <- function(sequences, ids, seqtype, sfile, inhouse )      #
## ---------------------------------------------------------------- ##
#' Multiple Sequence Alignment
#' @description Aligns multiple protein, DNA or CDS sequences.
#' @usage msa(sequences, ids = names(sequences), seqtype = "prot", sfile = FALSE, inhouse = FALSE)
#' @param sequences vector containing the sequences as strings.
#' @param seqtype it should be either "prot" of "dna" or "cds" (see details).
#' @param ids character vector containing the sequences' ids.
#' @param sfile if different to FALSE, then it should be a string indicating the path to save a fasta alignment file.
#' @param inhouse logical, if TRUE the in-house MUSCLE software is used. It must be installed on your system and in the search path for executables.
#' @details If seqtype is set to "cds" the sequences must not contain stop codons and they will be translated using the standard code. Afterward, the amino acid alignment will be used to lead the codon alignment.
#' @return Returns a list of four elements. The first one ($seq) provides the sequences analyzed, the second element ($id) returns the identifiers, the third element ($aln) provides the alignment in fasta format and the fourth element ($ali) gives the alignment in matrix format.
#' @examples msa(sequences = c("APGW", "AGWC", "CWGA"),ids = c("a", "b", "c"))
#' @importFrom bio3d seqbind
#' @importFrom bio3d seqaln
#' @importFrom bio3d write.fasta
#' @importFrom seqinr translate
#' @export

msa <- function (sequences, ids = names(sequences), seqtype = "prot", sfile = FALSE, inhouse = FALSE){

  tr <- function(seq, genetic_code = 1){
    seq <- gsub(" ", "", seq)
    seq <- strsplit(seq, "")[[1]]
    output <- paste(translate(seq, numcode = genetic_code), collapse = "")

    return(output)
  }

  if (length(sequences) < 2) {
    stop("At least two sequences are required!")
  } else if (length(sequences) != length(ids)) {
    stop("The number of sequences and sequences' ids doesn't match!")
  }
  if (inhouse) {
    if (seqtype == "cds"){
      dnaSeq <- sequences
      cod <- strsplit(gsub("(.{3})", "\\1 ", dnaSeq), split = " ")
      sequences <- unlist(lapply(sequences, function(x) tr(x)))
    }
    seqs <- lapply(sequences, function(x) strsplit(x, split = "")[[1]])
    sqs <- bio3d::seqbind(seqs[[1]], seqs[[2]], blank = "-")
    c <- 2
    while (c < length(sequences)) {
      c <- c + 1
      sqs <- bio3d::seqbind(sqs, seqs[[c]], blank = "-")
    }
    aln <- bio3d::seqaln(sqs, id = ids, exefile = "muscle")
    aln$seq <- sequences
    if (seqtype == "cds"){
      aln$cod <- aln$ali
      for (i in 1:length(cod)){
        contador <- 1
        for (j in 1:ncol(aln$ali)){
          if (aln$ali[i,j] != "-"){
            aln$cod[i,j] <- cod[[i]][contador]
            contador <- contador + 1
          } else {
            aln$cod[i,j] <- "---"
          }
        }
      }
    }
    if (sfile != FALSE & seqtype == "cds") {
      for (i in 1:length(ids)){
        t <- paste(">", ids[i], sep = "")
        cat(t, "\n", file = sfile, append = TRUE)
        tt <- paste(aln$cod[i,], collapse = "")
        cat(tt, "\n", file = sfile, append = TRUE)
      }
    } else if (sfile != FALSE & seqtype != "cds") {
      bio3d::write.fasta(aln, file = sfile)
    }
    if (file.exists("aln.fa")) {
      system("rm aln.fa")
    }
    return(aln)
  }  else {
    ## --- Using the Biostring and muscle R packages
    # seq <- Biostrings::AAStringSet(sequences)

    if (requireNamespace('Biostrings', quietly = TRUE)){
      if (seqtype == "prot"){
        seq <- Biostrings::AAStringSet(sequences)
      } else if (seqtype == "dna"){
        seq <- Biostrings::DNAStringSet(sequences)
      }
    } else {
      stop("You must install the package Biostrings in order to use this function")
    }

    if (requireNamespace('muscle', quietly = TRUE)){
      aln1 <- muscle::muscle(seq)
    } else {
      stop("You must install the package muscle in order to use this function")
    }
    aln <- list()
    aln$seq <- sequences
    aln$ids <- ids
    aln$aln <- as.character(aln1)
    l <- sapply(aln$aln, function(x) strsplit(x, split = ""))
    aln$ali <- matrix(unlist(l), nrow = length(sequences),
                      byrow = TRUE)
    if (sfile != FALSE) {
      for (i in 1:length(aln$aln)) {
        t <- paste(">", aln$ids[i], sep = "")
        cat(t, file = sfile, append = TRUE)
        if (seqtype == "cds"){
          tt <- paste("\n", aln$cod[i], "\n", sep = "")
        } else {
          tt <- paste("\n", aln$aln[i], "\n", sep = "")
        }
        cat(tt, file = sfile, append = TRUE)
      }
    }
    return(aln)
  }
}

## ---------------------------------------------------------------- ##
#                     mltree <- function()                           #
## ---------------------------------------------------------------- ##
#' Build Up a ML Tree
#' @description Given an alignment builds an ML tree.
#' @usage mltree(msa, df = TRUE, gapl = TRUE, model = "WAG")
#' @param msa input alignment.
#' @param df logical. When TRUE msa should be a dataframe, when FALSE msa should be a string giving the path to a fasta file containing the alignment.
#' @param gapl logical, when TRUE a gapless alignment is used.
#' @param model allows to choose an amino acid models (see the function phangorn::as.pml)
#' @details The function makes a NJ tree and then improvove it using an optimization procedure based on ML.
#' @return a ML optimized tree (and parameters)
#' @examples
#' a <- msa(sequences=c("RAPGT", "KMPGT", "ESGGT"), ids = letters[1:3])$ali
#' rownames(a) <- letters[1:3]
#' tr <- mltree(a)$tree
#' @seealso gapless_msa
#' @importFrom ape nj
#' @importFrom ape dist.aa
#' @importFrom phangorn as.phyDat
#' @importFrom phangorn pml
#' @importFrom phangorn optim.pml
#' @export

mltree <- function(msa, df = TRUE, gapl = TRUE, model = "WAG"){
  if (df == TRUE){
    aln <- msa
  } else {
    aln <- bio3d::read.fasta(msa)$ali
  }
  if (gapl == TRUE){
    data <- gapless_msa(aln)
  } else {
    data <- aln
  }

  tre.ini <- nj(ape::dist.aa(data))
  fit.ini <- pml(tre.ini, as.phyDat(as.matrix(data), type = "AA"), model = model)
  fit <- optim.pml(fit.ini, model = model)
  return(fit)
}

## ---------------------------------------------------------------- ##
#                   gapless_msa <- function()                        #
## ---------------------------------------------------------------- ##
#' Remove Gaps in a MSA
#' @description Removes gaps in a given msa.
#' @usage gapless_msa(msa, seqtype = 'AA', df = TRUE, sfile = FALSE)
#' @param msa input alignment.
#' @param seqtype the nature of the sequences: 'DNA' or 'AA'.
#' @param df logical. When TRUE msa should be a matrix, when FALSE msa should be a string giving the path to a fasta file containing the alignment.
#' @param sfile if different to FALSE, then it should be a string indicating the path to save a fasta alignment file.
#' @details It should be noted that this function does not carry out the alignment itself.
#' @return an alignment without gaps in form of matrix or a file containing such an alignment in fasta format.
#' @examples gapless_msa(msa(sequences = c("APGW", "AGWC", "CWGA"),ids = c("a", "b", "c"))$ali)
#' @seealso msa
#' @importFrom seqinr read.fasta
#' @export

gapless_msa <- function(msa, seqtype = "AA", df = TRUE, sfile = FALSE){
  if (df == FALSE){
    msa <- read.fasta(msa, seqtype = seqtype)
    M <- matrix(rep(NA, length(msa)*length(msa[[1]])), nrow = length(msa))
    rownames(M) <- attributes(msa)$names
    for (i in 1:length(msa)){
      M[i,] <- as.character(msa[[i]])
    }
  } else {
    M <- as.matrix(msa)
  }
  gapless <- c()
  for (i in 1:ncol(M)){
    if (! "-" %in% M[,i]){
      gapless <- c(gapless, i)
    }
  }
  output <- M[, gapless]

  if (sfile != FALSE) {
    for (i in 1:nrow(output)) {
      t <- paste(">", rownames(output)[i], sep = "")
      cat(t, file = sfile, append = TRUE)
      tt <- paste(output[i,], collapse = "")
      tt <- paste("\n", tt, "\n", sep = "")
      cat(tt, file = sfile, append = TRUE)
    }
  }

  return(as.data.frame(output))
}

## --------------------------------------- ##
##                   madRoot               ##
## --------------------------------------- ##
#' Find The Root of a Phylogenetic Tree Using MAD Method
#' @description Finds the root of an unrooted phylogenetic tree by minimizing the relative deviation from the molecular clock.
#' @usage madRoot(tree, output_mode = 'phylo')
#' @param tree unrooted tree string in newick format or a tree object of class 'phylo'.
#' @param output_mode amount of information to return. If 'phylo' (default) only the rooted tree is returned. If 'stats' also a structure with the ambiguity index, clock cv, the minimum ancestor deviation and the number of roots. If 'full' also an unrooted tree object, the index of the root branch, the branch ancestor deviations and a rooted tree object.
#' @details This function is a slight modification of the code provided by Tria et al at https://www.mikrobio.uni-kiel.de/de/ag-dagan/ressourcen.
#' @return a rooted tree and supplementary information if required.
#' @author Tria, F. D. K., Landan, G. and Dagan, T.
#' @examples
#' a <- msa(sequences=c("RAPGT", "KMPGT", "ESGGT"), ids = letters[1:3])$ali
#' rownames(a) <- letters[1:3]
#' tr <- mltree(a)$tree
#' rtr <- madRoot(tr)
#' @references Tria, F. D. K., Landan, G. and Dagan, T. Nat. Ecol. Evol. 1, 0193 (2017).
#' @importFrom ape read.tree
#' @importFrom ape is.binary
#' @importFrom ape is.rooted
#' @importFrom ape dist.nodes
#' @importFrom ape write.tree
#' @importFrom ape unroot
#' @importFrom ape multi2di
#' @importFrom ape drop.tip
#' @importFrom phytools reroot
#' @importFrom stats sd
#' @export

madRoot <- function(tree, output_mode = "phylo"){
  unrooted_newick <- tree
  # "Dependencies: 'ape' and 'phytools'","","Version: 1.1, 03-May-2017",sep="\n"))
  t <- NA
  if(inherits(unrooted_newick, "phylo")){
     t <- unrooted_newick
  }
  else{
    t <- read.tree(text = unrooted_newick)
  }
  if(is.rooted(t)){
    t <- unroot(t)
  }
  # t$node.label<-NULL #To allow parsing when identical OTUs are present
  if(!is.binary(t)){
    warning("Input tree is not binary! Internal multifurcations will be converted to branches of length zero and identical OTUs will be collapsed!")
    t <- multi2di(t)
  }
  tf <- t$edge.length < 0
  if(any(tf)){
    warning("Input tree contains negative branch lengths. They will be converted to zeros!")
    t$edge.length[tf] <- 0
  }

  notu <- length(t$tip.label)
  nbranch <- dim(t$edge)[1]
  npairs <- notu*(notu-1)/2
  nodeids <- 1:(nbranch+1)
  otuids <- 1:notu
  dis <- dist.nodes(t) # phenetic distance. All nodes
  sdis <- dis[1:notu,1:notu] # phenetic distance. otus only

  ## --- Start recursion to collapse identical OTUs, if present.
  ii <- which(sdis == 0, arr.ind = TRUE)
  k <- which(ii[,1] != ii[,2])
  if(length(k)){
    r <- ii[k[1],1]
    c <- ii[k[1],2]
    vv <- c(paste('@#',t$tip.label[r],'@#',sep=""),paste('(',t$tip.label[r],':0,',t$tip.label[c],':0)',sep=""))
    st <- drop.tip(t,c)
    st$tip.label[st$tip.label==t$tip.label[r]] <- vv[1]
    res <- madRoot(st,output_mode) # mad -> madRoot
    if(is.list(res)){
      res[[1]] <- sub(vv[1],vv[2],res[[1]])
    }
    else{
      res <- sub(vv[1],vv[2],res)
    }
    return(res) # create the list 'res' to return the results
  }
  ## --- End of recursion
  t2 <- t
  t2$edge.length <- rep(1,nbranch)
  disbr <- dist.nodes(t2) # split distance. All nodes
  sdisbr <- disbr[1:notu,1:notu] # split distance. otus only
  rho <- vector(mode = "numeric",length = nbranch) # Store position of the optimized root nodes (branch order as in the input tree)
  bad <- vector(mode = "numeric",length = nbranch) # Store branch ancestor deviations (branch order as in the input tree)
  i2p <- matrix(nrow = nbranch+1, ncol = notu)
  for (br in 1:nbranch){
    # collect the deviations associated with straddling otu pairs
    dij <- t$edge.length[br]
    if(dij == 0){
      rho[br] <- NA
      bad[br] <- NA
      next
    }
    rbca <- numeric(npairs)
    i <- t$edge[br,1]
    j <- t$edge[br,2]
    sp <- dis[1:notu,i]<dis[1:notu,j] # otu split for 'br'
    dbc <- matrix(sdis[sp,!sp],nrow=sum(sp),ncol=sum(!sp))
    dbi <- replicate(dim(dbc)[2],dis[(1:notu)[sp],i])

    rho[br] <- sum((dbc-2*dbi)*dbc^-2)/(2*dij*sum(dbc^-2)) # optimized root node relative to 'i' node
    rho[br] <- min(max(0,rho[br]),1)
    dab <- dbi+(dij*rho[br])
    ndab <- length(dab)
    rbca[1:ndab] <- as.vector(2*dab/dbc-1)
    ## --- Collect the remaining deviations (non-traversing otus)
    bcsp <- rbind(sp,!sp)
    ij <- c(i,j)
    counter <- ndab
    for (w in c(1,2)){
      if(sum(bcsp[w,])>=2){
        disbrw <- disbr[,ij[w]]
        pairids <- otuids[bcsp[w,]]
        for (z in pairids){
          i2p[,z] <- disbr[z,]+disbrw==disbrw[z]
        }
        for (z in 1:(length(pairids)-1)){
          p1 <- pairids[z]
          disp1 <- dis[p1,]
          pan <- nodeids[i2p[,p1]]
          for (y in (z+1):length(pairids)){
            p2 <- pairids[y]
            pan1 <- pan[i2p[pan,p2]]
            an <- pan1[which.max(disbrw[pan1])]
            counter <- counter+1
            rbca[counter] <- 2*disp1[an]/disp1[p2]-1
          }
        }
      }
    }
    if(length(rbca)!=npairs){
      stop("Unexpected number of pairs. Report this error to ftria@ifam.uni-kiel.de")
    }
    bad[br] <- sqrt(mean(rbca^2)) # branch ancestor deviation
  }
  ## --- Select the branch with the minum ancestor deviation and calculate the root ambiguity index
  jj <- sort(bad,index.return = TRUE)
  tf <- bad == jj$x[1]
  tf[is.na(tf)]<-FALSE
  nroots <- sum(tf)
  if (nroots>1){
    warning("More than one possible root position. Multiple newick strings printed")
  }
  madr <- which(tf) # Index of the mad root branch(es)
  rai <- jj$x[1]/jj$x[2] # Root ambiguity index
  badr <- bad[tf] # Branch ancestor deviations value for the root(s)
  ## --- Root the tree object, calculate the clock CV and retrieve the newick string
  rt <- vector(mode = "list", nroots) # Rooted tree object
  ccv <- vector(mode = "numeric",nroots) # Clock CV
  rooted_newick <- vector(mode = "character", nroots)
  for (i in 1:length(madr)){
    pp <- rho[madr[i]]*t$edge.length[madr[i]]
    nn <- t$edge[madr[i],]
    # rt[[i]] <- reroot(t,nn[2],pos = pp)
    rt[[i]] <- reroot(t,nn[2],position = pp)
    rooted_newick[i] <- write.tree(rt[[i]])
    dd <- dis[1:notu,nn]
    sp <- dd[,1]<dd[,2]
    otu2root <- vector(mode="numeric",notu)
    otu2root[sp] <- dd[sp,1] + pp
    otu2root[!sp] <- dd[!sp,1] - pp
    ccv[i] <- 100*sd(otu2root)/mean(otu2root)
  }
  rooted_newick <- sub(')Root;',');', rooted_newick)

  ## --- Output the result(s)
  rtr <- read.tree(text = rooted_newick)
  if(output_mode == "phylo")
  {
    return(rtr)
  }
  else if (output_mode == "stats"){
    root_stats <- data.frame(ambiguity_index = rai,
                             clock_cv = ccv,
                             ancestor_deviation = badr,
                             n_roots = nroots)
    return(list(rtr, root_stats))
  }
  else if(output_mode == 'full'){ #Rooted newick, stats, unrooted tree object, index of the branch root, ancestor deviations, rooted tree object
    root_stats <- data.frame(ambiguity_index = rai,
                             clock_cv=ccv,
                             ancestor_deviation=badr,
                             n_roots=nroots)
    return(list(rooted_newick,root_stats,t,madr,bad,rtr))
  }
  else{
    return(rtr)
  }
}

