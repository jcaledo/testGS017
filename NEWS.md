# v0.1.8
The function msa() has been updated to call the executables of either MUSCLE (v3.8.31) or Clustal Omega (v1.2.4) directly from R. Further details can be found in the documentation for the msa() function. However, users are free to use alternative methods to perform the required alignments. A new vignette titled "Performing Sequence Alignment in R" has been added, providing guidance on this topic.

# v0.1.7
A new function, orthology(), has been added. This function, which takes as input a file containing the species tree and the gene tree provided by the user, builds an orthology network graph for genes and species beyond those pre-selected in orthGS. Thus, this new version extends the functionality of the package beyond the GS family. A new vignette titled "Orthology beyond GS" has been included to illustrate the use of the new function.

# v0.1.6
Vignettes modified.

The user is given the choice to used the standalone version of MUSCLE or, alternatively, to utilize the R packages 'muscle'. Further details can be found in the documentation of the function msa().
