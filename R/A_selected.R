#' Adjacency Matrix for Orthology Graph
#'
#' 155 x 155 square matrix (155 GS proteins from 45 seed plant species)
#'
#' @format
#' A matrix with 155 rows and 155 columns
#'
#' @source
#' It has been generated using the function orthG::mapTrees()
#' and the reconciliation output file 'selected'.
#' Verbigracia: orthG::mapTrees('./inst/extdata/selected')
#' The reconciliation was carried out using RANGER-DTL with parameters D = 1, T = 10 and L = 1.
"A_selected"
