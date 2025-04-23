library(orthGS)

## ---------------------------------------------- ##
#                 Testing orthG                    #
## ---------------------------------------------- ##
test_that("orthG() works properly with UniProt",{

  skip_on_cran()
  skip_on_travis()

  a <- orthG(set = "Pp")
  b <- orthG(set = c("Pp", "Psy", "Psm", "Ap"))
  c <- orthG(set = "all")

  expect_is(a, 'list')
  expect_is(a[[1]], 'matrix')
  expect_is(a[[2]], 'igraph')
  expect_equal(nrow(a[[1]]), 3)

  expect_is(b, 'list')
  expect_is(b[[1]], 'matrix')
  expect_is(b[[2]], 'igraph')
  expect_equal(nrow(b[[1]]), 11)
  expect_equal(ncol(b[[1]]), 11)

  expect_is(c, 'list')
  expect_is(c[[1]], 'matrix')
  expect_is(c[[2]], 'igraph')
  expect_equal(nrow(c[[1]]), 155)
  expect_equal(ncol(c[[1]]), 155)

})

## ---------------------------------------------- ##
#                 Testing orthP                    #
## ---------------------------------------------- ##
test_that("orthP() works properly with UniProt",{

  skip_on_cran()
  skip_on_travis()

  a <- orthP(phylo_id = "Zm_GS1b_1", set = c("Arabidopsis thaliana", "Oryza sativa"))
  b <- orthP(phylo_id = "Zm_GS1b_1", set = "Ath")
  c <- orthP(phylo_id = "Zm_GS1b_1", set = "all")

  expect_is(a, 'list')
  expect_equal(length(a), 3)
  expect_is(a[[2]], 'character')
  expect_is(a[[3]], 'character')
  expect_true("Zm_GS1b_1" %in% a[[3]])

  expect_is(b, 'character')

  expect_is(c, 'list')
  expect_equal(length(c), 3)
  expect_is(c[[2]], 'character')
  expect_is(c[[3]], 'character')
  expect_true("Zm_GS1b_1" %in% a[[3]])

})

## ---------------------------------------------- ##
#               Testing getseqGS                   #
## ---------------------------------------------- ##
test_that("getseqGS() works properly with UniProt",{

  skip_on_cran()
  skip_on_travis()

  a <- getseqGS(phylo_id = "Pp_GS1b_2", molecule = "Prot")
  b <- getseqGS(phylo_id = "Pp_GS1b_2", molecule = "CDS")

  expect_is(a, 'character')
  expect_equal(nchar(a), 357)

  expect_is(b, 'character')
  expect_equal(nchar(b), 1074)

})

## ---------------------------------------------- ##
#                Testing subsetGS                  #
## ---------------------------------------------- ##
test_that("subsetGS() works properly with UniProt",{

  skip_on_cran()
  skip_on_travis()

  a <- subsetGS(sp = c("Arabidopsis thaliana", "Oryza sativa"))
  b <- subsetGS(sp = "Ath")
  # c <- subsetGS(sp = c("Arabidopsis thaliana", "Oryza sativa", "Homo sapiens"))
  # d <- subsetGS(sp = c("Mus musculus", "Arabidopsis thaliana", "Oryza sativa", "Homo sapiens"))

  expect_is(a, "data.frame")
  expect_equal(dim(a), c(10,23))
  expect_is(b, "data.frame")
  expect_equal(dim(b), c(6,23))
  expect_warning(subsetGS(sp = c("Arabidopsis thaliana", "Oryza sativa", "Homo sapiens")),
                 "The following species has not been found in our database:   Homo sapiens")
  expect_warning(subsetGS(sp = c("Mus musculus", "Arabidopsis thaliana", "Oryza sativa", "Homo sapiens")),
                 "The following species have not been found in our database:   Mus musculus , Homo sapiens")

})


## ---------------------------------------------- ##
#                Testing speciesGS                 #
## ---------------------------------------------- ##
test_that("subsetGS() works properly with UniProt",{

  skip_on_cran()
  skip_on_travis()

  a <- speciesGS(c("Pinus pinaster", "Ath"))
  b <- speciesGS("Pinus pinaster")
  c <- speciesGS("Ath")

  expect_is(a, "data.frame")
  expect_equal(dim(a), c(2,2))
  expect_equal(a[1,2], "Pp")
  expect_equal(a[2,2], "Arabidopsis thaliana")

  expect_is(b, "data.frame")
  expect_equal(dim(b), c(1,2))
  expect_equal(b[1,2], "Pp")

  expect_is(c, "data.frame")
  expect_equal(dim(c), c(1,2))
  expect_equal(c[1,2], "Arabidopsis thaliana")
})

## ---------------------------------------------- ##
#                Testing orthology                 #
## ---------------------------------------------- ##
test_that("orthology() works properly with UniProt",{

  skip_on_cran()
  skip_on_travis()

  path <- system.file("extdata", "input.trees", package = "orthGS")
  a <- orthology(trees = path, plot = FALSE, saverec = FALSE)
  b <- orthology(trees = path, plot = FALSE, saverec = "./vistoynovisto")

  expect_is(a, "list")
  expect_equal(length(a), 4)
  expect_is(a[[1]], "phylo")
  expect_is(a[[2]], "data.frame")
  expect_is(a[[3]], "matrix")
  expect_is(a[[4]], "igraph")

  expect_is(b, "list")
  expect_equal(length(b), 4)
  expect_is(b[[1]], "phylo")
  expect_is(b[[2]], "data.frame")
  expect_is(b[[3]], "matrix")
  expect_is(b[[4]], "igraph")
  expect_true(file.exists("./vistoynovisto"))
  file.remove("./vistoynovisto")
})
