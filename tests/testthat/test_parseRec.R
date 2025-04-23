library(orthGS)

## ---------------------------------------------- ##
#                 Testing coltips                  #
## ---------------------------------------------- ##
test_that("coltips() works properly",{

  skip_on_cran()
  skip_on_travis()
  tr <- ape::read.tree(text = "(((Pp_GS1a, Atr_GS1a), Pp_GS1b_1), whatever);")
  a <- coltips(tr)

    expect_is(a, 'character')
    expect_equal(length(a), 4)
    expect_true(a[4] == 'purple')
})


## ---------------------------------------------- ##
#                 Testing mapTrees                 #
## ---------------------------------------------- ##
test_that("mapTrees() works properly with UniProt",{

  skip_on_cran()
  skip_on_travis()

  a <- mapTrees(path2rec = "./data_t/representatives1")

  expect_is(a, 'list')
  expect_equal(length(a), 3)
  expect_is(a[[1]], 'phylo')
  expect_is(a[[2]], 'data.frame')
  expect_equal(dim(a[[2]]), c(19,4))
  expect_is(a[[3]], 'matrix')
  expect_equal(dim(a[[3]]), c(20,20))

})


## ---------------------------------------------- ##

