library(PPCKO)
context("1d domain case")



test_that(" in the 1d domain case KO without CV and hps check work", {
  
  data("data_1d", package = "PPCKO")
  
  expect_equal(length(
    PPCKO::KO_check_hps( X = data_1d )), 1)
  
  expect_equal(length(
    PPCKO::PPC_KO( X = data_1d )), 17)
})



test_that(" in the 1d domain case KO with CV for regularization parameter works", {
  
  data("data_1d", package = "PPCKO")
  alpha_vec <- c(1e-3,1e-2,1e-1,1,1e1,1e2)
  
  expect_equal(length(
    PPCKO::PPC_KO( X = data_1d,
                   id_CV = "CV_alpha",
                   alpha_vec = alpha_vec,
                   min_size_ts = 90,
                   max_size_ts = 92,
                   err_ret = 0)), 17)
  
  expect_equal(length(
    PPCKO::PPC_KO( X = data_1d,
                   id_CV = "CV_alpha",
                   alpha_vec = alpha_vec,
                   min_size_ts = 90,
                   max_size_ts = 92,
                   err_ret = 1)), 18)
})



test_that(" in the 1d domain case KO with CV for number of PPCs works", {
  
  data("data_1d", package = "PPCKO")
  k_vec     <- c(1,2,3,4)
  
  expect_equal(length(
    PPCKO::PPC_KO( X = data_1d,
                   id_CV = "CV_k",
                   k_vec = k_vec,
                   min_size_ts = 90,
                   max_size_ts = 92,
                   err_ret = 0)), 17)
  
  expect_equal(length(
    PPCKO::PPC_KO( X = data_1d,
                   id_CV = "CV_k",
                   k_vec = k_vec,
                   min_size_ts = 90,
                   max_size_ts = 92,
                   err_ret = 1)), 18)
})



test_that(" in the 1d domain case KO with CV for boht regularization parameter and number of PPCs works", {
  
  data("data_1d", package = "PPCKO")
  alpha_vec <- c(1e-3,1e-2,1e-1,1,1e1,1e2)
  k_vec     <- c(1,2,3,4)
  
  expect_equal(length(
    PPCKO::PPC_KO( X = data_1d,
                   id_CV = "CV",
                   alpha_vec = alpha_vec,
                   k_vec = k_vec,
                   min_size_ts = 90,
                   max_size_ts = 92,
                   err_ret = 0)), 17)
  
  expect_equal(length(
    PPCKO::PPC_KO( X = data_1d,
                   id_CV = "CV",
                   alpha_vec = alpha_vec,
                   k_vec = k_vec,
                   min_size_ts = 90,
                   max_size_ts = 92,
                   err_ret = 1)), 18)
})