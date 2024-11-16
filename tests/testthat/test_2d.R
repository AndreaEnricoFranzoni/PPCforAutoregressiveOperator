library(PPCKO)
context("2d domain case")

test_that(" in the 2d domain case KO without CV and hps check work", {
  
  data("data_2d", package = "PPCKO")
  alpha_vec <- c(1e-3,1e-2,1e-1,1,1e1,1e2)
  k_vec     <- c(1,2,3,4)
  
  left_ex_x1  <- 0
  right_ex_x1 <- 1
  left_ex_x2  <- 0
  right_ex_x2 <- 1
  dim_grid_x1 <- 10
  dim_grid_x2 <- 10
  x1.grid = seq(from=left_ex_x1, to=right_ex_x1, length=dim_grid_x1)
  x2.grid = seq(from=left_ex_x2, to=right_ex_x2, length=dim_grid_x2)
  
  
  x_t = PPCKO.local2::data_2d_wrapper_from_list(data_2d)
  expect_equal(nrow(x_t), 100)
  expect_equal(ncol(x_t), 20)
  
  
  expect_equal(length(
    PPCKO.local2::KO_check_hps_2d( X = x_t,
                                   dim_x1 = dim_grid_x1,
                                   dim_x2 = dim_grid_x2)), 1)
  
  expect_equal(length(
    PPCKO.local2::PPC_KO_2d( X = x_t )), 17)
})



test_that(" in the 2d domain case KO with CV for regularization parameter works", {
  
  data("data_2d", package = "PPCKO")
  left_ex_x1  <- 0
  right_ex_x1 <- 1
  left_ex_x2  <- 0
  right_ex_x2 <- 1
  dim_grid_x1 <- 10
  dim_grid_x2 <- 10
  x1.grid = seq(from=left_ex_x1, to=right_ex_x1, length=dim_grid_x1)
  x2.grid = seq(from=left_ex_x2, to=right_ex_x2, length=dim_grid_x2)
  
  x_t = PPCKO.local2::data_2d_wrapper_from_list(data_2d)
  expect_equal(nrow(x_t), 100)
  expect_equal(ncol(x_t), 20)
  
  alpha_vec <- c(1e-3,1e-2,1e-1,1,1e1,1e2)
  
  
  expect_equal(length(
    PPCKO.local2::PPC_KO_2d( X = x_t,
                             id_CV = "CV_alpha",
                             alpha_vec = alpha_vec,
                             min_size_ts = 10,
                             max_size_ts = 12,
                             err_ret = 0)), 17)
  
  expect_equal(length(
    PPCKO.local2::PPC_KO_2d( X = x_t,
                             id_CV = "CV_alpha",
                             alpha_vec = alpha_vec,
                             min_size_ts = 10,
                             max_size_ts = 12,
                             err_ret = 1)), 18)
})



test_that(" in the 2d domain case KO with CV for the number of PPCs works", {
  
  data("data_2d", package = "PPCKO")
  left_ex_x1  <- 0
  right_ex_x1 <- 1
  left_ex_x2  <- 0
  right_ex_x2 <- 1
  dim_grid_x1 <- 10
  dim_grid_x2 <- 10
  x1.grid = seq(from=left_ex_x1, to=right_ex_x1, length=dim_grid_x1)
  x2.grid = seq(from=left_ex_x2, to=right_ex_x2, length=dim_grid_x2)
  
  x_t = PPCKO.local2::data_2d_wrapper_from_list(data_2d)
  expect_equal(nrow(x_t), 100)
  expect_equal(ncol(x_t), 20)
  
  k_vec <- c(1,2,3,4)
  
  
  expect_equal(length(
    PPCKO.local2::PPC_KO_2d( X = x_t,
                             id_CV = "CV_k",
                             k_vec = k_vec,
                             min_size_ts = 10,
                             max_size_ts = 12,
                             err_ret = 0)), 17)
  
  expect_equal(length(
    PPCKO.local2::PPC_KO_2d( X = x_t,
                             id_CV = "CV_k",
                             k_vec = k_vec,
                             min_size_ts = 10,
                             max_size_ts = 12,
                             err_ret = 1)), 18)
})



test_that(" in the 2d domain case KO with CV for both reg param and number of PPCs works", {
  
  data("data_2d", package = "PPCKO")
  left_ex_x1  <- 0
  right_ex_x1 <- 1
  left_ex_x2  <- 0
  right_ex_x2 <- 1
  dim_grid_x1 <- 10
  dim_grid_x2 <- 10
  x1.grid = seq(from=left_ex_x1, to=right_ex_x1, length=dim_grid_x1)
  x2.grid = seq(from=left_ex_x2, to=right_ex_x2, length=dim_grid_x2)
  
  x_t = PPCKO.local2::data_2d_wrapper_from_list(data_2d)
  expect_equal(nrow(x_t), 100)
  expect_equal(ncol(x_t), 20)
  
  alpha_vec <- c(1e-3,1e-2,1e-1,1,1e1,1e2)
  k_vec <- c(1,2,3,4)
  
  
  expect_equal(length(
    PPCKO.local2::PPC_KO_2d( X = x_t,
                             id_CV = "CV",
                             alpha_vec = alpha_vec,
                             k_vec = k_vec,
                             min_size_ts = 10,
                             max_size_ts = 12,
                             err_ret = 0)), 17)
  
  expect_equal(length(
    PPCKO.local2::PPC_KO_2d( X = x_t,
                             id_CV = "CV_k",
                             alpha_vec = alpha_vec,
                             k_vec = k_vec,
                             min_size_ts = 10,
                             max_size_ts = 12,
                             err_ret = 1)), 18)
})