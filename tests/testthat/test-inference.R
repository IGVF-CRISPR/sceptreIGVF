test_that("Gasperini inference with default args works", {
  data("mudata_inference_gasperini")
  mudata_out <- inference_sceptre(mudata_inference_gasperini)
  pairs_to_test <- MultiAssayExperiment::metadata(mudata_inference_gasperini)$pairs_to_test |>
    as.data.frame()
  # test that metadata(mudata_out)$test_results exists
  expect_true("test_results" %in% names(MultiAssayExperiment::metadata(mudata_out)))
  test_results <- MultiAssayExperiment::metadata(mudata_out)$test_results
  # test that test_results is a data frame
  expect_true(is.data.frame(test_results))
  # test that test_results has the same number of rows as pairs_to_test
  expect_equal(nrow(test_results), nrow(pairs_to_test))
  # test that test_results has correct column names
  expect_equal(names(test_results), c(names(pairs_to_test), "p_value", "log2_fc"))
  # test that median p_value for positive_control pairs is small
  expect_true(test_results |>
                dplyr::filter(pair_type == "positive_control") |>
                dplyr::summarize(median(p_value) < 1e-10) |>
                dplyr::pull())
})

test_that("inference with extra args works", {
  data("mudata_inference_gasperini")
  mudata_out <- inference_sceptre(
    mudata_inference_gasperini,
    side = "left",
    formula_object = formula(~ prep_batch + log(response_n_nonzero) + log(response_n_umis))
  )
  pairs_to_test <- MultiAssayExperiment::metadata(mudata_inference_gasperini)$pairs_to_test |>
    as.data.frame()
  # test that metadata(mudata_out)$test_results exists
  expect_true("test_results" %in% names(MultiAssayExperiment::metadata(mudata_out)))
  test_results <- MultiAssayExperiment::metadata(mudata_out)$test_results
  # test that test_results is a data frame
  expect_true(is.data.frame(test_results))
  # test that test_results has the same number of rows as pairs_to_test
  expect_equal(nrow(test_results), nrow(pairs_to_test))
  # test that test_results has correct column names
  expect_equal(names(test_results), c(names(pairs_to_test), "p_value", "log2_fc"))
  # test that median p_value for positive_control pairs is small
  expect_true(test_results |>
                dplyr::filter(pair_type == "positive_control") |>
                dplyr::summarize(median(p_value) < 1e-10) |>
                dplyr::pull())
})

test_that("Papalexi inference with default args works", {
  data("mudata_inference_papalexi")
  mudata_out <- inference_sceptre(mudata_inference_papalexi)
  pairs_to_test <- MultiAssayExperiment::metadata(mudata_inference_papalexi)$pairs_to_test |>
    as.data.frame()
  # test that metadata(mudata_out)$test_results exists
  expect_true("test_results" %in% names(MultiAssayExperiment::metadata(mudata_out)))
  test_results <- MultiAssayExperiment::metadata(mudata_out)$test_results
  # test that test_results is a data frame
  expect_true(is.data.frame(test_results))
  # test that test_results has the same number of rows as pairs_to_test
  expect_equal(nrow(test_results), nrow(pairs_to_test))
  # test that test_results has correct column names
  expect_equal(names(test_results), c(names(pairs_to_test), "p_value", "log2_fc"))
  # test that median p_value for positive_control pairs is small
  expect_true(test_results |>
                dplyr::filter(pair_type == "positive_control") |>
                dplyr::summarize(median(p_value) < 1e-10) |>
                dplyr::pull())
})

test_that("Papalexi inference with extra args works", {
  data("mudata_inference_papalexi")
  mudata_out <- inference_sceptre(
    mudata_inference_papalexi,
    formula_object = formula(~ replicate + log(response_n_nonzero) + log(response_n_umis))
  )
  pairs_to_test <- MultiAssayExperiment::metadata(mudata_inference_papalexi)$pairs_to_test |>
    as.data.frame()
  # test that metadata(mudata_out)$test_results exists
  expect_true("test_results" %in% names(MultiAssayExperiment::metadata(mudata_out)))
  test_results <- MultiAssayExperiment::metadata(mudata_out)$test_results
  # test that test_results is a data frame
  expect_true(is.data.frame(test_results))
  # test that test_results has the same number of rows as pairs_to_test
  expect_equal(nrow(test_results), nrow(pairs_to_test))
  # test that test_results has correct column names
  expect_equal(names(test_results), c(names(pairs_to_test), "p_value", "log2_fc"))
  # test that median p_value for positive_control pairs is small
  expect_true(test_results |>
                dplyr::filter(pair_type == "positive_control") |>
                dplyr::summarize(median(p_value) < 1e-10) |>
                dplyr::pull())
})
