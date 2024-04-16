test_that("Gasperini guide assignment followed by inference works", {
  data("mudata_guide_assignment_gasperini")
  data("mudata_inference_gasperini")
  pairs_to_test <- MultiAssayExperiment::metadata(mudata_inference_gasperini)$pairs_to_test |>
    as.data.frame()
  mudata_inference_gasperini <- assign_grnas_sceptre(mudata_guide_assignment_gasperini)
  MultiAssayExperiment::metadata(mudata_inference_gasperini)$pairs_to_test <- pairs_to_test
  mudata_out <- inference_sceptre(mudata_inference_gasperini)
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

test_that("Papalexi guide assignment followed by inference works", {
  data("mudata_guide_assignment_papalexi")
  data("mudata_inference_papalexi")
  pairs_to_test <- MultiAssayExperiment::metadata(mudata_inference_papalexi)$pairs_to_test |>
    as.data.frame()
  mudata_inference_papalexi <- assign_grnas_sceptre(mudata_guide_assignment_papalexi)
  MultiAssayExperiment::metadata(mudata_inference_papalexi)$pairs_to_test <- pairs_to_test
  mudata_out <- inference_sceptre(mudata_inference_papalexi)
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
