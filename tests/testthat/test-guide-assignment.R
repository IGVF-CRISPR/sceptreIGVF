test_that("Gasperini guide assignment works", {
  data("mudata_guide_assignment_gasperini")
  mudata_out <- assign_grnas_sceptre(mudata_guide_assignment_gasperini)
  # test that mudata_out still has two experiments with names "gene" and "guide"
  expect_equal(names(mudata_out), c("gene", "guide"))
  # test that mudata_out[['guide']] has two assays, one of which is named "guide_assignment"
  assay_names <- SummarizedExperiment::assayNames(mudata_out[['guide']])
  expect_equal(length(assay_names), 2)
  expect_true("guide_assignment" %in% assay_names)
  # test that the "guide_assignment" has the same dimensions as the guide counts matrix
  guide_counts_assay_index <- which(assay_names != "guide_assignment")
  guide_counts_matrix <- SummarizedExperiment::assay(mudata_out[['guide']], guide_counts_assay_index)
  guide_assignment_matrix <- SummarizedExperiment::assay(mudata_out[['guide']], "guide_assignment")
  expect_equal(dim(guide_assignment_matrix), dim(guide_counts_matrix))
  # test that all entries of the guide assignment matrix are zero or one
  expect_true(all(guide_assignment_matrix == 0 | guide_assignment_matrix == 1))
})

test_that("Papalexi guide assignment works", {
  data("mudata_guide_assignment_papalexi")
  mudata_out <- assign_grnas_sceptre(mudata_guide_assignment_papalexi)
  # test that mudata_out still has two experiments with names "gene" and "guide"
  expect_equal(names(mudata_out), c("gene", "guide"))
  # test that mudata_out[['guide']] has two assays, one of which is named "guide_assignment"
  assay_names <- SummarizedExperiment::assayNames(mudata_out[['guide']])
  expect_equal(length(assay_names), 2)
  expect_true("guide_assignment" %in% assay_names)
  # test that the "guide_assignment" has the same dimensions as the guide counts matrix
  guide_counts_assay_index <- which(assay_names != "guide_assignment")
  guide_counts_matrix <- SummarizedExperiment::assay(mudata_out[['guide']], guide_counts_assay_index)
  guide_assignment_matrix <- SummarizedExperiment::assay(mudata_out[['guide']], "guide_assignment")
  expect_equal(dim(guide_assignment_matrix), dim(guide_counts_matrix))
  # test that all entries of the guide assignment matrix are zero or one
  expect_true(all(guide_assignment_matrix == 0 | guide_assignment_matrix == 1))
})
