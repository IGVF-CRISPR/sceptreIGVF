test_that("threshold-based gRNA assignment works", {
  # set gRNA threshold
  threshold <- 5
  # load sample MuData
  data(sample_mudata)
  # assign gRNAs via assign_grnas_sceptre
  mudata <- assign_grnas_sceptre(sample_mudata, method = "thresholding", threshold = threshold)
  package_grna_assignment_matrix <- mudata[["grna_assignment"]] |>
    SingleCellExperiment::counts() |>
    as.matrix() > 0
  # assign gRNAs by directly thresholding gRNA matrix
  grna_matrix <- sample_mudata[["guides"]]@assays@data@listData[[1]] |> as.matrix()
  correct_grna_assignment_matrix <- grna_matrix >= threshold
  # check that the gRNA assignment matrix is correct
  expect_identical(
    package_grna_assignment_matrix,
    correct_grna_assignment_matrix
  )
})
