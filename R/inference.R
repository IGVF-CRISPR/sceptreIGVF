#' SCEPTRE inference
#'
#' @param mudata A MuData object
#'
#' @return A MuData object with sceptre inference results
#' @export
#'
#' @examples
#' library(sceptreIGVF)
#' # load sample MuData
#' data(mudata_inference)
#' # inference
#' mudata_out <- inference_sceptre(mudata_inference)
#' mudata_out
inference_sceptre <- function(mudata) {
  # convert MuData object to sceptre object
  sceptre_object <- convert_mudata_to_sceptre_object(mudata)

  pairs_to_test <- MultiAssayExperiment::metadata(mudata)$pairs_to_test |>
    as.data.frame()

  discovery_pairs <- pairs_to_test |>
    dplyr::rename(
      grna_target = intended_target_name,
      response_id = gene_id
    )

  # set analysis parameters
  sceptre_object <- sceptre_object |>
    sceptre::set_analysis_parameters(
      discovery_pairs = discovery_pairs,
      formula_object = formula(~ log(response_n_nonzero) + log(response_n_umis))
    )

  # assign gRNAs and run default QC
  sceptre_object <- sceptre_object |>
    sceptre::assign_grnas(method = "thresholding", threshold = 1) |>
    sceptre::run_qc()

  # run discovery analysis
  sceptre_object <- sceptre_object |>
    sceptre::run_discovery_analysis()

  # get results
  discovery_results <- sceptre_object |>
    sceptre::get_result(analysis = "run_discovery_analysis") |>
    dplyr::select(response_id, grna_target, p_value, log_2_fold_change) |>
    dplyr::rename(gene_id = response_id,
                  intended_target_name = grna_target,
                  log2_fc = log_2_fold_change)

  test_results <- pairs_to_test |>
    dplyr::left_join(discovery_results, by = c("intended_target_name", "gene_id"))

  # add results to MuData
  metadata(mudata)$test_results <- test_results

  # return MuData
  return(mudata)
}


#' Extra gRNA assignment matrix from `sceptre` object
#'
#' @param sceptre_object A sceptre_object
#'
#' @return A sparse logical matrix of gRNA assignments
#' @export
extract_grna_assignment_matrix <- function(sceptre_object){
  grna_assignments <- sceptre_object@initial_grna_assignment_list
  if(length(grna_assignments) == 0){
    return(NULL)
  }
  grna_names <- rownames(sceptre_object@grna_matrix)
  cell_barcodes <- colnames(sceptre_object@grna_matrix)
  grna_assignments <- grna_assignments[grna_names]
  j <- unlist(grna_assignments)
  p <- c(0, cumsum(sapply(grna_assignments, length)))
  num_rows <- length(grna_assignments)
  num_cols <- sceptre_object@grna_matrix |> ncol()
  grna_assignment_matrix <- Matrix::sparseMatrix(
    j = j,
    p = p,
    dims = c(num_rows, num_cols),
    dimnames = list(grna_names, cell_barcodes),
    x = rep(1, length(j))
  )
  return(grna_assignment_matrix)
}
