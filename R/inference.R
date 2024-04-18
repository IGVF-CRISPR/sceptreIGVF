#' SCEPTRE inference
#'
#' @param mudata A MuData object
#' @param ... Additional arguments passed to sceptre::set_analysis_parameters();
#'  see `?sceptre::set_analysis_parameters` for details.
#'
#' @return A MuData object with sceptre inference results
#' @export
#'
#' @examples
#' library(sceptreIGVF)
#' # load sample MuData
#' data(mudata_inference_gasperini)
#' # inference
#' mudata_out <- inference_sceptre(mudata_inference_gasperini)
#' mudata_out
inference_sceptre <- function(mudata, ...) {
  # convert MuData object to sceptre object
  sceptre_object <- convert_mudata_to_sceptre_object(mudata)

  # extract set of discovery pairs to test
  pairs_to_test <- MultiAssayExperiment::metadata(mudata)$pairs_to_test |>
    as.data.frame()
  discovery_pairs <- pairs_to_test |>
    dplyr::rename(
      grna_target = intended_target_name,
      response_id = gene_id
    )

  # assemble arguments to set_analysis_parameters()
  args_list <- list(...)
  if("discovery_pairs" %in% names(args_list)){
    warning("The `discovery_pairs` argument is ignored. The `discovery_pairs` are set from the `pairs_to_test` metadata.")
  }
  args_list[["discovery_pairs"]] <- discovery_pairs
  if (!"formula_object" %in% names(args_list) | formula_object == "default") {
    args_list$formula_object <- stats::formula(~ log(response_n_nonzero) + log(response_n_umis))
  }
  args_list$sceptre_object <- sceptre_object

  # set analysis parameters
  sceptre_object <- do.call(sceptre::set_analysis_parameters, args_list)

  # extract gRNA assignment and turn off QC
  sceptre_object <- sceptre_object |>
    sceptre::assign_grnas(method = "thresholding", threshold = 1) |>
    sceptre::run_qc(n_nonzero_trt_thresh = 0L,
                    n_nonzero_cntrl_thresh = 0L,
                    p_mito_threshold = 1)

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

  # add results to MuData
  test_results <- pairs_to_test |>
    dplyr::left_join(discovery_results, by = c("intended_target_name", "gene_id"))
  MultiAssayExperiment::metadata(mudata)$test_results <- test_results

  # return MuData
  return(mudata)
}
