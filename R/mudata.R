#' Convert IGVF CRISPR MuData to sceptre_object
#'
#' @param mudata IGVF CRISPR MuData object
#' @param remove_collinear_covariates Logical flag indicating whether to remove
#' multicollinear covariates. If true, checks whether multicollinearity exists
#' among covariates, and if so, removes all covariates.
#'
#' @return A sceptre_object
#' @export
#'
#' @examples
#' library(sceptreIGVF)
#' # load sample MuData for guide assignment
#' data(mudata_guide_assignment_gasperini)
#' # convert MuData object to sceptre object
#' sceptre_object_guide_assignment <- convert_mudata_to_sceptre_object(
#'   mudata_guide_assignment_gasperini
#' )
#' sceptre_object_guide_assignment
#' # load sample MuData for inference
#' data(mudata_inference_gasperini)
#' # convert MuData object to sceptre object
#' sceptre_object_inference <- convert_mudata_to_sceptre_object(
#'   mudata_inference_gasperini
#' )
#' sceptre_object_inference
convert_mudata_to_sceptre_object <- function(mudata, remove_collinear_covariates = FALSE){
  # extract information from MuData
  moi <- MultiAssayExperiment::metadata(mudata[['guide']])$moi
  if(is.null(SummarizedExperiment::assayNames(mudata[['gene']]))){
    SummarizedExperiment::assayNames(mudata[['gene']]) <- 'counts'
  } else{
    SummarizedExperiment::assayNames(mudata[['gene']])[[1]] <- 'counts'
  }
  if(is.null(SummarizedExperiment::assayNames(mudata[['guide']]))){
    SummarizedExperiment::assayNames(mudata[['guide']]) <- 'counts'
  } else{
    SummarizedExperiment::assayNames(mudata[['guide']])[[1]] <- 'counts'
  }

  scRNA_data <- mudata@ExperimentList$gene
  guides_data <- mudata@ExperimentList$guide
  response_matrix <- scRNA_data@assays@data@listData[["counts"]]

  if(!is.null(SummarizedExperiment::colData(mudata))){
    covariates <- SummarizedExperiment::colData(mudata) |>
      as.data.frame()
    if(remove_collinear_covariates){
      model_matrix <- stats::model.matrix(object = ~ ., data = covariates)
      multicollinear <- Matrix::rankMatrix(model_matrix) < ncol(model_matrix)
      if(multicollinear){
        print("Removing multicollinear covariates")
        extra_covariates <- data.frame()
      } else{
        extra_covariates <- covariates
      }
    } else{
      extra_covariates <- covariates
    }
  } else{
    extra_covariates <- data.frame()
  }

  # if guide assignments not present, then extract guide counts
  if(length(guides_data@assays@data@listData) == 1){
    grna_matrix <- guides_data@assays@data@listData[["counts"]]
    # otherwise, extract guide assignments
  } else{
    grna_matrix <- guides_data@assays@data@listData[["guide_assignment"]]
  }

  grna_ids <- rownames(SingleCellExperiment::rowData(mudata[['guide']]))
  rownames(grna_matrix) <- grna_ids

  gene_ids <- rownames(SingleCellExperiment::rowData(mudata[['gene']]))
  rownames(response_matrix) <- gene_ids
  grna_target_data_frame <- SingleCellExperiment::rowData(mudata[['guide']]) |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "grna_id") |>
    dplyr::rename(grna_target = intended_target_name) |>
    dplyr::select(grna_id, grna_target)

  # assemble information into sceptre object
  sceptre_object <- sceptre::import_data(
    response_matrix = response_matrix,
    grna_matrix = grna_matrix,
    grna_target_data_frame = grna_target_data_frame,
    moi = moi,
    extra_covariates = extra_covariates
  )

  # return sceptre object
  return(sceptre_object)
}

#' Convert sceptre_object to MuData
#'
#' @param sceptre_object A sceptre_object
#'
#' @return A MuData object
#' @export
sceptre_object_to_mudata <- function(sceptre_object){
  # 1. Extract objects and subset to the set of cells in use
  cells_in_use <- sceptre_object@cells_in_use
  if(length(cells_in_use) == 0){
    cells_in_use <- 1:ncol(sceptre_object@response_matrix)
  }
  response_matrix <- sceptre_object@response_matrix[, cells_in_use]
  grna_matrix <- sceptre_object@grna_matrix[, cells_in_use]
  grna_assignment_matrix <- sceptre_object |>
    sceptre::get_grna_assignments() |>
    methods::as("dsparseMatrix")
  if(!is.null(grna_assignment_matrix)){
   grna_assignment_matrix <- grna_assignment_matrix[, cells_in_use]
  }
  covariate_df <- sceptre_object@covariate_data_frame[cells_in_use, ]
  grna_target_data_frame <- sceptre_object@grna_target_data_frame
  positive_control_pairs <- sceptre_object@positive_control_pairs
  discovery_pairs <- sceptre_object@discovery_pairs
  moi <- if(sceptre_object@low_moi) "low" else "high"
  gene_names <- sceptre_object@response_names

  # 2. Extract batch info, if present
  sample_df <- covariate_df |>
    dplyr::select(-dplyr::any_of(c("grna_n_nonzero", "grna_n_umis",
                  "response_n_nonzero", "response_n_umis", "response_p_mito")))

  # 3. Extra gRNA and gene information
  grna_ids <- rownames(grna_matrix)
  grna_rowdata <- grna_target_data_frame |>
    dplyr::arrange(match(grna_id, grna_ids)) |>
    tibble::column_to_rownames(var = "grna_id") |>
    dplyr::mutate(targeting = ifelse(grna_target != "non-targeting", "TRUE", "FALSE")) |>
    dplyr::relocate(targeting) |>
    dplyr::rename(intended_target_name = grna_target)
  if("chr" %in% colnames(grna_rowdata)){
    grna_rowdata <- grna_rowdata |>
      dplyr::mutate(chr = ifelse(is.na(chr), "", chr)) |>
      dplyr::rename(intended_target_chr = chr)
  }
  if("start" %in% colnames(grna_rowdata)){
    grna_rowdata <- grna_rowdata |>
      dplyr::mutate(start = ifelse(is.na(start), -9, start)) |>
      dplyr::rename(intended_target_start = start)
  }
  if("end" %in% colnames(grna_rowdata)){
    grna_rowdata <- grna_rowdata |>
      dplyr::mutate(end = ifelse(is.na(end), -9, end)) |>
      dplyr::rename(intended_target_end = end)
  }
  if("grna_group" %in% colnames(grna_rowdata)){
    grna_rowdata <- grna_rowdata |>
      dplyr::select(-grna_group)
  }
  grna_coldata <- covariate_df |>
    dplyr::select(grna_n_nonzero, grna_n_umis) |>
    dplyr::rename(num_expressed_guides = grna_n_nonzero, total_guide_umis = grna_n_umis)

  gene_coldata <- covariate_df |>
    dplyr::select(response_n_nonzero, response_n_umis) |>
    dplyr::rename(num_expressed_genes = response_n_nonzero, total_gene_umis = response_n_umis)

  # 4. Prepare for conversion to MuData
  response_matrix <- methods::as(response_matrix, "CsparseMatrix")
  grna_matrix <- methods::as(grna_matrix, "CsparseMatrix")
  if(!is.null(colnames(response_matrix))){
    cell_ids <- colnames(response_matrix)
  } else if(!is.null(colnames(grna_matrix))){
    cell_ids <- colnames(grna_matrix)
  } else if(!is.null(rownames(covariate_df))){
    cell_ids <- rownames(covariate_df)
  } else{
    cell_ids <- paste0("cell_", 1:ncol(response_matrix))
  }

  colnames(response_matrix) <- cell_ids
  colnames(grna_matrix) <- cell_ids
  rownames(sample_df) <- cell_ids
  if(!is.null(grna_assignment_matrix)){
    colnames(grna_assignment_matrix) <- cell_ids
  }
  pairs <- rbind(
    positive_control_pairs |> dplyr::mutate(pair_type = "positive_control"),
    discovery_pairs |> dplyr::mutate(pair_type = "discovery")
  )
  metadata <- list()

  metadata[["pairs_to_test"]] <- pairs |>
    dplyr::select(grna_target, response_id, pair_type) |>
    dplyr::rename(intended_target_name = grna_target, gene_id = response_id) |>
    MultiAssayExperiment::DataFrame()

  metadata[["test_results"]] <- pairs |>
    dplyr::select(grna_target, response_id, pair_type) |>
    dplyr::left_join(rbind(sceptre_object@power_result |>
                      dplyr::select(grna_target, response_id, p_value, log_2_fold_change),
                    sceptre_object@discovery_result |>
                      dplyr::select(grna_target, response_id, p_value, log_2_fold_change)),
              by = c("grna_target", "response_id")) |>
    dplyr::rename(intended_target_name = grna_target,
                  gene_id = response_id,
                  log2_fc = log_2_fold_change) |>
    MultiAssayExperiment::DataFrame()

  # 5. Create MuData object
  gene_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = response_matrix),
    colData = gene_coldata,
  )
  grna_assays <- list(counts = grna_matrix)
  if (!is.null(grna_assignment_matrix)) {
    grna_assays[["guide_assignment"]] <- grna_assignment_matrix
  }
  grna_se <- SummarizedExperiment::SummarizedExperiment(
    assays = grna_assays,
    rowData = grna_rowdata,
    colData = grna_coldata,
    metadata = list(moi = moi)
  )
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(gene = gene_se, guide = grna_se),
    colData = sample_df,
    metadata = metadata
  )

  return(mae)
}

#' Convert sceptre_object a set of MuData objects to be used as inputs and outputs
#' for both gRNA assignment and inference
#'
#' @param sceptre_object Sceptre object
#' @param num_discovery_pairs Number of discovery pairs to sample
#' @param gene_info Gene information data frame (optional)
#' @param guide_capture_method Guide capture method (optional)
#'
#' @return List of MuData objects
#' @export
sceptre_object_to_mudata_inputs_outputs <- function(
    sceptre_object,
    num_discovery_pairs,
    gene_info = NULL,
    guide_capture_method = NULL
  ){
  positive_control_pairs_2 <- sceptre_object |>
    sceptre::get_result("run_power_check") |>
    stats::na.omit() |>
    dplyr::select(response_id, grna_target)

  discovery_results <- sceptre_object |>
    sceptre::get_result("run_discovery_analysis") |>
    stats::na.omit()
  num_significant <- discovery_results |>
    dplyr::summarize(sum(significant)) |>
    dplyr::pull()
  num_non_significant <- discovery_results |>
    dplyr::summarize(sum(!significant)) |>
    dplyr::pull()
  num_significant_to_keep <- min(num_significant, round(num_discovery_pairs/2))
  num_non_significant_to_keep <- min(num_non_significant,
                                     num_discovery_pairs - num_significant_to_keep)
  discovery_pairs_2 <- rbind(
    discovery_results |>
      dplyr::filter(significant) |>
      dplyr::slice_sample(n = num_significant_to_keep),
    discovery_results |>
      dplyr::filter(!significant) |>
      dplyr::slice_sample(n = num_non_significant_to_keep)
  ) |>
    dplyr::select(response_id, grna_target)

  sceptre_object_2 <- sceptre_object |>
    sceptre::set_analysis_parameters(
      discovery_pairs = discovery_pairs_2,
      positive_control_pairs = positive_control_pairs_2,
      formula = sceptre_object@formula_object,
      side = c("left", "both", "right")[sceptre_object@side_code + 2],
    ) |>
    sceptre::assign_grnas(parallel = TRUE) |>
    sceptre::run_qc() |>
    sceptre::run_power_check(parallel = TRUE) |>
    sceptre::run_discovery_analysis(parallel = TRUE)

  mae_inference_output <- sceptre_object_to_mudata(sceptre_object_2)
  if(!is.null(gene_info)){
    SummarizedExperiment::rowData(mae_inference_output[["gene"]]) <- gene_info
  }
  if(!is.null(guide_capture_method)){
    MultiAssayExperiment::metadata(mae_inference_output[["guide"]])$capture_method <- guide_capture_method
  }

  mae_inference_input <- mae_inference_output
  MultiAssayExperiment::metadata(mae_inference_input)$test_results <- NULL

  mae_guide_assignment_output <- mae_inference_input
  MultiAssayExperiment::metadata(mae_guide_assignment_output)$pairs_to_test <- NULL

  mae_guide_assignment_input <- mae_guide_assignment_output
  guide_assays_list <- SummarizedExperiment::assays(mae_guide_assignment_input[['guide']])
  guide_assays_list$guide_assignment <- NULL
  SummarizedExperiment::assays(mae_guide_assignment_input[['guide']]) <- guide_assays_list

  mae_inference_output_minimal <- mae_inference_output
  MultiAssayExperiment::colData(mae_inference_output_minimal) <- MultiAssayExperiment::DataFrame(
    row.names = rownames(MultiAssayExperiment::colData(mae_inference_output_minimal))
  )
  SummarizedExperiment::rowData(mae_inference_output_minimal[['gene']]) <- NULL
  SummarizedExperiment::colData(mae_inference_output_minimal[['gene']]) <- NULL
  SummarizedExperiment::rowData(mae_inference_output_minimal[['guide']]) <-
    SummarizedExperiment::rowData(mae_inference_output_minimal[['guide']])[,c("targeting", "intended_target_name")]
  SummarizedExperiment::colData(mae_inference_output_minimal[['guide']]) <- NULL
  MultiAssayExperiment::metadata(mae_inference_output_minimal)$pairs_to_test <-
    MultiAssayExperiment::metadata(mae_inference_output_minimal)$pairs_to_test[,c("intended_target_name", "gene_id")]
  MultiAssayExperiment::metadata(mae_inference_output_minimal)$test_results <-
    MultiAssayExperiment::metadata(mae_inference_output_minimal)$test_results[,c("intended_target_name", "gene_id", "p_value")]

  mae_inference_input_minimal <- mae_inference_output_minimal
  MultiAssayExperiment::metadata(mae_inference_input_minimal)$test_results <- NULL

  mae_guide_assignment_output_minimal <- mae_inference_input_minimal
  MultiAssayExperiment::metadata(mae_guide_assignment_output_minimal)$pairs_to_test <- NULL

  mae_guide_assignment_input_minimal <- mae_guide_assignment_output_minimal
  guide_assays_list <- SummarizedExperiment::assays(mae_guide_assignment_input_minimal[['guide']])
  guide_assays_list$guide_assignment <- NULL
  SummarizedExperiment::assays(mae_guide_assignment_input_minimal[['guide']]) <- guide_assays_list

  return(list(
    inference_output = mae_inference_output,
    inference_input = mae_inference_input,
    guide_assignment_output = mae_guide_assignment_output,
    guide_assignment_input = mae_guide_assignment_input,
    inference_output_minimal = mae_inference_output_minimal,
    inference_input_minimal = mae_inference_input_minimal,
    guide_assignment_output_minimal = mae_guide_assignment_output_minimal,
    guide_assignment_input_minimal = mae_guide_assignment_input_minimal)
  )
}


#' Save a list of MuData objects to disk
#'
#' @param mudata_list List of MuData objects
#' @param path Path to save files to
#' @param prefix Prefix to add to file names
#'
#' @return NULL
#' @export
save_mudata_list <- function(mudata_list, path = ".", prefix = ""){
  # create subdirectories of path for guide assignment and inference
  guide_assignment_path <- file.path(path, "guide_assignment")
  inference_path <- file.path(path, "inference")
  paths <- c(guide_assignment_path, inference_path)
  for(p in paths){
    if(!dir.exists(p)){
      dir.create(p, recursive = TRUE)
    }
  }
  # save guide assignment and inference MuData objects to disk
  for(mudata_name in names(mudata_list)){
    # determine subfolder to save file in
    if(grepl("guide_assignment", mudata_name)){
      path <- guide_assignment_path
    } else if(grepl("inference", mudata_name)){
      path <- inference_path
    } else {
      stop("mudata_name must contain either 'guide_assignment' or 'inference'")
    }
    MuData::writeH5MU(object = mudata_list[[mudata_name]],
                      file = file.path(path, paste0(prefix, mudata_name, ".h5mu")))
  }
}
