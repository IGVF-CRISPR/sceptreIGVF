#' Convert IGVF CRISPR MuData to sceptre_object
#'
#' @param mudata IGVF CRISPR MuData object
#' @param moi A string indicating the multiplicity of infection ("high" or "low")
#' @param include_covariates A logical indicating whether to include the
#' covariates from colData in the sceptre object as extra covariates
#'
#' @return A sceptre_object
#' @export
#'
#' @examples
#' library(sceptreIGVF)
#' # load sample MuData
#' data(sample_mudata)
#' # convert MuData object to sceptre object
#' sceptre_object <- convert_mudata_to_sceptre_object(sample_mudata)
#' sceptre_object
convert_mudata_to_sceptre_object <- function(mudata, moi = "high", include_covariates = TRUE){
  # extract information from MuData
  scRNA_data <- mudata@ExperimentList$scRNA
  guides_data <- mudata@ExperimentList$guides
  response_matrix <- scRNA_data@assays@data@listData[[1]]
  grna_matrix <- guides_data@assays@data@listData[[1]]
  if(include_covariates){
    extra_covariates <- scRNA_data@colData |>
      as.data.frame()
  } else{
    extra_covariates <- data.frame()
  }
  response_names <- scRNA_data@rowRanges@elementMetadata$feature_name |>
    as.character()
  grna_target_data_frame <- guides_data@rowRanges@elementMetadata |>
    as.data.frame() |>
    dplyr::rename(grna_id = feature_name, grna_target = target_elements,
                  chr = guide_chr, start = guide_start, end = guide_end) |>
    dplyr::select(grna_id, grna_target, chr, start, end)

  # assemble information into sceptre object
  sceptre_object <- sceptre::import_data(
    response_matrix = response_matrix,
    grna_matrix = grna_matrix,
    grna_target_data_frame = grna_target_data_frame,
    moi = moi,
    extra_covariates = extra_covariates,
    response_names = response_names
  )

  # return sceptre object
  return(sceptre_object)
}

#' Add a matrix as a new experiment to existing MuData object
#'
#' @param new_matrix A matrix to be added to the MuData object
#' @param mudata An existing MuData object
#' @param experiment_name A string indicating the name of the new experiment
#'
#' @return A MuData object with the new experiment added
#' @export
#'
#' @examples
#' library(sceptreIGVF)
#' library(MultiAssayExperiment)
#' # load sample MuData
#' data(sample_mudata)
#' # create new matrix
#' cell_barcodes <- sample_mudata[["scRNA"]] |> colnames()
#' num_cells <- length(cell_barcodes)
#' num_features <- 100
#' new_matrix <- matrix(
#'   rnorm(num_features * num_cells),
#'   nrow = num_features, ncol = num_cells,
#'   dimnames = list(NULL, cell_barcodes)
#' )
#' # add new matrix to MuData object
#' mudata <- add_matrix_to_mudata(new_matrix, sample_mudata, "new_experiment")
#' mudata
add_matrix_to_mudata <- function(new_matrix, mudata, experiment_name){
  # convert matrix to SingleCellExperiment
  grna_assignment_sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = new_matrix)
  )
  # add SingleCellExperiment to MuData
  cell_barcodes <- colnames(new_matrix)
  mudata@ExperimentList[[experiment_name]] <- grna_assignment_sce
  mudata@sampleMap <- rbind(
    MultiAssayExperiment::sampleMap(mudata),
    data.frame(
      assay = rep(experiment_name, length(cell_barcodes)),
      primary = cell_barcodes,
      colname = cell_barcodes
    )
  )
  # return MuData
  return(mudata)
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
  grna_assignment_matrix <- extract_grna_assignment_matrix(sceptre_object)
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
  batch_cols <- grep("rep|batch", names(covariate_df), ignore.case = TRUE)
  # Check the number of matching columns and act accordingly
  if (length(batch_cols) == 0) {
    # No matching columns, create a DataFrame with all ones
    sample_df <- MultiAssayExperiment::DataFrame(batch = rep(1, nrow(covariate_df)))
  } else if (length(batch_cols) == 1) {
    # One matching column, create a DataFrame with its contents
    sample_df <- MultiAssayExperiment::DataFrame(
      batch = covariate_df[[batch_cols]]
    )
  } else {
    # More than one matching column, throw an error
    stop("Error: More than one column found containing 'rep' or 'batch'")
  }

  # 3. Extra gRNA and gene information
  grna_ids <- rownames(grna_matrix)
  grna_rowdata <- grna_target_data_frame |>
    dplyr::arrange(match(grna_id, grna_ids)) |>
    tibble::column_to_rownames(var = "grna_id") |>
    dplyr::mutate(targeting = grna_target != "non-targeting") |>
    dplyr::relocate(targeting) |>
    dplyr::rename(guide_target = grna_target)
  if("chr" %in% colnames(grna_rowdata)){
    grna_rowdata <- grna_rowdata |>
      dplyr::mutate(chr = ifelse(is.na(chr), "", chr)) |>
      dplyr::rename(guide_chr = chr)
  }
  if("start" %in% colnames(grna_rowdata)){
    grna_rowdata <- grna_rowdata |>
      dplyr::mutate(start = ifelse(is.na(start), -9, start)) |>
      dplyr::rename(guide_start = start)
  }
  if("end" %in% colnames(grna_rowdata)){
    grna_rowdata <- grna_rowdata |>
      dplyr::mutate(end = ifelse(is.na(end), -9, end)) |>
      dplyr::rename(guide_end = end)
  }
  if("grna_group" %in% colnames(grna_rowdata)){
    grna_rowdata <- grna_rowdata |>
      dplyr::select(-grna_group)
  }
  if(sum(!is.na(gene_names)) > 0){
    gene_rowdata <- data.frame(symbol = gene_names)
  } else{
    gene_rowdata <- NULL
  }

  # 4. Prepare for conversion to MuData
  response_matrix <- methods::as(response_matrix, "CsparseMatrix")
  grna_matrix <- methods::as(grna_matrix, "CsparseMatrix")
  if(is.null(colnames(response_matrix))){
    cell_ids <- paste0("cell_", 1:ncol(response_matrix))
  } else{
    cell_ids <- colnames(response_matrix)
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
  metadata <- list(moi = moi)
  if(nrow(pairs) > 0){
    metadata[["inference_results"]] = pairs |>
      dplyr::select(grna_target, response_id, pair_type) |>
      dplyr::rename(gene_id = response_id) |>
      dplyr::mutate(p_value = -9, log_2_FC = -9) |>
      MultiAssayExperiment::DataFrame()
  }

  # 5. Create MuData object
  gene_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = response_matrix),
    rowData = gene_rowdata
  )
  grna_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = grna_matrix),
    rowData = grna_rowdata
  )
  if (!is.null(grna_assignment_matrix)) {
    grna_assignment_se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = grna_assignment_matrix),
    )
  } else {
    grna_assignment_se <- NULL
  }
  experiment_list <- list(gene = gene_se, guide = grna_se)
  if (!is.null(grna_assignment_se)) {
    experiment_list[["guide_assignment"]] <- grna_assignment_se
  }
  mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = experiment_list,
    colData = sample_df,
    metadata = metadata
  )

  return(mae)
}

