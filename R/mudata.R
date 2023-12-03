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
