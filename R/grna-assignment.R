#' Assign gRNAs to cells via sceptre's mixture method
#'
#' @param mudata A MuData object
#'
#' @return A MuData object with sceptre gRNA assignments added to the `guide` modality
#' as a new assay called "guide_assignment"
#' @export
#'
#' @examples
#' library(sceptreIGVF)
#' # load sample MuData
#' data(mudata_guide_assignment_gasperini)
#' # assign gRNAs
#' mudata_out <- assign_grnas_sceptre(mudata_guide_assignment_gasperini)
#' mudata_out
assign_grnas_sceptre <- function(mudata) {
  # convert MuData object to sceptre object, removing multicollinear covariates
  sceptre_object <- convert_mudata_to_sceptre_object(
    mudata,
    remove_collinear_covariates = TRUE
  )

  # set analysis parameters
  sceptre_object <- sceptre_object |>
   sceptre::set_analysis_parameters(
     discovery_pairs = data.frame(
       grna_target = character(0),
       response_id = character(0)
     )
   )

  # assign gRNAs
  sceptre_object <- sceptre_object |> sceptre::assign_grnas(method = "mixture")

  # extract sparse logical matrix of gRNA assignments
  grna_assignment_matrix <- extract_grna_assignment_matrix(sceptre_object)
  colnames(grna_assignment_matrix) <- colnames(MultiAssayExperiment::assay(mudata[['guide']]))

  # add gRNA assignment matrix to MuData
  SummarizedExperiment::assays(mudata[['guide']])[['guide_assignment']] <- grna_assignment_matrix

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
