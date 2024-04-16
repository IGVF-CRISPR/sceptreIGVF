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
  grna_assignment_matrix <- sceptre_object |>
    sceptre::get_grna_assignments() |>
    methods::as("dsparseMatrix")
  colnames(grna_assignment_matrix) <- colnames(MultiAssayExperiment::assay(mudata[['guide']]))

  # add gRNA assignment matrix to MuData
  SummarizedExperiment::assays(mudata[['guide']])[['guide_assignment']] <- grna_assignment_matrix

  # return MuData
  return(mudata)
}
