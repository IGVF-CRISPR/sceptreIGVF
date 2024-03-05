library(MuData)
mudata_guide_assignment <- readH5MU("~/code/research/sceptre-igvf/data/guide_assignment/gasperini_guide_assignment_input.h5mu")
usethis::use_data(mudata_guide_assignment, overwrite = TRUE)

mudata_inference <- readH5MU("~/code/research/sceptre-igvf/data/inference/gasperini_inference_input.h5mu")
usethis::use_data(mudata_inference, overwrite = TRUE)
