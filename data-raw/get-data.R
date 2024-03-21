library(MuData)

mudata_guide_assignment_gasperini <- readH5MU("~/code/research/sceptre-igvf/data/guide_assignment/gasperini_guide_assignment_input.h5mu")
usethis::use_data(mudata_guide_assignment_gasperini, overwrite = TRUE)

mudata_inference_gasperini <- readH5MU("~/code/research/sceptre-igvf/data/inference/gasperini_inference_input.h5mu")
usethis::use_data(mudata_inference_gasperini, overwrite = TRUE)

mudata_guide_assignment_papalexi <- readH5MU("~/code/research/sceptre-igvf/data/guide_assignment/papalexi_guide_assignment_input.h5mu")
usethis::use_data(mudata_guide_assignment_papalexi, overwrite = TRUE)

mudata_inference_papalexi <- readH5MU("~/code/research/sceptre-igvf/data/inference/papalexi_inference_input.h5mu")
usethis::use_data(mudata_inference_papalexi, overwrite = TRUE)


