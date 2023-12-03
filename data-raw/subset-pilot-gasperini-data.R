library(MuData)

# load MuData
mudata <- readH5MU("data-raw/Gasperini_2019_sample_pilot.h5mu")

# subset MuData
cell_barcodes <- mudata[["scRNA"]] |> colnames()
cell_barcodes_subset <- cell_barcodes[1:500]
sample_mudata <- mudata[, cell_barcodes_subset]

# write MuData to package
usethis::use_data(sample_mudata)
