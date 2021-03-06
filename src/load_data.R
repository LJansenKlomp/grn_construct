# script for opening raw GEO data and saving this data
# in a .mat file, R file and .csv file

library(dplyr)
library(Seurat)
library(R.matlab)

data_dir <- "C:/Data/Code/scimap_code/grn_construct/data"

data_list <- list.dirs(path = data_dir, full.names = TRUE, recursive = FALSE)
data_names <- basename(data_list)

for (k in 1:length(data_names)) {
    # open data from directory
    data_loc <- data_list[k]
    day <- substring(data_names[k], 2)

    # read and normalize data
    c.data <- Read10X(data.dir = data_loc)
    c <- CreateSeuratObject(counts = c.data, project = "ipsc")

    c <- NormalizeData(c, normalization.method = "LogNormalize",
     scale.factor = 10000)

    # extract feature names from data
    all.genes <- rownames(c)
    genes <- list(all.genes)[[1]]

    # save normalized counts and feature names
    writeMat(paste(data_dir, "/ncounts_D", day, sep = ""),
    mat = c[["RNA"]]@data)
    mat <- c[["RNA"]]@data
    save(mat, file = paste(data_dir, "/R_ncounts_D", day, sep = ""))

    write.csv(genes, paste(data_dir, "/feature_nm_D", day, sep = ""))
}