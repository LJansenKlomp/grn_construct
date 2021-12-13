# script for different stages of differentiation using
# the STRING database as input for which genes need to be included

library(GENIE3)
library(data.table)
library(R.matlab)

# parameters to be set by user
cell_tp <- "chondrogenic"
target_genes <- c("COL2A1", "ACAN")
day <- 42

localpath <- "C:/Data/phdCode/"

# load and organize data
print("Loading data...")

genes <- read.csv(paste(localpath, "data/feature_nm_D1", sep = ""))[, 2]
load(paste(localpath, "data/R_ncounts_D", as.character(day), sep = ""))

# obtain numbers of target genes as in scRNA-seq data
node_nms <- rep()
for (k in 1:length(target_genes)) {

  # included genes per target gene taken from degree .tsv file as
  # downloaded from STRING
  fl_nm <- paste(localpath, "data/", target_genes[k], "_inter.tsv", sep = "")
  nd_dat <- as.data.frame(fread(fl_nm))

  for (j in 1:dim(nd_dat)[1]) {

    # if the current gene is not already a 
    # member of the list of added nodes, add it
    nm <- nd_dat[j, 1]
    ind <- match(nm, genes)
    if (!ind %in% node_nms) {
      node_nms <- append(node_nms, ind)
    }
  }
}

# for ease of reference, sort the obtained indices
node_nms <- sort(node_nms)

# construct network using GENIE3
mat <- as.matrix(mat)
mat <- mat[node_nms, ]
gene_set <- genes[node_nms]

print("Data loaded, constructing network...")
net <- GENIE3(
  mat,
  regulators = NULL,
  targets = NULL,
  treeMethod = "RF",
  K = "sqrt",
  nTrees = 1000,
  nCores = 1,
  returnMatrix = TRUE,
  verbose = FALSE
)

# write results (network and names of features) to .mat files
writeMat(paste(localpath, "data/", cell_tp, "_GENIEnet.mat", sep = ""), x = net)

writeMat(paste(localpath, "data/", cell_tp, "_GENIEgenes.mat", sep = ""),
 g = gene_set)

print("Done")