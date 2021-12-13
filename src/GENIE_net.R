library(GENIE3)
library(data.table)
library(R.matlab)

cell_tp <- "chondrogenic"
day <- 42
localpath <- "C:/Data/phdCode/"
target_genes <- c("COL2A1", "ACAN")
genes <- read.csv(paste(localpath, "data/feature_nm_D1", sep = ""))[, 2]

print("Loading data...")
load(paste(localpath, "data/R_ncounts_D", as.character(day), sep = ""))

node_nms <- rep()
for (k in 1:length(target_genes)) {
  fl_nm <- paste(localpath, "data/", target_genes[k], "_inter.tsv", sep = "")
  nd_dat <- as.data.frame(fread(fl_nm))

  for (j in 1:dim(nd_dat)[1]) {
    nm <- nd_dat[j, 1]
    ind <- match(nm, genes)
    if (!ind %in% node_nms) {
      node_nms <- append(node_nms, ind)
    }
  }
}

node_nms <- sort(node_nms)

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
writeMat(paste(localpath, "data/", cell_tp, "_GENIEnet.mat", sep = ""), x = net)

writeMat(paste(localpath, "data/", cell_tp, "_GENIEgenes.mat", sep = ""), g = gene_set)

print("Done")