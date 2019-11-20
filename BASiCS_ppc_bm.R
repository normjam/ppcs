library("Seurat")
## This requires a modified version of basics without a batch requirement
# library("devtools")
# install_github("catavallejos/BASiCS", ref = "batches")
library("BASiCS")
library("ggplot2")
library("rhdf5")

options(stringsAsFactors = FALSE)


bm_genes <- read.csv("genes_1000_bm.cite.csv", header = FALSE)[, 2]
bm <- readRDS("bm.cite.rds")

table(bm$celltype)


bm <- bm[, bm$celltype == "CD4 Naive"]
table(rowSums(bm@assays$RNA@counts[bm_genes, ]) != 0)
counts <- as.matrix(bm@assays$RNA@counts)[bm_genes, ]
counts <- counts[rowSums(counts) != 0, ]

sce <- SingleCellExperiment(assays = list(counts = counts))

fit <- BASiCS_MCMC(
  sce,
  N = 20000,
  Thin = 10,
  Burn = 10000,
  Regression = FALSE,
  WithSpikes = FALSE
)
saveRDS(fit, "cd4naive_mcmc.rds")

zero_genes <- setdiff(bm_genes, rownames(counts))
zero_matrix <- matrix(
  0,
  nrow = length(zero_genes),
  ncol = ncol(counts),
  dimnames = list(zero_genes, colnames(counts))
)
mats <- lapply(
  sample(1000, 50),
  function(i) {
    c <- counts(BASiCS_Draw(fit,  rep("1", ncol(fit@parameters[["nu"]])), i))
    rownames(c) <- colnames(fit@parameters[["mu"]])
    colnames(c) <- colnames(counts)
    c <- rbind(c, zero_matrix)
    c[bm_genes, ]
  }
)



a <- array(dim = c(ncol(counts), 1000, length(mats)))
for (i in 1:length(mats)) {
  a[, , i] <- t(mats[[i]])
}

h5write(a, "basics_bm_samples.h5", "x")
