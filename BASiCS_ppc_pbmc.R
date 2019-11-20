library("Seurat")
## This requires a modified version of basics without a batch requirement
# library("devtools")
# install_github("catavallejos/BASiCS", ref = "batches")
library("BASiCS")
library("ggplot2")
library("rhdf5")

options(stringsAsFactors = FALSE)


pbmc_genes <- read.csv("genes_1000_pbmc3k.csv", header = FALSE)[, 2]
pbmc <- readRDS("pbmc3k.rds")

# table(pbmc$celltype)

pbmc <- pbmc[, pbmc$celltype == "CD4 Memory"]
# table(rowSums(pbmc@assays$RNA@counts[pbmc_genes, ]) != 0)
counts <- as.matrix(pbmc@assays$RNA@counts)[pbmc_genes, ]
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

saveRDS(fit, "cd4mem_mcmc.rds")


zero_genes <- setdiff(pbmc_genes, rownames(counts))
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
    c[pbmc_genes, ]
  }
)

a <- array(dim = c(ncol(counts), 1000, length(mats)))
for (i in 1:length(mats)) {
  a[, , i] <- t(mats[[i]])
}

h5write(a, "basics_pbmc_samples.h5", "x")


# dropout <- function(counts) rowMeans(counts == 0)
# mean <- function(counts) rowMeans(counts)
# complexity <- function(counts) colSums(counts != 0)
# libsize <- function(counts) colSums(counts)

# funs <- c(dropout, mean, complexity, libsize)
# d <- lapply(seq_along(mats),
#   function(ind) {
#     data.frame(
#       ind = ind,
#       dropout = dropout(mats[[ind]]),
#       mean = mean(mats[[ind]])
#     )
#   }
# )

# df <- do.call(rbind, d)
# df <- rbind(df, 
#   data.frame(
#     ind = "true",
#     dropout = dropout(counts),
#     mean = mean(counts)
#   )
# )
# ggplot(df, aes(x = ind, y = dropout)) +
#   geom_violin()

# ggplot(df, aes(x = ind, y = mean)) +
#   geom_violin()




# d <- lapply(seq_along(mats),
#   function(ind) {
#     data.frame(
#       ind = ind,
#       complexity = complexity(mats[[ind]]),
#       libsize = libsize(mats[[ind]])
#     )
#   }
# )

# df <- do.call(rbind, d)
# df <- rbind(df, 
#   data.frame(
#     ind = "true",
#     complexity = complexity(counts),
#     libsize = libsize(counts)
#   )
# )
# ggplot(df, aes(x = ind, y = complexity)) +
#   geom_violin()

# ggplot(df, aes(x = ind, y = libsize)) +
#   geom_violin()

