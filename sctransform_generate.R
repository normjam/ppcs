# Learn a GLM with regularized NB error distribution, then sample data
# Use the sctransform package with commit hash 3a389a7 (that's where the generate function was added)
# https://github.com/ChristophH/sctransform/commit/3a389a7ab3fd624f9ebf87c48b3a3e2bbbeb8e36

library('sctransform')
library('Seurat')
library('rhdf5')

# some of the vst steps can use multiple cores, set number here
options(future.fork.enable = TRUE)
options(future.globals.maxSize = 10 * 1024 ^ 3)
future::plan(strategy = 'multicore', workers = 4)

# number of samples
N <- 50

pars <- list(dataset_name = 'pbmc3k',
             genes_file = '~/Projects/normjam/genes_1000_pbmc3k.csv',
             seurat_obj_file = '~/Projects/normjam/pbmc3k.rds')

pars <- list(dataset_name = 'bone_marrow_cd4_naive',
             genes_file = '~/Projects/normjam/genes_1000_bm.cite.csv',
             seurat_obj_file = '~/Projects/normjam/bm.cite.rds',
             use_celltype = 'CD4 Naive')

genes <- read.csv(pars$genes_file, header = FALSE, stringsAsFactors = FALSE)$V2
s <- readRDS(pars$seurat_obj_file)
if (!is.null(pars$use_celltype)) {
  cm <- s[['RNA']]@counts[, s@meta.data$celltype == pars$use_celltype]
} else {
  cm <- s[['RNA']]@counts
}

rm(s)

set.seed(68294280)
# learn the model below
vst_out <- vst(cm, return_cell_attr = TRUE, res_clip_range = c(-Inf, Inf))

# now sample data N times
res <- list()
for (i in 1:N) {
  message('iteration\t', i)
  res[[i]] <- array(0, dim = c(ncol(cm), length(genes)), dimnames = list(colnames(cm), genes))
  generated_data <- generate(vst_out, genes = genes)
  genes_use <- genes[genes %in% rownames(generated_data)]
  res[[i]][, genes_use] <- t(generated_data[genes_use, ])
}

# combine results and write to h5 file
res <- array(unlist(res), dim = c(ncol(cm), length(genes), N), dimnames = list(colnames(cm), genes, NULL))

h5_file <- sprintf('~/Projects/normjam/%s_generated_sctransform.h5', pars$dataset_name)
unlink(x = h5_file)
h5write(obj = res, file = h5_file, name = 'x', createnewfile = TRUE)
