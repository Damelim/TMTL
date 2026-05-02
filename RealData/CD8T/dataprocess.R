##### FINAL: GSE164378 PBMC CITE-seq (RNA 3P + ADT 3P) #####
library(Matrix)
library(data.table)
library(Seurat)

## -------------------------------------------------
## 0. Paths
## -------------------------------------------------
raw_dir <- "../../../GSE164378_RAW"

rna_matrix <- "GSM5008740_RNA_5P-matrix.mtx.gz"
rna_feat   <- "GSM5008740_RNA_5P-features.tsv.gz"
rna_bar    <- "GSM5008740_RNA_5P-barcodes.tsv.gz"

adt_matrix <- "GSM5008741_ADT_5P-matrix.mtx.gz"
adt_feat   <- "GSM5008741_ADT_5P-features.tsv.gz"
adt_bar    <- "GSM5008741_ADT_5P-barcodes.tsv.gz"

meta_file  <- "GSE164378_sc.meta.data_5P.csv.gz"
hvg_n <- 500

## 1. Load metadata
meta <- fread(file.path(raw_dir, meta_file))
barcode_col <- "V1"
donor_col   <- "donor"

## 2. Load RNA
rna_counts <- readMM(file.path(raw_dir, rna_matrix))
rna_feats  <- fread(file.path(raw_dir, rna_feat), header = FALSE)[[2]]
rna_cells  <- fread(file.path(raw_dir, rna_bar), header = FALSE)[[1]]

rownames(rna_counts) <- rna_feats
colnames(rna_counts) <- rna_cells

dim(meta)
colnames(meta)

table(meta$V1) # cell names
hist(meta$nCount_ADT)
hist(meta$nCount_RNA)
hist(meta$nFeature_ADT)
hist(meta$nFeature_RNA)
table(meta$orig.ident)
table(meta$lane)
table(meta$donor)
library(dplyr)

table(meta$time) # 0, 2, 7 ??
table(meta$celltype.l1) ; table(meta$celltype.l1) %>% sum
# table(meta$celltype.l2) ; table(meta$celltype.l2) %>% sum
# table(meta$celltype.l3) ; table(meta$celltype.l3) %>% sum

table(meta$Phase) #  G1 G2M S 

table(meta$Batch) # Batch1 Batch2 

## 3. Load ADT
adt_counts <- readMM(file.path(raw_dir, adt_matrix))
adt_feats  <- fread(file.path(raw_dir, adt_feat), header = FALSE)[[2]]
adt_cells  <- fread(file.path(raw_dir, adt_bar), header = FALSE)[[1]]

rownames(adt_counts) <- adt_feats
colnames(adt_counts) <- adt_cells

## 4. Define RECORDS = common cells
common_cells <- Reduce(
  intersect,
  list(colnames(rna_counts), colnames(adt_counts), meta[[barcode_col]])
)

subset_cells <- meta$celltype.l1 == "CD8 T"
rna_counts <- rna_counts[, subset_cells]
adt_counts <- adt_counts[, subset_cells]
meta <- meta[subset_cells, ]

stopifnot(
  identical(colnames(rna_counts), colnames(adt_counts)),
  identical(colnames(rna_counts), meta[[barcode_col]])
)
cat("Cells used:", ncol(rna_counts), "\n")









##### hvg selection --> log1p --> scale


gene_names <- rownames(rna_counts)
keep_genes <- !grepl("^MT-|^RPL|^RPS|^HBA|^HBB", gene_names)
rna_counts_filt <- rna_counts[keep_genes, ]
rm(rna_counts)

seu <- CreateSeuratObject(counts = rna_counts_filt)
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = hvg_n, verbose = FALSE)
hvg <- VariableFeatures(seu)



X = t(rna_counts_filt[hvg,]) # BEFORE transformation
apply(X, 2, var) %>% sort(decreasing =T) %>% round(2) 

XX = log1p(X) ; XX = scale(XX)

cormat = cor(XX)

# library(corrplot)
# corrplot(
#   cormat,
#   method = "color",
#   col = colorRampPalette(c("#053061", "#2166AC", "#F7F7F7", "orange", "#67001F"))(200),
#   tl.cex = 0.2,      
#   tl.pos = "n",   
#   cl.cex = 0.8
# )
thr <- 0.8  
high_cor_pairs <- which(abs(cormat) > thr & upper.tri(cormat), arr.ind = TRUE)
nrow(high_cor_pairs)


gene_names <- colnames(XX)
hvg_rank <- setNames(seq_along(hvg), hvg)
to_drop <- c()

for (i in seq_len(nrow(high_cor_pairs))) {
  g1 <- gene_names[high_cor_pairs[i, 1]]
  g2 <- gene_names[high_cor_pairs[i, 2]]
  
  if (!(g1 %in% to_drop) && !(g2 %in% to_drop)) {
    if (hvg_rank[g1] > hvg_rank[g2]) {
      to_drop <- c(to_drop, g1)
    } else {
      to_drop <- c(to_drop, g2)
    }
  }
}
length(to_drop)






#### drop and re-transform ####

X_pruned = X[,!(colnames(X) %in% to_drop)]

XX = X_pruned %>% log1p %>% scale





## 6. ADT → CLR (no scaling)
clr_transform <- function(x) {
  log1p(x / exp(rowMeans(log1p(x + 1e-8))))
}

Y <- t(adt_counts)

zeroprop = function(x){mean(x==0)}
erase_tasks = which((apply(Y,2,zeroprop) > 0.8)) %>% names
erase_tasks
Y = Y[,!(colnames(Y) %in% erase_tasks)]

YY <- clr_transform(Y)

rm(adt_counts)

## 7. Sanity check  
stopifnot(nrow(XX) == nrow(YY), identical(rownames(XX), rownames(YY)))

cat("Final shapes:\n")
cat("X:", dim(XX)[1], "cells x", dim(XX)[2], "genes\n")
cat("Y:", dim(YY)[1], "cells x", dim(YY)[2], "proteins\n")

## 8. Donor-wise split & save

donors <- unique(meta[[donor_col]])
for (i in 1:length(donors)) {
  idx <- which(meta[[donor_col]] == donors[i])
  X <- XX[idx, , drop = FALSE] ; print(dim(X))
  Y <- YY[idx, , drop = FALSE]
  cat("Saving donor", donors[i], "- cells:", nrow(X), "\n")
  save(X, Y, file = paste("clr_data_donor_nofiltered_1000hvg", i, ".Rdata", sep = ""))
}
cat("ALL DONE. X/Y records perfectly aligned.\n")