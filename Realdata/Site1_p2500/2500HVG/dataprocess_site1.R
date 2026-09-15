##### GSE194122 (NeurIPS-2021 CITE-seq, BMMC) — T cells only, domain = T-cell subtype #####
##### preprocessing pipeline identical to the GSE164378 script                       #####
library(zellkonverter)
library(SingleCellExperiment)
library(Matrix)
library(Seurat)
library(dplyr)

## -------------------------------------------------
## 0. Paths and settings
## -------------------------------------------------
h5ad_file <- "~/Library/CloudStorage/GoogleDrive-96limtotoro@gmail.com/My Drive/research_multitask/dataset/RNA_protein/GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad"

hvg_n      <- 2500
det_thr    <- 0.3
min_cells  <- 100
zero_thr_Y <- 0.8     # 단백질: 0 비율이 이 값 초과면 제거 (Y만 사용)
cor_thr    <- 0.8     # 공선성 pruning: |r| > cor_thr 인 쌍에서 하나 제거
cor_XY_thr <- 0.2     # 어떤 유전자와도 |cor| < 이 값이면 그 task 제거

subtype_col <- "cell_type"      # <- domain (T-cell subtype)
donorid_col <- "DonorID"
exclude_donors <- character(0)  # <- c("15078") 하면 제외

## 도메인 간 차이가 donor/batch 효과와 교란되지 않도록 단일 site로 제한
restrict_col   <- "Site"        # 안 쓰려면 NA
restrict_value <- "site1"

## T-cell lineage pattern (ILC 등은 제외)
t_pattern <- "^(CD4\\+ T|CD8\\+ T|T reg|Treg|MAIT|gdT|dnT|T prog)"

## version-safe accessor for the normalized matrix
get_data <- function(obj){
  if (utils::packageVersion("Seurat") >= "5.0.0")
    SeuratObject::LayerData(obj, assay = "RNA", layer = "data")
  else
    Seurat::GetAssayData(obj, slot = "data")
}

clr_transform <- function(x) {
  log1p(x / exp(rowMeans(log1p(x + 1e-8))))   # geometric mean per cell
}

## -------------------------------------------------
## 1. Load and split RNA / ADT
## -------------------------------------------------
sce <- readH5AD(h5ad_file)
stopifnot(subtype_col %in% colnames(colData(sce)),
          donorid_col %in% colnames(colData(sce)))

is_gex <- rowData(sce)$feature_types == "GEX"
is_adt <- rowData(sce)$feature_types == "ADT"
cat("GEX:", sum(is_gex), " ADT:", sum(is_adt), "\n")

rna_counts <- as(assay(sce[is_gex, ], "counts"), "dgCMatrix")   # genes x cells
adt_counts <- as(assay(sce[is_adt, ], "counts"), "dgCMatrix")   # proteins x cells

meta <- as.data.frame(colData(sce))
meta$barcode <- colnames(sce)
stopifnot(identical(colnames(rna_counts), meta$barcode),
          identical(colnames(adt_counts), meta$barcode))
rm(sce); gc()

cat("\n=== Site x DonorID ===\n");   print(table(meta$Site, meta$DonorID))
cat("\n=== cells per Site ===\n");    print(sort(table(meta$Site), decreasing = TRUE))
cat("\n=== cells per DonorID ===\n"); print(sort(table(meta$DonorID), decreasing = TRUE))
cat("\n=== cell type x Site ===\n");  print(table(meta[[subtype_col]], meta$Site))
cat("\nall cell types present:\n")
print(sort(table(meta[[subtype_col]]), decreasing = TRUE))

## -------------------------------------------------
## 1b. Exclude specified donors
## -------------------------------------------------
if (length(exclude_donors) > 0) {
  cat("\ncells before donor exclusion:", ncol(rna_counts), "\n")
  valid_cells <- !(as.character(meta[[donorid_col]]) %in% exclude_donors)
  cat("excluding donor(s):", paste(exclude_donors, collapse = ", "),
      "->", sum(!valid_cells), "cells removed\n")
  rna_counts <- rna_counts[, valid_cells]
  adt_counts <- adt_counts[, valid_cells]
  meta       <- meta[valid_cells, , drop = FALSE]
  cat("cells after  donor exclusion:", ncol(rna_counts), "\n")
}

## -------------------------------------------------
## 1c. Restrict to the T-cell lineage
## -------------------------------------------------
ct <- as.character(meta[[subtype_col]])
t_cells <- grepl(t_pattern, ct)
cat("\nT-cell subtypes matched:\n")
print(sort(table(ct[t_cells]), decreasing = TRUE))
cat("\nNON-T cell types dropped:\n")
print(sort(unique(ct[!t_cells])))
cat("\nT cells:", sum(t_cells), "of", length(t_cells), "\n")
stopifnot(sum(t_cells) > 0)

rna_counts <- rna_counts[, t_cells]
adt_counts <- adt_counts[, t_cells]
meta       <- meta[t_cells, , drop = FALSE]

## -------------------------------------------------
## 1d. Restrict to a single site (or donor)
##     -> heterogeneity across domains reflects cell type, not batch/donor
## -------------------------------------------------
if (!is.na(restrict_col)) {
  stopifnot(restrict_col %in% colnames(meta))
  cat("\ncells before", restrict_col, "restriction:", ncol(rna_counts), "\n")
  keep <- as.character(meta[[restrict_col]]) == restrict_value
  stopifnot(sum(keep) > 0)
  rna_counts <- rna_counts[, keep]
  adt_counts <- adt_counts[, keep]
  meta       <- meta[keep, , drop = FALSE]
  cat("restricted to", restrict_col, "=", restrict_value, "->", ncol(rna_counts), "cells\n")
  cat("cell types remaining:\n")
  print(sort(table(droplevels(factor(meta[[subtype_col]]))), decreasing = TRUE))
}

## -------------------------------------------------
## 2. Drop subtypes with too few cells
## -------------------------------------------------
tab <- sort(table(droplevels(factor(meta[[subtype_col]]))), decreasing = TRUE)
cat("\ncells per subtype (before min_cells filter):\n"); print(tab)
cat("\nsorted sizes:", paste(sort(as.numeric(tab)), collapse = ", "), "\n")
cat("min_cells =", min_cells, "-> keeps", sum(tab >= min_cells), "subtypes\n")

keep_types   <- names(tab)[tab >= min_cells]
subset_cells <- as.character(meta[[subtype_col]]) %in% keep_types
rna_counts <- rna_counts[, subset_cells]
adt_counts <- adt_counts[, subset_cells]
meta       <- meta[subset_cells, , drop = FALSE]
n_input <- ncol(rna_counts)
cat("cells retained:", n_input, "\n")

## -------------------------------------------------
## 3. Gene filtering (detection rate)
## -------------------------------------------------
gene_names <- rownames(rna_counts)
keep_genes <- gene_names   # !grepl("^MT-|^RPL|^RPS|^HBA|^HBB", gene_names)
rna_counts_filt <- rna_counts[keep_genes, ]
rm(rna_counts); gc()

det_rate <- Matrix::rowMeans(rna_counts_filt > 0)
cat("genes before detection filter:", nrow(rna_counts_filt), "\n")
cat("detection rate summary:\n"); print(summary(det_rate))
rna_counts_filt <- rna_counts_filt[det_rate >= det_thr, ]
cat("genes after  detection filter:", nrow(rna_counts_filt), "\n")

## -------------------------------------------------
## 4. Seurat object; drop extreme rows (1% / 99%)
## -------------------------------------------------
seu <- CreateSeuratObject(counts = rna_counts_filt)
seu$nCount_ADT <- Matrix::colSums(adt_counts)

qr <- quantile(seu$nCount_RNA, c(0.01, 0.99))
qa <- quantile(seu$nCount_ADT, c(0.01, 0.99))
cat("nCount_RNA cutoffs:", qr, "\n")
cat("nCount_ADT cutoffs:", qa, "\n")

seu <- subset(seu,
              subset = nCount_RNA >= qr[1] & nCount_RNA <= qr[2] &
                nCount_ADT >= qa[1] & nCount_ADT <= qa[2])

keep_cells <- colnames(seu)
cat("rows kept:", length(keep_cells), "of", n_input, "\n")
adt_counts <- adt_counts[, keep_cells]
meta       <- meta[match(keep_cells, meta$barcode), , drop = FALSE]
stopifnot(identical(colnames(seu), colnames(adt_counts)),
          identical(colnames(seu), meta$barcode))
cat("cells per subtype (after row filter):\n")
print(sort(table(droplevels(factor(meta[[subtype_col]]))), decreasing = TRUE))

## -------------------------------------------------
## 5. LogNormalize -> HVG -> scale
##    X_ij -> log(1 + 1e4 * X_ij / sum_j X_ij)
## -------------------------------------------------
seu <- NormalizeData(seu, normalization.method = "LogNormalize",
                     scale.factor = 1e4, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst",
                            nfeatures = hvg_n, verbose = FALSE)
hvg <- VariableFeatures(seu)
cat("HVGs selected:", length(hvg), "\n")

X  <- as.matrix(Matrix::t(get_data(seu)[hvg, ]))   # cells x genes, already log-normalized
XX <- scale(X)

## -------------------------------------------------
## 6. Collinearity pruning
## -------------------------------------------------
cormat <- cor(XX)
high_cor_pairs <- which(abs(cormat) > cor_thr & upper.tri(cormat), arr.ind = TRUE)
cat("highly correlated pairs:", nrow(high_cor_pairs), "\n")

gene_names <- colnames(XX)
hvg_rank <- setNames(seq_along(hvg), hvg)
to_drop <- c()
for (i in seq_len(nrow(high_cor_pairs))) {
  g1 <- gene_names[high_cor_pairs[i, 1]]
  g2 <- gene_names[high_cor_pairs[i, 2]]
  if (!(g1 %in% to_drop) && !(g2 %in% to_drop)) {
    if (hvg_rank[g1] > hvg_rank[g2]) to_drop <- c(to_drop, g1) else to_drop <- c(to_drop, g2)
  }
}
cat("dropped:", length(to_drop), "\n")

#### drop and re-standardize (do NOT log again) ####
X_pruned <- X[, !(colnames(X) %in% to_drop), drop = FALSE]
XX <- scale(X_pruned)
cat("final gene detection rate:\n")
print(summary(Matrix::rowMeans(get_data(seu)[colnames(XX), ] > 0)))

## -------------------------------------------------
## 7. ADT -> zero-proportion filter -> CLR
## -------------------------------------------------
Y <- as.matrix(Matrix::t(adt_counts))
zeroprop <- function(x){ mean(x == 0) }
zp <- apply(Y, 2, zeroprop)
cat("\nY zero-proportion quantiles:\n"); print(round(quantile(zp, c(0,.25,.5,.75,.9,1)), 3))
erase_tasks <- names(which(zp > zero_thr_Y))
cat("proteins dropped (", length(erase_tasks), "):\n"); print(erase_tasks)
Y  <- Y[, !(colnames(Y) %in% erase_tasks), drop = FALSE]
YY <- clr_transform(Y)
rm(adt_counts); gc()

## -------------------------------------------------
## 7b. Drop tasks with no linear association to any predictor
##     (max_j |cor(X_j, Y_k)| < cor_XY_thr)
## -------------------------------------------------
cxy  <- abs(cor(XX, YY))                 # p x K
cmax <- apply(cxy, 2, max)
cat("\nmax |cor(X, Y_k)| quantiles:\n")
print(round(quantile(cmax, c(0,.1,.25,.5,.75,.9,1)), 3))
drop_tasks <- names(which(cmax < cor_XY_thr))
cat("tasks dropped by cor_XY_thr =", cor_XY_thr, ":", length(drop_tasks),
    "of", ncol(YY), "\n")
if (length(drop_tasks)) print(drop_tasks)
YY <- YY[, cmax >= cor_XY_thr, drop = FALSE]
cat("K after task screening:", ncol(YY), "\n")
stopifnot(ncol(YY) >= 2)

## -------------------------------------------------
## 8. Sanity check
## -------------------------------------------------
stopifnot(nrow(XX) == nrow(YY), identical(rownames(XX), rownames(YY)))
cat("\nFinal shapes:\n")
cat("X:", dim(XX)[1], "cells x", dim(XX)[2], "genes\n")
cat("Y:", dim(YY)[1], "cells x", dim(YY)[2], "proteins\n")

## 8b. Diagnostic scatterplots on the matrices used for fitting
par(mfrow = c(6,6), mar = c(2,2,1,1))
for(j in 1:6){
  for(i in 1:6){
    plot(XX[,i], YY[,j], pch = 16, cex = 0.3, col = rgb(0,0,0,0.15),
         xlab = "", ylab = "")
    lines(lowess(XX[,i], YY[,j]), col = "red", lwd = 1.5)
  }
}
par(mfrow = c(1,1))

## -------------------------------------------------
## 9. Subtype-wise split & save
## -------------------------------------------------
domains <- sort(unique(as.character(meta[[subtype_col]])))
cat("domains (T-cell subtype):", paste(domains, collapse = ", "), "\n")
cat("number of domains:", length(domains), "\n")
write.csv(data.frame(index = seq_along(domains), celltype = domains),
          "domain_index.csv", row.names = FALSE)

for (i in seq_along(domains)) {
  idx <- which(as.character(meta[[subtype_col]]) == domains[i])
  X <- XX[idx, , drop = FALSE]
  Y <- YY[idx, , drop = FALSE]
  print(dim(X)); print(dim(Y))
  cat("Saving domain", domains[i], "- cells:", nrow(X), "\n")
  save(X, Y, file = paste("clr_data_donor_nofiltered_2500hvg", i, ".Rdata", sep = ""))
}
cat("ALL DONE. X/Y records perfectly aligned.\n")
