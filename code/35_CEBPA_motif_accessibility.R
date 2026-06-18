# CEBPA ATAC motif accessibility analysis
# Quantifies per-cell CEBPA motif accessibility in Mid vs Term Hofbauer cells
# Used for Figure 2G
library(Signac); library(Seurat)
library(MotifDb); library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings); library(parallel)

# Load ATAC object
obj <- readRDS("results/Hofbauer_ATAC_mid_term.rds")
peaks <- granges(obj)

# Filter to standard chromosomes
std_chrs <- paste0("chr", c(1:22, "X", "Y"))
keep <- as.character(seqnames(peaks)) %in% std_chrs
peaks_filt <- peaks[keep]
seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, peaks_filt)

# Get CEBPA PWM from MotifDb (human)
cebpa_motifs <- query(MotifDb, c("CEBPA", "hsapiens"))
pwm <- cebpa_motifs[[1]]
max_score <- sum(apply(pwm, 2, max))
min_score <- 0.85 * max_score

# Parallel scan for CEBPA motifs
seq_list <- as.list(seqs)
hit_counts <- unlist(mclapply(seq_list, function(s) {
  tryCatch(countPWM(pwm, s, min.score = min_score), error = function(e) 0)
}, mc.cores = 8))

# Map hits back to original peak indices
cebpa_hit_idx <- which(keep)[which(hit_counts > 0)]

# Per-cell CEBPA motif accessibility
cebpa_mat <- obj@assays$ATAC@counts[cebpa_hit_idx, , drop = FALSE]
obj$CEBPA_motif_access <- log1p(colSums(cebpa_mat) / obj$nCount_ATAC * 10000)

# Statistics
wt <- wilcox.test(CEBPA_motif_access ~ group, data = obj@meta.data)

# Save
saveRDS(obj, "results/Hofbauer_ATAC_CEBPA_motif.rds")
