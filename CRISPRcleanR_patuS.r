# Instalar y cargar paquetes necesarios
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
library(devtools)
library(readxl)
devtools::install_github("francescojm/CRISPRcleanR")
library(CRISPRcleanR)

# --- 1. Leer y preparar la anotación ---
anno_raw <- read_excel("TFM_12O/TFM_previo/A32170 1843168 LentiArray Epigenetics Glycerol Target Plan.xlsx")
LibraryNorm <- data.frame(
  CODE      = as.character(anno_raw$crispr_id),
  GENES     = sub("-Human$", "", anno_raw$gene_name),
  EXONE     = sapply(strsplit(anno_raw$exon, ", "), function(x) {
    exons_clean <- gsub("exon", "", x)
    exon_mas_frecuente <- names(sort(table(exons_clean), decreasing = TRUE))[1]
    paste0("ex", exon_mas_frecuente)
  }),
  CHRM      = gsub("chr", "", anno_raw$chromosome),
  STRAND    = gsub("[()]", "", anno_raw$chr_direction),
  STARTpos  = anno_raw$chr_start,
  ENDpos    = anno_raw$chr_stop,
  seq       = anno_raw$`Target Sequence`
)
rownames(LibraryNorm) <- LibraryNorm$CODE

# --- 2. Leer y preparar la matriz de counts ---
df_counts <- read.table("csv2/PatuS_X.tsv", header=TRUE, sep="\t", check.names=FALSE)
rownames(df_counts) <- df_counts$sample; df_counts$sample <- NULL
counts <- as.data.frame(t(df_counts))
counts$sgRNA <- rownames(counts)
counts$sgRNA <- as.character(counts$sgRNA)
counts$gene <- LibraryNorm$GENES[match(counts$sgRNA, LibraryNorm$CODE)]
counts <- counts[, c("sgRNA", "gene", setdiff(colnames(counts), c("sgRNA", "gene")))]

# --- 3. Definir controles y reordenar columnas ---
control_names <- c(
  "PatuS-7d-1", "PatuS-7d-2", "PatuS-7d-3", "PatuS-7d-4", "PatuS-7d-5",
  "PatuS-14d-1", "PatuS-14d-2", "PatuS-14d-3", "PatuS-14d-4", "PatuS-14d-5",
  "PatuS-21d-1", "PatuS-21d-2", "PatuS-21d-3", "PatuS-21d-4", "PatuS-21d-5",
  "pEPIGEN1", "pEPIGEN2"
)
other_names <- setdiff(colnames(counts), c("sgRNA", "gene", control_names))
counts <- counts[, c("sgRNA", "gene", control_names, other_names)]
numcontrols <- length(control_names)

# --- 4. Leer y preparar matriz CNV ---
cnv <- read.table("TFM_12O/CNV_correccion/PatuS/PatuS_CNV.txt", header=TRUE, sep="\t")
colnames(cnv)[1:2] <- c("GENES", "LogFC")
cnv_filtrado <- cnv[cnv$GENES %in% LibraryNorm$GENES, ]
sgRNA_CNV <- merge(LibraryNorm[, c("CODE", "GENES")], cnv_filtrado, by = "GENES")
sgRNAs_faltantes <- LibraryNorm[!(LibraryNorm$CODE %in% sgRNA_CNV$CODE), c("GENES", "CODE")]
sgRNAs_faltantes$LogFC <- 0.00
sgRNA_CNV_fill <- rbind(sgRNA_CNV, sgRNAs_faltantes)
orden_samples <- c(
  control_names,
  "PatuS-Cas9-7d-1", "PatuS-Cas9-7d-2", "PatuS-Cas9-7d-3", "PatuS-Cas9-7d-4", "PatuS-Cas9-7d-5",
  "PatuS-Cas9-14d-1", "PatuS-Cas9-14d-2", "PatuS-Cas9-14d-3", "PatuS-Cas9-14d-4", "PatuS-Cas9-14d-5",
  "PatuS-Cas9-21d-1", "PatuS-Cas9-21d-2", "PatuS-Cas9-21d-3", "PatuS-Cas9-21d-4", "PatuS-Cas9-21d-5"
)
CNV_logFCs <- data.frame(
  CODE = sgRNA_CNV_fill$CODE,
  GENES = sgRNA_CNV_fill$GENES
)
for (rep in orden_samples) CNV_logFCs[[rep]] <- sgRNA_CNV_fill$LogFC
rownames(CNV_logFCs) <- as.character(CNV_logFCs$CODE)
CNV_logFCs <- CNV_logFCs[match(LibraryNorm$CODE, rownames(CNV_logFCs)), ]
colnames(CNV_logFCs)[1:2] <- c("sgRNA", "gene")
CNV_logFCs$sgRNA <- as.character(CNV_logFCs$sgRNA)
rownames(CNV_logFCs) <- CNV_logFCs$sgRNA

# --- 5. Corrección CNV con CRISPRcleanR ---
stopifnot(all(rownames(CNV_logFCs) %in% rownames(LibraryNorm)))
stopifnot(!any(is.na(LibraryNorm[rownames(CNV_logFCs), "CHRM"])))
mappedLogFCs <- ccr.logFCs2chromPos(
  foldchanges = CNV_logFCs,
  libraryAnnotation = LibraryNorm
)
CNVLibrary_norm_patuS <- ccr.GWclean(
  gwSortedFCs = mappedLogFCs,
  label = "CNVLibrary_norm_patuS",
  display = TRUE,
  saveTO = "TFM_12O/CNV_correccion/PatuS/"
)
corrected_counts <- ccr.correctCounts(
  CL = "PatuS",
  normalised_counts = counts,
  libraryAnnotation = LibraryNorm,
  correctedFCs_and_segments = CNVLibrary_norm_patuS,
  OutDir = "TFM_12O/CNV_correccion/PatuS/",
  ncontrols = numcontrols
)

# --- 6. Guardar resultados ---
write.table(corrected_counts, file = "csv2/PatuS_CNVcorrected.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, dec = ".")