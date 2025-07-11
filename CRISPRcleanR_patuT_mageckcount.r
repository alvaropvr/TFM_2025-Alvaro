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
LibraryNorm$CODE <- as.character(LibraryNorm$CODE)
rownames(LibraryNorm) <- LibraryNorm$CODE

# --- 2. Leer y preparar la matriz de counts ---
df_counts <- read.table("TFM_12O/02.Counts/mageckcount/all_samples_PatuT.count.txt", header=TRUE, sep="\t", check.names=FALSE)
cols_to_convert <- setdiff(colnames(df_counts), c("sgRNA", "Gene"))
df_counts[cols_to_convert] <- lapply(df_counts[cols_to_convert], function(x) as.numeric(as.character(x)))
counts <- df_counts
counts$sgRNA <- as.character(counts$sgRNA)
rownames(counts) <- counts$sgRNA

# --- 3. Definir controles y reordenar columnas ---
control_names <- c(
  "PatuT-7d-1", "PatuT-7d-2", "PatuT-7d-3", "PatuT-7d-4", "PatuT-7d-5",
  "PatuT-14d-1", "PatuT-14d-2", "PatuT-14d-3", "PatuT-14d-4", "PatuT-14d-5",
  "PatuT-21d-1", "PatuT-21d-2", "PatuT-21d-3", "PatuT-21d-4", "PatuT-21d-5",
  "pEPIGEN1", "pEPIGEN2"
)
other_names <- setdiff(colnames(counts), c("sgRNA", "Gene", control_names))
counts <- counts[, c("sgRNA", "Gene", control_names, other_names)]
colnames(counts)[1:2] <- c("sgRNA", "gene")
numcontrols <- length(control_names)

# --- 4. Leer y preparar matriz CNV ---
cnv <- read.table("TFM_12O/CNV_correccion/PatuT/PatuT_CNV.txt", header=TRUE, sep="\t")
colnames(cnv)[1:2] <- c("GENES", "LogFC")
cnv$GENES <- as.character(cnv$GENES)
LibraryNorm$GENES <- as.character(LibraryNorm$GENES)

# Merge robusto: asegura que todos los sgRNAs de LibraryNorm están presentes y bien alineados
sgRNA_CNV <- merge(LibraryNorm[, c("CODE", "GENES")], cnv, by = "GENES", all.x = TRUE, sort = FALSE)
sgRNA_CNV$LogFC[is.na(sgRNA_CNV$LogFC)] <- 0.00
sgRNA_CNV$CODE <- as.character(sgRNA_CNV$CODE)

# Asegura que el orden es el mismo que LibraryNorm
sgRNA_CNV <- sgRNA_CNV[match(LibraryNorm$CODE, sgRNA_CNV$CODE), ]
stopifnot(all(sgRNA_CNV$CODE == LibraryNorm$CODE))
stopifnot(!any(is.na(sgRNA_CNV$CODE)))

# Construir CNV_logFCs con los mismos sgRNAs y orden que counts
orden_samples <- c(
  control_names,
  "PatuT-Cas9-7d-1", "PatuT-Cas9-7d-2", "PatuT-Cas9-7d-3", "PatuT-Cas9-7d-4", "PatuT-Cas9-7d-5",
  "PatuT-Cas9-14d-1", "PatuT-Cas9-14d-2", "PatuT-Cas9-14d-3", "PatuT-Cas9-14d-4", "PatuT-Cas9-14d-5",
  "PatuT-Cas9-21d-1", "PatuT-Cas9-21d-2", "PatuT-Cas9-21d-3", "PatuT-Cas9-21d-4", "PatuT-Cas9-21d-5"
)
CNV_logFCs <- data.frame(
  sgRNA = counts$sgRNA,
  gene = counts$gene
)
for (rep in orden_samples) CNV_logFCs[[rep]] <- sgRNA_CNV$LogFC[match(counts$sgRNA, sgRNA_CNV$CODE)]
CNV_logFCs$sgRNA <- as.character(CNV_logFCs$sgRNA)
rownames(CNV_logFCs) <- CNV_logFCs$sgRNA

# --- 5. Corrección CNV con CRISPRcleanR ---
stopifnot(nrow(CNV_logFCs) == nrow(counts))
stopifnot(all(CNV_logFCs$sgRNA == counts$sgRNA))
stopifnot(all(CNV_logFCs$sgRNA == LibraryNorm$CODE))
stopifnot(!any(is.na(CNV_logFCs$sgRNA)))
stopifnot(!any(is.na(rownames(CNV_logFCs))))

mappedLogFCs <- ccr.logFCs2chromPos(
  foldchanges = CNV_logFCs,
  libraryAnnotation = LibraryNorm
)
CNVLibrary_norm_patuT <- ccr.GWclean(
  gwSortedFCs = mappedLogFCs,
  label = "CNVLibrary_norm_patuT_mageck",
  display = TRUE,
  saveTO = "TFM_12O/CNV_correccion/PatuT/"
)

corrected_counts <- ccr.correctCounts(
  CL = "PatuT_mageck",
  normalised_counts = counts,
  libraryAnnotation = LibraryNorm,
  correctedFCs_and_segments = CNVLibrary_norm_patuT,
  OutDir = "TFM_12O/CNV_correccion/PatuT/",
  ncontrols = numcontrols
)


# --- 6. Guardar resultados ---
write.table(corrected_counts, file = "csv2/PatuT_mageck_CNVcorrected.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, dec = ".")
