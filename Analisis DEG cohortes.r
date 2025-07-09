# Instalar BiocManager primero si no está disponible
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Instalar paquetes de Bioconductor
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}
if (!requireNamespace("enrichplot", quietly = TRUE)) {
  BiocManager::install("enrichplot")
}

# Instalar paquetes de CRAN
if (!requireNamespace("msigdbr", quietly = TRUE)) {
  install.packages("msigdbr")
}

library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)

library(arrow)
library(tidyverse)
library(limma)
library(ggplot2)
if (!require(ggrepel)) {
  install.packages("ggrepel")
  library(ggrepel)
}

# Leer datos
data <- read_feather("../Cohortes/Combat_Data/COMBAT_with5datasets/Combat5_Final_AUseq_Log2TPMplus1.feather")

# Columnas clínicas reales
clin_cols <- c("ID_Paciente", "Tiempo", "Supervivencia", "Edad", "Sexo", "Grado agrupado", "Estadío agrupado", "Subtipo molecular", "Dataset")

# Extrae solo la matriz de expresión
gene_expr <- data %>% select(-all_of(clin_cols))

# Crear clasificación High/Low por mediana de PRMT7
prmt7_expr <- gene_expr$PRMT7
grupo <- ifelse(prmt7_expr > median(prmt7_expr, na.rm=TRUE), "High", "Low")

# Elimina PRMT7 de la matriz de expresión para el DEG
gene_expr <- gene_expr %>% select(-PRMT7)

# Transponer para limma: genes en filas, muestras en columnas
expr_mat <- as.matrix(t(gene_expr))

# Crear diseño
grupo <- factor(grupo)
design <- model.matrix(~ 0 + grupo)
colnames(design) <- levels(grupo)

# DEG con limma
fit <- lmFit(expr_mat, design)
contrast <- makeContrasts(High - Low, levels = design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

# Resultados
res <- topTable(fit2, adjust = "fdr", number = Inf)

# Volcano plot para PRMT7 - estilo exacto como referencia
res$gene <- rownames(res)
res$log10FDR <- -log10(res$adj.P.Val)

# Sistema de colores basado en significancia y magnitud del fold change
res$color <- ifelse(res$adj.P.Val < 0.05 & res$logFC >= 1, "red",
                   ifelse(res$adj.P.Val < 0.05 & res$logFC > 0 & res$logFC < 1, "lightpink",
                          ifelse(res$adj.P.Val < 0.05 & res$logFC <= -1, "blue",
                                 ifelse(res$adj.P.Val < 0.05 & res$logFC < 0 & res$logFC > -1, "lightblue", "grey"))))

# Etiquetas de genes removidas

# Contar genes significativos para las anotaciones (usando FDR < 0.05)
down_sig <- sum(res$adj.P.Val < 0.05 & res$logFC < 0, na.rm=TRUE)
down_2fold <- sum(res$adj.P.Val < 0.05 & res$logFC < -1, na.rm=TRUE)
up_sig <- sum(res$adj.P.Val < 0.05 & res$logFC > 0, na.rm=TRUE)
up_2fold <- sum(res$adj.P.Val < 0.05 & res$logFC > 1, na.rm=TRUE)

p_prmt7 <- ggplot(res, aes(x = logFC, y = -log10(adj.P.Val), color = color)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 0.5) +
  scale_color_manual(
    values = c("grey" = "#D3D3D3", "lightblue" = "#ADD8E6", "blue" = "#0000FF", 
               "lightpink" = "#FFB6C1", "red" = "#FF0000"),
    labels = c("NS", "Down (NS)", "Down (Sig)", "Up (NS)", "Up (Sig)")) +
  labs(title = "Volcano Plot: Low vs High PRMT7 Expression", 
       x = expression(log[2]~fold~change), 
       y = expression(-log[10]~(FDR)), 
       color = "Expression") +  
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        legend.position = "top",
        legend.box.spacing = unit(0.2, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", size = 0.8),
        axis.ticks = element_line(color = "black", size = 0.6)) +
  scale_x_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 0.5)) +
  scale_y_continuous(limits = c(0, 2.5), 
                     breaks = c(0, 1.0, 1.5, 2.0, 2.5),
                     labels = c("0", "1.0", "1.5", "2.0", "2.5")) +
  annotate("text", x = -1.5, y = 2.2, 
           label = paste0("Downregulated\n", down_sig, " genes\n(", down_2fold, " > 2-fold)"), 
           color = "black", size = 3.5, fontface = "bold",
           hjust = 0.5, vjust = 1) +
  annotate("text", x = 1.5, y = 2.2, 
           label = paste0("Upregulated\n", up_sig, " genes\n(", up_2fold, " > 2-fold)"), 
           color = "black", size = 3.5, fontface = "bold",
           hjust = 0.5, vjust = 1)

print(p_prmt7)

# Guardar el plot de PRMT7
ggsave("../Gráficos e imágenes/volcano_PRMT7.png", plot = p_prmt7, width = 10, height = 8, dpi = 300)
cat("Plot de PRMT7 guardado como volcano_PRMT7.png\n")

# ===== ANÁLISIS GSEA PARA PRMT7 =====

# 1. Obtener gene sets desde MSigDB
cat("Obteniendo gene sets...\n")

# Hallmarks de cáncer (H)
hallmarks <- msigdbr(species = "Homo sapiens", category = "H")

# Preparar lista rankeada de genes para GSEA (todos los genes)
gene_list_prmt7 <- res$logFC
names(gene_list_prmt7) <- rownames(res)
gene_list_prmt7 <- sort(gene_list_prmt7, decreasing = TRUE)

# Diagnóstico de datos para GSEA
cat("=== DIAGNÓSTICO PRMT7 ===\n")
cat("GSEA PRMT7 - Total de genes rankeados:", length(gene_list_prmt7), "\n")
cat("GSEA PRMT7 - Rango logFC:", round(min(gene_list_prmt7), 3), "a", round(max(gene_list_prmt7), 3), "\n")
cat("Genes con logFC > |1|:", sum(abs(gene_list_prmt7) > 1), "\n")
cat("Genes con logFC > |0.5|:", sum(abs(gene_list_prmt7) > 0.5), "\n")
cat("Distribución logFC - Q1:", round(quantile(gene_list_prmt7, 0.25), 3), 
    "Mediana:", round(median(gene_list_prmt7), 3), 
    "Q3:", round(quantile(gene_list_prmt7, 0.75), 3), "\n")


# Ejecutar GSEA con Cancer Hallmarks
cat("Ejecutando GSEA con Cancer Hallmarks (con semilla fija)...\n")
gsea_hallmarks <- GSEA(geneList = gene_list_prmt7,
                      TERM2GENE = hallmarks[,c("gs_name", "gene_symbol")],
                      minGSSize = 15,
                      maxGSSize = 500,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      verbose = FALSE,
                      seed = TRUE)

cat("Hallmarks enriquecidos en GSEA (PRMT7):", nrow(gsea_hallmarks@result), "\n")

if(nrow(gsea_hallmarks@result) > 0) {
  # Mostrar todos los términos significativos
  print(gsea_hallmarks@result[,c("Description", "setSize", "enrichmentScore", "NES", "p.adjust")])
  
  # Crear gráfica GSEA dotplot
  gsea_dot_prmt7 <- dotplot(gsea_hallmarks, showCategory = nrow(gsea_hallmarks@result), split = ".sign") +
    facet_grid(.~.sign) +
    labs(title = "GSEA Analysis - PRMT7 (Cancer Hallmarks)",
         x = "Gene Ratio") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 10),
          legend.title = element_text(size = 10, face = "bold"),
          strip.text = element_text(size = 10, face = "bold"))
  
  print(gsea_dot_prmt7)
  
  # Guardar GSEA dotplot
  ggsave("../Gráficos e imágenes/gsea_hallmarks_PRMT7.png", plot = gsea_dot_prmt7, 
         width = 14, height = 8, dpi = 300)
  cat("GSEA hallmarks de PRMT7 guardado como gsea_hallmarks_PRMT7.png\n")

 # Crear enrichment plots para los términos más significativos
num_plots <- min(2, nrow(gsea_hallmarks@result))  # Máximo 2 plots

for(i in 1:num_plots) {
  pathway_id <- gsea_hallmarks@result$ID[i]
  pathway_desc <- gsea_hallmarks@result$Description[i]
  pathway_nes <- round(gsea_hallmarks@result$NES[i], 2)
  pathway_fdr <- round(gsea_hallmarks@result$p.adjust[i], 3)
  
  enrich_plot <- gseaplot2(gsea_hallmarks, geneSetID = pathway_id, 
                          title = paste0("PRMT7 - ", pathway_desc, 
                                       "\nNES: ", pathway_nes, ", FDR: ", pathway_fdr),
                          color = "red", base_size = 11)
  
  print(enrich_plot)
  
  # Guardar cada plot con nombre diferente
  filename <- paste0("../Gráficos e imágenes/gsea_enrichment_PRMT7_hallmarks_top", i, ".png")
  ggsave(filename, plot = enrich_plot, width = 12, height = 8, dpi = 300)
  cat(paste("GSEA enrichment plot", i, "de PRMT7 guardado como gsea_enrichment_PRMT7_hallmarks_top", i, ".png\n"))
}

} else {
  cat("No se encontraron hallmarks enriquecidos en GSEA para PRMT7\n")
}

cat("\n=== Análisis GSEA para PRMT7 completado ===\n\n")

# ===== ANÁLISIS DEG PARA EHMT2 =====

# Extraer matriz de expresión (sin columnas clínicas)
gene_expr_ehmt2 <- data %>% select(-all_of(clin_cols))

# Crear clasificación High/Low por mediana de EHMT2
ehmt2_expr <- gene_expr_ehmt2$EHMT2
grupo_ehmt2 <- ifelse(ehmt2_expr > median(ehmt2_expr, na.rm=TRUE), "High", "Low")

# Elimina EHMT2 de la matriz de expresión para el DEG
gene_expr_ehmt2 <- gene_expr_ehmt2 %>% select(-EHMT2)

# Transponer para limma: genes en filas, muestras en columnas
expr_mat_ehmt2 <- as.matrix(t(gene_expr_ehmt2))

# Crear diseño
grupo_ehmt2 <- factor(grupo_ehmt2)
design_ehmt2 <- model.matrix(~ 0 + grupo_ehmt2)
colnames(design_ehmt2) <- levels(grupo_ehmt2)

# DEG con limma
fit_ehmt2 <- lmFit(expr_mat_ehmt2, design_ehmt2)
contrast_ehmt2 <- makeContrasts(High - Low, levels = design_ehmt2)
fit2_ehmt2 <- contrasts.fit(fit_ehmt2, contrast_ehmt2)
fit2_ehmt2 <- eBayes(fit2_ehmt2)

# Resultados
res_ehmt2 <- topTable(fit2_ehmt2, adjust = "fdr", number = Inf)

# Volcano plot para EHMT2 - estilo mejorado con etiquetas
res_ehmt2$gene <- rownames(res_ehmt2)
res_ehmt2$log10FDR <- -log10(res_ehmt2$adj.P.Val)

# Sistema de colores basado en significancia y magnitud del fold change
res_ehmt2$color <- ifelse(res_ehmt2$adj.P.Val < 0.05 & res_ehmt2$logFC >= 1, "red",
                   ifelse(res_ehmt2$adj.P.Val < 0.05 & res_ehmt2$logFC > 0 & res_ehmt2$logFC < 1, "lightpink",
                          ifelse(res_ehmt2$adj.P.Val < 0.05 & res_ehmt2$logFC <= -1, "blue",
                                 ifelse(res_ehmt2$adj.P.Val < 0.05 & res_ehmt2$logFC < 0 & res_ehmt2$logFC > -1, "lightblue", "grey"))))

# Etiquetas de genes removidas

# Contar genes significativos para las anotaciones (usando FDR < 0.05)
down_sig_ehmt2 <- sum(res_ehmt2$adj.P.Val < 0.05 & res_ehmt2$logFC < 0, na.rm=TRUE)
down_2fold_ehmt2 <- sum(res_ehmt2$adj.P.Val < 0.05 & res_ehmt2$logFC < -1, na.rm=TRUE)
up_sig_ehmt2 <- sum(res_ehmt2$adj.P.Val < 0.05 & res_ehmt2$logFC > 0, na.rm=TRUE)
up_2fold_ehmt2 <- sum(res_ehmt2$adj.P.Val < 0.05 & res_ehmt2$logFC > 1, na.rm=TRUE)

p_ehmt2 <- ggplot(res_ehmt2, aes(x = logFC, y = -log10(adj.P.Val), color = color)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size = 0.5) +
  scale_color_manual(
    values = c("grey" = "#D3D3D3", "lightblue" = "#ADD8E6", "blue" = "#0000FF", 
               "lightpink" = "#FFB6C1", "red" = "#FF0000"),
    labels = c("NS", "Down (NS)", "Down (Sig)", "Up (NS)", "Up (Sig)")) +
  labs(title = "Volcano Plot: Low vs High EHMT2 Expression", 
       x = expression(log[2]~fold~change), 
       y = expression(-log[10]~(FDR)), 
       color = "Expression") +  
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        legend.position = "top",
        legend.box.spacing = unit(0.2, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", size = 0.8),
        axis.ticks = element_line(color = "black", size = 0.6)) +
  scale_x_continuous(limits = c(-2, 2), breaks = seq(-2, 2, 0.5)) +
  scale_y_continuous(limits = c(0, 2.5), 
                     breaks = c(0, 1.0, 1.5, 2.0, 2.5),
                     labels = c("0", "1.0", "1.5", "2.0", "2.5")) +
  annotate("text", x = -1.5, y = 2.2, 
           label = paste0("Downregulated\n", down_sig_ehmt2, " genes\n(", down_2fold_ehmt2, " > 2-fold)"), 
           color = "black", size = 3.5, fontface = "bold",
           hjust = 0.5, vjust = 1) +
  annotate("text", x = 1.5, y = 2.2, 
           label = paste0("Upregulated\n", up_sig_ehmt2, " genes\n(", up_2fold_ehmt2, " > 2-fold)"), 
           color = "black", size = 3.5, fontface = "bold",
           hjust = 0.5, vjust = 1)

print(p_ehmt2)

# Guardar el plot de EHMT2
ggsave("../Gráficos e imágenes/volcano_EHMT2.png", plot = p_ehmt2, width = 10, height = 8, dpi = 300)
cat("Plot de EHMT2 guardado como volcano_EHMT2.png\n")

# ===== ANÁLISIS GSEA PARA EHMT2 CON HALLMARKS =====

# Preparar lista rankeada de genes para GSEA (todos los genes)
# Rankear por logFC (todos los genes, no solo significativos)
gene_list_ehmt2 <- res_ehmt2$logFC
names(gene_list_ehmt2) <- rownames(res_ehmt2)
gene_list_ehmt2 <- sort(gene_list_ehmt2, decreasing = TRUE)

cat("GSEA EHMT2 - Total de genes rankeados:", length(gene_list_ehmt2), "\n")
cat("GSEA EHMT2 - Rango logFC:", round(min(gene_list_ehmt2), 3), "a", round(max(gene_list_ehmt2), 3), "\n")

# GSEA con Hallmarks de cáncer
gsea_ehmt2 <- GSEA(geneList = gene_list_ehmt2,
                   TERM2GENE = hallmarks[,c("gs_name", "gene_symbol")],
                   minGSSize = 15,
                   maxGSSize = 500,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   verbose = FALSE)

cat("Hallmarks enriquecidos en GSEA (EHMT2):", nrow(gsea_ehmt2@result), "\n")

if(nrow(gsea_ehmt2@result) > 0) {
  # Mostrar todos los términos significativos
  print(gsea_ehmt2@result[,c("Description", "setSize", "enrichmentScore", "NES", "p.adjust")])
  
  # Crear gráfica GSEA dotplot
  gsea_dot_ehmt2 <- dotplot(gsea_ehmt2, showCategory = nrow(gsea_ehmt2@result), split = ".sign") +
    facet_grid(.~.sign) +
    labs(title = "GSEA Analysis - EHMT2 (Cancer Hallmarks)",
         x = "Gene Ratio") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 10),
          legend.title = element_text(size = 10, face = "bold"),
          strip.text = element_text(size = 10, face = "bold"))
  
  print(gsea_dot_ehmt2)
  
  # Guardar GSEA dotplot
  ggsave("../Gráficos e imágenes/gsea_hallmarks_EHMT2.png", plot = gsea_dot_ehmt2, 
         width = 14, height = 8, dpi = 300)
  cat("GSEA hallmarks de EHMT2 guardado como gsea_hallmarks_EHMT2.png\n")

 # Crear enrichment plots para los términos más significativos
num_plots <- min(2, nrow(gsea_ehmt2@result))  # Máximo 2 plots

for(i in 1:num_plots) {
  pathway_id <- gsea_ehmt2@result$ID[i]
  pathway_desc <- gsea_ehmt2@result$Description[i]
  pathway_nes <- round(gsea_ehmt2@result$NES[i], 2)
  pathway_fdr <- round(gsea_ehmt2@result$p.adjust[i], 3)
  
  enrich_plot <- gseaplot2(gsea_ehmt2, geneSetID = pathway_id, 
                          title = paste0("EHMT2 - ", pathway_desc, 
                                       "\nNES: ", pathway_nes, ", FDR: ", pathway_fdr),
                          color = "red", base_size = 11)
  
  print(enrich_plot)
  
  # Guardar cada plot con nombre diferente
  filename <- paste0("../Gráficos e imágenes/gsea_enrichment_EHMT2_hallmarks_top", i, ".png")
  ggsave(filename, plot = enrich_plot, width = 12, height = 8, dpi = 300)
  cat(paste("GSEA enrichment plot", i, "de EHMT2 guardado como gsea_enrichment_EHMT2_hallmarks_top", i, ".png\n"))
}

} else {
  cat("No se encontraron hallmarks enriquecidos en GSEA para EHMT2\n")
}

cat("\n=== Análisis GSEA para EHMT2 completado ===\n")