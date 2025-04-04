######
# Script : Visualizacion grafica de los resultados de DEG
# Author: Emiliano Ferro Rodriguez, Jorge Alfredo Suazo Victoria, Sofia Gamiño Estrada
# Date: 02/04/2025
# Description: El siguiente script nos permite realiza el Analisis de Terminos GO
# a partir de los datos provenientes del Analisis de DEG
# Primero correr el script "DEG_analysis.R"
# Usage: Correr las lineas en un nodo de prueba en el cluster.
# Arguments:
#   - Input: 
#       - dds_TypesnCondition.RData (dds), 
#       - vst_TypesnCondition.RData (vsdata) 
#       - archivos de salida de DEG en formato CSV (DE_GCvsControl.csv) 
#   - Output: Volcano plot y Heatmap
#######


# Limpiar el entorno de trabajo
rm(list = ls())

# Cargar librerías necesarias
library(dplyr)
library(pheatmap)
library(ggplot2)
library(DESeq2) 

# Definir directorio de figuras
figdir <- '/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo1/results/Analisis_diferencial/figures/'

# Cargar datos procesados previamente en DEG_analysis.R
load("/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo1/results/Analisis_diferencial/dds_TypesnCondition.RData")  # Variable "dds"
load("/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo1/results/Analisis_diferencial/vst_TypesnCondition.RData")  # Variable "vsdata"

# Cargar resultados de análisis diferencial
DE_GCvsControl <- read.csv("/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo1/results/Analisis_diferencial/DE_GCvsControl.csv")

df <- as.data.frame(DE_GCvsControl)

# Clasificar genes según expresión diferencial
# Se consideran "Up-regulated" los genes con log2FoldChange >= 2 y padj < 0.05
# Se consideran "Down-regulated" los genes con log2FoldChange <= -2 y padj < 0.05
# Los demás se clasifican como "Unchanged"
df <- df %>% 
  mutate(Expression = case_when(log2FoldChange >= 2 & padj < 0.05 ~ "Up-regulated",
                                log2FoldChange <= -2 & padj < 0.05 ~ "Down-regulated",
                                TRUE ~ "Unchanged"))

# Generar Volcano Plot
ggplot(df, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(color = Expression), size = 0.7) +
  labs(title = "Tratamiento vs Control") +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"p-adj")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  ggsave(filename = paste0(figdir, "VolcanoPlotGCvsControl.png"), width = 10, height = 6, dpi = 300)

# Seleccionar los 20 genes más significativos para el heatmap
topGenes <- head(order(DE_GCvsControl$padj), 20)

# Crear dataframe con información de muestras
sample_info <- data.frame(
  SampleType = colData(dds)$Condition,  # Categorías: Lesional, Nonlesional, Healthy
  Group = colData(dds)$Type  # Categorías: GC o Control
)
rownames(sample_info) <- colnames(dds)  # Asegurar coincidencia con la matriz de expresión

# Definir colores para anotaciones en el heatmap
annotation_colors <- list(
  SampleType = c(Lesional = "red", Nonlesional = "orange", Healthy = "green"),
  Group = c(GC = "red", Control = "blue")
)

# Convertir IDs de Ensembl a símbolos de genes
library(org.Hs.eg.db)
library(AnnotationDbi)

ensembl_ids <- rownames(vsdata)
gene_symbols <- mapIds(org.Hs.eg.db, keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
rownames(vsdata) <- gene_symbols  # Reemplazar nombres en la matriz

# Generar heatmap de los genes más significativos
png(file = paste0(figdir, "Heatmap_vsd_topgenes.png"), width = 3000, height = 2000, res = 300)
pheatmap(
  assay(vsdata)[topGenes,], 
  cluster_rows = FALSE, 
  cluster_cols = TRUE, 
  show_rownames = TRUE,
  annotation_col = sample_info,
  annotation_colors = annotation_colors
)
dev.off()

# ---- Generación de heatmap de coeficientes log2FoldChange ----
betas <- coef(dds)
colnames(betas)

# Crear matriz con coeficientes de los genes más significativos
mat <- betas[topGenes, -1] 

# Aplicar umbral de log2FoldChange para limitar valores extremos
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr

# Guardar el heatmap en archivo PNG
png(file = paste0(figdir, "Heatmap_log2FoldChange_topgenes.png"), width = 3000, height = 2000, res = 300)
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101), cluster_col=FALSE)
dev.off()
