######
# Script : Analisis de expresion diferencial
# Author: Emiliano Ferro Rodriguez, Jorge Alfredo Suazo Victoria, Sofia Gamiño Estrada
# Date: 02/04/2025
# Description: El siguiente script nos permite realiza el Analisis de expresion Diferencial
# a partir de los datos provenientes del alineamiento Star-Salmon del pipeline de nf-core a R,
# Primero Obtener el Resultado del Pipeline: nf-core/rna-seq/3.12.0
# Usage: Correr las lineas en un nodo de prueba en el cluster.
# Arguments:
#   - Input: Cargar el objeto DESeqDataSet previamente guardado
#   - Output: DEG
#######

rm(list=ls())  # Limpiar el entorno de trabajo

# Cargar paquetes necesarios
library(DESeq2)   # Para análisis de expresión diferencial
library(dplyr)    # Para manipulación de datos
tlibrary(ggplot2) # Para visualización de datos
library(ggExtra)  # Para agregar gráficos marginales
tlibrary(DOSE)    # Para análisis de enriquecimiento funcional

# Cargar el objeto DESeqDataSet previamente guardado
load("/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo1/DEG_output/deseq2_qc/deseq2.dds.RData")

# Definir directorios de salida
outdir <- "/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo1/results/"  # Directorio para guardar resultados
figdir <- '/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo1/results/figures/'  # Directorio para guardar figuras

# Cargar la metadata actualizada
metadata <- read.table("/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo1/datas/metadata.tsv", header=TRUE, sep="\t", row.names=1)

# Convertir las columnas en factores para el análisis
tmetadata$Condition <- as.factor(metadata$Condition)
metadata$Type <- as.factor(metadata$Type)
metadata$Individuals <- as.factor(metadata$Individuals)

# Extraer nombres de muestras y factores de escala del objeto DESeqDataSet
sample <- colData(dds)$sample
sizeFactor <- colData(dds)$sizeFactor

# Filtrar metadata para que coincida con las muestras en el objeto DESeqDataSet
metadata <- metadata[sample,]

# Actualizar colData del objeto DESeqDataSet con la nueva metadata
colData(dds) <- DataFrame(metadata)

# Definir el diseño experimental para el análisis de expresión diferencial
design(dds) <- ~Condition * Type

# Ejecutar el análisis de expresión diferencial con DESeq2
dds <- DESeq(dds)

# Transformación de variancia estabilizada (VST) para normalización
tvsdata <- vst(dds, blind = FALSE) 

# Guardar el objeto actualizado con los resultados de DESeq2
save(metadata, dds, file = paste0(outdir, 'dds_TypesnCondition.RData'))

# Obtener matriz de expresión normalizada
vs_mat <- assay(vsdata)

# Realizar análisis de componentes principales (PCA) sobre la matriz normalizada
pca <- prcomp(t(vs_mat))

# Convertir los resultados de PCA en un dataframe
pca_df <- as.data.frame(pca$x)
pca_df$Condition <- colData(dds)$Condition  # Agregar metadata al PCA
pca_df$Type <- colData(dds)$Type

# Guardar la figura del PCA en un archivo PNG
png(file = paste0(figdir, "PCA_vsdata.png"), width = 3000, height = 2000, res = 300)

# Crear gráfico PCA con ggplot2
p1 <- ggplot(pca_df, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 3, aes(shape = Type)) +  # Agregar puntos con diferente forma según el tipo
  theme_minimal() +  # Aplicar tema minimalista
  labs(title = "PCA de datos VST", x = "PC1", y = "PC2") +  # Etiquetas del gráfico
  theme(legend.position = "right")  # Posición de la leyenda

# Agregar densidad en los ejes del gráfico
p1 <- ggMarginal(p1, type = "density", groupColour = TRUE, groupFill = TRUE)

# Mostrar gráfico y cerrar dispositivo gráfico
print(p1)
dev.off()

# Ver los nombres de los resultados generados
resultsNames(dds)

# Guardar los datos transformados por VST
save(metadata, vsdata, file = paste0(outdir, 'vst_TypesnCondition.RData'))

# Comparaciones de expresión diferencial y guardado de resultados en archivos CSV

# Comparación: Lesional vs Healthy
res_LesvsHea <- results(dds, name = "Condition_Lesional_vs_Healthy")
summary(res_LesvsHea)
write.csv(res_LesvsHea, file=paste0(outdir, 'DE_LesionalvsHealthy.csv'))

# Comparación: Nonlesional vs Healthy
res_NlesvsHea <- results(dds, name = "Condition_Nonlesional_vs_Healthy")
summary(res_NlesvsHea)
write.csv(res_NlesvsHea, file=paste0(outdir, 'DE_NonLesionalvsHealthy.csv'))

# Comparación: GC vs Control
res_GCvsCtrl <- results(dds, name = "Type_GC_vs_Control")
summary(res_GCvsCtrl)
write.csv(res_GCvsCtrl, file=paste0(outdir, 'DE_GCvsControl.csv'))

# Comparación: Efecto de la condición Lesional en GC
res_Lesional.GC <- results(dds, name = "ConditionLesional.TypeGC")
summary(res_Lesional.GC)
write.csv(res_Lesional.GC, file=paste0(outdir, 'DE_LesionalGC.csv'))

# Comparación: Efecto de la condición Nonlesional en GC
res_Nonlesional.GC <- results(dds, name = "ConditionNonlesional.TypeGC")
summary(res_Nonlesional.GC)
write.csv(res_Nonlesional.GC, file=paste0(outdir, 'DE_NonlesionalGC.csv'))
