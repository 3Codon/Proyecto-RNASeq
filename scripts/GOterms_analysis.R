#######
# Script : Analisis de terminos GO y KEGG
# Author: Emiliano Ferro Rodriguez, Jorge Alfredo Suazo Victoria, Sofia Gamiño Estrada
# Date: 02/04/2024
# Description: El siguiente script nos permite realiza determinacion funcional
# a partir de los datos provenientes del alineamiento Star-Salmon del pipeline de nf-core a R,
# Primero Obtener el Resultado del Pipeline: nf-core/rna-seq/3.12.0
# Usage: Correr las lineas en un nodo de prueba en el cluster.
# Arguments:
#   - Input: metadata.csv, cargar el objeto DESeqDataSet previamente guardado
#   - Output: Matriz de cuentas (CSV y RData)
#######

# --- Cargar paquetes necesarios ---
library(gprofiler2)
library(enrichplot)
library(DOSE)
library(clusterProfiler)
library(ggplot2)
library(tidyverse)
library(dplyr)

# --- Definir directorios de entrada y salida ---
indir <- "/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo1/results/Analisis_diferencial/"
outdir <- "/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo1/results/Analisis_diferencial/"
figdir <- '/mnt/atgc-d1/bioinfoII/rnaseq/BioProject_2025/Equipo1/results/Analisis_diferencial/figures/'

# --- Definir bases de datos para análisis de términos GO ---
sources_db <- c("GO:BP", "KEGG", "REAC", "TF", "MIRNA", "CORUM", "HP", "HPA", "WP")

# --- Seleccionar archivos CSV disponibles ---
files <- dir(indir, pattern = "^DE_(.+)\\.csv$") 

# --- Cargar y procesar el primer archivo ---
plot_name <- gsub("^DE_(.+)\\.csv$", "\\1", files[1]) # Extraer nombre del archivo
df <- read.csv(file = paste0(indir, files[1]), row.names = 'X') # Cargar datos

# --- Clasificar genes según expresión diferencial ---
abslogFC <- 2 # Umbral para considerar cambios significativos
df <- df %>% 
  dplyr::mutate(Expression = case_when(
    log2FoldChange >= abslogFC & padj < 0.05 ~ "Up-regulated",
    log2FoldChange <= -abslogFC & padj < 0.05 ~ "Down-regulated",
    TRUE ~ "Unchanged")) 

# --- Obtener nombres de los genes diferencialmente expresados ---
# Genes up-regulated
up_genes <- df %>% filter(Expression == 'Up-regulated') %>% arrange(padj, desc(abs(log2FoldChange)))
up_genes <- rownames(up_genes) 

# Convertir identificadores de ENSEMBL a símbolos de genes
library(org.Hs.eg.db)
library(AnnotationDbi)
gene_symbols <- mapIds(org.Hs.eg.db, keys = up_genes, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
up_genes <- gene_symbols

# Genes down-regulated
down_genes <- df %>% filter(Expression == 'Down-regulated') %>% arrange(padj, desc(abs(log2FoldChange)))
down_genes <- rownames(down_genes) 

gene_symbols <- mapIds(org.Hs.eg.db, keys = down_genes, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
down_genes <- gene_symbols

# --- Análisis de enriquecimiento funcional con g:Profiler ---
multi_gp <- gost(list("Upregulated" = up_genes, "Downregulated" = down_genes), 
                 correction_method = "fdr", user_threshold = 0.05,
                 multi_query = FALSE, ordered_query = TRUE, 
                 sources = sources_db, 
                 evcodes = TRUE,  # Incluir genes asociados a cada término
                 organism = 'hsapiens') 

# --- Definir colores para cada categoría ---
Category_colors <- data.frame(
  category = c("GO:BP", "GO:CC", "GO:MF", "KEGG", 'REAC', 'TF', 'MIRNA', 'HPA', 'CORUM', 'HP', 'WP'), 
  label = c('Biological Process', 'Cellular Component', 'Molecular Function',  "KEGG",
            'REAC', 'TF', 'MIRNA', 'HPA', 'CORUM', 'HP', 'WP'),
  colors =  c('#FF9900', '#109618','#DC3912', '#DD4477', '#3366CC','#5574A6', '#22AA99', '#6633CC', '#66AA00', '#990099', '#0099C6'))

# --- Generar y guardar gráfico Manhattan ---
gostp1 <- gostplot(multi_gp, interactive = FALSE)
ggsave(paste0(figdir, "ManhattanGO_", plot_name, ".png"), plot = gostp1, dpi = 300)

# --- Convertir resultados a dataframe ---
gost_query <- as.data.frame(multi_gp$result)
bar_data <- data.frame("term" = as.factor(gost_query$term_name), "condition" = gost_query$query, 
                       "count" = gost_query$term_size, "p.adjust" = gost_query$p_value, 
                       'category' = as.factor(gost_query$source), "go_id" = as.factor(gost_query$term_id),
                       'geneNames' = gost_query$intersection)

# --- Procesar genes down-regulated ---
bar_data_down <- subset(bar_data, condition == 'Downregulated')
bar_data_down <- head(bar_data_down[order(bar_data_down$p.adjust),], 40)
bar_data_down$p.val <- round(-log10(bar_data_down$p.adjust), 2)
bar_data_down$num <- seq(1:nrow(bar_data_down)) 
bar_data_down <- left_join(bar_data_down, Category_colors, by= "category")
save(bar_data_down, file = paste0(outdir, "DOWN_GO_", plot_name, ".RData"))

# --- Graficar términos enriquecidos para genes down-regulated ---
g.down <- ggplot(bar_data_down, aes(p.val, reorder(term, -num), fill = category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = p.val), color = "black", hjust = 0, size = 2.2) +
  labs(x = "-log10(p-value)", y = NULL) +
  scale_fill_manual(name='Category', labels = unique(bar_data_down$label), values = unique(bar_data_down$colors)) +
  theme_classic()
ggsave(paste0(figdir, "barplotDOWN_GO_", plot_name, ".png"), plot = g.down, dpi = 600, width = 10, height = 5)

# --- Procesar genes up-regulated ---
bar_data_up <- subset(bar_data, condition == 'Upregulated')
bar_data_up <- head(bar_data_up[order(bar_data_up$p.adjust),], 40)
bar_data_up$p.val <- round(-log10(bar_data_up$p.adjust), 2)
bar_data_up$num <- seq(1:nrow(bar_data_up)) 
bar_data_up <- left_join(bar_data_up, Category_colors, by= "category")
save(bar_data_up, file = paste0(outdir, "UP_GO_", plot_name, ".RData"))

# --- Graficar términos enriquecidos para genes up-regulated ---
g.up <- ggplot(bar_data_up, aes(p.val, reorder(term, -num), fill = category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = p.val), color = "black", hjust = 0, size = 2.2) +
  labs(x = "-log10(p-value)", y = NULL) +
  scale_fill_manual(name='Category', labels = unique(bar_data_up$label), values = unique(bar_data_up$colors)) +
  theme_classic()
ggsave(paste0(figdir, "barplotUP_GO_", plot_name, ".png"), plot = g.up, dpi = 600, width = 10, height = 5)

kegg_gene_list <- df$log2FoldChange
names(kegg_gene_list) <- rownames(df)

library(org.Hs.eg.db)
library(AnnotationDbi)

entrez_ids <- mapIds(org.Hs.eg.db, keys = names(kegg_gene_list),
                     column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")

kegg_gene_list <- kegg_gene_list[!is.na(entrez_ids)]
names(kegg_gene_list) <- entrez_ids[!is.na(entrez_ids)]

library(pathview)

kegg_organism <- "hsa"
pathway_id <- "hsa04657"

# Generar el gráfico nativo KEGG (PNG)
pv_out <- pathview(gene.data = kegg_gene_list,
                   pathway.id = pathway_id,
                   species = kegg_organism)

# Generar una versión alternativa en PDF (opcional)
pv_out_pdf <- pathview(gene.data = kegg_gene_list,
                       pathway.id = pathway_id,
                       species = kegg_organism,
                       kegg.native = FALSE)
