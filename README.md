# Expresión Diferencial y Análisis Funcional de Pacientes con Psoriasis en Tratamiento con Glucocorticoides

- Emiliano Ferro Rodriguez
- Jorge Alfredo Suazo Victoria
- Sofia Gamiño Estrada

### Fecha de elaboracion: 01-04-2025

## Material Utilizado

|Descripcion|Informacion|
|-----------|-----------|
|Bioproject|[PRJNA494527](https://www.ebi.ac.uk/ena/browser/view/PRJNA494527)|
|Especie|*Homo Sapiens*|
|Tipo de biblioteca|sinlge-end|
|Metodo|RNA-Total|
|Numero de transcriptomas|34|
|Numero de replicas| 17 Replicas biologicas, con una replica tecnica por cada una ( Control y Firma génica inducida por glucocorticoides en la piel humana)|
|Secuenciador Empleado|Ilumina NextSeq500|
|Profundidad de secuenciacion de cada transcriptoma|12M a 40M|
|Tamaño de las lecturas|75bp|
|Articulo Cientifico|Sarkar MK, Kaplan N, Tsoi LC, Xing X et al. Endogenous Glucocorticoid Deficiency in Psoriasis Promotes Inflammation and Abnormal Differentiation. J Invest Dermatol 2017 Jul;137(7):1474-1483. PMID: 28259685 Los datos se pueden descargar desde NCBI o usando ENA.|




- [nf-core/rna-seq/3.12.0](https://nf-co.re/rnaseq/3.12.0/)
- [nextflow/23.10.0](https://www.nextflow.io/)
- [singularity/3.7.0](https://docs.sylabs.io/guides/3.7/user-guide/)
- [Supercomputo LAVIS (Cluster DNA)](https://lavis.unam.mx/servicios/)
- [Rstudio/Posit](https://posit.co/download/rstudio-desktop/)
### Paquetes de R utilizados

```
library(gprofiler2)
library(enrichplot)
library(clusterProfiler)
library(ggplot2)
library(tidyverse)
library(pheatmap)
library(DESeq2)
library(dplyr)
library(ggExtra)      
library(DOSE)
library(pathview)
library(BiocManager)
library(org.Hs.eg.db) 
library(AnnotationDbi)
library(tree)
```
# Pipeline
![](images/Diagrama_de_flujo_RNASeq.png)


## Elección de datos

Para elegir los datos primero se piensa en el organismo con el que se va a tratar en este caso escogimos al *Homo sapiens* (Humano)

Se seleccionó el artículo "[Endogenous Glucocorticoid Deficiency in Psoriasis Promotes Inflammation and Abnormal Differentiation.](https://europepmc.org/article/PMC/5545780)" en el cual se habla del desconocimiento que se tiene de la biologia molecular sobre la terapia con glucocorticoides en pacientes con psoriasis.

## nf-core/rnaseq

Utilizamos el pipeline nf-core/rnaseq basado en Nextflow, una herramienta reconocida por su capacidad de ejecutar flujos de trabajo escalables y reproducibles. Este enfoque es ideal para análisis complejos, ya que reduce la dependencia de configuraciones manuales y facilita la trazabilidad del proceso.

El pipeline requiere tan solo un archivo CSV (samplesheet.csv) que especifica los nombres de las muestras y la ubicación de los archivos FASTQ, simplificando la integración de datos experimentales. Además, se incorpora el genoma de referencia hg38 de Ensembl, lo que proporciona una base robusta para el alineamiento y análisis transcriptómico.

