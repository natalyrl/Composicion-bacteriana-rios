#1.1 Descargar las secuencias. En este estudio se utilizaron seis muestras por río:
#Para el río Bogotá: https://www.ebi.ac.uk/ena/browser/view/PRJEB22915
#Para el río Pinheiros: https://zenodo.org/records/1172783

#1.2 Instalar y cargar DADA2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.16")
library(BiocManager)
library(dada2)
library(Rcpp)

#1.3 Pipeline DADA2 (análisis de calidad y corte) para las secuencias del río Bogotá
rutaCol <- "/Users/Ariana/Downloads/Colombia"
list.files(rutaCol)

fnFs_Col <- sort(list.files(rutaCol, pattern="_1.fastq", full.names = TRUE))
fnRs_Col <- sort(list.files(rutaCol, pattern="_2.fastq", full.names = TRUE))

sample.names_Col <- sapply(strsplit(basename(fnFs_Col), "_"), `[`, 1)

plotQualityProfile(fnFs_Col[1:2])
plotQualityProfile(fnFs_Col[3:4])
plotQualityProfile(fnFs_Col[5:6])

plotQualityProfile(fnRs_Col[1:2])
plotQualityProfile(fnRs_Col[3:4])
plotQualityProfile(fnRs_Col[5:6])

filtFs_Col <- file.path(rutaCol, "filtrados_Col", paste0(sample.names_Col, "_1_filt.fastq.gz"))
filtRs_Col <- file.path(rutaCol, "filtrados_Col", paste0(sample.names_Col, "_2_filt.fastq.gz"))
names(filtFs_Col) <- sample.names_Col
names(filtRs_Col) <- sample.names_Col

out_Col <- filterAndTrim(fnFs_Col, filtFs_Col, fnRs_Col, filtRs_Col, truncLen=c(150,150),
                       maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                       compress=TRUE, multithread=FALSE)
head(out_Col)

#1.4 Pipeline DADA2 (análisis de calidad y corte) para las secuencias del río Pinheiros
rutaBrasil<- "/Users/Ariana/Downloads/Brasil"
list.files(rutaBrasil)

fnFs_Brasil <- sort(list.files(rutaBrasil, pattern="_1.fastq", full.names = TRUE))
fnRs_Brasil <- sort(list.files(rutaBrasil, pattern="_2.fastq", full.names = TRUE))

sample.names_Brasil <- sapply(strsplit(basename(fnFs_Brasil), "_"), `[`, 1)

plotQualityProfile(fnFs_Brasil[1:2])
plotQualityProfile(fnFs_Brasil[3:4])
plotQualityProfile(fnFs_Brasil[5:6])

plotQualityProfile(fnRs_Brasil[1:2])
plotQualityProfile(fnRs_Brasil[3:4])
plotQualityProfile(fnRs_Brasil[5:6])

filtFs_Brasil <- file.path(rutaBrasil, "filtrados_Brasil", paste0(sample.names_Brasil, "_1_filt.fastq.gz"))
filtRs_Brasil <- file.path(rutaBrasil, "filtrados_Brasil", paste0(sample.names_Brasil, "_2_filt.fastq.gz"))
names(filtFs_Brasil) <- sample.names_Brasil
names(filtRs_Brasil) <- sample.names_Brasil

out_Brasil <- filterAndTrim(fnFs_Brasil, filtFs_Brasil, fnRs_Brasil, filtRs_Brasil, truncLen=c(0,250),
                       maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                       compress=TRUE, multithread=FALSE)
head(out_Brasil)

#Las carpetas "filtrados_Col" y "filtrados_Brasil" serán el input en el paso 2. Ambas carpetas se pasaron del computador al clúster usando FileZilla
