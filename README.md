# Comparación de la composición bacteriana de dos ríos urbanos: Río Bogotá (Colombia) y Río Pinheiros (Brasil)

**Autores:** Ariana Delgadillo Pérez*, Nataly Rodríguez Lugo y Sergio A. Sánchez León 

Escuela de Ciencias e Ingeniería, Universidad del Rosario

*Email: ariana.delgadillo@urosario.edu.co

Este repositorio contiene los scripts, datos y documentación del análisis de comunidades bacterianas en dos ríos urbanos (Pinheiros y Bogotá), basado en secuencias
del gen 16S rRNA.

## Contenido del repositorio

- Scripts de: (1) DADA2 para evaluar la calidad de las secuencias y cortarlas, (2) el resto del pipeline de DADA2 en el clúster para obtener dos tablas de conteos luego de la asignación taxonómica,
y (3) los análisis de diversidad alfa, beta y abundancia diferencial en R. Disponibles en [scripts](https://github.com/natalyrl/Composicion-bacteriana-rios/tree/main/scripts).
- Gráficos de los resultados de los análisis. Disponibles en [imágenes](https://github.com/natalyrl/Composicion-bacteriana-rios/tree/main/imagenes).

## Descripción del análisis

Los datos provienen de bases de datos públicas:

- Río Bogotá: [ENA - PRJEB22915](https://www.ebi.ac.uk/ena/browser/view/PRJEB22915)
- Río Pinheiros: [Zenodo Dataset](https://zenodo.org/records/1172783)

 *Nota*: Se tomaron solo seis pares de secuencias para cada río.

Los análisis incluyen:

1. **Calidad de secuencias y corte** con paquetes `BiocManager` y `DADA2`.
2. **Procesamiento de secuencias** con DADA2 desde el clúster.
3. **Análisis de diversidad alfa, beta y comparación diferencial** con paquetes `phyloseq`, `ANCOMBC2` y `vegan`:
   - Índices de Shannon y Simpson
   - Evaluación de géneros diferencialmente abundantes
   - Wilcoxon
   - PCoA basado en el índice de Bray-Curtis
   - PERMANOVA
     
  ![resultados](https://github.com/natalyrl/Composicion-bacteriana-rios/blob/main/imagenes/Resultados.png)
