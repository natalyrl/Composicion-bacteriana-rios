#El input son las dos tablas de conteos por géneros

#3.1 Instalar y cargar paquetes
if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")
library(phyloseq)
BiocManager::install("microbiome")
library(microbiome)
BiocManager::install("ANCOMBC")
library(ANCOMBC)
library(dplyr)
library(readr)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(readxl)
library(ggsci)
library(reshape2)
library(RColorBrewer)
library(vegan)

#3.2 Combinar las dos tablas por géneros
colombia <- read_csv("/Users/Ariana/Downloads/tabla_conteos_por_genero_Col.csv")
brasil   <- read_csv("/Users/Ariana/Downloads/tabla_conteos_por_genero_Bra.csv")

#Unir por 'Genus', manteniendo todos los géneros
colnames(brasil)
colnames(colombia)
unida <- full_join(colombia, brasil, by = "Genus")

#Reemplazar NAs con ceros
unida[is.na(unida)] <- 0

#Guardar tabla combinada
write_csv(unida, "tabla_conteos_colombia_brasil.csv")

#3.3 Crear el objeto phyloseq
otu<- read_excel("/Users/Ariana/Downloads/tabla_colombia_brasil.xlsx", sheet = "OTU")
tax<- read_excel("/Users/Ariana/Downloads/tabla_colombia_brasil.xlsx", sheet = "Taxa")
samples <- read_excel("/Users/Ariana/Downloads/datos_rios.xlsx")

rownames(otu) <- otu$ASV
otu <- otu %>% select(-ASV)

rownames(tax) <- tax$ASV
tax<- tax %>% select(-ASV)

rownames(samples) <- samples$Muestra
samples <- samples %>% select(-Muestra)

otu <- as.matrix(otu)
tax <- as.matrix(tax)

#Transformación a objetos phyloseq
OTUGEN = otu_table(otu, taxa_are_rows = TRUE)
TAXGEN = tax_table(tax)
samples = sample_data(samples)
samples_name <- c("Uno_B","Dos_B","Tres_B","Cuatro_B","Cinco_B","Seis_B","Uno_P","Dos_P",
                  "Tres_P", "Cuatro_P", "Cinco_P", "Seis_P")
sample_names(samples)=samples_name
Generos <- phyloseq(OTUGEN, TAXGEN, samples)
Generos

#verificar
sample_names(Generos)
rank_names(Generos)
sample_variables(Generos)

#3.4 Escoger el top 10 de géneros con más conteos
top10_taxa <- names(sort(taxa_sums(Generos), decreasing = TRUE))[1:10]
Generos_top10 <- prune_taxa(top10_taxa, Generos)
Generos_rel <- transform_sample_counts(Generos_top10, function(x) x / sum(x))

#3.5 Graficar las abundancias relativas del top 10
Abundancias<-plot_bar(Generos_rel, fill = "Genero") +
  facet_grid(~Rio, scales = "free_x") +
  scale_fill_brewer(palette = "Spectral")+
  theme_bw() +
  labs(title = "a",
       x = "Muestra",
       y = "Abundancia relativa")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
Abundancias

#3.6 Diversidad alfa: índices de Shannon y Simpson
samples$Muestra <- rownames(samples)
sample_data(Generos)$Muestra <- samples$Muestra

Alfa<-plot_richness(Generos, x="Muestra", color="Rio", measures=c("Shannon", "Simpson"))+
  geom_point(size=4, alpha=0.7)+
  theme_bw()+
  labs(x = "Muestra",
       y = "Medida de Diversidad Alfa",
       fill = "Río",
       title = "b")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
Alfa

#3.7 Diversidad beta: índice de Bray-Curtis
Generos_ab <- transform_sample_counts(Generos, function(x) x / sum(x))
bray_dist <- distance(Generos_ab, method = "bray")
PCoA <- ordinate(Generos_ab, method = "PCoA", distance = bray_dist)

Beta <- plot_ordination(Generos_ab, PCoA, color = "Rio") +
  geom_point(size = 3) +
  stat_ellipse(type = "t", linetype = "solid") +
  theme_minimal()+
  labs(fill="Río",
       title = "c")
Beta

#3.8 Pruebas estadísticas

##3.8.1 Normalidad
unida_larga <- unida %>%
  pivot_longer(cols = -Genus, names_to = "Muestra", values_to = "Abundancia")
normalidad <- unida_larga %>%
  group_by(Genus) %>%
  summarise(
    p_shapiro = shapiro.test(Abundancia)$p.value
  ) %>%
  mutate(
    normal = ifelse(p_shapiro >= 0.05, "Sí", "No")
  ) %>%
  arrange(p_shapiro)
normalidad
#No son normales

##3.8.2 ANCOMBC
levels(sample_data(Generos)$Rio)
unique(sample_data(Generos)$Rio)
ANCOMBC_rios<- ancombc2(
  data        = Generos,           
  assay_name  = "counts",     
  fix_formula = "Rio",
  group       = "Rio",
  p_adj_method= "BH",
  prv_cut     = 0.10,
  lib_cut     = 0,
  struc_zero  = TRUE,
  neg_lb      = TRUE,
  alpha       = 0.05,
  pairwise    = TRUE,
  global      = TRUE,
  verbose     = TRUE)
head(ANCOMBC_rios$res)
View(ANCOMBC_rios$res) #Solo significativo en el río Pinheiros
sigRios <- subset(ANCOMBC_rios$res, diff_RioRío_Pinheiros)
View(sigRios)
View(Generos@tax_table)

##3.8.3 Comparar índices de diversidad alfa entre grupos
wilcox_shannon <- wilcox.test(Shannon ~ Rio, data = diversidad)
wilcox_simpson <- wilcox.test(Simpson ~ Rio, data = diversidad)
wilcox_shannon
wilcox_simpson

##3.8.4 Comparar índice de Bray-Curtis entre grupos
beta <- as(sample_data(Generos_ab), "data.frame")
permanova <- adonis2(bray_dist ~ Rio, data = beta, permutations = 999)
print(permanova)
