# Libraries
library(ggplot2)
library(DESeq2)
library(magrittr)      
library(tidyverse)
library(openxlsx)
library(dplyr)
library(clusterProfiler)
library(gprofiler2)
library(openxlsx)
library(writexlsx)
library(org.Mm.eg.db)
library(WGCNA)

# Count Data
count_data_SM539 <- read.table("SM539_countTable.txt", header = TRUE, row.names = 1)
count_data_SM539 <- as.matrix(count_data_SM539)
mode(count_data_SM539) <- "numeric"
counts.filtered <- count_data_SM539[rowSums(count_data_SM539) > 0, ]

dim(counts.filtered)

head(counts.filtered)

colnames(counts.filtered) <- sub("^X", "", colnames(counts.filtered))

colnames(counts.filtered)

# SampleInfo

sampleinfo=read.table("SM539_sampleinfo.txt", sep="\t", header=T, row=1)
dim(sampleinfo)
head(sampleinfo)


# OBJETO DESEQ2 De los datos SM539 Completos Condition * Month

sampleinfo$Condition <- as.factor(sampleinfo$Condition)
sampleinfo$Months <- as.factor(sampleinfo$Months)
dds <- DESeqDataSetFromMatrix(countData = counts.filtered ,
                                    colData = sampleinfo ,
                                    design = ~ Condition * Months)
dds <- DESeq(dds)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

normalized_counts <- normalized_counts[rowSums(normalized_counts != 0) > 0, ]
head(normalized_counts)
transposed_data <- t(normalized_counts)

dim(transposed_data)

# Filtrado de genes con poca varianza

# Calcular la varianza de cada gen
gene_variances <- apply(transposed_data, 2, var)

# Filtrar genes con varianza mayor a un umbral (ajusta según necesidad)
var_threshold <- quantile(gene_variances, 0.5) # Usamos el percentil 50 como umbral
filtered_data <- transposed_data[, gene_variances > var_threshold]

# Verificar dimensiones después del filtrado
dim(filtered_data)

# Verificar calidad de las muestras y los genes
gsg <- goodSamplesGenes(filtered_data, verbose = 3)

# Si hay problemas, eliminar genes o muestras defectuosas
if (!gsg$allOK) {
  filtered_data <- filtered_data[gsg$goodSamples, gsg$goodGenes]
}

# Verificar nuevamente las dimensiones
dim(filtered_data)

# Agrupación Clustering

# Verificar calidad de las muestras y los genes
gsg <- goodSamplesGenes(filtered_data, verbose = 3)

# Si hay problemas, eliminar genes o muestras defectuosas
if (!gsg$allOK) {
  filtered_data <- filtered_data[gsg$goodSamples, gsg$goodGenes]
}

# Verificar nuevamente las dimensiones
dim(filtered_data)

sampleinfo$Condition <- as.factor(sampleinfo$Condition)
sampleinfo$Months <- as.factor(sampleinfo$Months)
sampleinfo$Sex <- as.factor(sampleinfo$Sex)

# Convertir las variables de factores en números
sampleinfo$Condition_num <- as.numeric(factor(sampleinfo$Condition))
sampleinfo$Sex_num <- as.numeric(factor(sampleinfo$Sex))
sampleinfo$Months_num <- as.numeric(factor(sampleinfo$Months))

# Crear el dataframe numérico con las columnas relevantes
traitData_numeric <- sampleinfo[, c("Condition_num", "Sex_num", "Months_num")]

# Usar colorRampPalette para generar un gradiente más amplio
myColorPalette <- colorRampPalette(c("pink", "white", "blue"))

# Asignar los colores a los datos
traitColors <- numbers2colors(traitData_numeric, 
                              colors = myColorPalette(100),  # Usamos 100 colores en el gradiente
                              signed = FALSE)

## Dendograma y Trait HeatMap
### Los rosas son:
#### - Control
#### - Hembras
#### - 15 meses

### Los azules son:

#### - Repro
#### - Machos
#### - 7 meses

# Graficar el dendrograma de las muestras con los colores de las características fenotípicas
plotDendroAndColors(sampleTree, traitColors, 
                    groupLabels = colnames(traitData_numeric), 
                    main = "Sample dendrogram and trait heatmap")

## Construcción de la red WGCNA
powers <- c(1:20)  
sft <- pickSoftThreshold(filtered_data, powerVector = powers, verbose = 5)

# Visualizar los resultados
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")
abline(h = 0.9, col = "blue")  # Indicador de R^2 > 0.9

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")

softPower <- 15  
adjacency <- adjacency(filtered_data, power = softPower)

# Calcular la matriz de topología
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM  # Diferencia para obtener la dissimilaridad

## Dendograma Por Módulos
# Agrupamiento jerárquico
geneTree <- hclust(as.dist(dissTOM), method = "average")

# Visualizar el árbol de agrupamiento
plot(geneTree, xlab = "", sub = "", main = "Clustering of Genes", labels = FALSE)

# Detectar módulos de coexpresión
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 0, pamRespectsDendro = FALSE,
                             minClusterSize = 50)

# Convertir a colores para visualizar
moduleColors <- labels2colors(dynamicMods)
table(moduleColors)

# Visualizar los módulos detectados
plotDendroAndColors(geneTree, moduleColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Calcular las eigengenes de los módulos
MEList = moduleEigengenes(expression, colors = moduleColors)
MEs = MEList$eigengenes

# Verifica las primeras filas de las eigengenes
head(MEs)
# Calcular la disimilitud de las eigengenes
MEDiss = 1 - cor(MEs)

# Verifica las primeras filas de la matriz de disimilitud
head(MEDiss)

# Realizar clustering jerárquico sobre las eigengenes
METree = hclust(as.dist(MEDiss), method = "average")

# Graficar el dendrograma del clustering de eigengenes
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.25  # Ajusta este valor según tus datos

# Identificar los módulos después del corte
moduleColors_clust = cutree(METree, h = MEDissThres)

# Verifica los módulos obtenidos
table(moduleColors_clust)

MEs0 <- moduleEigengenes(filtered_data, moduleColors)$eigengenes
MEs0 <- orderMEs(MEs0)
MEs0$treatment = row.names(MEs0)

mME <- MEs0 %>%
  pivot_longer(-treatment) %>%
  filter(!is.na(value)) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = unique(name))  
  )

# Crear el gráfico
ggplot(mME, aes(x = treatment, y = name, fill = value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1, 1)
  ) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "Relaciones Módulo-Rasgo", y = "Módulos", fill = "Correlación")

library(ggplot2)
# Crear un gráfico de burbujas para visualizar la relación entre módulos y rasgos
ggplot(mME, aes(x = treatment, y = name, size = value)) +
  geom_point(alpha = 0.5) +
  scale_size_continuous(range = c(2, 10)) +
  theme_bw() +
  labs(title = "Relación Módulo-Rasgo", x = "Rasgo", y = "Módulo")

## Associating Modules and Phenotypes

# Definir el número de genes y muestras
nGenes = ncol(filtered_data) 
nSamples = nrow(filtered_data)

moduleTraitCor = cor(MEs, traitData_numeric, use = "p")  
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)  

# Mostrar las correlaciones y sus p-valores
textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)

par(mar = c(3, 5, 3, 3))  

# Mostrar las correlaciones y p-valores en un gráfico de calor
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(traitData_numeric),  # Etiquetas de las variables fenotípicas
               yLabels = names(MEs),  # Etiquetas de los módulos
               ySymbols = names(MEs),  # Nombres de los módulos
               colorLabels = FALSE,
               colors = blueWhiteRed(50),  # Colores de la escala
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,  # Ajustar el tamaño de las etiquetas dentro del gráfico de calor
               cex.lab = 0.5,   # Tamaño de las etiquetas de los módulos (yLabels)
               zlim = c(-1, 1),  # Ajuste de la escala de colores entre -1 y 1
               main = paste("Module-trait relationships"))

## Enrichment Analysis

# Condition
magenta_genes <- colnames(filtered_data)[moduleColors == "magenta"]
yellow_genes <- colnames(filtered_data)[moduleColors == "yellow"]
brown_genes <- colnames(filtered_data)[moduleColors == "brown"]

# Month UpRegulated
green_genes <- colnames(filtered_data)[moduleColors == "green"]
turquoise_genes <- colnames(filtered_data)[moduleColors == "turquoise"]
greenyellow_genes <- colnames(filtered_data)[moduleColors == "greenyellow"]

# Month DownRegulated

cyan_genes <- colnames(filtered_data)[moduleColors == "cyan"]
salmon_genes <- colnames(filtered_data)[moduleColors == "salmon"]
tan_genes <- colnames(filtered_data)[moduleColors == "tan"]

length(magenta_genes)
length(yellow_genes)
length(brown_genes)
length(green_genes)
length(turquoise_genes)
length(greenyellow_genes)
length(cyan_genes)
length(salmon_genes)
length(tan_genes)

library(clusterProfiler)
library(org.Mm.eg.db)
# Convertir a ENTREZID para cada módulo
magenta_genes_ids <- bitr(magenta_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
yellow_genes_ids <- bitr(yellow_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
brown_genes_ids <- bitr(brown_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

green_genes_ids <- bitr(green_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
turquoise_genes_ids <- bitr(turquoise_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
greenyellow_genes_ids <- bitr(greenyellow_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

cyan_genes_ids <- bitr(cyan_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
salmon_genes_ids <- bitr(salmon_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
tan_genes_ids <- bitr(tan_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

### Análisis GO

library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(writexl)
library(openxlsx)
# Función para agregar nombres ENS y SYMBOL a los resultados GO
add_ens_and_symbol_to_go <- function(go_results) {
  if (!is.null(go_results)) {
    # Obtener el mapeo ENTREZID -> ENSG y ENTREZID -> SYMBOL
    entrez_to_info <- AnnotationDbi::select(org.Mm.eg.db, 
                                            keys = unique(unlist(strsplit(go_results@result$geneID, "/"))), 
                                            keytype = "ENTREZID", 
                                            columns = c("ENSEMBL", "SYMBOL"))
    
    # Agregar columnas con los nombres ENS y SYMBOL
    go_results@result$ENSEMBL <- sapply(go_results@result$geneID, function(genes) {
      entrez_ids <- unlist(strsplit(genes, "/"))
      ens_ids <- entrez_to_info$ENSEMBL[match(entrez_ids, entrez_to_info$ENTREZID)]
      paste(ens_ids, collapse = ";")
    })
    
    go_results@result$SYMBOL <- sapply(go_results@result$geneID, function(genes) {
      entrez_ids <- unlist(strsplit(genes, "/"))
      symbols <- entrez_to_info$SYMBOL[match(entrez_ids, entrez_to_info$ENTREZID)]
      paste(symbols, collapse = ";")
    })
  }
  return(go_results)
}

# Modificar el análisis para usar la nueva función
go_magenta <- enrichGO(gene = magenta_genes_ids$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
go_magenta <- add_ens_and_symbol_to_go(go_magenta)

go_yellow <- enrichGO(gene = yellow_genes_ids$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
go_yellow <- add_ens_and_symbol_to_go(go_yellow)

go_brown <- enrichGO(gene = brown_genes_ids$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
go_brown <- add_ens_and_symbol_to_go(go_brown)

go_green <- enrichGO(gene = green_genes_ids$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
go_green <- add_ens_and_symbol_to_go(go_green)

go_turquoise <- enrichGO(gene = turquoise_genes_ids$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
go_turquoise <- add_ens_and_symbol_to_go(go_turquoise)

go_greenyellow <- enrichGO(gene = greenyellow_genes_ids$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
go_greenyellow <- add_ens_and_symbol_to_go(go_greenyellow)

go_cyan <- enrichGO(gene = cyan_genes_ids$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
go_cyan <- add_ens_and_symbol_to_go(go_cyan)

go_salmon <- enrichGO(gene = salmon_genes_ids$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
go_salmon <- add_ens_and_symbol_to_go(go_salmon)

go_tan <- enrichGO(gene = tan_genes_ids$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
go_tan <- add_ens_and_symbol_to_go(go_tan)

# Combinar los resultados de enriquecimiento en una sola tabla
enrichment_table <- rbind(
  data.frame(module = "magenta", go_magenta@result),
  data.frame(module = "yellow", go_yellow@result),
  data.frame(module = "brown", go_brown@result),
  data.frame(module = "green", go_green@result),
  data.frame(module = "turquoise", go_turquoise@result),
  data.frame(module = "greenyellow", go_greenyellow@result),
  data.frame(module = "cyan", go_cyan@result),
  data.frame(module = "salmon", go_salmon@result),
  data.frame(module = "tan", go_tan@result)
)

# Guardar la tabla completa como archivo Excel
write.xlsx(enrichment_table, file = "GOEnrichmentTable.xlsx", rowNames = FALSE)

# Filtrar los 5 mejores términos por módulo
best_terms_table <- enrichment_table %>%
  group_by(module) %>%
  slice_max(n = 5, order_by = -p.adjust)

# Guardar los mejores términos en otro archivo Excel
write.xlsx(best_terms_table, file = "Best_GO_Terms.xlsx", rowNames = FALSE)

library(ggplot2)
# Visualización con dotplot para los mejores términos enriquecidos y títulos personalizados
dotplot(go_magenta, showCategory = 10) + ggtitle("Módulo Magenta")
dotplot(go_yellow, showCategory = 10) + ggtitle("Módulo Yellow")
dotplot(go_brown, showCategory = 10) + ggtitle("Módulo Brown")
dotplot(go_green, showCategory = 10) + ggtitle("Módulo Green")
dotplot(go_turquoise, showCategory = 10) + ggtitle("Módulo Turquoise")
dotplot(go_greenyellow, showCategory = 10) + ggtitle("Módulo GreenYellow")
dotplot(go_cyan, showCategory = 10) + ggtitle("Módulo Cyan")
dotplot(go_salmon, showCategory = 10) + ggtitle("Módulo Salmon")
dotplot(go_tan, showCategory = 10) + ggtitle("Módulo Tan")

### Análisis KEGG

# Cargar paquetes necesarios
library(clusterProfiler)
library(openxlsx)

# Realizar el análisis KEGG para cada módulo
kegg_magenta <- enrichKEGG(gene = magenta_genes_ids$ENTREZID, organism = 'mmu', pAdjustMethod = "BH", qvalueCutoff = 0.05)
kegg_yellow <- enrichKEGG(gene = yellow_genes_ids$ENTREZID, organism = 'mmu', pAdjustMethod = "BH", qvalueCutoff = 0.05)
kegg_brown <- enrichKEGG(gene = brown_genes_ids$ENTREZID, organism = 'mmu', pAdjustMethod = "BH", qvalueCutoff = 0.05)

kegg_green <- enrichKEGG(gene = green_genes_ids$ENTREZID, organism = 'mmu', pAdjustMethod = "BH", qvalueCutoff = 0.05)
kegg_turquoise <- enrichKEGG(gene = turquoise_genes_ids$ENTREZID, organism = 'mmu', pAdjustMethod = "BH", qvalueCutoff = 0.05)
kegg_greenyellow <- enrichKEGG(gene = greenyellow_genes_ids$ENTREZID, organism = 'mmu', pAdjustMethod = "BH", qvalueCutoff = 0.05)

kegg_cyan <- enrichKEGG(gene = cyan_genes_ids$ENTREZID, organism = 'mmu', pAdjustMethod = "BH", qvalueCutoff = 0.05)
kegg_salmon <- enrichKEGG(gene = salmon_genes_ids$ENTREZID, organism = 'mmu', pAdjustMethod = "BH", qvalueCutoff = 0.05)
kegg_tan <- enrichKEGG(gene = tan_genes_ids$ENTREZID, organism = 'mmu', pAdjustMethod = "BH", qvalueCutoff = 0.05)

# Cargar paquete necesario para mapeo
library(org.Mm.eg.db)

# Función para agregar columnas ENS y SYMBOL a los resultados KEGG
add_symbol_and_ens_column <- function(kegg_results) {
  if (!is.null(kegg_results)) {
    # Obtener el mapeo ENTREZID -> ENSG y SYMBOL
    entrez_to_ens_symbol <- AnnotationDbi::select(org.Mm.eg.db, 
                                                  keys = unique(unlist(strsplit(kegg_results@result$geneID, "/"))), 
                                                  keytype = "ENTREZID", 
                                                  columns = c("ENSEMBL", "SYMBOL"))
    
    # Agregar columnas con nombres ENS y SYMBOL
    kegg_results@result$ENSEMBL <- sapply(kegg_results@result$geneID, function(genes) {
      entrez_ids <- unlist(strsplit(genes, "/"))
      ens_ids <- entrez_to_ens_symbol$ENSEMBL[match(entrez_ids, entrez_to_ens_symbol$ENTREZID)]
      paste(ens_ids, collapse = ";")
    })
    
    kegg_results@result$SYMBOL <- sapply(kegg_results@result$geneID, function(genes) {
      entrez_ids <- unlist(strsplit(genes, "/"))
      symbol_ids <- entrez_to_ens_symbol$SYMBOL[match(entrez_ids, entrez_to_ens_symbol$ENTREZID)]
      paste(symbol_ids, collapse = ";")
    })
  }
  return(kegg_results)
}

# Aplicar la función para agregar las columnas ENS y SYMBOL a cada módulo
kegg_magenta <- add_symbol_and_ens_column(kegg_magenta)
kegg_yellow <- add_symbol_and_ens_column(kegg_yellow)
kegg_brown <- add_symbol_and_ens_column(kegg_brown)
kegg_green <- add_symbol_and_ens_column(kegg_green)
kegg_turquoise <- add_symbol_and_ens_column(kegg_turquoise)
kegg_greenyellow <- add_symbol_and_ens_column(kegg_greenyellow)
kegg_cyan <- add_symbol_and_ens_column(kegg_cyan)
kegg_salmon <- add_symbol_and_ens_column(kegg_salmon)
kegg_tan <- add_symbol_and_ens_column(kegg_tan)

# Crear un archivo Excel con una hoja por módulo incluyendo ENS y SYMBOL
wb <- createWorkbook()

# Función para añadir una hoja al workbook
add_kegg_sheet <- function(wb, sheet_name, kegg_results) {
  if (!is.null(kegg_results)) {
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, kegg_results@result)
  }
}

# Añadir hojas para cada módulo
add_kegg_sheet(wb, "Magenta", kegg_magenta)
add_kegg_sheet(wb, "Yellow", kegg_yellow)
add_kegg_sheet(wb, "Brown", kegg_brown)
add_kegg_sheet(wb, "Green", kegg_green)
add_kegg_sheet(wb, "Turquoise", kegg_turquoise)
add_kegg_sheet(wb, "GreenYellow", kegg_greenyellow)
add_kegg_sheet(wb, "Cyan", kegg_cyan)
add_kegg_sheet(wb, "Salmon", kegg_salmon)
add_kegg_sheet(wb, "Tan", kegg_tan)

# Guardar el archivo Excel con ENS y SYMBOL
saveWorkbook(wb, file = "KEGG_Enrichment.xlsx", overwrite = TRUE)

# Cargar paquetes necesarios


# Filtrar los resultados KEGG significativos
filter_significant_kegg <- function(kegg_results) {
  if (!is.null(kegg_results)) {
    significant_results <- kegg_results@result[kegg_results@result$p.adjust <= 0.05, ]
    return(significant_results)
  }
  return(NULL)
}

# Función para agregar la columna SYMBOL a los resultados KEGG significativos
add_symbol_to_significant <- function(kegg_results) {
  if (!is.null(kegg_results)) {
    significant_results <- filter_significant_kegg(kegg_results)
    
    if (!is.null(significant_results)) {
      # Obtener los ENTREZID de los resultados significativos
      entrez_ids <- unique(unlist(strsplit(as.character(significant_results$geneID), "/")))
      entrez_ids <- as.character(entrez_ids)  # Asegurarse de que los ENTREZID sean caracteres
      
      # Obtener el mapeo ENTREZID -> SYMBOL
      entrez_to_symbol <- AnnotationDbi::select(org.Mm.eg.db, 
                                                keys = entrez_ids, 
                                                keytype = "ENTREZID", 
                                                columns = "SYMBOL")
      
      # Agregar la columna SYMBOL a los resultados
      significant_results$SYMBOL <- sapply(significant_results$geneID, function(genes) {
        entrez_ids <- unlist(strsplit(as.character(genes), "/"))
        symbols <- entrez_to_symbol$SYMBOL[match(entrez_ids, entrez_to_symbol$ENTREZID)]
        paste(symbols, collapse = ";")  # Combinar los símbolos por cada geneID
      })
      
      return(significant_results)
    }
  }
  return(NULL)
}

# Crear un archivo Excel con una hoja por módulo para KEGG significativos con SYMBOL
wb_significant <- createWorkbook()

# Función para añadir una hoja con KEGG significativos (y SYMBOL) al workbook
add_significant_kegg_sheet <- function(wb, sheet_name, kegg_results) {
  significant_results <- add_symbol_to_significant(kegg_results)
  if (!is.null(significant_results)) {
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, significant_results)
  }
}

# Añadir hojas con KEGG significativos para cada módulo
add_significant_kegg_sheet(wb_significant, "Magenta", kegg_magenta)
add_significant_kegg_sheet(wb_significant, "Yellow", kegg_yellow)
add_significant_kegg_sheet(wb_significant, "Brown", kegg_brown)
add_significant_kegg_sheet(wb_significant, "Green", kegg_green)
add_significant_kegg_sheet(wb_significant, "Turquoise", kegg_turquoise)
add_significant_kegg_sheet(wb_significant, "GreenYellow", kegg_greenyellow)
add_significant_kegg_sheet(wb_significant, "Cyan", kegg_cyan)
add_significant_kegg_sheet(wb_significant, "Salmon", kegg_salmon)
add_significant_kegg_sheet(wb_significant, "Tan", kegg_tan)

# Guardar el archivo Excel con KEGG significativos y SYMBOL
saveWorkbook(wb_significant, file = "Significant_KEGG.xlsx", overwrite = TRUE)


### Análisis GProfiler


# Crear una función para realizar el análisis de g:Profiler
gprofiler_analysis <- function(gene_list, module_name) {
  # Realizar el análisis
  gost_res <- gost(query = gene_list, organism = "mmusculus", significant = TRUE)
  
  # Si hay resultados, devolverlos como un data frame
  if (!is.null(gost_res$result)) {
    gost_res$result$module <- module_name  # Añadir la información del módulo
    gost_res$result$ENSEMBL <- paste(gene_list, collapse = ",")  # Añadir ENS usados
    return(gost_res$result)
  } else {
    return(NULL)
  }
}

# Realizar análisis por cada módulo
gprofiler_magenta <- gprofiler_analysis(magenta_genes, "Magenta")
gprofiler_yellow <- gprofiler_analysis(yellow_genes, "Yellow")
gprofiler_brown <- gprofiler_analysis(brown_genes, "Brown")
gprofiler_green <- gprofiler_analysis(green_genes, "Green")
gprofiler_turquoise <- gprofiler_analysis(turquoise_genes, "Turquoise")
gprofiler_greenyellow <- gprofiler_analysis(greenyellow_genes, "GreenYellow")
gprofiler_cyan <- gprofiler_analysis(cyan_genes, "Cyan")
gprofiler_salmon <- gprofiler_analysis(salmon_genes, "Salmon")
gprofiler_tan <- gprofiler_analysis(tan_genes, "Tan")

# Combinar todos los resultados en una única tabla
gprofiler_all <- do.call(rbind, list(
  gprofiler_magenta, gprofiler_yellow, gprofiler_brown,
  gprofiler_green, gprofiler_turquoise, gprofiler_greenyellow,
  gprofiler_cyan, gprofiler_salmon, gprofiler_tan
))

# Guardar los resultados en un archivo Excel
write.xlsx(gprofiler_all, file = "gProfiler_Enrichment_AllModules.xlsx", rowNames = FALSE)
