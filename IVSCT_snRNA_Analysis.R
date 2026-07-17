pacman::p_load(
  tidyverse,
  Seurat,
  clustree,
  pheatmap,
  ggsci,
  plyr,
  reshape2,
  SingleR,
  tidyr,
  dplyr,
  DoubletFinder,
  SeuratDisk,
  symphony,
  data.table,
  matrixStats,
  Matrix,
  harmony,
  ggplot2,
  ggthemes,
  ggrastr,
  RColorBrewer,
  patchwork,
  presto,
  ggrepel,
  speckle,
  limma,
  scProportionTest,
  CellChat,
  ggpubr,
  decontX
)

set.seed(0422)
options(future.globals.maxSize = 20 * 1000 * 1024 ^ 2)

fromList <- function(input) {
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    elements %in% x
  }))
  data <- matrix(data, ncol = length(input), byrow = FALSE)
  data <- data.frame(data)
  names(data) <- names(input)
  rownames(data) <- elements
  return(data)
}
#drawdata
base_dir <- "../snRNA1024/"
files <- list.files(base_dir,  full.names = TRUE)
file_name <- list.files(base_dir,  full.names = F)
obj = list()
for (i in 1:6) {
  counts <-
    Read10X_h5(paste0(files[i], '/filtered_feature_bc_matrix.h5'))
  colnames(counts) = paste0(file_name[i], colnames(counts))
  obj[[file_name[i]]] = CreateSeuratObject(
    counts,
    min.cells = 3,
    min.features = 200,
    project = file_name[[i]]
  )
  obj[[file_name[i]]]$sample = file_name[[i]]
}

scRNA = merge(obj[[1]], obj[-1])
scRNA
scRNA <-
  PercentageFeatureSet(object = scRNA,
                       pattern = "^mt-",
                       col.name = "percent.mt")
scRNA <-
  PercentageFeatureSet(object = scRNA,
                       pattern = "^Hb[ab]",
                       col.name = "percent.hb")
scRNA <-
  PercentageFeatureSet(object = scRNA,
                       pattern = "^Rp[sl]",
                       col.name = "percent.rb")
scRNA = subset(
  scRNA,
  nFeature_RNA > 200 &
    nFeature_RNA < 7000 & percent.mt < 5 &
    percent.hb < 5 & percent.rb < 15
)
scRNA  = SCTransform(scRNA, verbose = TRUE)
scRNA = RunPCA(scRNA, verbose = FALSE)
scRNA <- RunHarmony(scRNA, 'sample', plot_convergence = T)
scRNA = RunUMAP(scRNA,
                dims = 1:30,
                reduction = 'harmony',
                verbose = FALSE)
scRNA_test = FindNeighbors(scRNA, reduction = 'harmony', verbose = FALSE)
scRNA_test = FindClusters(scRNA_test, resolution = c(seq(0, 1.6, 0.1)), verbose = FALSE)
clustree(scRNA_test, prefix = "SCT_snn_res.")
scRNA = FindNeighbors(scRNA,  reduction = 'harmony', verbose = FALSE)
scRNA = FindClusters(scRNA, resolution = 0.2, verbose = FALSE)
DimPlot(scRNA, label = T)
gene <-
  c(
    'Astrocyte' = 'Aldh1l1',
    'Astrocyte' = 'Aqp4',
    'Astrocyte' = 'Gpc5',
    'Astrocyte' = 'Gja1',
    'Astrocyte' = 'Gfap',
    'Microglia' = 'C1qa',
    'Microglia' = 'Cd74',
    'Microglia' = 'Csf1r',
    'Microglia' = 'P2ry12',
    'Microglia' = 'Tmem119',
    'Microglia' = 'Cx3cr1',
    'Microglia' = 'Hexb',
    'Inh' = 'Gad1',
    'Inh' = 'Gad2',
    'Inh' = 'Slc32a1',
    'Exc' = 'Nrgn',
    'Exc' = 'Slc17a7',
    'Exc' = 'Neurod6',
    'GC' = 'Prox1',
    'Endothelial' = 'Pecam1',
    'Endothelial' = 'Cdh5',
    'Endothelial' = 'Ebf1',
    'Endothelial' = 'Flt1',
    'Endothelial' = 'Cldn5',
    'Oligo' = 'Mbp',
    'Oligo' = 'Mobp',
    'Oligo' = 'Plp1',
    'Oligo' = 'Olig1',
    'OPC' = 'Vcan'
  )
scRNA$annotation <- plyr::mapvalues(
  scRNA$seurat_clusters,
  c(0, 3, 4, 7, 8, 2, 5, 10, 11, 1, 6,
    12, 9, 13, 14),
  c(
    rep('Exc', 4),
    rep('GC', 1),
    rep('Inh', 3),
    rep('Astrocyte', 1),
    rep('Oligo', 2),
    rep('OPC', 1),
    rep('Microglia', 1),
    rep('Doublets', 1),
    rep('Endothelial', 1)
  )
)
scRNA = subset(scRNA, annotation != 'Doublets')

#Decount
counts <- scRNA@assays$RNA@counts
decontX_results <- decontX(counts)
scRNA$Contamination = decontX_results$contamination
scRNA$Con = ifelse(scRNA$Contamination < 0.2, 'No', 'Yes')
table(scRNA$annotation, scRNA$Con)
scRNA = subset(scRNA, Con == 'No')
#Doubletfinder
scRNA_list <- list()
sample_names <- unique(scRNA$sample)
for (sample_id in sample_names) {
  sc_one <- subset(scRNA, sample == sample_id)
  sweep.res <- paramSweep_v3(sc_one, PCs = 1:30, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  homotypic.prop <- modelHomotypic(sc_one$annotation)
  nExp_poi <- round(ncol(sc_one) * 8 * 1e-6)
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  sc_one <-
    doubletFinder_v3(
      sc_one,
      PCs = 1:30,
      pN = 0.25,
      pK = mpK,
      nExp = nExp_poi,
      reuse.pANN = FALSE,
      sct = TRUE
    )
  col_names <- colnames(sc_one@meta.data)
  pANN_col <- col_names[grep("pANN", col_names)]
  class_col <- col_names[grep("DF.classifications", col_names)]
  sc_one <-
    doubletFinder_v3(
      sc_one,
      PCs = 1:30,
      pN = 0.25,
      pK = mpK,
      nExp = nExp_poi.adj,
      reuse.pANN = pANN_col,
      sct = TRUE
    )
  new_class_col <-
    colnames(sc_one@meta.data)[grep("DF.classifications", colnames(sc_one@meta.data))][2]
  sc_one$Doublets <- "Singlet"
  sc_one$Doublets[sc_one[[class_col]] == "Doublet"] <-
    "Doublet-High Confidience"
  sc_one$Doublets[sc_one[[class_col]] == "Doublet" &
                    sc_one[[new_class_col]] == "Singlet"] <- "Doublet-Low Confidience"
  scRNA_list[[sample_id]] <- sc_one
}
merged_obj <- merge(scRNA_list[[1]], y = scRNA_list[-1])
scRNA$Doublets <- merged_obj$Doublets
scRNA <- subset(scRNA, Doublets %in% c('Singlet'))

#ratio
pro = sc_utils(scRNA)
prop_test <- permutation_test(
  pro,
  cluster_identity = "annotation",
  sample_1 = "MAPT",
  sample_2 = "IVSCT",
  sample_identity = "orig.ident"
)

#DEG
DefaultAssay(scRNA) = 'RNA'
scRNA = NormalizeData(scRNA, assay = 'RNA')
target_groups <-
  c("Global", "Exc", "Inh", "Oligo", "GC", "Astrocyte", "Microglia")
DEG.list <- list()
for (group_name in target_groups) {
  if (group_name == "Global") {
    sub_obj <- scRNA
  } else {
    sub_obj <- subset(scRNA, annotation == group_name)
  }
  de_res <-
    presto::wilcoxauc(
      sub_obj,
      group_by = 'orig.ident',
      seurat_assay = 'RNA',
      assay = 'data'
    )
  DEG.list[[group_name]] = de_res %>% filter(group == "IVSCT")
}

##subcluster
Inh = subset(scRNA, annotation == 'Inh')
DefaultAssay(Inh) = 'RNA'
Inh = NormalizeData(Inh, assay = 'RNA')
Inh  = SCTransform(Inh, verbose = TRUE)
Inh = RunPCA(Inh, verbose = FALSE)
Inh <- RunHarmony(Inh, 'sample', plot_convergence = T)
Inh = RunUMAP(Inh,
              dims = 1:10,
              reduction = 'harmony',
              verbose = FALSE)
Inh = FindNeighbors(Inh,
                    reduction = 'harmony',
                    verbose = FALSE,
                    dims = 1:10)
Inh = FindClusters(Inh, resolution = 0.1, verbose = FALSE)
Inh = subset(Inh, seurat_clusters %in% c(0, 1, 2, 3, 4, 5, 6))
Inh  = SCTransform(Inh, verbose = TRUE, variable.features.n = 2000)
Inh = RunPCA(Inh, verbose = FALSE)
Inh <- RunHarmony(Inh, 'sample', plot_convergence = T)
Inh = RunUMAP(Inh,
              dims = 1:10,
              reduction = 'harmony',
              verbose = FALSE)
Inh = FindNeighbors(Inh,  reduction = 'harmony', verbose = FALSE)
Inh = FindClusters(Inh, resolution = 0.1, verbose = FALSE)
pro = sc_utils(Inh)
prop_test <- permutation_test(
  pro,
  cluster_identity = "seurat_clusters",
  sample_1 = "MAPT",
  sample_2 = "IVSCT",
  sample_identity = "orig.ident"
)

Exc = subset(scRNA, annotation == 'Exc')
Exc  = SCTransform(Exc, verbose = TRUE)
Exc = RunPCA(Exc, verbose = FALSE)
Exc <- RunHarmony(Exc, 'sample', plot_convergence = T)
Exc = RunUMAP(Exc,
              dims = 1:20,
              reduction = 'harmony',
              verbose = FALSE)
Exc = FindNeighbors(Exc,  reduction = 'harmony', verbose = FALSE)
Exc = FindClusters(Exc, resolution = 0.1, verbose = FALSE)
pro = sc_utils(Exc)
prop_test <- permutation_test(
  pro,
  cluster_identity = "seurat_clusters",
  sample_1 = "MAPT",
  sample_2 = "IVSCT",
  sample_identity = "orig.ident"
)


##DG snRNA human deg
DG.counts <- readRDS("DG_granule_counts.rds")
DG.meta <- readRDS("DG_granule_metadata.rds")
DG.H <- CreateSeuratObject(counts = DG.counts, meta.data = DG.meta)
DG.H$Braak_stage = ifelse(DG.H$Braak > 3, 'High', 'Low')
DG.H <- NormalizeData(DG.H, assay = "RNA")
Idents(DG.H) = DG.H$Braak_stage
DG.H$MAST_Group = ifelse(DG.H$Annotation_sub == 'DG granule cells', 'GC', DG.H$Annotation)
deg = list()
for (i in c('Global', 'Ast', 'End', 'Exc', 'GC', 'Inh', 'Mic', 'Oli', 'OPC')) {
  if (i == 'Global') {
    sub <- DG.H
  } else {
    sub <- subset(DG.H, MAST_Group == i)
  }
  deg[[i]] = FindMarkers(
    sub,
    ident.1 = 'High',
    ident.2 = 'Low',
    test.use = "MAST",
    latent.vars = c('PMI', 'Age', 'Gender'),
    assay = "RNA" ,
    logfc.threshold = 0.01
  )
}
#PSP bulk deg
PSP.deg = read.csv('Mayo.PSP_vs_Con.deg.unselect.csv')
PSP.deg = PSP.deg[PSP.deg$adj.P.Val < 0.05, ]
PSP.deg = PSP.deg[abs(PSP.deg$logFC) > log2(1.5), ]
##PSP snRNA deg
PSP.deg = read.csv('./PMID39648200DEG.csv')
PSP.deg.neuron = PSP.deg[PSP.deg$cluster %in% c('Exc', 'Inh'), ]
PSP.deg.neuron = PSP.deg.neuron[!is.na(PSP.deg.neuron$DEG.ID), ]
PSP.deg.neuron = PSP.deg.neuron[PSP.deg.neuron$p_val_adj < 0.05,]
Inh1 = subset(plot_data, celltype.cluster == 'Inh1')
clusters <- unique(Inh1$celltype.cluster)
deg_results_list <- list()
DefaultAssay(plot_data_sub) <- "RNA"
for (cluster in clusters) {
  cell_subset <-
    subset(plot_data_sub, subset = celltype.cluster == cluster)
  table_groups <- table(cell_subset$disease.status)
  if (!all(c("PSP", "Control") %in% names(table_groups)) ||
      any(table_groups < 3)) {
    next
  }
  Idents(cell_subset) <- "disease.status"
  deg_res <- FindMarkers(
    object = cell_subset,
    ident.1 = "PSP",
    ident.2 = "Control",
    test.use = "MAST",
    assay = 'RNA',
    slot = 'counts',
    latent.vars = c("PMI", "percent.ribo", "percent.mt", "age"),
    logfc.threshold = 0,
    min.pct = 0.1
  )
  deg_res$gene <- rownames(deg_res)
  deg_res$cluster <- cluster
  deg_results_list[[cluster]] <- deg_res
}
all_psp_vs_control_degs <- do.call(rbind, deg_results_list)

#GSEA snRNA human AD
human.deg = readRDS('HPC.braak.deg.rds')
for (i in 1:9) {
  human.deg[[i]]$gene = rownames(human.deg[[i]])
}
species_name <- "Homo sapiens"
m_t2g_kegg <-
  msigdbr(species = species_name,
          category = "C2",
          subcategory = "CP:KEGG") %>%
  dplyr::select(gs_name, gene_symbol)
run_gsea_from_df <- function(deg_df, group_name) {
  message("Analyzing GSEA for: ", group_name)
  logfc_col <-
    intersect(colnames(deg_df), c("avg_log2FC", "avg_logFC", "logFC"))[1]
  if (is.na(logfc_col)) {
    warning("No logFC column found for ", group_name)
    return(NULL)
  }
  if ("gene" %in% colnames(deg_df)) {
    gene_names <- deg_df$gene
  } else {
    gene_names <- rownames(deg_df)
  }
  geneList <- deg_df[[logfc_col]]
  names(geneList) <- gene_names
  geneList <- geneList[!is.na(geneList)]
  geneList <- geneList + rnorm(length(geneList), sd = 1e-10)
  geneList <- sort(geneList, decreasing = TRUE)
  gsea_res <- GSEA(
    geneList,
    TERM2GENE = m_t2g_kegg,
    pvalueCutoff = 0.05,
    minGSSize = 10,
    maxGSSize = 500,
    verbose = FALSE
  )
  return(as.data.frame(gsea_res))
}
all_results_human_kegg <- list()
for (grp in names(human.deg)) {
  res <- tryCatch({
    run_gsea_from_df(human.deg[[grp]], grp)
  }, error = function(e) {
    message("Error in ", grp, ": ", e$message)
    NULL
  })
  if (!is.null(res) && nrow(res) > 0) {
    all_results_human_kegg[[grp]] <- res
  }
}

####GSEA psp human bulk
psp.deg = read.csv("Mayo.PSP_vs_Con.deg.unselect.csv")
psp.deg = psp.deg[!duplicated(psp.deg$GeneName), ]
psp.deg = psp.deg[!is.na(psp.deg$GeneName), ]
geneList <- psp.deg[, 'logFC']
names(geneList) <- psp.deg$GeneName
geneList <- geneList[!is.na(geneList)]
geneList <- geneList + rnorm(length(geneList), sd = 1e-10)
geneList <- sort(geneList, decreasing = TRUE)
species_name <- "Homo sapiens"
m_t2g_kegg <-
  msigdbr(species = species_name,
          category = "C2",
          subcategory = "CP:KEGG") %>%
  dplyr::select(gs_name, gene_symbol)
gsea_res <- GSEA(
  geneList,
  TERM2GENE = m_t2g_kegg,
  pvalueCutoff = 1,
  minGSSize = 10,
  maxGSSize = 500,
  verbose = FALSE
)

m_t2g_go <-
  msigdbr(species = species_name,
          category = "C5",
          subcategory = "GO:BP") %>%
  dplyr::select(gs_name, gene_symbol)

gsea_res_go <- GSEA(
  geneList,
  TERM2GENE = m_t2g_go,
  pvalueCutoff = 1,
  minGSSize = 10,
  maxGSSize = 500,
  verbose = FALSE
)
gsea_res_go@result$Description <-
  gsub("Gobp ", "", gsea_res_go@result$Description) %>%
  gsub("_", " ", .) %>% str_to_title()

##GSEA PSP snRNA
deg_df <- PSP.deg.Inh1[!is.na(PSP.deg.Inh1$gene.symbol),]
gene_map <-
  bitr(
    deg_df$gene.symbol,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
deg_df <-
  merge(deg_df, gene_map, by.x = "gene.symbol", by.y = "SYMBOL")
deg_df <- deg_df[order(abs(deg_df$avg_log2FC), decreasing = TRUE),]
deg_df <- deg_df[!duplicated(deg_df$ENTREZID),]
gene_list <- deg_df$avg_log2FC
names(gene_list) <- deg_df$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)
gsea_results <- list()
gse_go <- gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  keyType = "ENTREZID",
  minGSSize = 10,
  maxGSSize  = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)
if (!is.null(gse_go) && nrow(gse_go) > 0) {
  gse_go <-
    setReadable(gse_go, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
}
gsea_results[["GO"]] <- gse_go
gse_kegg <-
  gseKEGG(
    geneList  = gene_list,
    organism     = "hsa",
    minGSSize    = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = FALSE
  )
if (!is.null(gse_kegg) && nrow(gse_kegg) > 0) {
  gse_kegg <-
    setReadable(gse_kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
}
gsea_results[["KEGG"]] <- gse_kegg
gse_reactome <-
  gsePathway(
    geneList = gene_list,
    organism = "human",
    minGSSize    = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = FALSE
  )

if (!is.null(gse_reactome) && nrow(gse_reactome) > 0) {
  gse_reactome <-
    setReadable(gse_reactome, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
}
gsea_results[["Reactome"]] <- gse_reactome

##AUCELL
species_name <- 'Mus musculus'
reac_df <-
  msigdbr(species = species_name,
          category = "C2",
          subcategory = "CP:REACTOME")
gene_sets <- reac_df %>% split(f = .$gs_name)
names(gene_sets) <- sapply(names(gene_sets), function(x) {
  x %>% gsub("KEGG_|REACTOME_", "", .) %>% gsub("_", " ", .) %>%
    tolower() %>% tools::toTitleCase() %>% gsub(" ", "_", .)
})
gene_sets <- lapply(gene_sets, function(x)
  x$gene_symbol)
cells_rankings <-
  AUCell_buildRankings(GetAssayData(scRNA, slot = "counts"), plotStats = FALSE)
cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings)
auc_matrix <- getAUC(cells_AUC)
for (pathway in rownames(auc_matrix)) {
  scRNA[[pathway]] <- auc_matrix[pathway,]
}

species_name <- "Homo sapiens"
reac_df <-
  msigdbr(species = species_name,
          category = "C2",
          subcategory = "CP:REACTOME")
gene_sets <- reac_df %>% split(f = .$gs_name)
names(gene_sets) <- sapply(names(gene_sets), function(x) {
  x %>% gsub("KEGG_|REACTOME_", "", .) %>% gsub("_", " ", .) %>%
    tolower() %>% tools::toTitleCase() %>% gsub(" ", "_", .)
})
gene_sets <- lapply(gene_sets, function(x)
  x$gene_symbol)


##pySCENIC
seu_subset <- subset(scRNA, annotation  %in% c('Exc', 'GC', 'Inh'))
seu_subset$annotation = droplevels(seu_subset$annotation)
loom_file <- "Exc_Inh_GC_for_pyscenic.loom"
DefaultAssay(seu_subset) = 'RNA'
Seu2Loom(seu_subset, filename = loom_file , layers = 'RNA')
loom_file <- "./snRNA1024/Final/Inh_pyscenic.loom"
DefaultAssay(Inh) = 'RNA'
Seu2Loom(Inh, filename = loom_file , layers = 'RNA')

#loom <- open_loom("aucell_output.inh.loom")
loom <- open_loom("aucell_output.loom")
regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name = 'RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)
close_loom(loom)
sub_regulonAUC <-
  regulonAUC[, match(colnames(seu_subset), colnames(regulonAUC))]
seu_subset
seu_subset$Group = paste0(seu_subset$orig.ident, '.', seu_subset$annotation)
Groups <-
  data.frame(row.names = colnames(seu_subset), Group = seu_subset$Group)
selectedResolution <- "Group"
cellsPerGroup <-
  split(rownames(Groups), Groups[, selectedResolution])
sub_regulonAUC <-
  sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)), ]
regulonActivity_byGroup <-
  sapply(cellsPerGroup, function(cells)
    rowMeans(getAUC(sub_regulonAUC)[, cells]))
regulonActivity_byGroup_Scaled <-
  t(scale(
    t(regulonActivity_byGroup),
    center = T,
    scale = T
  ))
regulonActivity_byGroup_Scaled = na.omit(regulonActivity_byGroup_Scaled)
rss <-
  calcRSS(AUC = getAUC(sub_regulonAUC), cellAnnotation = Groups[colnames(sub_regulonAUC), selectedResolution])
rss = regulonActivity_byGroup_Scaled

Inh3.deg <- c(
  "Rasgrf1",
  "Mapt",
  "Atp2a2",
  "Plekha5",
  "Gm42439",
  "Hsf5",
  "E4f1",
  "Slc15a2",
  "Kmt2e",
  "mt-Co1",
  "Cplane1",
  "Luc7l3",
  "Gphn",
  "mt-Co2",
  "Luc7l",
  "Rif1",
  "Ankrd12",
  "Cmss1",
  "Lars2",
  "ENSMUSG00000095041"
)
regu = rownames(sub_regulonAUC)
regulons[[regu]][regulons[[regu]] %in% Inh3.deg]

#hdWGCNA
theme_set(theme_cowplot())
scRNA@active.assay <- "RNA"
scRNA$seurat_group = scRNA$annotation
Idents(scRNA) <- scRNA$seurat_group
scRNA <-
  SetupForWGCNA(
    scRNA,
    gene_select = "fraction",
    fraction = 0.05,
    wgcna_name = 'GC'
  )
scRNA <-
  MetacellsByGroups(
    scRNA = scRNA,
    group.by = c("seurat_group", 'sample'),
    k = 25,
    max_shared = 10,
    reduction = 'harmony',
    ident.group = 'seurat_group'
  )
scRNA <- NormalizeMetacells(scRNA)
scRNA <-
  SetDatExpr(
    scRNA,
    group_name = 'GC',
    group.by = 'seurat_group',
    assay = 'RNA',
    slot = 'data'
  )
scRNA <- TestSoftPowers(scRNA, networkType = 'signed')
scRNA <- ConstructNetwork(scRNA, overwrite_tom = T, tom_name = 'GC')
TOM <- GetTOM(scRNA)
scRNA <- ScaleData(scRNA, features = VariableFeatures(scRNA))
scRNA <- ModuleEigengenes(scRNA, group.by.vars = 'sample')
hMEs <- GetMEs(scRNA)
MEs <- GetMEs(scRNA, harmonized = FALSE)
scRNA <-
  ModuleConnectivity(scRNA, group.by = 'seurat_group', group_name = 'GC')
scRNA <- ResetModuleNames(scRNA, new_name = paste0('GC', '-M'))
modules <- GetModules(scRNA) %>% subset(module != 'grey')
hub_df <- GetHubGenes(scRNA, n_hubs = 25)
head(hub_df)
MEs <- GetMEs(scRNA, harmonized = TRUE)
modules <- GetModules(scRNA)
mods <- levels(modules$module)
mods <- mods[mods != 'grey']
scRNA@meta.data <- cbind(scRNA@meta.data, MEs)
hme = scRNA@meta.data
colnames(hme) = gsub('-', '.', colnames(hme))

stat_test <-
  hme %>% group_by(annotation) %>% wilcox_test(GC.M1 ~ orig.ident) %>%
  adjust_pvalue(method = "BH") %>%   add_significance("p.adj") %>% add_xy_position(x = "orig.ident")
stat_test <-
  hme %>% group_by(annotation) %>% wilcox_test(GC.M1 ~ orig.ident) %>%
  adjust_pvalue(method = "bonferroni") %>%  add_significance("p.adj") %>%
  filter(p.adj < 0.05) %>% add_xy_position(x = "orig.ident")

#cellchat
MAPT <- subset(scRNA, orig.ident == 'MAPT')
IVSCT <- subset(scRNA, orig.ident == 'IVSCT')
data_list <- list(MAPT, IVSCT)
names(data_list) <- c('MAPT', 'IVSCT')
CellChatDB <- CellChatDB.mouse
cellchat_list <- list()
for (i in 1:2) {
  Seuratobject = data_list[[i]]
  Seuratobject = createCellChat(Seuratobject@assays$SCT@data,
                                meta = Seuratobject@meta.data,
                                group.by = 'cc_group')
  Seuratobject@DB <- CellChatDB
  Seuratobject <- subsetData(Seuratobject)
  Seuratobject <- identifyOverExpressedGenes(Seuratobject)
  Seuratobject <- identifyOverExpressedInteractions(Seuratobject)
  Seuratobject <- projectData(Seuratobject, adj = PPI.mouse)
  Seuratobject <- computeCommunProb(Seuratobject)
  Seuratobject <- filterCommunication(Seuratobject, min.cells = 10)
  Seuratobject <- computeCommunProbPathway(Seuratobject)
  Seuratobject <- aggregateNet(Seuratobject)
  cellchat_list[[i]] <-  Seuratobject
}
