### DA-seq applied on immune cells from melanoma patients
# Author: Jun Zhao (jun.zhao@yale.edu)

# Original paper: https://www.sciencedirect.com/science/article/pii/S0092867418313941

library(Seurat)
library(reshape2)

#! Please set working directory to your local github repo of DA-seq
source("./DA-seq.R")
dir.create("./plots/")



## Set some parameters

# Seurat version, default 2, can also be 3
Sversion <- 2

# set the python path to use
python <- "/usr/bin/python"

# set GPU number to use
GPU <- ""



## Data preprocessing

# load data
download.file(
  "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120575/suppl/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz",
  "./data/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz"
)
data_exp <- read.table(
  "./data/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz",
  sep = "\t", header = F, row.names = 1, stringsAsFactors = F, skip = 2
)
data_exp <- data_exp[,-16292]

data_colInfo <- read.table(
  "./data/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz",
  sep = "\t", header = F, stringsAsFactors = F, nrows = 2
)
data_colInfo <- data_colInfo[,-1]

colnames(data_exp) <- data_colInfo[1,]

data_exp <- Matrix(as.matrix(data_exp), sparse = T)


# patient info
download.file(
  "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120575/suppl/GSE120575_patient_ID_single_cells.txt.gz",
  "./data/GSE120575_patient_ID_single_cells.txt.gz"
)
patient_info <- read.table(
  "./data/GSE120575_patient_ID_single_cells.txt.gz", 
  sep = "\t", header = T, stringsAsFactors = F, skip = 19, nrows = 16291
)
mean(colnames(data_exp) == patient_info$title)
rownames(patient_info) <- patient_info$title


# cluster info
cluster_info <- read.table(
  "./data/cluster_info", sep = "\t", header = T, stringsAsFactors = F
)
rownames(cluster_info) <- cluster_info$Cell.Name


# patient lesion info
lesion_info <- read.table(
  "./data/patient_lesion", sep = "\t", header = T, stringsAsFactors = F
)


# Seurat
if(Sversion == 2){
  data_S <- CreateSeuratObject(
    raw.data = data_exp, project = "melanoma.immune"
  )
} else if(Sversion == 3){
  data_S <- CreateSeuratObject(
    counts = data_exp, project = "melanoma.immune"
  )
} else {
  warning("Please set Sversion to 2 or 3!")
}


# set metadata for each cell
data_S@meta.data$condition <- patient_info[colnames(data_exp), "characteristics..response"]
data_S@meta.data$lesion <- patient_info[colnames(data_exp), 
                                        "characteristics..patinet.ID..Pre.baseline..Post..on.treatment."]
data_S@meta.data$cluster <- cluster_info$Cluster.number

data_S <- ScaleData(data_S)


# calculate gene variance to set variable genes
gene_var <- apply(data_exp, 1, var)
if(Sversion == 2){
  data_S@var.genes <- names(gene_var)[gene_var > 6]
} else if(Sversion == 3){
  data_S@assays$RNA@var.features <- names(gene_var)[gene_var > 6]
} else {
  warning("Please set Sversion to 2 or 3!")
}


# dimension reduction
if(Sversion == 2){
  data_S <- RunPCA(data_S, pcs.compute = 10, do.print = F)
  data_S <- RunTSNE(data_S, dims.use = 1:10, seed.use = 1)
  TSNEPlot(data_S, group.by = "condition", pt.size = 0.5)
  TSNEPlot(data_S, group.by = "cluster", do.label = T, pt.size = 0.5)
  pca_embedding <- data_S@dr$pca@cell.embeddings
  tsne_embedding <- data_S@dr$tsne@cell.embeddings
} else if(Sversion == 3){
  data_S <- RunPCA(data_S, npcs = 10, assay = "RNA", features = VariableFeatures(data_S), verbose = F)
  data_S <- RunTSNE(data_S, dims = 1:10, seed.use = 1)
  TSNEPlot(data_S, group.by = "condition", pt.size = 0.5)
  TSNEPlot(data_S, group.by = "cluster", label = T, pt.size = 0.5)
  pca_embedding <- data_S@reductions$pca@cell.embeddings
  tsne_embedding <- data_S@reductions$tsne@cell.embeddings
} else {
  warning("Please set Sversion to 2 or 3!")
}



## DA-seq

# set k values to use
k.vector <- seq(50, 500, 50)


# sample labels for responders and non-responders
labels_res <- lesion_info[lesion_info$Response.status..R.responder..NR.non.responder == "R","Sample.name"]
labels_nonres <- lesion_info[lesion_info$Response.status..R.responder..NR.non.responder == "NR","Sample.name"]


# find top DA cells
da_cells <- getDAcells(
  X = pca_embedding, 
  cell.labels = data_S@meta.data$lesion, 
  labels.1 = labels_res, 
  labels.2 = labels_nonres, 
  k.vector = k.vector, 
  k.folds = 10, n.runs = 10, pred.thres = c(0.075,0.925), 
  plot.embedding = tsne_embedding, 
  python.use = python, 
  source.code = "./DA_logit.py", GPU = GPU
)

da_cells <- updateDAcells(
  X = da_cells, pred.thres = c(0.075,0.925), 
  do.plot = T, plot.embedding = tsne_embedding
)

da_cells$pred.plot
da_cells$da.cells.plot

ggsave(da_cells$pred.plot + xlab("tSNE_1") + ylab("tSNE_2"), 
       filename = "./plots/melanoma_prediction.pdf", width = 5, height = 4.5)
ggsave(da_cells$da.cells.plot + xlab("tSNE_1") + ylab("tSNE_2"),
       filename = "./plots/melanoma_DA_cells.pdf", width = 4.5, height = 4.5)



# cluster DA cells to get DA regions
da_regions <- getDAregion(
  X = pca_embedding, 
  cell.idx = da_cells$da.cell.idx, 
  k = 5, alpha = 0.25, iter.max = 30, 
  cell.labels = data_S@meta.data$lesion,
  labels.1 = labels_res, 
  labels.2 = labels_nonres, 
  plot.embedding = tsne_embedding
)
da_regions$da.region.plot
da_regions$DA.stat

ggsave(da_regions$da.region.plot + xlab("tSNE_1") + ylab("tSNE_2"),
       filename = "./plots/melanoma_DA_regions.pdf", width = 5, height = 4.5)
write.table(da_regions$DA.stat, file = "./plots/melanoma_DA_region_stats.txt", sep = "\t", 
            col.names = T, row.names = F, quote = F)


n.da <- length(unique(da_regions$cluster.res)) - 1


# STG marker finder for DA regions
STG.markers <- STGmarkerFinder(
  X = as.matrix(data_S@data), 
  cell.idx = da_cells$da.cell.idx, 
  da.region.label = da_regions$cluster.res, 
  lambda = 1.5, n.runs = 5,  
  python.use = "/path/to/your/python", 
  source.code = "./DA_STG.py"
)

head(STG.markers$da.markers[["1"]])


# get markers as well as the model
STG.genes.model <- STGmarkerFinder(
  X = as.matrix(data_S@data), 
  cell.idx = da_cells$da.cell.idx, 
  da.region.label = da_regions$cluster.res, 
  lambda = 1.2, n.runs = 5, return.model = T, 
  python.use = "/path/to/your/python", 
  source.code = "./DA_STG.py"
)

# plot linear predictions from the model for DA region 5
gg <- plotCellScore(
  X = tsne_embedding, 
  score = STG.genes.model$model[["5"]][["pred"]], 
  cell.col = viridis_pal()(10)
)
ggsave(gg + xlab("tSNE_1") + ylab("tSNE_2"), 
       filename = "./plots/melanoma_DA5_STGlinearPred.pdf", width = 5, height = 4.5)



## Interpret DA regions

# plot some markers
marker_genes <- c(
  "CD19","MS4A1","IGHM","CD79A",
  "VCAM1","LAG3","CD27","CD38",
  "CD14","CSF3R","VCAN","LYZ",
  "IL7R","TCF7","CD8A","CCL5",
  "CCR7","LEF1","SELL"
)

# add "da" slot
data_S@meta.data$da <- 0
data_S@meta.data$da[da_cells$da.cell.idx] <- da_regions$cluster.res


# add STG information
STG.marker.info <- do.call(rbind, lapply(STG.genes.model$da.markers, function(x,inputgenes){
  as.numeric(inputgenes %in% x$gene)
}, inputgenes = rev(marker_genes)))
STG.marker.info <- rbind(0, STG.marker.info)
colnames(STG.marker.info) <- rev(marker_genes)
rownames(STG.marker.info) <- c(1:(n.da+1))
STG.marker.info[STG.marker.info == 0] <- NA

STG.marker.info.m <- melt(STG.marker.info)
STG.marker.info.m <- STG.marker.info.m[-which(is.na(STG.marker.info.m$value)),]
STG.marker.info.m$value <- 10 * STG.marker.info.m$value


# create dot plot with Seurat
if(Sversion == 2){
  STG.marker.info.m$value <- STG.marker.info.m$value / 100
  ggdot <- DotPlot(
    data_S, genes.plot = marker_genes, cols.use = c("gray","red"), group.by = "da", 
    x.lab.rot = T, do.return = T
  )
} else if(Sversion == 3){
  ggdot <- DotPlot(data_S_new, features = marker_genes, cols = c("gray","red"), group.by = "da") + 
    RotatedAxis()
} else {
  warning("Please set Sversion to 2 or 3!")
}

gg <- ggdot + geom_point(data = STG.marker.info.m, aes(x = Var2, y = Var1, size = value))
ggsave(plot = gg, filename = "./plots/melanoma_DotPlot.pdf", width = 10, height = 5)



