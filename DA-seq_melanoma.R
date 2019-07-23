### DA-seq applied on immune cells from melanoma patients
# Author: Jun Zhao (jun.zhao@yale.edu)

# Original paper: https://www.sciencedirect.com/science/article/pii/S0092867418313941

library(Seurat) # Seurat version 2.3.0

#! Please set working directory to your local github repo of DA-seq
source("./DA-seq.R")



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
data_S <- CreateSeuratObject(
  raw.data = data_exp, project = "melanoma.immune"
)


# set metadata for each cell
data_S@meta.data$condition <- patient_info[data_S@cell.names, "characteristics..response"]
data_S@meta.data$lesion <- patient_info[data_S@cell.names, "characteristics..patinet.ID..Pre.baseline..Post..on.treatment."]
data_S@meta.data$cluster <- cluster_info$Cluster.number

data_S <- ScaleData(data_S)


# calculate gene variance to set variable genes
gene_var <- apply(data_exp, 1, var)
data_S@var.genes <- names(gene_var)[gene_var > 6]


# dimension reduction
data_S <- RunPCA(data_S, pcs.compute = 10, do.print = F)

data_S <- RunTSNE(data_S, dims.use = 1:10)
TSNEPlot(data_S, group.by = "condition", pt.size = 0.5)
TSNEPlot(data_S, group.by = "cluster", do.label = T, pt.size = 0.5)

pca_embedding <- data_S@dr$pca@cell.embeddings
tsne_embedding <- data_S@dr$tsne@cell.embeddings



## DA-seq

# calculate diffusion coordinates
library(diffusionMap)
diffuse.res <- diffuse(D = dist(pca_embedding), neigen = 20)


# set k values to use
k.vector <- seq(50, 500, 50)


# sample labels for responders and non-responders
labels_res <- lesion_info[lesion_info$Response.status..R.responder..NR.non.responder == "R","Sample.name"]
labels_nonres <- lesion_info[lesion_info$Response.status..R.responder..NR.non.responder == "NR","Sample.name"]


# find top DA cells
da_cells <- getDAcells(
  X = diffuse.res$X, do.diffuse = F, 
  cell.labels = data_S@meta.data$lesion, 
  labels.1 = labels_res, 
  labels.2 = labels_nonres, 
  k.vector = k.vector, plot.embedding = tsne_embedding
)
da_cells$da.cells.plot


# cluster DA cells to get DA regions
da_regions <- getDAregion(
  X = pca_embedding, cell.idx = da_cells$da.cell.idx, k = 5, alpha = 0.1, 
  cell.labels = data_S@meta.data$lesion,
  labels.1 = labels_res, 
  labels.2 = labels_nonres, 
  plot.embedding = tsne_embedding
)
da_regions$da.region.plot
da_regions$DA.stat



## Marker detection

# set ident in Seurat
n.da <- length(unique(da_regions$cluster.res)) - 1
data_S@meta.data$da <- 0
data_S@meta.data$da[da_cells$da.cell.idx] <- da_regions$cluster.res
data_S <- SetAllIdent(data_S, id = "da")
TSNEPlot(data_S, pt.size = 0.5)


# identify markers for each DA region
da.markers <- list()
for(i in 1:n.da){
  da.markers[[i]] <- FindMarkers(data_S, ident.1 = i, only.pos = T, min.pct = 0.1, min.diff.pct = 0.09)
  da.markers[[i]]$pct.diff <- da.markers[[i]]$pct.1 - da.markers[[i]]$pct.2
}


# plot some markers
marker_genes <- c(
  "CD19","MS4A1","IGHM","CD79A",
  "CD14","CD33","CSF3R","AIF1",
  "VCAM1","LAG3","CD27","CD38",
  "CCR7","LEF1","SELL","ACTN1",
  "IL7R","TCF7","CD8A","CCL5"
)
DotPlot(data_S, genes.plot = marker_genes, cols.use = c("gray","red"), group.by = "da", x.lab.rot = T)


