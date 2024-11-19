##########
### Load files into Seurat
##########

require(Seurat) # Seurat version 4 not 5 
library(tidyverse)
require(stringr)
library(mapscvi)
seed=12345

# adjust to you local copy of files from cellranger or GEO
# male : GSM8615373
path_sample_male = "/beegfs/v0/bruening_group/CCG/2021-12-maggie-scseq/TX176_A006850169/cellranger/SID160887/outs/filtered_feature_bc_matrix/"
# female: GSM8615374
path_sample_female = "/beegfs/v0/bruening_group/CCG/2021-12-maggie-scseq/TX176_A006850169/cellranger/SID160889/outs/filtered_feature_bc_matrix/"

# adjust to downloaded version of HypoMap: https://doi.org/10.17863/CAM.87955
hypomap_file = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/hypoMap_publication/hypoMap.rds"

# load male
male_wt_sample_matrix = Seurat::Read10X(path_sample_male)
sample_name="SID160887"
colnames(male_wt_sample_matrix) =  paste0(sample_name,"_",colnames(male_wt_sample_matrix))
male_wt_sample_meta = data.frame(Cell_ID = colnames(male_wt_sample_matrix),sample = sample_name,sex = "male")
rownames(male_wt_sample_meta) = colnames(male_wt_sample_matrix)
male_wt_sample_seurat = CreateSeuratObject(male_wt_sample_matrix,meta.data = male_wt_sample_meta,project = "male_mch")
# load female
female_wt_sample_matrix = Seurat::Read10X(path_sample_female)
sample_name="SID160889"
colnames(female_wt_sample_matrix) =  paste0(sample_name,"_",colnames(female_wt_sample_matrix))
female_wt_sample_meta = data.frame(Cell_ID = colnames(female_wt_sample_matrix),sample = sample_name,sex = "female")
rownames(female_wt_sample_meta) = colnames(female_wt_sample_matrix)
female_wt_sample_seurat = CreateSeuratObject(female_wt_sample_matrix,meta.data = female_wt_sample_meta,project = "female_mch")

# Merge
mch_seurat = merge(male_wt_sample_seurat,y = female_wt_sample_seurat)
mch_seurat # 12442 cells

hist(mch_seurat@meta.data$nCount_RNA[mch_seurat@meta.data$nCount_RNA<2000],breaks = 30)

##########
### basic QC
##########

# filtering
mch_seurat[["percent_mt"]]<- Seurat::PercentageFeatureSet(mch_seurat, pattern = "^mt-")
mch_seurat_process = subset(mch_seurat,subset = 
                      nCount_RNA >= 800 &
                      nFeature_RNA >= 500 &
                      percent_mt <= 5)
mch_seurat_process # 11242 cells

##########
### basic processing
##########

# I rerun processing
mch_seurat_process <- Seurat::NormalizeData(object = mch_seurat_process,  verbose = F)
mch_seurat_process <- Seurat::FindVariableFeatures(object = mch_seurat_process, selection.method = "vst", nfeatures = 2000)
mch_seurat_process <- Seurat::ScaleData(object = mch_seurat_process, verbose = TRUE)
mch_seurat_process <- Seurat::RunPCA(object = mch_seurat_process,npcs = 50, seed.use = seed)
mch_seurat_process <- Seurat::RunUMAP(object = mch_seurat_process,  dims = 1:50, n.neighbors = 20, verbose = F, seed.use = seed)

FeaturePlot(mch_seurat_process,features = c("nCount_RNA","nFeature_RNA","percent_mt","Pmch"),ncol =4,order = TRUE)
DimPlot(mch_seurat_process,group.by = "sex")
FeaturePlot(mch_seurat_process,features="Csf1r",order = TRUE)

##########
### map onto hypomap
##########

### map onto full map
hypoMap = readRDS(hypomap_file)

mch_seurat_mapped = mapscvi::map_new_seurat_hypoMap(query_seurat_object = mch_seurat_process,reference_seurat = hypoMap,label_col = "C286_named",max_epochs=20)

inlcude_clusters = mch_seurat_mapped@meta.data %>% dplyr::group_by(predicted) %>% dplyr::count() %>% dplyr::filter(n >= 10)
mch_seurat_mapped@meta.data$predicted_clean = mch_seurat_mapped@meta.data$predicted
mch_seurat_mapped@meta.data$predicted_clean[! mch_seurat_mapped@meta.data$predicted_clean %in% inlcude_clusters$predicted] =NA

DimPlot(mch_seurat_mapped,reduction = "umap",group.by = "predicted_clean",label = TRUE)+NoLegend()

##########
### Remove cells
##########

# initial clustering
mch_seurat_mapped = FindNeighbors(mch_seurat_mapped,reduction = "pca",k.param = 20)
mch_seurat_mapped = FindClusters(mch_seurat_mapped,resolution = 1)
DimPlot(mch_seurat_mapped,reduction = "umap",group.by = "seurat_clusters",label = TRUE)+NoLegend()

cluster_classification = mch_seurat_mapped@meta.data %>%
  dplyr::group_by(predicted,seurat_clusters) %>% 
  dplyr::count() %>% dplyr::group_by(seurat_clusters) %>% dplyr::slice_max(n,n = 1)
cluster_qc =  mch_seurat_mapped@meta.data %>%
  dplyr::group_by(seurat_clusters) %>% 
  dplyr::summarise(mean_nCount_RNA = mean(nCount_RNA),mean_nFeature_RNA = mean(nFeature_RNA), mean_percent_mt = mean(percent_mt) )
cluster_classification = dplyr::left_join(cluster_classification,cluster_qc,by="seurat_clusters")

remove_clusters = cluster_classification$seurat_clusters[!cluster_classification$predicted %in% c("C286-84: Pmch.GLU-7","C286-2: Otx1.GLU-1") | cluster_classification$mean_nFeature_RNA < 1500]
cells_rm = rownames(mch_seurat_mapped@meta.data)[mch_seurat_mapped@meta.data$seurat_clusters %in% remove_clusters]
cells_keep = rownames(mch_seurat_mapped@meta.data)[!mch_seurat_mapped@meta.data$seurat_clusters %in% remove_clusters]
DimPlot(mch_seurat_mapped,reduction = "umap",group.by = "seurat_clusters",cells.highlight = cells_rm,label = TRUE)+NoLegend()
DimPlot(mch_seurat_mapped,reduction = "umap",group.by = "seurat_clusters",cells.highlight = cells_keep,label = TRUE)+NoLegend()

mch_seurat_subset = subset(mch_seurat_mapped,cells = cells_keep)

##########
### Re-run subset obejct
##########

# I rerun processing
mch_seurat_subset <- Seurat::NormalizeData(object = mch_seurat_subset,  verbose = F)
mch_seurat_subset <- Seurat::FindVariableFeatures(object = mch_seurat_subset, selection.method = "vst", nfeatures = 500)
mch_seurat_subset <- Seurat::ScaleData(object = mch_seurat_subset, verbose = TRUE)
mch_seurat_subset <- Seurat::RunPCA(object = mch_seurat_subset,npcs = 50, seed.use = seed)
mch_seurat_subset <- Seurat::RunUMAP(object = mch_seurat_subset,  dims = 1:50, n.neighbors = 20, verbose = F, seed.use = seed)
# mch clustering
mch_seurat_subset = FindNeighbors(mch_seurat_subset,reduction = "pca",k.param = 20)
mch_seurat_subset = FindClusters(mch_seurat_subset,resolution = 0.7)

DimPlot(mch_seurat_subset,reduction = "umap",group.by = "seurat_clusters",label = TRUE)+NoLegend()

FeaturePlot(mch_seurat_subset,features=c("Hcrtr1","Hcrtr2"),order = TRUE)

##########
### Update cluster numbers
##########

cluster_relabel_df = data.frame(RNA_snn_res.0.7 = unique(mch_seurat_subset@meta.data$RNA_snn_res.0.7))
cluster_relabel_df
cluster_relabel_df$mch_clusters = c(5,1,8,0,10,7,9,6,4,3,2)

temp_meta = dplyr::left_join(mch_seurat_subset@meta.data,cluster_relabel_df,by="RNA_snn_res.0.7")
rownames(temp_meta) = temp_meta$Cell_ID

mch_seurat_subset@meta.data = temp_meta

DimPlot(mch_seurat_subset,reduction = "umap",group.by = "mch_clusters",label = TRUE)+NoLegend()
FeaturePlot(mch_seurat_subset,features=c("Sim1"),order = TRUE)

##########
#### Save

data_dir ="/beegfs/scratch/bruening_scratch/lsteuernagel/data/2021-12-maggie-analysis/"
dir.create(data_dir)
saveRDS(mch_seurat_subset,paste0(data_dir,"mch_seurat_clean.rds"))





