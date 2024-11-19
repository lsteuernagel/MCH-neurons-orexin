
##########
### Load object as used in anlyses & app
##########

library(Seurat)
library(tidyverse)
library(ggplot2)

# use the output of the prcoessing script or the object from GEO and adjust the file path
data_dir ="/beegfs/scratch/bruening_scratch/lsteuernagel/data/2021-12-maggie-analysis/"
mch_seurat = readRDS(paste0(data_dir,"mch_seurat_clean.rds"))

# define result dir
figure_dir = "figures/"
dir.create(figure_dir,showWarnings = F)

##########
### MCH DATASET general overview
##########

general_pt_size = 0.7
general_label_size = 7

## extra colors:
short_palette = as.character(palette.colors(palette = "Okabe-Ito"))
short_palette = short_palette[!short_palette %in% c("#999999","#000000")]
getOkabeItoPalette = colorRampPalette(short_palette)

cluster_umap = DimPlot(mch_seurat,group.by = "mch_clusters",cols=getOkabeItoPalette(length(unique(mch_seurat@meta.data$mch_clusters))),
                       label = TRUE,label.size = general_label_size,repel = F,pt.size = general_pt_size)+ggtitle(NULL)
cluster_umap

cluster_umap_nolabel = DimPlot(mch_seurat,group.by = "mch_clusters",cols=getOkabeItoPalette(length(unique(mch_seurat@meta.data$mch_clusters))),label = F)+ggtitle(NULL)
cluster_umap_nolabel
## SAVE 
ggsave(filename = paste0(figure_dir,"FIG_A_UMAP_Clusters",".pdf"),plot = cluster_umap, "pdf",dpi=300,width=240,height = 220,units="mm")#
ggsave(filename = paste0(figure_dir,"FIG_A_UMAP_Clusters_nolabel",".pdf"),plot = cluster_umap_nolabel, "pdf",dpi=300,width=240,height = 220,units="mm")#

##########
### Run markers between all clusters
##########

Idents(mch_seurat) = "mch_clusters"
markers_mch_clusters = FindAllMarkers(mch_seurat,min.pct = 0.1,min.diff.pct = 0.05,logfc.threshold = 0.3)
markers_mch_clusters$pct_diff = markers_mch_clusters$pct.1 - markers_mch_clusters$pct.2

## table for supplement:
markers_mch_clusters_filtered_print = markers_mch_clusters %>% dplyr::select(gene,cluster,avg_log2FC,pct.1,pct.2,pct_diff,p_val,p_val_adj) %>%
  dplyr::filter((pct_diff >= 0.1 & avg_log2FC >= 0.333 ) & p_val_adj < 0.001) %>%
  dplyr::arrange(cluster,desc(avg_log2FC))
data.table::fwrite(markers_mch_clusters_filtered_print,paste0(figure_dir,"TABLE_FIG_A_markers_all_clusters",".txt"),sep="\t")
WriteXLS::WriteXLS(x = markers_mch_clusters_filtered_print,ExcelFileName = paste0(figure_dir,"TABLE_FIG_A_markers_all_clusters",".xlsx"),SheetNames = "markers_mch_clusters")

##########
### Orexin receptor Average expression and pct -- Figure
##########

mch_seurat@meta.data$Hcrtr2_expr = FetchData(mch_seurat,"Hcrtr2")[,1]
nrow(mch_seurat@meta.data[mch_seurat@meta.data$Hcrtr2_expr > 0,])

a1 = t(Seurat::AverageExpression(mch_seurat,group.by = "mch_clusters",features = c("Hcrtr1","Hcrtr2"))$RNA) %>% as.data.frame()
a2 = scUtils::gene_pct_cluster(mch_seurat,col_name ="mch_clusters",genes = c("Hcrtr1","Hcrtr2") )

orx1_color = "#FEEB11"
orx2_color = "#FD7170"
both_color="#A900AD"

a1b = a1 %>%
  dplyr::mutate(cluster = rownames(a1), ratio = log2(Hcrtr1 / Hcrtr2)) %>%
  tidyr::gather(key="gene",value="AvgExpr",-cluster,-ratio) %>%
  dplyr::arrange(desc(ratio))
a1b$cluster = factor(a1b$cluster,levels = unique(a1b$cluster))
orxr_barplot = ggplot(a1b,aes(x=cluster,y=AvgExpr,group=gene,fill=gene))+
  geom_col(position="dodge",color="grey70")+
  cowplot::theme_cowplot()+
  theme(text = element_text(size=20),axis.text.x = element_text(size=18,angle=90),axis.title.x = element_blank())+
  scale_fill_manual(values = c("Hcrtr1"=orx1_color,"Hcrtr2"=orx2_color))
#geom_rect()
orxr_barplot

## SAVE 
ggsave(filename = paste0(figure_dir,"FIG_C_ORXR_barplot",".pdf"),plot = orxr_barplot, "pdf",dpi=300,width=300,height = 200,units="mm")#

# then define:
mch_seurat@meta.data$orx_types = "ORX1+/ORX2-"
#mch_seurat@meta.data$orx_types[mch_seurat@meta.data$mch_clusters %in% c("0_Otx1/Parpbp/Tacr3/Cartpt","1_Otx1/Parpbp/Tacr3/Cartpt","12_Otx1/Parpbp/Tacr3/Sim1")] = "ORX1+/ORX2+"
#mch_seurat@meta.data$orx_types[mch_seurat@meta.data$mch_clusters %in% c("10_Otx1/Cpne8/Ano1","8_Otx1/Cpne8/Reln")] = "ORX1-/ORX2+"
# version with just 1 shared cluster
mch_seurat@meta.data$orx_types[mch_seurat@meta.data$mch_clusters %in% c(1,2)] = "ORX1+/ORX2+"
mch_seurat@meta.data$orx_types[mch_seurat@meta.data$mch_clusters %in% c(9,10)] = "ORX1-/ORX2+"

p_colored = DimPlot(mch_seurat,group.by = "orx_types",label = TRUE,pt.size = general_pt_size,label.size = 5)+scale_color_manual(values = c("ORX1+/ORX2+"=both_color,"ORX1+/ORX2-"=orx1_color,"ORX1-/ORX2+"=orx2_color))
plist1 = FeaturePlot(mch_seurat,features =  c("Hcrtr1","Hcrtr2"),pt.size = general_pt_size,order = TRUE,keep.scale = "all",combine = F)
maxval = max(sapply(plist1,function(x) return(max(x$data[,4]))))
hcrtr_umaps = cowplot::plot_grid(plotlist = list(
  hcrtr1 = plist1[[1]]+NoAxes()+scale_color_gradient(low = "lightgrey", high = "blue",limits =  c(0,maxval)),
  hcrtr2 = plist1[[2]]+NoAxes()+scale_color_gradient(low = "lightgrey", high = "blue",limits =  c(0,maxval)),
  p_colored+NoAxes()+ggtitle("ORXR Types")
),nrow = 1)
hcrtr_umaps

## SAVE 
ggsave(filename = paste0(figure_dir,"FIG_B_HCRTR_UMAPs",".pdf"),plot = hcrtr_umaps, "pdf",dpi=300,width=520,height = 150,units="mm")#

# or blend
FeaturePlot(mch_seurat,features =  c("Hcrtr1","Hcrtr2"),pt.size = general_pt_size,blend = TRUE,order = TRUE,cols = c("blue",orx2_color),keep.scale = "all")

# calculate actual co-expression
mch_seurat@meta.data$Hcrtr1_Hcrtr2 = scUtils::CalculateMultScore(mch_seurat,features = c("Hcrtr1","Hcrtr2"))
FeaturePlot(mch_seurat,features =  c("Hcrtr1_Hcrtr2"),pt.size = general_pt_size,order = TRUE,cols = c("grey80",both_color),keep.scale = "all")

##########
### Supplementary plot for Cartpt/Tacr3
##########

cartpt_tacr3_plot = FeaturePlot(mch_seurat,features =  c("Cartpt","Tacr3"),pt.size = general_pt_size,blend = TRUE,order = TRUE,cols = c("darkred","darkblue"),blend.threshold = 0,keep.scale = "all",combine = F)
cartpt_tacr3_plot[[3]]+NoLegend()
cartpt_tacr3_plot[[4]]

ggsave(filename = paste0(figure_dir,"SUP_B_TACR3_CARTPT_UMAP",".pdf"),plot = cartpt_tacr3_plot[[3]]+NoLegend(), "pdf",dpi=300,width=300,height = 300,units="mm")#
ggsave(filename = paste0(figure_dir,"SUP_B_TACR3_CARTPT_legend",".pdf"),plot = cartpt_tacr3_plot[[4]], "pdf",dpi=300,width=300,height = 300,units="mm")#

plist2 = FeaturePlot(mch_seurat,features =  c("Cartpt","Tacr3"),pt.size = general_pt_size,order = TRUE,keep.scale = "all",combine = F)
maxval = max(sapply(plist2,function(x) return(max(x$data[,4]))))
cartptTacr3_umaps = cowplot::plot_grid(plotlist = list(
  hcrtr1 = plist2[[1]]+NoAxes()+scale_color_gradient(low = "lightgrey", high = "blue",limits =  c(0,maxval)),
  hcrtr2 = plist2[[2]]+NoAxes()+scale_color_gradient(low = "lightgrey", high = "blue",limits =  c(0,maxval))
),nrow = 1)
cartptTacr3_umaps

ggsave(filename = paste0(figure_dir,"SUP_B_TACR3_CARTPT_sideByside",".pdf"),plot = cartptTacr3_umaps, "pdf",dpi=430,width=350,height = 150,units="mm")#

##########
### Supplementary plot for Slc17a6 / Slc32a1
##########

# percentages
mch_seurat@meta.data$Slc17a6 = FetchData(mch_seurat,vars = "Slc17a6")[,1]
mch_seurat@meta.data$Slc32a1 = FetchData(mch_seurat,vars = "Slc32a1")[,1]

mch_seurat@meta.data$Slc17a6[mch_seurat@meta.data$Slc17a6 > 0 ] = 1
mch_seurat@meta.data$Slc32a1[mch_seurat@meta.data$Slc32a1 > 0 ] = 1
sum(mch_seurat@meta.data$Slc17a6) / nrow(mch_seurat@meta.data) *100 # 34.62 % Slc17a6
sum(mch_seurat@meta.data$Slc32a1) / nrow(mch_seurat@meta.data) *100 # 1.24 % Slc32a1

plist2 = FeaturePlot(mch_seurat,features =  c("Slc17a6","Slc32a1"),pt.size = general_pt_size,order = TRUE,keep.scale = "all",combine = F)
maxval = max(sapply(plist2,function(x) return(max(x$data[,4]))))
vgatglut_umaps = cowplot::plot_grid(plotlist = list(
  hcrtr1 = plist2[[1]]+NoAxes()+scale_color_gradient(low = "lightgrey", high = "blue",limits =  c(0,maxval)),
  hcrtr2 = plist2[[2]]+NoAxes()+scale_color_gradient(low = "lightgrey", high = "blue",limits =  c(0,maxval))
),nrow = 1)
vgatglut_umaps

ggsave(filename = paste0(figure_dir,"SUP_B_vGLUT_VGAT_sideByside",".pdf"),plot = vgatglut_umaps, "pdf",dpi=430,width=350,height = 150,units="mm")#

##########
### Additional plot with co expression
##########

# not included currently

a1 = t(Seurat::AverageExpression(mch_seurat,group.by = "mch_clusters",features = c("Hcrtr1","Hcrtr2","Hcrtr1_Hcrtr2"))$RNA) %>% as.data.frame()
Hcrtr1_Hcrtr2_column = mch_seurat@meta.data %>% dplyr::mutate(Hcrtr1_Hcrtr2_expm1 = expm1(Hcrtr1_Hcrtr2)) %>%
  dplyr::group_by(mch_clusters) %>% 
  dplyr::summarise(avg = (mean(Hcrtr1_Hcrtr2_expm1))) %>% as.data.frame() # not using log1p because avgexpression does not seem to use it either !?
colnames(Hcrtr1_Hcrtr2_column)[2] = "Hcrtr1_Hcrtr2"
Hcrtr1_Hcrtr2_column = Hcrtr1_Hcrtr2_column[match(rownames(a1),Hcrtr1_Hcrtr2_column$mch_clusters),]
a1 = cbind(a1,Hcrtr1_Hcrtr2 = Hcrtr1_Hcrtr2_column[,2])
a2 = scUtils::gene_pct_cluster(mch_seurat,col_name ="mch_clusters",genes = c("Hcrtr1","Hcrtr2") )

orx1_color = "#FEEB11"
orx2_color = "#FD7170"
both_color="#A900AD"

a1b = a1 %>%
  dplyr::mutate(cluster = rownames(a1), ratio = log2(Hcrtr1 / Hcrtr2)) %>%
  tidyr::gather(key="gene",value="AvgExpr",-cluster,-ratio) %>%
  dplyr::arrange(desc(ratio))
a1b$cluster = factor(a1b$cluster,levels = unique(a1b$cluster))
orxr_barplot = ggplot(a1b,aes(x=cluster,y=AvgExpr,group=gene,fill=gene))+
  geom_col(position="dodge",color="grey70")+
  cowplot::theme_cowplot()+
  theme(text = element_text(size=20),axis.text.x = element_text(size=15,angle=90),axis.title.x = element_blank())+
  scale_fill_manual(values = c("Hcrtr1"=orx1_color,"Hcrtr2"=orx2_color,"Hcrtr1_Hcrtr2"=both_color))
#geom_rect()
orxr_barplot

##########
### Male femal distribution
##########

# how many of these neurons are male female
orx_type_counts = mch_seurat@meta.data %>%
  dplyr::group_by(sex) %>% dplyr::add_count(name="n_sex") %>%
  dplyr::group_by(orx_types,sex) %>% dplyr::add_count(name="n_sex_type") %>%
  dplyr::distinct(orx_types,sex,n_sex,n_sex_type) %>%
  dplyr::mutate(pct = n_sex_type / n_sex * 100) %>%
  dplyr::arrange(desc(pct))
orx_type_counts

orx_type_counts2 = orx_type_counts %>% 
  dplyr::group_by(sex) %>%
  mutate(ypos = cumsum(pct)-0.5*pct, 
         orx_types = factor(orx_types, levels=rev(unique(orx_types[order(-(pct))], ordered=TRUE))))

# orx_type_counts$pct_hack = orx_type_counts$pct
# orx_type_counts$pct_hack[orx_type_counts$orx_types=="ORX1+/ORX2+"] = -1*orx_type_counts$pct_hack[orx_type_counts$orx_types=="ORX1+/ORX2+"]

# pie chart with sex distribution
sex_dist_piechart = ggplot(orx_type_counts2, aes(x="", y=pct, fill=orx_types)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  geom_text(aes(y=ypos,label=paste0(round(pct,2),"%")),size=5) + 
  facet_wrap(facets = ~sex)+
  scale_fill_manual(values = c("ORX1+/ORX2+"=both_color,"ORX1+/ORX2-"=orx1_color,"ORX1-/ORX2+"=orx2_color))+
  cowplot::theme_cowplot()+NoAxes()+
  theme(text = element_text(size=15))
# theme(text = element_text(size=20),axis.text.x = element_text(size=15,angle=90),axis.title.x = element_blank())+

sex_dist_piechart

## SAVE 
ggsave(filename = paste0(figure_dir,"SUP_A_SexDist_Piechart",".pdf"),plot = sex_dist_piechart, "pdf",dpi=200,width=300,height = 150,units="mm")#

##########
### Markers and violin plot (general --> (B))
##########

Idents(mch_seurat) = "orx_types"
markers_orxr_types = FindMarkers(mch_seurat,ident.1 = "ORX1+/ORX2-",ident.2 = "ORX1-/ORX2+",min.pct = 0.1,min.diff.pct = 0.05,logfc.threshold = 0.3)
markers_orxr_types$gene = rownames(markers_orxr_types)
markers_orxr_types$pct_diff = markers_orxr_types$pct.1 - markers_orxr_types$pct.2

markers_orxr_types_filtered = markers_orxr_types %>%
  dplyr::filter((pct_diff > 0 & avg_log2FC>0) | (pct_diff < 0 & avg_log2FC < 0) & p_val_adj < 0.001) %>%
  dplyr::arrange(desc(avg_log2FC))

n_top = 10
Idents(mch_seurat) = "orx_types"
v1 = VlnPlot(mch_seurat,features = c(head(markers_orxr_types_filtered$gene,n = n_top),tail(markers_orxr_types_filtered$gene,n=n_top)),stack = TRUE,fill.by = "ident")+
  scale_fill_manual(values =c("ORX1+/ORX2+"=both_color,"ORX1+/ORX2-"=orx1_color,"ORX1-/ORX2+"=orx2_color))+ylab(NULL)
v1

## SAVE 
ggsave(filename = paste0(figure_dir,"FIG_D_Violin_top10Markers_all",".pdf"),plot = v1, "pdf",dpi=200,width=450,height = 150,units="mm")#



##########
### Re-run markers between the three orxr2+ subpopulations And show those in a violinplot similar to (B)
##########

# define clusters
orxr2_clusters = c(1,2,9,10)
cells = mch_seurat@meta.data$Cell_ID[mch_seurat@meta.data$mch_clusters %in% orxr2_clusters]
mch_seurat_orxr2 = subset(mch_seurat,cells = cells)

mch_seurat_orxr2@meta.data$mch_clusters = factor(mch_seurat_orxr2@meta.data$mch_clusters,levels=rev(c(1,2,9,10)))

# run marker genes
Idents(mch_seurat_orxr2) = "mch_clusters"
markers_orxr2_clusters = FindAllMarkers(mch_seurat_orxr2,min.pct = 0.1,min.diff.pct = 0.05,logfc.threshold = 0.3)
markers_orxr2_clusters$pct_diff = markers_orxr2_clusters$pct.1 - markers_orxr2_clusters$pct.2


markers_orxr2_clusters_filtered = markers_orxr2_clusters %>%
  dplyr::filter((abs(pct_diff) > 0.2 & (avg_log2FC ) > 0.5) & p_val_adj < 0.001 & pct.2 < 0.5) %>%
  dplyr::arrange(desc(avg_log2FC))

n_top = 12
top_markers_per_cluster = markers_orxr2_clusters_filtered %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = desc(p_val_adj),n = n_top,with_ties = F)

# Idents(mch_seurat_orxr2) = "mch_clusters"
# v2 = VlnPlot(mch_seurat_orxr2,features = unique(c(top_markers_per_cluster_up$gene,top_markers_per_cluster_down$gene)),stack = TRUE,fill.by = "ident")+
#   scale_fill_manual(values =plot_data_color_vec)+
#   ylab(NULL)
# v2

# get color vector from figure A!!
g <- ggplot_build(cluster_umap)
plot_data = cbind(cluster_umap$data,g$data[[1]])
plot_data_color = plot_data %>% dplyr::distinct(mch_clusters,colour)
plot_data_color_vec = plot_data_color$colour
names(plot_data_color_vec) = plot_data_color$mch_clusters

Idents(mch_seurat) = "mch_clusters"
v2 = VlnPlot(mch_seurat_orxr2,features = unique(c(top_markers_per_cluster$gene)),stack = TRUE,fill.by = "ident")+
  scale_fill_manual(values =plot_data_color_vec)+
  ylab(NULL)
v2

## SAVE 
ggsave(filename = paste0(figure_dir,"FIG_D_Violin_top12Markers_orxr2_clusters",".pdf"),plot = v2, "pdf",dpi=200,width=450,height = 150,units="mm")#


# heatmap
markers_orxr2_clusters_heat = markers_orxr2_clusters_filtered %>%
  dplyr::arrange(cluster,desc(avg_log2FC))
heatmap_markers = Seurat::DoHeatmap(mch_seurat_orxr2,features = unique(markers_orxr2_clusters_heat$gene),group.by = "mch_clusters",slot = "data",group.colors = plot_data_color_vec)+
  scale_fill_viridis_c()+theme(axis.text.y = element_blank())

ggsave(filename = paste0(figure_dir,"FIG_D_Heatmap_Markers_orxr2_clusters",".pdf"),plot = heatmap_markers, "pdf",dpi=300,width=300,height = 250,units="mm")#


## table for supplement:
markers_orxr2_clusters_filtered_print = markers_orxr2_clusters %>% dplyr::select(gene,cluster,avg_log2FC,pct.1,pct.2,pct_diff,p_val,p_val_adj) %>%
  dplyr::filter((pct_diff >= 0.1 & avg_log2FC >= 0.5 ) & p_val_adj < 0.001) %>%
  dplyr::arrange(cluster,desc(avg_log2FC))
data.table::fwrite(markers_orxr2_clusters_filtered_print,paste0(figure_dir,"TABLE_FIG_D_markers_orxr2_clusters",".txt"),sep="\t")
WriteXLS::WriteXLS(x = markers_orxr2_clusters_filtered_print,ExcelFileName = paste0(figure_dir,"TABLE_FIG_D_markers_orxr2_clusters",".xlsx"),SheetNames = "Markers_ORXR2_clusters")

##########
### Gi and Gq analysis
##########

# have a look at gi and gq associated genes

## get expressed genes
all_genes = data.frame(gene = rownames(mch_seurat_orxr2@assays$RNA@data), expression_sum = rowSums(mch_seurat_orxr2@assays$RNA@data))
expressed_genes = rowSums(mch_seurat_orxr2@assays$RNA@data) 
expressed_genes = expressed_genes[expressed_genes > 0]

# I manually make a list based on: https://en.wikipedia.org/wiki/Heterotrimeric_G_protein
# adding RGSs based on https://www.cell.com/cell/pdf/S0092-8674(20)31140-5.pdf
gai_pathway = stringr::str_to_title(c("GNAO1", "GNAI1", "GNAI2", "GNAI3", "GNAT1", "GNAT2", "GNAT3", "GNAZ", "Pde6a","Pde6b","Pde6c","Pde6d","Pde6g","Pde6h","Adcy1","Adcy5","Adcy6"))
gai_pathway_with_rgs = stringr::str_to_title(c(gai_pathway,c("RGS14","RGS12","RGS10","RGS18","RGS3","RGS4","RGS21","RGS16","RGS1","RGS8","RGS20","RGS17","RGS19","RGS9","RGS11","RGS7","RGS6","KCNJ3","KCNJ6","KCNJ9")))
gaqs_pathway = stringr::str_to_title(c("GNAQ", "GNA11", "GNA14", "GNA15","GNAS","GNAL", "Plce1","Plcb1","Plcb2","Plcb3","Plcb4","Adcy2","Adcy4","Adcy9"))
gaqs_pathway_with_rgs = stringr::str_to_title(c(gaqs_pathway,c("RGS13","RGS1","RGS5","RGS8","RGS19","RGS16","RGS18","RGS17","RGS2")))

v3 = VlnPlot(mch_seurat_orxr2,features = gai_pathway_with_rgs[gai_pathway_with_rgs %in% names(expressed_genes)],stack = TRUE,fill.by = "ident")+
  scale_fill_manual(values =plot_data_color_vec)+
  ylab(NULL)+ggtitle("gai_pathway")
v3
ggsave(filename = paste0(figure_dir,"SUP_C_gai_pathway_violin",".pdf"),plot = v3, "pdf",dpi=200,width=550,height = 150,units="mm")#


v4 = VlnPlot(mch_seurat_orxr2,features = gaqs_pathway_with_rgs[gaqs_pathway_with_rgs %in% names(expressed_genes)],stack = TRUE,fill.by = "ident")+
  scale_fill_manual(values =plot_data_color_vec)+
  ylab(NULL)+ggtitle("gaqs_pathway")
v4
ggsave(filename = paste0(figure_dir,"SUP_C_gaqs_pathway_violin",".pdf"),plot = v4, "pdf",dpi=200,width=450,height = 150,units="mm")#

# gi,gq marker gene stats:
Idents(mch_seurat_orxr2) = "mch_clusters"
markers_gai_pathway = FindAllMarkers(mch_seurat_orxr2,features = gai_pathway_with_rgs[gai_pathway_with_rgs %in% all_genes$gene],min.pct = 0,logfc.threshold = 0)
markers_gaqs_pathway = FindAllMarkers(mch_seurat_orxr2,features = gaqs_pathway_with_rgs[gaqs_pathway_with_rgs %in% all_genes$gene],min.pct = 0,logfc.threshold = 0)

data.table::fwrite(markers_gai_pathway,paste0(figure_dir,"TABLE_SUP_C_gai_between_orxr2_clusters",".txt"),sep="\t")
data.table::fwrite(markers_gaqs_pathway,paste0(figure_dir,"TABLE_SUP_C_gaqs_between_orxr2_clusters",".txt"),sep="\t")


WriteXLS::WriteXLS(x = list(markers_gai_pathway,markers_gaqs_pathway),
                   ExcelFileName = paste0(figure_dir,"TABLE_SUP_C_ga_between_orxr2_clusters",".xlsx"),
                   SheetNames = c("gai_between_orxr2_clusters","aqs_between_orxr2_clusters"))



