library(CellTrek)
library(tictoc)
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
source("funcs.R")

## data use ----
# load("data/TLS_self/ST.rdata")
# data<-readRDS("data/TLS_self/data.fine.map.rds")
# 
# st_data<-S2
# st_data<-UpdateSeuratObject(st_data)
# 
# sc_data<-subset(data,orig.ident=="S2")
# sc_data<-UpdateSeuratObject(sc_data)
# 
# sc_data@assays$SCT<-NULL
# sc_data@assays$integrated<-NULL
# DefaultAssay(sc_data)<-"RNA"
# 
# saveRDS(sc_data,"data/TLS_self/sc_data.rds")
# 
# st_data@assays$SCT<-NULL
# saveRDS(st_data,"data/TLS_self/st_data.rds")

sc_data<-readRDS("data/sc_data.rds")
st_data<-readRDS("data/st_data.rds")


## Run ----
set.seed(42)

celltrek.int <- self.traint(st_data=st_data, sc_data=sc_data, sc_assay='RNA', cell_names='celltype_fine')

# CellTrek::celltrek
celltrek.out <- self.celltrek(st_sc_int=celltrek.int, int_assay='traint', sc_data=sc_data, sc_assay = 'RNA', 
                              reduction='pca', intp=T, intp_pnt=5000, intp_lin=F, nPCs=30, ntree=1000, 
                              dist_thresh=0.55, top_spot=5, spot_n=5, repel_r=20, repel_iter=20, keep_model=T)$celltrek


# CellTrek::celltrek_vis(celltrek.out@meta.data %>% dplyr::select(coord_x, coord_y, celltype_fine:id_new),
#                        celltrek.out@images$slice1@image, celltrek.out@images$slice1@scale.factors$lowres)

# p0<- DimPlot(celltrek.int, group.by = "type") 
# p1<-SpatialDimPlot(celltrek.out,group.by = "celltype",crop = F,pt.size.factor = 0.9)
# p2<-SpatialDimPlot(celltrek.out,group.by = "celltype_fine",crop = F)
# p3<-SpatialDimPlot(subset(celltrek.out,celltype=="Immune"),group.by = "celltype_fine",crop = F,pt.size.factor = 1.2)
# ggsave(cowplot::plot_grid(p0,p1,p3,p2,ncol = 2),filename = paste0("out/celltrek.RNA.celltype2.pdf"),width = 20,height = 20)


# print(celltrek.out)

# LayerData(celltrek.out,layer = "counts",assay = "RNA")[1:20,1:20]%>% pheatmap::pheatmap()

# celltrek.out@images$slice1@coordinates %>% plot()
