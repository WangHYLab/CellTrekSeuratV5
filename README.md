## Fixes to Tool [CellTrek](https://github.com/navinlabcode/CellTrek) to enable running in SeuratV5

### 1. Install
```R
devtools::install_github("WangHYLab/CellTrekSeuratV5")
```
### 2. Use demo
```R
library(CellTrek)
library(Seurat)
library(CellTrekSeuratV5)
data("sc_data")
data("sc_data")
celltrek.int <- CellTrekSeuratV5::self.traint(  st_data,
                                                sc_data,
                                                st_assay = "Spatial",
                                                sc_assay = "RNA",
                                                norm = "LogNormalize",
                                                nfeatures = 2000,
                                                cell_names = "celltype_fine")

celltrek.out <- CellTrekSeuratV5::self.celltrek(st_sc_int=celltrek.int, int_assay='traint', sc_data=sc_data, sc_assay = 'RNA', 
                                                reduction='pca', intp=T, intp_pnt=5000, intp_lin=F, nPCs=30, ntree=1000, 
                                                dist_thresh=0.55, top_spot=5, spot_n=5, repel_r=20, repel_iter=20, keep_model=T)$celltrek
```
