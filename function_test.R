library(Seurat)
library(tidyverse)

bp_genes = read.delim("broad_peak_anno_genes.txt", header = FALSE)
bp_dmv_genes = read.delim("broad_peak_anno_dmv_genes.txt", header = FALSE)
top95_genes = read.delim("top95_genes.txt", header = FALSE)

sample_sheet = read.table("cerebellum_cell_metadata.tsv", sep = "\t", header = TRUE)
files = list.files("~/Documents/data/")
files = files[-match("No_purkinje", files)]

seurat_processing = function(folder, sample_sheet) {
  matrix = Read10X(paste0("~/Documents/data/",folder))
  seu_object = CreateSeuratObject(raw.data = matrix, min.cells = 3, min.genes = 200, 
                                  project = folder)
  seu_object_meta = sample_sheet %>% filter(sample_id==folder)
  seu_object_barcode = seu_object_meta %>% select(cell_type)
  rownames(seu_object_barcode) = substr(seu_object_meta[,"sample_name"], 6, 19)
  seu_object = AddMetaData(object = seu_object,
                           metadata = seu_object_barcode,
                           col.name = 'cell_type')
  seu_object = subset(x = seu_object, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mito < 0.05)
  seu_object =  FilterCells(object = seu_object, subset.names = "nGene", 
                            low.thresholds = 300, high.thresholds = 5000)
  seu_object = NormalizeData(object = seu_object , 
                             normalization.method = "LogNormalize", 
                             scale.factor = 10000)
  seu_object_purk = SubsetData(object = seu_object, 
                               cells.use = seu_object@meta.data$cell_type=="Purkinje")
  seu_object_purk_df = as.matrix(seu_object_purk@data)
  seu_object_purk_df = as.data.frame(seu_object_purk_df) %>% rownames_to_column()
  # calculate vari
  bp_dmv_genes_var = e14a_purk_df %>% filter(rowname %in% bp_dmv_genes$V1) %>%
    select_if(is.numeric) %>% apply(1,var)
  
  bp_genes_var = e14a_purk_df %>% filter(rowname %in% bp_genes$V1) %>%
    select_if(is.numeric) %>% apply(1,var)
  
  top95_genes_var = e14a_purk_df %>% filter(rowname %in% top95_genes$V1) %>%
    select_if(is.numeric) %>% apply(1,var)
  return(seu_object_purk_df)
}

test = seurat_processing("E14B", sample_sheet)


