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
  mito.genes = grep(pattern = "^mt-", x = rownames(x = seu_object@data), value = TRUE)
  percent.mito = Matrix::colSums(seu_object@raw.data[mito.genes, ])/Matrix::colSums(seu_object@raw.data)
  seu_object = AddMetaData(object = seu_object, metadata = percent.mito, col.name = "percent.mito")
  seu_object = FilterCells(object = seu_object, subset.names = c("nGene", "percent.mito"),
                           low.thresholds = c(300, -Inf), high.thresholds = c(5000, 0.05))
  seu_object = NormalizeData(object = seu_object , 
                             normalization.method = "LogNormalize", 
                             scale.factor = 10000)
  return(seu_object)
}

timepoint_variance = function(folder, sample_sheet) {
  seu_object = seurat_processing(folder, sample_sheet)
  
  # select only purkinje cells
  seu_object_purk = SubsetData(object = seu_object, 
                               cells.use = seu_object@meta.data$cell_type=="Purkinje")
  seu_object_purk_df = as.matrix(seu_object_purk@data)
  seu_object_purk_df = as.data.frame(seu_object_purk_df) %>% rownames_to_column()
  
  # calculate variance
  bp_dmv_genes_var = seu_object_purk_df %>% filter(rowname %in% bp_dmv_genes$V1) %>%
    select_if(is.numeric) %>% apply(1,var)
  
  bp_genes_var = seu_object_purk_df %>% filter(rowname %in% bp_genes$V1) %>%
    select_if(is.numeric) %>% apply(1,var)
  
  top95_genes_var = seu_object_purk_df %>% filter(rowname %in% top95_genes$V1) %>%
    select_if(is.numeric) %>% apply(1,var)
  
  bound = bind_rows(
    tibble(val = bp_dmv_genes_var, var = 'bp_dmv'),
    tibble(val = bp_genes_var, var = 'bp'),
    tibble(val = top95_genes_var, var = 'top95')) %>% add_column(timepoint = folder)
  
  return(bound)
}

variances = map_dfr(files, ~timepoint_variance(., sample_sheet))

ggplot(variances) + geom_boxplot(aes(var, val)) + facet_wrap(~timepoint)
