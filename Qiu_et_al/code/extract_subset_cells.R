
## load libraries
library(Matrix)

## specify where you want to read the input and save the output
## The input folder should contain the mtx and other folders.
input_path = "path/to/input/folder/"
output_path = "path/to/output/folder/"

## define functions
ExtractDataByExp <- function(df_cell, df_gene, exp_id) {
  experiment_id = paste0("run_", exp_id) 
  print(experiment_id)
  gene_count_tmp = Matrix::readMM(paste0(input_path, "mtx/gene_count.", 
                                         experiment_id, 
                                         ".mtx.gz"))
  df_gene_tmp = read.csv(paste0(input_path, "mtx/gene_annotation.csv.gz"),
                         row.names = 1, as.is = T)
  df_cell_tmp = read.csv(paste0(input_path, "mtx/cell_annotation.",
                                experiment_id, ".csv.gz"),
                         row.names = 1, as.is = T)
  rownames(gene_count_tmp) = as.vector(rownames(df_gene_tmp))
  colnames(gene_count_tmp) = as.vector(rownames(df_cell_tmp))
  
  gene_count = gene_count_tmp[intersect(rownames(gene_count_tmp), 
                                        as.vector(df_gene$gene_ID)), 
                              intersect(colnames(gene_count_tmp), 
                                        as.vector(df_cell$cell_id)), 
                              drop=FALSE]
  return(gene_count)
}

## load cell annotations
pd_all = readRDS(paste0(input_path, "other/df_cell_graph.rds"))
rownames(pd_all) = as.vector(pd_all$cell_id)

## keep a subset of cell annotations
example_i = "liver"
celltype_include = c("Gut",
                     "Hepatocytes")
pd_sub = pd_all[pd_all$celltype_new %in% celltype_include, ]

## read gene annotations
mouse_gene = read.table(paste0(input_path, "other/mouse.v12.geneID.txt"), header = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)

## keep a subset of gene annotations
mouse_gene_sub = mouse_gene[(mouse_gene$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & 
                              mouse_gene$chr %in% paste0("chr", c(1:19, "M")),]

## keep the selected subset of scTx data (time consuming step)
gene_count = do.call('cbind', lapply(c(4, 13, 14, 15, 16, "17_sub1", "17_sub2", 18, 19, 20,
                                       21, 22, 23, 24, 25, 26, 27, 28), function(i) {
                                         ExtractDataByExp(pd_sub, mouse_gene_sub, i)
                                       }))

## write output to files
Matrix::writeMM(t(gene_count), paste0(output_path, example_i, ".gene_count.mtx"))
write.csv(pd_sub, paste0(output_path, example_i, ".df_cell.csv"))
write.csv(mouse_gene_sub, paste0(output_path, "df_gene.csv"))
