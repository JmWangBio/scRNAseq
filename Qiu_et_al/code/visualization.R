
## load libraries
library(dplyr)
library(Matrix)
library(ggplot2)

## load intermediate output
example_i = "liver"
work_data_path = "path/to/output/from/subset_scRNAseq_wf.py"

pd_x = read.csv(paste0(work_data_path, example_i, "_adata_scale.obs.csv"), 
                header=T, row.names=1, as.is=T)

nn_matrix = read.csv(paste0(work_data_path, example_i, "_adata_scale.kNN_15.csv"), 
                     as.is=T, header=F)
nn_matrix = as.matrix(nn_matrix)
nn_matrix = nn_matrix + 1 ### python and R using different start index

x = data.frame(i = rep(1:nrow(nn_matrix), ncol(nn_matrix)),
               j = c(nn_matrix), stringsAsFactors = FALSE)

dat = Matrix::sparseMatrix(i = as.numeric(as.vector(x$i)),
                           j = as.numeric(as.vector(x$j)),
                           x = 1)

dat_t = t(dat) + dat
x = data.frame(summary(dat_t))
x = x[x$x == 2 & x$i > x$j,]
x$x = NULL   ### x saves the MNN pairs

y = data.frame(i = 1:nrow(pd_x),
               j = 1:nrow(pd_x),
               meta_group = as.vector(pd_x$meta_group), 
               stringsAsFactors = FALSE)

# check mapping between celltype ID and celltype names
pd_x %>% dplyr::select(celltype_new, meta_group) %>% distinct()

x_rev = x
x_rev$i = x$j
x_rev$j = x$i
x_rev = rbind(x, x_rev)
dat = x_rev %>% left_join(y %>% select(i, meta_group), by = "i") %>%
  left_join(y %>% select(j, meta_group), by = "j") %>%
  filter(meta_group.x != meta_group.y)
dat = dat[dat$meta_group.x == "G_M8" & dat$meta_group.y == "G_M9",]  # G_M8: gut; G_M9: hepatocytes

group = rep("other", nrow(pd_x))
group[c(1:nrow(pd_x)) %in% as.vector(dat$i)] = "group_1"  # gut cells bearing highest resemblance to hepatocytes
group[c(1:nrow(pd_x)) %in% as.vector(dat$j)] = "group_2"  # hepatocytes bearing highest resemblance to gut cells
pd_x$group = as.vector(group)

color_plate = c("#b69340", "#8d70c9")
names(color_plate) = paste0("group_", 1:2)

## UMAP with gut and hepatocytes bearing the highest resemblance to the other cell type highlighted
p1 = ggplot() +
  geom_point(data = pd_x, aes(x = UMAP_2d_1, y = UMAP_2d_2), size=0.6, color = "grey80") +
  geom_point(data = pd_x[pd_x$group != "other",], aes(x = UMAP_2d_1, y = UMAP_2d_2), size=1) +
  geom_point(data = pd_x[pd_x$group != "other",], aes(x = UMAP_2d_1, y = UMAP_2d_2, color = group), size=0.5) +
  theme_void() +
  scale_color_manual(values=color_plate) +
  theme(legend.position="none")

## load cell-level gene counts
work_data_path = "path/to/output/from/extract_subset_cells.R"
gene_count <- Matrix::readMM(paste0(work_data_path, example_i, ".gene_count.mtx"))

## read cell-level metadata
pd_sub <- read.csv(paste0(work_data_path, example_i, ".df_cell.csv"), row.names = 1, as.is = T)

## read gene metadata
mouse_gene_sub <- read.csv(paste0(work_data_path, "df_gene.csv"), row.names = 1, as.is = T)

## assign row and column names
colnames(gene_count) <- rownames(mouse_gene_sub)
rownames(gene_count) <- rownames(pd_sub)

## UMAP highlighting the expression of a gene (e.g., Gata6, Alb)
gene_sym <- "Gata6"

gene_id = mouse_gene_sub[mouse_gene_sub$gene_short_name == gene_sym, "gene_ID"]
gene_count_subset <- gene_count[, gene_id]

p2 = ggplot() +
  geom_point(data = pd_x %>%
               mutate(gene = log10(as.numeric(gene_count_subset))), 
             aes(x = UMAP_2d_1, y = UMAP_2d_2, color = gene), size=0.6) +
  theme_void() + 
  labs(color = "log10(exp)") +
  theme(legend.position="bottom")
