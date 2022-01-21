library(Seurat)
library(infercnv)
setwd("/media/inspur/AS2150G2/LDK/sclc_metastasis/")
getwd()

# load data
load("data_out/2.SCLC.filter.rdata")
ls()

table(SCLC.filter$celltype, SCLC.filter$orig.ident)
table(Idents(SCLC.filter))

subset(SCLC.filter,idents = "tumor")->SCLC.tumor
SCLC.tumor
table(SCLC.tumor$orig.ident)

### data process-counts_matrix
counts_matrix = GetAssayData(SCLC.filter, slot="counts")
head(counts_matrix)[,1:5]
dim(counts_matrix)
# write.table(counts_matrix,"infercnv_data/infercnv_count.txt",sep = "\t",quote = FALSE)

table(SCLC.filter@meta.data[,"celltype"])
table(SCLC.filter@meta.data[,"orig.ident"])

### data process-cell annotation
as.vector(SCLC.filter@meta.data[,"celltype"])-> cell_annotation_tmp
cell_annotation_tmp[which(cell_annotation_tmp=="tumor")]<-"malignant"
paste(cell_annotation_tmp,SCLC.filter@meta.data[,"orig.ident"],sep="_")-> cell_annotation_tmp2
grep("CD8T|NK|CD4T|Macrophages|Neutrophils|Megakaryocytes|B cells|Monocytes|Plasma|DC", cell_annotation_tmp2)->idx2
cell_annotation_tmp2[idx2]<-"normal"
cell_annotation_tmp2-> cell_annotation
cell_annotation <- data.frame(V1=rownames(SCLC.filter@meta.data), V2=cell_annotation)

head(cell_annotation)
dim(cell_annotation)
table(cell_annotation$V2)

identical(cell_annotation$V1, colnames(counts_matrix))


# write.table(cell_annotation,"infercnv_data/cell_annotation.txt",sep="\t",col.names = F,row.names=F,quote=F)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file="infercnv_data/cell_annotation.txt",
                                    delim="\t",
                                    gene_order_file="data_input/05_gene_order_file_only_autosome.txt",
                                    ref_group_names=c("normal"))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="infercnv_res/output_dir_i3",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T,HMM_type="i3")



## infercnv2--running in tmux session sclc_2 
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="infercnv_res/output_dir_i6",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T,HMM_type="i6",plot_steps=F,analysis_mode="subclusters"
)

