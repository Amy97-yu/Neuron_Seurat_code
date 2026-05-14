# Authored By LiuHs
# 注意必须要安装R 4.2.0
# 之后载入修改后的monocle2 版本;主要修改了Beam.R 等函数
# 不要用devtools::load_all("./pkg/monocle_2.26.0/monocle/") # 载入修改好的monocle
# 用install.packages 安装
# 注意需要安装 Rtools42 for Windows  不然会出现编译错误
# R版本、monocle修改后的包和Rtools42已经附录在pkg文件夹，如需售后帮忙安装可以联系客服


install.packages("./pkg/monocle_2.26.0.tar/monocle/", repos = NULL, type = "source")
#install.packages("Seurat")

#devtools::load_all("./pkg/monocle_2.26.0/monocle/") # 不要用这种方式载入

# 后续可视化部分如果报错了，请(1) 核对包是否安装好 (2) 是否载入分析好的rds (3) 环境是 R 4.2 还是 4.3
################################################################################
# 开始分析
# monocle2 可以进行降维聚类和ordering
# 有无监督的区别在于，如何选择ordering genes
# 半监督在于基于注释好的celltype选择marker基因；
# 无监督，自动聚类后选择maker；或者选择高变基因
# 如果客户已经Seurat中提供了注释的celltype 信息
# 在这里使用的为半监督聚类
# setwd("./NY081_monocle2/")
rm(list = ls())
gc()
options(stringsAsFactors = F)
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(optparse)
  library(monocle)
  library(dplyr)
})
source("./Data/source.R")
option_list <- list(
  make_option(c("-m","--method"),help = "the method using for selecting ordering genes,celltype/cluster/disp",default = "celltype"),
  make_option(c("-r", "--rds"), help = " ", default = "./sce_npc_v4.rds"),
  make_option(c("-o","--out"),help = "out dir",default = "Result_20250921"),
  make_option(c("-c","--core"),help = "Core",default = 1),
  make_option(c("-v","--reverse"),help = "Reverse the biological process ",default = NULL), # NULL FALSE
  make_option(c("-g","--genes"),help = "the gene list",default = "./Data/6.30.gene.txt")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

sce <- readRDS(opt$rds)

if (!dir.exists(opt$out)) {
  dir.create(opt$out)
}
setwd(opt$out)

#data <- as(as.matrix(sce@assays$RNA@counts), 'sparseMatrix')
#data <- as(as.matrix(GetAssayData(sce,assay = "RNA",slot = "counts")), 'sparseMatrix')
#data <- readRDS("../Data/neural stem cell.data_exp.rds")
#1. 从Seurat正确提取数据
#data <- GetAssayData(sce, slot = "counts")

# 2. 检查数据维度
#stopifnot(dim(data)[1] > 0 & dim(data)[2] > 0)

# 3. 验证数据类型
#if (class(data)[1] != "dgCMatrix") {
  #cat("转换为稀疏矩阵...\n")
  #data <- Matrix(as.matrix(data), sparse = TRUE)
#}

#data <- as(as.matrix(data), 'sparseMatrix')
#genes <- readRDS("../Data/neural stem cell.genes.rds")
#pd <- new('AnnotatedDataFrame', data = sce@meta.data)
#fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
#fData <- data.frame(gene_short_name = genes, row.names = genes)
#fd <- new('AnnotatedDataFrame', data = fData)

#cds <- newCellDataSet(data,
                     # phenoData = pd,
                      #featureData = fd,
                      #expressionFamily = negbinomial.size())#从seurat到monocle的数据转换

data <- readRDS("../Data/25.6.30_Neural Progenitor Cells.ZHUANHUAN_V4.rds")

# 2. 获取count矩阵 (保留原始稀疏格式)
counts <- GetAssayData(sce, slot = "counts")
stopifnot(class(counts) == "dgCMatrix")  # 确认稀疏格式

# 3. 创建统一基因标识
gene_names <- rownames(counts)  # 直接从数据获取基因名
genes <- readRDS("../Data/neural stem cell.genes.rds")

# 4. 对齐基因标识（核心修正）
if(!identical(gene_names, genes)) {
  # 策略：使用表达矩阵基因标识覆盖外部基因列表
  warning("基因标识不一致，使用表达矩阵行名对齐")
  genes <- gene_names  # 强制统一
}

# 5. 构建元数据对象（关键步骤）
## 细胞元数据
pd <- new('AnnotatedDataFrame', data = sce@meta.data)

## 基因元数据
fData <- data.frame(
  gene_short_name = genes,  # 必须包含此列
  row.names = genes         # 行名必须与表达矩阵相同
)
fd <- new('AnnotatedDataFrame', data = fData)

# 6. 创建CellDataSet（直接使用稀疏矩阵）
cds <- newCellDataSet(
  counts,  # 原始稀疏矩阵
  phenoData = pd,
  featureData = fd,
  expressionFamily = negbinomial.size()
)

#expressionFamily 需要根据输入数据做出相应的选择
#negbinomial.size()与negbinomial()都是适用于UMI count的输入，而前者计算速度更快，后者运行的更加准确
#tobit() ，适用于截断正态分布(truncated normal distributions)的数据,如FPKM、TPM，monocle处理这类数据会先取log
#gaussianff()，看名字就知道，处理高斯分布的数据，通常是log后的FPKM、TPM，不推荐，因为你所normalization、log的方式未必与monocle相同
# 假设 seurat_obj 是 Seurat 对象
# 提取表达矩阵（使用原始计数）

print("############01_QC##############")
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds, cores = opt$core, relative_expr = TRUE)
cds <- detectGenes(cds, min_expr = 0.1)

cds_expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 0.01*ncol(cds)))
#cds_expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10))
cds <- cds[cds_expressed_genes,]

print("############02_Selecting ordering genes##############")
if (opt$method == "celltype"| opt$method == "NULL") {
  # using celltype
  print("Now using annotation celltype")
  diff_test_res <- differentialGeneTest(cds[cds_expressed_genes,],
                                        fullModelFormulaStr = "~celltype")
  ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
  cds <- setOrderingFilter(cds, ordering_genes)
  #Reversed Graph Embedding降维
  cds <- reduceDimension(cds, max_components = 2, method = "DDRTree")
  cds <- orderCells(cds,reverse = opt$reverse)
} else if (opt$method == "cluster") {
  # using clustering deno
  print("Now deno cluster")
  cds <- setOrderingFilter(cds, ordering_genes = cds_expressed_genes)
  # cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
  cds <- reduceDimension(cds,
                           max_components = 2,
                           norm_method = 'log',
                           num_dim = 3,
                           reduction_method = 'tSNE',
                           verbose = T)
  cds <- clusterCells(cds,
                        rho_threshold = 2,
                        delta_threshold = 4,
                        skip_rho_sigma = T)
  #plot_cell_clusters(cds, color_by = 'as.factor(Cluster)')|
   # plot_cell_clusters(cds, color_by = 'as.factor(celltype)')
  plot_cell_clusters(cds, color_by = 'as.factor(Cluster)')
  ggsave("./plot_cell_cluster.pdf",width = 8,height = 8)
  ggsave("./plot_cell_cluster.png",width = 8,height = 8)
  clustering_DEG_genes <- differentialGeneTest(cds[cds_expressed_genes,],
                         fullModelFormulaStr = '~Cluster',
                         cores = 10)

  cds_ordering_genes <- row.names(subset(clustering_DEG_genes, qval < 0.01))
  cds <- setOrderingFilter(cds,
                             ordering_genes = cds_ordering_genes)

  cds <- reduceDimension(cds, method = 'DDRTree')
  cds <- orderCells(cds,reverse = opt$reverse)
} else if (opt$method == "celltype") {
  # using method disp
  print("using disper genes for ordering")
  disp_table <- dispersionTable(cds)
  disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
  cds <- setOrderingFilter(cds, disp.genes)
  # Reversed Graph Embedding
  cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
  cds <- orderCells(cds,reverse = opt$reverse)
}
saveRDS(cds,"./2025.11.2NPC.rds")

########## 后续可视化部分记得需要反复变更环境，记得载入保存好的rds############
##########
#cds <- readRDS("./Result_0303/cds.rds")
cds <- readRDS("./2025.11.2NPC.rds")
print("############03_plotting_for_trajectory##############")
#plot_cell_trajectory(cds, color_by = "celltype")##Cluster轨迹分布图
#ggsave("./trajectory_by_celltype.pdf",width = 8,height = 8)
#ggsave("./trajectory_by_celltype.png",width = 8,height = 8)
plot_cell_trajectory(cds, color_by = "Pseudotime")##Pseudotime轨迹图
ggsave("./trajectory_by_Pseudotime.pdf",width = 10,height = 8)
ggsave("./trajectory_by_Pseudotime.png",width = 10,height = 8)
plot_cell_trajectory(cds, color_by = "State")##State轨迹图
ggsave("./trajectory_by_State.pdf",width = 10,height = 8)
ggsave("./trajectory_by_State.png",width = 10,height = 8)

plot_ordering_genes(cds)
ggsave("./ordering_genes_distribution.pdf",height = 10,width = 8)
ggsave("./ordering_genes_distribution.png",height = 10,width = 8)

plot_pc_variance_explained(cds, return_all = F)
ggsave("./plot_pc_variance_explained.pdf",height = 10,width = 8)
ggsave("./plot_pc_variance_explained.png",height = 10,width = 8)


plot_cell_trajectory(cds, color_by = "Pseudotime") + facet_wrap(~ celltype, nrow = 2)
ggsave("./trajectory_splitby_celltype.pdf",height = 10,width = 8)
ggsave("./trajectory_splitby_celltype.png",height = 10,width = 8)


plot_cell_trajectory(cds, color_by = "celltype") + facet_wrap(~ celltype, nrow = 2)
ggsave("./trajectory_by_celltype_splitby_celltype.pdf",height = 10,width = 8)
ggsave("./trajectory_by_celltype_splitby_celltype.png",height = 10,width = 8)

plot_cell_trajectory(cds, color_by = "celltype") + facet_wrap(~ group, nrow = 1)
ggsave("./trajectory_splitby_celltype.pdf")
ggsave("./trajectory_splitby_celltype.png")


plot_cell_trajectory(cds, color_by = "Pseudotime") + facet_wrap(~ group, nrow = 1)
ggsave("./trajectory_splitby_group.pdf")
ggsave("./trajectory_splitby_group.png")

print("############04_DEG_along_trajectory##############")
diff_test_res <- differentialGeneTest(cds,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
# Rik
#Rik_gene <- diff_test_res$gene_short_name[str_detect(diff_test_res$gene_short_name,"Rik")]
#index <- which(rownames(diff_test_res) %in% Rik_gene)
#diff_test_res <- diff_test_res[-index,]

diff_test_res[,c("gene_short_name", "pval", "qval")] %>% head()

sig.gene <- row.names(subset(diff_test_res, qval < 0.0001))
#plot_genes_in_pseudotime(cds[sig.gene[1:5]],
#                         color_by = 'celltype',ncol = 5)
#如果显著的太多了
## 只看自己需要的基因
#genes_tmp = c("Ercc5","Cd9")
#plot_genes_in_pseudotime(cds[genes_tmp],
#                         color_by = 'celltype',ncol = 2)

plot_genes_in_pseudotime(cds[c("Nes")],
                         color_by = 'CellType',ncol = 2)

plot_genes_in_pseudotime(cds[c("Gfap")],
                         color_by = 'celltype',ncol = 2)


plot_genes_in_pseudotime(cds[c("1810058I24Rik")],
                         color_by = 'group',ncol = 2)


plot_genes_in_pseudotime(cds[c("Pik3r1")],
                         color_by = 'group',ncol = 2)

plot_genes_in_pseudotime(cds[c("Akt1")],
                         color_by = 'group',ncol = 2)

plot_genes_in_pseudotime(cds[c("Nes","Gfap","Sox2","Pax6","Ascl1", "Nrxn1", "Map2","Gja1", "Aqp4", "Agt","Cspg4", "Olig2", "Vcan")],
                         color_by = 'celltype',ncol = 2)
plot_genes_in_pseudotime(cds[c("Dnm1l", "Fis1", "Mff", "Mief1")],
                         color_by = 'group',ncol = 1)
# 定义您的自定义基因列表
custom_genes <- c("Nes", "Gfap", "Sox2", "Pax6", "Ascl1", "Nrxn1", "Map2",
                  "Gja1", "Aqp4", "Agt", "Cspg4", "Olig2", "Vcan")
custom_genes <- c("Dnm1l", "Fis1", "Mff", "Mief1","")
plot_pseudotime_heatmap2(cds[custom_genes, ],
                         num_clusters = 3,
                         cores = 2,
                         show_rownames = TRUE,
                         return_heatmap = TRUE)





plot_genes_in_pseudotime(cds[c("Gfap")],
                         color_by = 'group',ncol = 2)



plot_genes_in_pseudotime(cds[c("Dnm1l","Fis1","Mff","Mfn1","Mfn2")],
                         color_by = 'group',ncol = 2)

plot_genes_in_pseudotime(cds[c("Sox2")],
                         color_by = 'group', ncol = 1) +
  theme(
    # 调整坐标轴标签大小
    axis.title = element_text(size = 14),
    # 调整坐标轴刻度文字大小
    axis.text = element_text(size = 12),
    # 调整图例文字大小
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    # 调整分面标题大小（每个基因的名称）
    strip.text = element_text(size = 13)
  )

plot_genes_in_pseudotime(cds[c("Fabp7")],
                         color_by = 'group',ncol = 2)



plot_genes_in_pseudotime(cds[sample(sig.gene,4,replace = FALSE)],
                         color_by = 'celltype',ncol = 2)
ggsave("./plot_genes_in_pseudotime.pdf",height = 10,width = 8)
ggsave("./plot_genes_in_pseudotime.png",height = 10,width = 8)

# 这里选择top1000，或者也可以先做功能分析 获得自己想要的功能再选基因
plot_pseudotime_heatmap(cds[sig.gene[1:1000],],
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = F)
ggsave("./plot_pseudotime_heatmap.pdf",height = 8,width = 9)
ggsave("./plot_pseudotime_heatmap.png",height = 8,width = 9)

print("############05_BEAM分析##############")
plot_cell_trajectory(cds,color_by = 'Pseudotime')

saveRDS(sig.gene,"./sig.gene.rds")
write.csv(sig.gene,"./sig.gene.csv") # v2.1
#cds <- readRDS("./cds.rds")


BEAM_res <- BEAM(cds,cores = 15,branch_point = 1)#查看分支点1两侧的基因表达变化
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
BEAM_gene <- row.names(subset(BEAM_res,qval < 0.01))
length(BEAM_gene)
## [1] 1537
plot_genes_branched_heatmap(cds[BEAM_gene[1:100],],
                            branch_point = 1, # 默认是以第一个分支
                            num_clusters = 2,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)


ggsave("./plot_genes_branched_heatmap.pdf",width = 10,height = 8)
ggsave("./plot_genes_branched_heatmap.png",width = 10,height = 8)

saveRDS(BEAM_res,"./BEAM_res.rds")
saveRDS(BEAM_gene,"./BEAM_gene.rds")
write.csv(BEAM_gene,"./BEAM_gene.csv") # v2.1
################################################################################
print("#################修饰拟时热图##########################################")
# 切换回去了
devtools::load_all("../pkg/ClusterGVis-main/")
#devtools::load_all("./pkg/ClusterGVis-main/ClusterGVis/")
# Note: please update your ComplexHeatmap to the latest version!
# install.packages("devtools")
# devtools::install_github("junjunlab/ClusterGVis")
library(ClusterGVis)
library(monocle)
library(widgetTools)
library(DynDoc)
# return plot
#setwd("./Result_min/")
cds <- readRDS("./Result_0303/cds.rds")
BEAM_res <- readRDS("./Result_0303/BEAM_res.rds")  # 这些sig_gene 自己选择
BEAM_gene <- readRDS("./Result_0303/BEAM_gene.rds")
sig.gene <- readRDS("./Result_0303/sig.gene.rds")
#################################################################
print("######################heatmap_pseudotime_heatmap#######################")

plot_pseudotime_heatmap2(cds[sig.gene[1:100],],
                         num_clusters = 4,
                         cores = 1,
                         show_rownames = T,
                         return_heatmap = T)
Mt.gene <- read.table("./Data/MT.txt")$V1
sig.gene_MT <- intersect(sig.gene,Mt.gene) # 既是差异基因又是线粒体基因


write.table(sig.gene_MT,"./Result_0624/sig.gene_MT.txt")

df <- plot_pseudotime_heatmap2(cds[sig.gene_MT,],
                               num_clusters = 4,
                               cores = 1)
saveRDS(df,"./Result_0624/df.rds")
visCluster(object = df,plot.type = "line",ncol = 2)
ggsave("./plot_pseudo_BEAM_viscluster.pdf")
ggsave("./plot_pseudo_BEAM_viscluster.png")

visCluster(object = df,plot.type = "heatmap")
ggsave("./plot_pseudo_BEAM_viscluster_heatmap.pdf")
ggsave("./plot_pseudo_BEAM_viscluster_heatmap.png")


gene = sample(df$wide.res$gene,20,replace = F)
visCluster(object = df,plot.type = "heatmap",
           markGenes = gene)
ggsave("plot_pseudo_BEAM_viscluster_heatmap_labels.pdf",height = 8,width = 6)
ggsave("plot_pseudo_BEAM_viscluster_heatmap_labels.png",height = 8,width = 6)

visCluster(object = df,plot.type = "both")
ggsave("plot_pseudo_BEAM_viscluster_heatmap_both.pdf",height = 8,width = 6)
ggsave("plot_pseudo_BEAM_viscluster_heatmap_both.png",height = 8,width = 6)

################################################################################
print("######################heatmap_BEAM_heatmap#######################")
#plot_genes_branched_heatmap2(cds[row.names(subset(BEAM_res,qval < 1e-5)),],
#                             branch_point = 1,
#                             num_clusters = 4,
#                             cores = 1,
#                             use_gene_short_name = T,
#                             show_rownames = F,
 #                            return_heatmap = T)

df <- plot_genes_branched_heatmap2(cds[row.names(subset(BEAM_res,qval < 1e-5)),],
                                   branch_point = 1,
                                   num_clusters = 4,
                                   cores = 1,
                                   use_gene_short_name = T,
                                   show_rownames = T)

visCluster(object = df,plot.type = "heatmap")
visCluster(object = df,plot.type = "heatmap",
           pseudotime_col = c("purple","yellow","green"))

pdf(file = "two-branch.pdf",height = 6,width = 7)
visCluster(object = df,plot.type = "both")
dev.off()

saveRDS(df,"./df_two_branched.rds")
#plot_multiple_branches_heatmap2(cds[row.names(BEAM_res)[1:100],],
#                                branches = c(1,2,4),
#                                num_clusters = 4,
#                                cores = 1,
#                                use_gene_short_name = T,
#                                show_rownames = T,
#                                return_heatmap = T)
### 接下来要进行功能分析，但是旧版本的org.Hs.eg.db;org.Mm.eg.db已经过期不适用了
### 所以需要切换为R最新版本
### 请在Tools > Global option 中切换R版本，并且重启Rstudio
# devtools::load_all("../pkg/ClusterGVis-main/")
library(org.Hs.eg.db)
library(org.Mm.eg.db)
df <- readRDS("./df_two_branched.rds")
enrich <- enrichCluster(object = df,
                        OrgDb = org.Mm.eg.db,
                        type = "BP",
                        organism = "mmu",
                        pvalueCutoff = 0.5,
                        topn = 5,
                        seed = 5201314)
markGenes = sample(unique(df$wide.res$gene),25,replace = F)

# PLOT
pdf('branch-enrich.pdf',height = 9,width = 16,onefile = F)
visCluster(object = df,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich,
           go.col = rep(jjAnno::useMyCol("calm",n = 4),each = 5),
           add.bar = T,
           line.side = "left")
dev.off()
################################################################################
# 绘制个性化基因
library(ClusterGVis)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
if (!is.null(opt$genes)) {
  #gl <- read.table("../Data/dif_gene.txt")$V1
  gl <- c("Fabp7", "Map2", "Nes", "Sox2","Olig1")
  #gl <- hsa2mmu(gl) # 如果是小鼠则需要同源转换，否则无需运行
  print(gl)
  cds_beam <- readRDS("./df_two_branched.rds")
  df_pesudotime <- readRDS("./df.rds")
  cds <- readRDS("./cds.rds")
  gene_beam <- readRDS("./BEAM_gene.rds")
  sig.gene <- readRDS("./sig.gene.rds")

  genes <- rownames(cds)[which(rownames(cds) %in% gl)]
  print(genes)
  plot_genes_in_pseudotime(cds[genes[1:4],],
                           color_by = 'celltype',ncol = 2)
  ggsave("./plot_genes_in_pseudotime_selected.pdf",height = 10,width = 8)
  ggsave("./plot_genes_in_pseudotime_selected.png",height = 10,width = 8)

  markGenes = gl[1:10] # 线粒体基因需要label
  visCluster(object = df_pesudotime,plot.type = "heatmap",
             markGenes = markGenes)
  ggsave("./plot_pseudotime_heatmap_selected.pdf",height = 8,width = 9)
  ggsave("./plot_pseudotime_heatmap_selected.png",height = 8,width = 9)

  enrich <- enrichCluster(object = cds_beam,
                          OrgDb = org.Mm.eg.db,
                          type = "BP",
                          organism = "mmu",
                          pvalueCutoff = 0.5,
                          topn = 5,
                          seed = 5201314)
  markGenes = gl[1:10]
  # PLOT

  pdf('branch-enrich_beam_selected_Genes.pdf',height = 9,width = 16,onefile = F)
  visCluster(object = cds_beam,
             plot.type = "both",
             column_names_rot = 45,
             show_row_dend = F,
             markGenes = markGenes,
             markGenes.side = "left",
             annoTerm.data = enrich,
             go.col = rep(jjAnno::useMyCol("calm",n = 4),each = 5),
             add.bar = T,
             line.side = "left")
  dev.off()
}
################################################################################
# 个性化作图 沿轨迹细胞变化
# ggplot个性化修饰monocle2结果作图
print("######################美化轨迹#########################################")
#提取数据=======================================================================
library(tidyverse)
cds <- readRDS("./cds.rds")
data_df <- t(reducedDimS(cds)) %>% as.data.frame() %>% #提取坐标
  select_(Component_1 = 1, Component_2 = 2) %>% #重命名
  rownames_to_column("cells") %>% #rownames命名
  mutate(pData(cds)$State) %>% #添加State
  mutate(pData(cds)$Pseudotime,
         pData(cds)$orig.ident,
         pData(cds)$celltype)#将这些需要作图的有用信息都添加上

colnames(data_df) <- c("cells","Component_1","Component_2","State",
                       "Pseudotime","orig.ident","celltype")
#==============================================================================
#轨迹数据提取---完全摘录于monocle包原函数
dp_mst <- minSpanningTree(cds)
reduced_dim_coords <- reducedDimK(cds)
ica_space_df <- Matrix::t(reduced_dim_coords) %>% as.data.frame() %>%
  select_(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>%
  mutate(sample_name = rownames(.), sample_state = rownames(.))

#构建一个做轨迹线图的数据
edge_df <- dp_mst %>% igraph::as_data_frame() %>%
  select_(source = "from", target = "to") %>%
  left_join(ica_space_df %>% select_(source = "sample_name",
                                     source_prin_graph_dim_1 = "prin_graph_dim_1",
                                     source_prin_graph_dim_2 = "prin_graph_dim_2"), by = "source") %>%
  left_join(ica_space_df %>% select_(target = "sample_name",
                                     target_prin_graph_dim_1 = "prin_graph_dim_1",
                                     target_prin_graph_dim_2 = "prin_graph_dim_2"), by = "target")

#==============================================================================
#计算细胞比例
data_df$orig.ident = data_df$celltype
Cellratio <- prop.table(table(data_df$State, data_df$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c('State',"orig.ident","Freq")
#==============================================================================
#ggplot作图
library(ggplot2)
library(tidydr)
#install.packages("tidydr")
library(ggforce)
#install.packages("ggforce")
library(ggrastr)
#BiocManager::install("ggrastr")
#install.packages("viridis")
library(viridis)
col <- c("#74d2e7","#48a9c5","#0085ad","#8db9ca","#4298b5","#005670","#00205b","#009f4d",
         "#efdf00","#fe5000","#da1884","#a51890","#0077c8","#008eaa")
g <- ggplot() +
  ggrastr::geom_point_rast(data = data_df, aes(x = Component_1,
                                      y = Component_2,
                                      color =Pseudotime)) + #散点图
  viridis::scale_color_viridis()+#密度色
  geom_segment(aes_string(x = "source_prin_graph_dim_1",
                          y = "source_prin_graph_dim_2",
                          xend = "target_prin_graph_dim_1",
                          yend = "target_prin_graph_dim_2"),
               linewidth = 1,
               linetype = "solid", na.rm = TRUE, data = edge_df)+#添加轨迹线
  tidydr::theme_dr(arrow = grid::arrow(length = unit(0, "inches")))+#坐标轴主题修改
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ggforce::geom_arc(arrow = arrow(length = unit(0.15, "inches"), #曲线箭头
                         type = "closed",angle=30),
           aes(x0=0,y0=-3,r=5, start=-0.4*pi, end=0.4*pi),lwd=1)+
  ggforce::geom_arc_bar(data=subset(Cellratio,State=='1'),stat = "pie",#添加饼图
               aes(x0=-15,y0=0,r0=0,r=2.5,amount=Freq,fill=orig.ident))+
  ggforce::geom_arc_bar(data=subset(Cellratio,State=='2'),stat = "pie",
               aes(x0=2,y0=9,r0=0,r=2.5,amount=Freq,fill=orig.ident))+
  ggforce::geom_arc_bar(data=subset(Cellratio,State=='3'),stat = "pie",
               aes(x0=10,y0=-8,r0=0,r=2.5,amount=Freq,fill=orig.ident))+
  scale_fill_manual(values = col)
g
ggsave("./plot_pseudo_branched_cellportiton.pdf",height = 10,width = 10)
ggsave("./plot_pseudo_branched_cellportiton.png",height = 10,width = 10)
################################################################################
print("#################细胞功能随pesudotime变化##############################")
library(dplyr)
library(monocle)
library(ggplot2)
library(Seurat)
library(reshape2)

genes <- c("Opa1","Inf2","Mff","Numb","Fis1","Mief1","Stat5b","Pgam5") # 修改为自己的即可
genes <- read.table("../Data/dif_gene.txt")$V1
genes = "Prdx6"
genes <- rownames(cds)[which(rownames(cds) %in% genes)]
print(genes)

cds_genes <- cds[genes,]
exp = cds_genes@assayData$exprs
#exp <- exprs(cds_genes)
exp <- as.data.frame(exp)
exp <- log2(exp + 1)
exp <- t(exp)

#将上述几个基因的拟时表达添加到monocle
pData(cds) = cbind(pData(cds), exp)

#提取作图数据，只需要游基因表达和拟时即可
data <- pData(cds)
colnames(data)
#选择需要的列即可，我这里的origin.ident就是分组
data <- data[,c("group","Pseudotime",genes)]


#data <- data[which(data$group == "D0" | data$group == "D7"),]
#ggplot作图

#使用分屏，应该就是文献种的办法
#首先将data宽数据转化为长数据
data_long_m<-melt(data, id.vars = c("group", "Pseudotime"), #需保留的不参与聚合的变量列名
                  measure.vars = genes,#选择需要转化的列
                  variable.name = c('gene'),#聚合变量的新列名
                  value.name = 'value')#聚合值的新列名
colnames(data_long_m)

ggplot(data_long_m, aes(x=Pseudotime, y=value, color=group))+
  geom_smooth(aes(fill= group))+ #平滑的填充
  xlab('pseudotime') +
  ylab('Relative expression') +
  facet_wrap(~gene, scales = "free_y")+ #分面，y轴用各自数据
  theme(axis.text = element_text(color = 'black',size = 12),
        axis.title = element_text(color = 'black',size = 14),
        strip.text = element_text(color = 'black',size = 14))+ #分面标题
  scale_color_manual(name=NULL, values = c("#efdf00","#a51890","blue","red"))+#修改颜色
  scale_fill_manual(name=NULL, values = c("#efdf00","#a51890","blue","red"))#修改颜色
ggsave("./4.18Genes_plot_along_pseudotime.pdf",width = 12,height = 9)
ggsave("./4.18Genes_plot_along_pseudotime.png",width = 12,height = 9)
################################################################################
print("#################个性化分析之通路活性分析##############################")
if (!require(homologene)) {
  install.packages("homologene")
  require(homologene)
}

dir = "../Data/Genelist/"
files <- dir(dir)
gene_name <- c()
for (i in files) {
  dir_files <- paste(dir,i,sep = "")
  genelist <- read.table(dir_files,header = FALSE)
  genelist <- genelist$V1

  genes <- homologene::homologene(genelist,inTax = 9606, outTax = 10090)
  genelist <- genes$`10090`

  genelist <- list(genelist)
  name <- strsplit(i,split = ".",fixed = T)[[1]][1]
  sce <- Seurat::AddModuleScore(object = sce,features = genelist,ctrl = 100,pool = rownames(sce),
                                k = F,nbin = 24,name = name)
  names(sce@meta.data)[which(names(sce@meta.data) == paste0(name,"1"))] <- name
  gene_name <- c(gene_name,name)
}

################################作图

#提取作图数据
data1 <- sce@meta.data
data2 <- pData(cds)

colnames(data1)
data1 <- data1[, c("group",gene_name)]
data2 <- data2[,c("Pseudotime","group")]
data_all <- cbind(data1,data2)

data1_long   <- melt(data_all, id.vars = c("group", "Pseudotime"), #需保留的不参与聚合的变量列名
                    measure.vars = 3:4,#选择需要转化的列
                    variable.name = c('pathway'),#聚合变量的新列名
                    value.name = 'value')#聚合值的新列名


data1_new <- data1_long
levels(data1_new$pathway) <- c(gene_name)


#作图和上述一样
ggplot(data1_new, aes(x=Pseudotime, y=value, color= group))+
  geom_smooth(aes(fill= group))+ #平滑的填充
  xlab('pseudotime') +
  ylab('Relative expression') +
  facet_wrap(~levels(pathway), scales = "free_y")+ #分面，y轴用各自数据
  theme(axis.text = element_text(color = 'black',size = 12),
        axis.title = element_text(color = 'black',size = 12),
        strip.text = element_text(color = 'black',size = 14))+ #分面标题
  scale_color_manual(name=NULL, values = c("#efdf00","#a51890","blue","red"))+#修改颜色
  scale_fill_manual(name=NULL, values = c("#efdf00","#a51890","blue","red"))#修改颜色

##如果有安装包的问题和售后问题咨询客服
# 直接使用 monocle2 的内置参数调整热图

