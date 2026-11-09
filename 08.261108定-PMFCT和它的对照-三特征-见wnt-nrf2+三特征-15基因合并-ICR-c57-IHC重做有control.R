#08.XHX构建小鼠以及改造,ICR小鼠-6组共21个样本
setwd("D:/cloudstudy-biotree zengjianming/XHXR/260408xhxmouse")
rm(list=ls())
options(stringsAsFactors = F)
library(dplyr)
library(stringr)
library(data.table)
library(readxl)
library(clusterProfiler)
library(GSEABase)
library(GSVA)
#读入小鼠数据
library(tidyverse)
library(data.table)
library(glue)
library(ggplot2)
library(RColorBrewer)
library(GSVA)
library(pheatmap)
library(clusterProfiler)
library(DESeq2)
library(rstatix)
library(ggpubr)
library(patchwork)
#读入tpm数据,最全
exp<- fread("../HBdataset/data/260305-论文/XHX跑/最新数据备用/292samples/328/260723-最全340samples-public+XHX-mousetoHUMAN-salmon_Symbol_tpm.txt",
            header = T,data.table = F)

row.names(exp) <- exp[, 1]
exp<-exp[,-1]
#过滤基因
exp = exp[apply(exp, 1, function(x) sum(x > 0) > 0.5*ncol(exp)), ]
nrow(exp)
#筛选出已发表组织，肿瘤
clinical<- read.csv("../HBdataset/data/260305-论文/XHX跑/最新数据备用/clinical-newGSVA3-最新newsubtype2-4rd-260626.csv",header = T)
#XHX构的10个模型和对照，减去PMFCT14-D28，FCT1008,YB1008，删除MYC-B-YAP-28
#PMFCT与它的对照们，包括公共数据，26-8-02重做-ICR
#clinical3<-clinical[c(40,42:43,6:8,27:29,32,34:35,102:107,165:167,170,172:173,204,276:278,282:283,233:235,301:304,309:311),]#ICR+两质粒+YB,XHX,172换176
#删除对照组102:107,
clinical3<-clinical[c(40,42:43,6:8,27:28,32,34:35,
                      78,81,83,276,102:107),]#ICR+两质粒+YB,XHX,172换176,加上赵斌的三组
clinical3<-clinical[c(34:35,336:337,339,325:327,170,172:173,165:167
                      ,168:169,174),]#c57小鼠,17个样本

rownames(clinical3) <- clinical3$sampleName


p = identical(colnames(exp),rownames(clinical3));p
if(!p) {
  s = intersect(colnames(exp),rownames(clinical3))
  exp = exp[,s]
  clinical3 = clinical3[s,]
}
gene.expr <-log2(exp+1)%>%as.matrix
#自定义基因集跑GSVA-未监督分类选用的+自定义----
load("D:/cloudstudy-biotree zengjianming/XHXR/HBdataset/data/260305-论文/gmt/geneSets2.RData")
#选择genesets
genesetspar <- geneSets2[c(
  "WONG_EMBRYONIC_STEM_CELL_CORE",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MYC_TARGETS_V2",
  "KEGG_ECM_RECEPTOR_INTERACTION",
  "GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION",
  "GOBP_MESENCHYME_DEVELOPMENT","GOBP_BILE_ACID_METABOLIC_PROCESS",
  "KEGG_DRUG_METABOLISM_CYTOCHROME_P450",
  "REACTOME_YAP1_AND_WWTR1_TAZ_STIMULATED_GENE_EXPRESSION",
  "ERYTHROPOIESIS","FETAL_LIVER","WNT_TARGETS")
  ]
#GSVA
gsvapar <- gsvaParam(gene.expr, genesetspar, maxDiff=TRUE,kcdf="Gaussian")

## 运行GSVA打分
gsva_es <- gsva(gsvapar)

#修改gsva_es行名
# 先提取当前的名称
original_names <- rownames(gsva_es)

# 批量替换：去掉前缀，并将下划线转为空格
short_names <- original_names %>%
  # 1. 去除常见前缀
  gsub("^(HALLMARK|KEGG|GOBP|REACTOME|AIZARANI_LIVER|WONG)_", "", .) %>% 
  # 2. 把剩下的下划线换成空格
  gsub("_", " ", .) %>% 
  # 3. 关键步骤：转为句首大写，其余小写
  str_to_title()

# 将简化后的名称重新赋给矩阵
rownames(gsva_es) <- short_names
class(gsva_es)

save(genesetspar,gsva_es,clinical,exp,gene.expr,file="PMFCT和它的对照们-GSVA-260901-有YB无control-ICR.RData")


#加载数据####
#load("已发表小鼠3个HB模型-XHX构共65samples-GSVA-260723-改造效果.RData")
#直接
gsva_es1<-gsva_es[-c(7,11:13),]#减去wnt和NRF2
#去除基因，用不scale的值
pathway.df <- data.frame(pathway = rownames(gsva_es1))
head(pathway.df)
unique(pathway.df$pathway)
# 1. 先把所有通路初始化为 "Fetal_liverGenes"
pathway.df$colour.group <- "NFE2L2"
# 2. 批量匹配定义 (利用 grepl 的 OR | 逻辑)---少
pathway.df$colour.group[grepl("STEM|G2M", pathway.df$pathway, ignore.case = T)] <- "PRO"
pathway.df$colour.group[grepl("MESENCHY|ECM|MATRIX|TGF|OSSIFICATION|MUSCLE|MYOGENESIS", pathway.df$pathway, ignore.case = T)] <- "MES"
pathway.df$colour.group[grepl("METABOLISM|BILE|DETOXIFICATION", pathway.df$pathway, ignore.case = T)] <- "HEP"
pathway.df$colour.group[grepl("MYC", pathway.df$pathway, ignore.case = T)] <- "MYC"
pathway.df$colour.group[grepl("YAP", pathway.df$pathway, ignore.case = T)] <- "YAP"

#最全注释####
sampleAnno <- setDF(clinical3)
sampleAnno$Pathology[sampleAnno$Pathology == ""] <- "Control"
sampleAnno$MYC[sampleAnno$MYC == ""] <- "No"
sampleAnno$CTNNB1.mut[sampleAnno$CTNNB1.mut == ""] <- "No"
sampleAnno$FGF19.CCND1[sampleAnno$FGF19.CCND1 == ""] <- "No"
sampleAnno$YAP1.mut[sampleAnno$YAP1.mut == ""] <- "No"
sampleAnno$NFE2L2.mut[sampleAnno$NFE2L2.mut == ""] <- "No"
unique(sampleAnno$GSE)
unique(sampleAnno$background)
color3 <- c(
  "ICR"      = "#B8C0C2",   # 淡灰蓝 (Cool Gray)
  "C57BL/6J" = "#3D405B"
)
unique(sampleAnno$treatment5)

color2 <- c(
  "MFCB"= "#EA4335", 
  "FCB"="#34A853","Control"="#4285F4",
 "MFC"               = "#c4a7e7",
 "MB"            = "#FFB5A7",  
  "YB"="#FBBC05"
) 


# 💡 建议在控制台运行下面这行，一秒预览34种颜色在热图里的堆叠色块效果
barplot(rep(1, 6), col = color2, names.arg = 1:6, las = 2, main = "34 Colors Preview")
unique(sampleAnno$group)
group <- c("Tumor"="#DE4599",
  "Normal" = "#A9ABBE"    
)
unique(sampleAnno$Pathology)
color5<-c("HB"="black",
          "HCC"="#D9D9D9",
           "Control"= "#9DEAF6FF")
sampleAnno$Pathology <- factor(
  sampleAnno$Pathology,
  levels = c("HB", "HCC","Control"
             ))#"Normal"


unique(sampleAnno$MYC)
MYC <- c(
  "No"  = "lightgrey",
  "Yes" = "black"
)
unique(sampleAnno$CTNNB1.mut)
CTNNB1.mut <- c(
  "No"        = "lightgrey",
 "S33/S37A"="#F9B982"
)
unique(sampleAnno$FGF19.CCND1)
FGF19.CCND1 <- c(
  "No"  = "lightgrey",
  "Yes" = "black"
)
unique(sampleAnno$YAP1.mut)
YAP1.mut = c(
  "No"    = "lightgrey",
  "S127A" = "black"
)

row_annotation <- data.frame(Pathway = pathway.df$colour.group)
rownames(row_annotation) <- pathway.df$pathway

unique(sampleAnno$treatment5)
sampleAnno$treatment5 <- factor(sampleAnno$treatment5, levels = c("MFCB" ,"MFC" ,"YB","FCB" ,"MB" , 
                                                                  "Control"))

sampleAnno$group <- factor(
  sampleAnno$group,
  levels = c("Tumor", "Normal"))
sampleAnno$GSE <- factor(
  sampleAnno$GSE,
  levels = c("XHX","PRJNA674008"))
sampleAnno$CTNNB1.mut <- factor(
  sampleAnno$CTNNB1.mut,
  levels = c("S33/S37A","No"))
sampleAnno$MYC <- factor(
  sampleAnno$MYC,
  levels = c("Yes","No"))
sampleAnno$FGF19.CCND1 <- factor(
  sampleAnno$FGF19.CCND1,
  levels = c("Yes","No"))
sampleAnno$background <- factor(
  sampleAnno$background,
  levels = c("ICR","C57BL/6J" ))
# 💡 确保加载了高级热图包-第二种作图方式####
require(corto)
library(ComplexHeatmap)
library(circlize)
desired_row_order <- c("PRO", "MES", "HEP","MYC","YAP")
pathway.df$colour.group<-factor(pathway.df$colour.group, levels = desired_row_order)

ha.r <- rowAnnotation(
  Pathway = pathway.df$colour.group, # 行注释的数据源
  col = list(
    Pathway = c(
      "PRO"  = "#962D2D",
      "MES"  = "#EFC000",
      "HEP"  = "#0073C2",
      "MYC" = "#8A4198","YAP"="#00A878"
    )
  ),
  show_annotation_name = FALSE,show_legend = FALSE,
  #annotation_name_gp = gpar(fontsize = 10),
  simple_anno_size = unit(1.5, "mm"))


#修改布局
#修改布局####
lgd_param <- list(title_gp = gpar(fontsize = 10, fontface = "bold"),
                  labels_gp = gpar(fontsize = 10),
                  grid_height = unit(3, "mm"),
                  grid_width = unit(3, "mm"))
lgd_param1 <- list(
  title_gp    = gpar(fontsize = 10, fontface = "bold"),
  labels_gp   = gpar(fontsize = 10),
  grid_height = unit(3, "mm"),
  grid_width  = unit(3, "mm"),
  ncol        = 2  # 设置图例项按 4 列横向排列
)
ha.b <- HeatmapAnnotation(
  df = data.frame(  
    type=sampleAnno$Pathology,
    treatment=sampleAnno$treatment5
  ),
  col = list(
    type=color5,
    treatment=color2 
  ),
  simple_anno_size = unit(3.5, "mm"),
  annotation_name_gp = gpar(fontsize = 10), # 设为 10 号加粗
  annotation_name_side = "right",na_col = "white",
  annotation_legend_param = list(
    type         = lgd_param,
    treatment   = lgd_param1
  )
)
col_factor <- factor(
  sampleAnno$treatment5, 
  levels = c("MFCB" ,"MFC" ,"YB",
             "FCB" ,"MB" , "Control"))
type_factor<-factor(sampleAnno$Pathology, 
levels = c("HB", "HCC","Control"))

# 2. 构建热图
set.seed(123578)
pdf("./results/260910/260901新-PMFCT和它的对照-21samples-GSVA新基因集-第二种热图10p-列名PMCT用赵斌的.pdf", width = 9,height = 5)
ht <-Heatmap(gsva_es1,
             col = colorRamp2(c(-0.8, 0, 0.8), c("blue", "white", "red")),
             column_order = order(type_factor,col_factor), 
             column_split = sampleAnno$Pathology,column_title = NULL,
             cluster_rows = FALSE,     column_gap = unit(1, "mm"),
             show_row_dend = FALSE,
             height = unit(nrow(gsva_es1) * 5, "mm"),
             width = unit(ncol(gsva_es1) * 3, "mm"),
             row_split = pathway.df$colour.group,
             row_title_side = "left",                 # 文字显示在最左侧（色块旁边）
             row_title_rot = 90,                      # 💡 旋转 90 度（如需不旋转横向显示，可改为 0）
             row_title_gp = gpar(
               fontsize = 10, 
               fontface = "bold", 
               col = "black"),
             row_names_gp = gpar(fontsize = 10),
             show_column_names = FALSE,
             row_dend_reorder = TRUE,
             row_gap = unit(1, "mm"),
             border = TRUE,
             name = "GSVA Score",heatmap_legend_param = list(
               title_gp = gpar(fontsize = 10, fontface = "bold"), # 标题字号（如 9 pt）
               labels_gp = gpar(fontsize = 10)                    # 刻度数字字号（如 8 pt）
             ),
             top_annotation = ha.b,
             left_annotation = ha.r
)

draw(ht, 
     heatmap_legend_side = "left",      # 主图例在右侧
     annotation_legend_side = "left",  # 注释图例全部移到底部
     merge_legend = TRUE
)
dev.off()
#鼠源基因表达情况####
exp.mouse<- fread("../HBdataset/data/260305-论文/XHX跑/最新数据备用/292samples/328/260723-最全340samples-public+XHX-mouse-salmon_Symbol_tpm.txt",
                  header = T,data.table = F)
row.names(exp.mouse) <- exp.mouse[, 1]
exp.mouse<-exp.mouse[,-1]
#用clinical3
p = identical(colnames(exp.mouse),rownames(clinical3));p
if(!p) {
  s = intersect(rownames(clinical3),colnames(exp.mouse))
  exp.mouse = exp.mouse[,s]
  clinical3 = clinical3[s,]
}

gene.expr.mouse <-log2(exp.mouse+1)%>%as.matrix
target_genes<-c("Axin2", "Lgr5","Nkd1", "Tnfrsf19", 
                "Gata1","Nfe2","Rhag","Slc4a1", "Trim58",   
                "Dlk1", "Hif3a", "Igsf1" , "Peg3","Trim71")#wnt-干性，红系,"Sfrp1","Apcdd1",
plot_matrix <- gene.expr.mouse[target_genes, ]
wnt_genes      <- c("Nkd1", "Lgr5",  "Axin2", "Tnfrsf19")#,"AXIN2","LGR5"
stemness_genes <- c("Dlk1", "Hif3a", "Trim71", "Peg3", "Igsf1")#"DUSP9",,"TRIM71"

# 剩下的基因全部自动归为红细胞相关
# setdiff 可以自动找出 target_genes 中扣除掉前两组后剩下的所有基因
erythroid_genes <- setdiff(target_genes, c(wnt_genes, stemness_genes))

# 3. 构建行注释数据框 (row_annotation)
# 核心逻辑：先创建一个空的 Data Frame，行名设为基因名
row_annotation <- data.frame(
  Gene_Class = rep(NA, length(target_genes)),
  row.names = target_genes
)

# 4. 批量填入分类标签
row_annotation[wnt_genes, "Gene_Class"]      <- "WNT"
row_annotation[stemness_genes, "Gene_Class"] <- "LP"
row_annotation[erythroid_genes, "Gene_Class"] <- "ERY"
gene_zscore <- t(scale(t(as.matrix(plot_matrix[, colnames(gsva_es1)]))))
#0-1范围
min_max_scale <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# 对表达矩阵的每一行（基因）应用该函数
gene_zscore <- t(apply(plot_matrix, 1, min_max_scale))
quantile(gene_zscore, c(0.01, 0.05, 0.95, 0.99))
row_annotation$Gene_Class <- factor(
  row_annotation$Gene_Class,
  levels = c("WNT", "LP", "ERY"))
ha.r2 <- rowAnnotation(
  Gene_Class = row_annotation$Gene_Class, # 行注释的数据源
  col = list(
    Gene_Class = c(
      "LP"  = "#fbb05b",   
      "WNT"  = "#7bc4e2",   
      "ERY"  = "#ed6ca4"   
    )
  ),
  show_annotation_name = FALSE,show_legend = FALSE,
  simple_anno_size = unit(1.5, "mm")# 是否在行注释上方显示“Pathway”这个标题
)

set.seed(123578)
pdf("./results/260922/260922新-PMFCT和它的对照-ICR-21samples-热图6P+15特征-3.pdf", width = 8, height = 5)
ht_gsva<-Heatmap(gsva_es1,
                 col = colorRamp2(c(-0.8, 0, 0.8),c("#377EB8", "white", "#E41A1C")),
                 column_order = order(type_factor,col_factor), 
                 column_split = sampleAnno$Pathology,column_title = NULL,
                 cluster_rows = FALSE,     column_gap = unit(1, "mm"),
                 show_row_dend = FALSE,
                 height = unit(nrow(gsva_es1) * 5, "mm"),
                 width = unit(ncol(gsva_es1) * 2.5, "mm"),
                 row_split = pathway.df$colour.group,
                 row_title_side = "left",                 # 文字显示在最左侧（色块旁边）
                 row_title_rot = 0,                      # 💡 旋转 90 度（如需不旋转横向显示，可改为 0）
                 row_title_gp = gpar(
                   fontsize = 10, 
                   fontface = "bold", 
                   col = "black"),
                 row_names_gp = gpar(fontsize = 10),
                 show_column_names = FALSE,
                 row_dend_reorder = TRUE,
                 row_gap = unit(1, "mm"),
                 border = TRUE,
                 name = "GSVA Score",heatmap_legend_param = list(
                   direction = "horizontal",title_gp = gpar(fontsize = 10, fontface = "bold"), 
                   labels_gp = gpar(fontsize = 10),legend_width   = unit(2, "cm"),legend_height = unit(1.5, "cm")),
                 top_annotation = ha.b,
                 left_annotation = ha.r
)

ht_gene <-Heatmap(gene_zscore,
                  col = colorRamp2(c( 0,0.25, 0.5,0.75,1),  c("#5687E7","#9BB8D5","#F5F7FA","#DD7B84","#BC3C48")),
                  column_order = order(type_factor,col_factor), 
                  column_split = sampleAnno$Pathology,column_title = NULL,
                  cluster_rows = FALSE,     column_gap = unit(1, "mm"),
                  show_row_names = FALSE, show_row_dend = FALSE,
                  height = unit(nrow(gene_zscore) * 1.5, "mm"),
                  width = unit(ncol(gsva_es1) * 2.5, "mm"),
                  row_split = row_annotation$Gene_Class,
                  row_title_side = "left",                 # 文字显示在最左侧（色块旁边）
                  row_title_rot = 0,                      # 💡 旋转 90 度（如需不旋转横向显示，可改为 0）
                  row_title_gp = gpar(
                    fontsize = 10, 
                    fontface = "bold", 
                    col = "black"),
                  row_names_gp = gpar(fontsize = 10,fontface = "italic"),
                  show_column_names = FALSE,
                  row_dend_reorder = TRUE,
                  row_gap = unit(1, "mm"),
                  border = TRUE,
                  name = "Z-score",heatmap_legend_param = list(
                    direction = "horizontal",
                    title_gp = gpar(fontsize = 10, fontface = "bold"), 
                    labels_gp = gpar(fontsize = 10) ,legend_width = unit(2, "cm"),legend_height = unit(1.5, "cm")) ,
                  left_annotation = ha.r2
)

ht_list <- ht_gsva %v% ht_gene
draw(ht_list,
     heatmap_legend_side = "left",      # 主图例在右侧
     annotation_legend_side = "left",  # 注释图例全部移到底部
     merge_legend = TRUE,
     ht_gap = unit(1, "mm")
)
dev.off()

##读入人源基因的tpm再转化为log2
exp1<-fread("../HBdataset/data/260305-论文/XHX跑/最新数据备用/292samples/328/260723-最全368samples-public+XHX-HUMANgene-salmon_Symbol_tpm.txt",
            header = T,data.table = F)
row.names(exp1)<-exp1$Symbol
exp1<-exp1[,-1]
colnames(exp1)[195]<-"PMFCT156.CL_1"
colnames(exp1)[198]<-"PMFCT156.CL_2"
#
p = identical(rownames(clinical3),colnames(exp1));p
if(!p) {
  s = intersect(colnames(exp1),rownames(clinical3))
  exp1 = exp1[,s]
  clinical3 = clinical3[s,]
}
gene.expr.human <-log2(exp1+1)%>%as.matrix
gene.expr.human<-gene.expr.human[c("CCND1","CTNNB1","FGF19","MYC","YAP1"),]
target_genes<-c("CCND1","CTNNB1","FGF19","MYC","YAP1")#wnt-干性，红系

# 3. 构建行注释数据框 (row_annotation)
# 核心逻辑：先创建一个空的 Data Frame，行名设为基因名
row_annotation3 <- data.frame(
  Gene_Class = rep("Exogenous", length(target_genes)),
  row.names = target_genes
)

gene_zscore2 <- t(scale(t(as.matrix(gene.expr.human[, colnames(gsva_es1)]))))

quantile(gene_zscore2, c(0.01, 0.05, 0.95, 0.99))
ha.r4 <- rowAnnotation(
  Gene_Class = row_annotation3$Gene_Class, # 行注释的数据源
  col = list(
    Gene_Class = c(
      "F/C"  = '#dcc6dc'
    )
  ),
  show_annotation_name = FALSE,show_legend = FALSE,
  simple_anno_size = unit(1.5, "mm")# 是否在行注释上方显示“Pathway”这个标题
)
set.seed(123578)
pdf("./results/260910/260910新-PMFCT和它的对照-ICR-21samples-热图6P+15特征-2-FC.pdf", width = 8, height = 6)
ht_gsva<-Heatmap(gsva_es1,
                 col = colorRamp2(c(-0.8, 0, 0.8),c("#377EB8", "white", "#E41A1C")),
                 column_order = order(type_factor,col_factor), 
                 column_split = sampleAnno$Pathology,column_title = NULL,
                 cluster_rows = FALSE,     column_gap = unit(1, "mm"),
                 show_row_dend = FALSE,
                 height = unit(nrow(gsva_es1) * 5, "mm"),
                 width = unit(ncol(gsva_es1) * 2.5, "mm"),
                 row_split = pathway.df$colour.group,
                 row_title_side = "left",                 # 文字显示在最左侧（色块旁边）
                 row_title_rot = 90,                      # 💡 旋转 90 度（如需不旋转横向显示，可改为 0）
                 row_title_gp = gpar(
                   fontsize = 10, 
                   fontface = "plain", 
                   col = "black"),
                 row_names_gp = gpar(fontsize = 10),
                 show_column_names = FALSE,
                 row_dend_reorder = TRUE,
                 row_gap = unit(1, "mm"),
                 border = TRUE,
                 name = "GSVA Enrichment",heatmap_legend_param = list(
                   direction = "horizontal",title_gp = gpar(fontsize = 10, fontface = "bold"), 
                   labels_gp = gpar(fontsize = 10),legend_width   = unit(2, "cm"),legend_height = unit(1.5, "cm")),
                 top_annotation = ha.b,
                 left_annotation = ha.r
)
ht_gene2 <-Heatmap(gene_zscore2,
                   col = colorRamp2(c(-2, 0, 3),  c("#440154","#21918C" , "#FDE725")),
                   column_order = order(type_factor,col_factor), 
                   column_split = sampleAnno$Pathology,column_title = NULL,
                   cluster_rows = FALSE,     column_gap = unit(1, "mm"),
                   show_row_names = FALSE, show_row_dend = FALSE,
                   height = unit(nrow(gene_zscore) * 1.5, "mm"),
                   width = unit(ncol(gsva_es1) * 2.5, "mm"),
                   row_split = row_annotation3$Gene_Class,
                   row_title_side = "left",                 # 文字显示在最左侧（色块旁边）
                   row_title_rot = 90,                      # 💡 旋转 90 度（如需不旋转横向显示，可改为 0）
                   row_title_gp = gpar(
                     fontsize = 10, 
                     fontface = "plain", 
                     col = "black"),
                   row_names_gp = gpar(fontsize = 10,fontface = "italic"),
                   show_column_names = FALSE,
                   row_dend_reorder = TRUE,
                   row_gap = unit(1, "mm"),
                   border = TRUE,
                   name = "Z_Score",heatmap_legend_param = list(
                     direction = "horizontal",
                     title_gp = gpar(fontsize = 10, fontface = "bold"), 
                     labels_gp = gpar(fontsize = 10) ,legend_height = unit(2, "cm")),  
                   left_annotation = ha.r4
)
ht_list <- ht_gsva %v% ht_gene2
ht_gene <-Heatmap(gene_zscore,
                  col = colorRamp2(c(-2, 0, 3),  c("#440154","#21918C" , "#FDE725")),
                  column_order = order(type_factor,col_factor), 
                  column_split = sampleAnno$Pathology,column_title = NULL,
                  cluster_rows = FALSE,     column_gap = unit(1, "mm"),
                  show_row_names = FALSE, show_row_dend = FALSE,
                  height = unit(nrow(gene_zscore) * 1.5, "mm"),
                  width = unit(ncol(gsva_es1) * 2.5, "mm"),
                  row_split = row_annotation$Gene_Class,
                  row_title_side = "left",                 # 文字显示在最左侧（色块旁边）
                  row_title_rot = 90,                      # 💡 旋转 90 度（如需不旋转横向显示，可改为 0）
                  row_title_gp = gpar(
                    fontsize = 10, 
                    fontface = "plain", 
                    col = "black"),
                  row_names_gp = gpar(fontsize = 10,fontface = "italic"),
                  show_column_names = FALSE,
                  row_dend_reorder = TRUE,
                  row_gap = unit(1, "mm"),
                  border = TRUE,
                  name = "Z_Score2",heatmap_legend_param = list(
                    direction = "horizontal",
                    title_gp = gpar(fontsize = 10, fontface = "bold"), 
                    labels_gp = gpar(fontsize = 10) ,legend_height = unit(2, "cm")),  
                  left_annotation = ha.r2
)



ht_list2 <- ht_list %v% ht_gene
draw(ht_list2,
     heatmap_legend_side = "left",      # 主图例在右侧
     annotation_legend_side = "left",  # 注释图例全部移到底部
     merge_legend = TRUE,
     ht_gap = unit(1, "mm")
)
dev.off()








#气泡图，人源基因画####
library(ggh4x)
library(dplyr)
library(tibble)
sample_info<-subset(sampleAnno[,c(1:2,6)])
df_long <- as.data.frame(t(gene.expr.human)) %>%
  rownames_to_column("sampleName") %>%
  left_join(sample_info, by = "sampleName") %>%
  pivot_longer(cols = -c(sampleName, GSE,treatment5), names_to = "gene", values_to = "expr")
df_long$treatment5 <- factor(df_long$treatment5, levels = c("MFCB" ,"MFC" ,"MB" , 
                                                         "FCB" ,"YB","Control" ))
unique(df_long$treatment5)  
colnames(df_long)[3]<-"group"
row_anno_df <- row_annotation %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene")
#对df_long的expr进行标准化
df_long1 <- df_long %>%
  group_by(gene) %>%                          # 按基因分组
  mutate(expr_scaled = as.vector(scale(expr))) %>% # 计算Z-score标准化后的值
  ungroup()
df_long1 <- df_long %>%
  group_by(gene)
bubble_data <- df_long1 %>%
  group_by(group, gene) %>%
  summarise(
    mean_expr = mean(expr, na.rm = TRUE), # 用于颜色
    pct_expressed = mean(expr > 0) * 100 ,
    .groups = "drop"
  )

# 2. 根据基因分类对基因进行排序，锁定绘图顺序
gene_order <- c("MYC","FGF19","CCND1","CTNNB1","YAP1")

bubble_data$gene <- factor(bubble_data$gene, levels = gene_order)


# 3. 锁定 5 个组的横坐标顺序
bubble_data$group <- factor(bubble_data$group, 
                            levels = c("MFCB" ,"MFC" ,
                                       "FCB" ,"YB","MB" , "Control"))

p_bubble <- ggplot(bubble_data, aes(x = group, y = gene)) +
  # 绘制气泡
  geom_point(aes(size = pct_expressed, fill = mean_expr), shape = 21, color = "grey20", stroke = 0.4) +
  
  # 颜色渐变
  scale_fill_gradient2(
    low = "#313695", mid = "#F7F7F7", high = "#A50026", 
    midpoint = mean(bubble_data$mean_expr, na.rm = TRUE),
    name = "expression"
  ) +
  
  # 气泡大小范围
  scale_size_continuous(range = c(1, 3), name = "Expressed\nSamples (%)") +
  
  # 横坐标标签
  scale_x_discrete(labels = c("MFCB" ,"MFC" ,
                              "FCB" ,"YB","MB" , "Control"
  )) + 
  
  
  # ── 标题与主题修饰 ──────────────────────────────────────────────────
  labs(
    title = NULL,
    x = NULL, 
    y = NULL
  ) +
  theme_bw(base_size = 8) +
  theme(
    plot.title         = element_text(face = "bold", size = 11, hjust = 0.5, margin = margin(b = 10)),
    axis.text.x        = element_text(size = 8, color = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y        = element_text(size = 8, color = "black", face = "italic"),
    
    # 分面位置与间距
    strip.placement    = "outside",           # 将分面标签置于坐标轴外侧
    strip.background   = element_blank(),    # 必须保留此项或删除原本的 element_blank()，使 strip_themed 生效
    panel.spacing      = unit(0.4, "lines"),  
    panel.border       = element_rect(color = "grey80", fill = NA, linewidth = 0.4), 
    
    # 图例与网格线
    legend.title       = element_text(size = 7.5, face = "bold"),
    legend.text        = element_text(size = 7),
    legend.key.size    = unit(0.35, "cm"),
    panel.grid.major   = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    plot.margin        = margin(10, 10, 10, 15)
  )

# 查看图形
print(p_bubble)

#倒置
my_colors <- c("#77AEBD","#ACCBD4","#F5F7FA","#FDCD9F","#EAAA6C")
range(bubble_data$mean_expr, na.rm = TRUE)

p_bubble <- ggplot(bubble_data, aes(x = gene, y = group)) +   # ← 交换 x/y 实现倒置
  
  # 1. 气泡
  geom_point(
    aes(size = pct_expressed, fill = mean_expr),
    shape = 21, color = "grey20", stroke = 0.4
  ) +
  
  # 2. 自定义 5 色渐变
  scale_fill_gradientn(
    colors = my_colors,
    values = scales::rescale(
      c(min(bubble_data$mean_expr, na.rm = TRUE),
        mean(bubble_data$mean_expr, na.rm = TRUE),
        max(bubble_data$mean_expr, na.rm = TRUE))
    ),
    name = "expression"
  )+
  
  # 3. 气泡大小
  scale_size_continuous(
    range = c(1, 3),
    name  = "Expressed\nSamples (%)"
  ) +
  
  # 4. 横坐标（现在是 group）
  scale_y_discrete(
    labels = c("MFCB", "MFC", "FCB", "YB", "MB", "Control")
  ) +
  
  # 5. 标签
  labs(title = NULL, x = NULL, y = NULL) +
  
  # 6. 主题
  theme_bw(base_size = 8) +
  theme(
    plot.title       = element_text(face = "bold", size = 11, hjust = 0.5,
                                    margin = margin(b = 10)),
    axis.text.x      = element_text(size = 8, color = "black", face = "italic",
                                    angle = 45, hjust = 1, vjust = 1),
    axis.text.y      = element_text(size = 8, color = "black"),
    
    strip.placement  = "outside",
    strip.background = element_blank(),
    panel.spacing    = unit(0.4, "lines"),
    panel.border     = element_rect(color = "grey80", fill = NA, linewidth = 0.4),
    
    legend.title     = element_text(size = 7.5, face = "bold"),
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(0.35, "cm"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(10, 10, 10, 15)
  )

print(p_bubble)



# PDF 保存
ggsave(
  filename = "./results/260922/260922新-PMFCT和它的对照-ICR-21samples-人源基因表达气泡图绝对值.pdf",
  plot     = p_bubble,
  width    = 3,   # 气泡图比较紧凑，不需要之前13.5那么宽了
  height   = 2,
  device   = cairo_pdf
)


#自定义三特征的GSVA特征箱线图####
#一个pdf，GSVA ，clinical2，treatment4
ssgsea_es<-gsva_es[c(11:13),]
plot_df <- as.data.frame(t(ssgsea_es))
plot_df$SampleID <- rownames(plot_df)
plot_df$Group <- sampleAnno$treatment5[match(plot_df$SampleID, rownames(sampleAnno))]

plot_df_long <- reshape2::melt(plot_df, 
                               id.vars = c("SampleID", "Group"), 
                               variable.name = "GeneSet", 
                               value.name = "GSVA_Score")
unique(plot_df_long$Group)
plot_df_long$Group <- factor(plot_df_long$Group, levels = c("MFCB" ,"MFC" ,"MB" , 
                                                            "FCB" ,"YB"))
df_long<-plot_df_long
colnames(df_long)[2]<-"group"
colnames(df_long)[3]<-"gene"
colnames(df_long)[4]<-"expr"
group_colors <- color2
# 分组顺序（影响x轴排列）
group_order <- c("MFCB" ,"MFC" ,"MB" , 
                 "FCB" ,"YB")

# 统计比较对（根据你的生物学问题按需增减）

my_comparisons <- list(
  c("MFCB", "MFC" ),
  c("MFCB", "MB"),
  c("MFCB", "FCB"),c("MFCB", "YB"),
  c("MFC", "MB"),c("FCB", "YB")
)

# ── 4. 统计检验（Wilcoxon，FDR校正）────────────────────────
stat_df <- df_long %>%
  mutate(
    expr = as.numeric(expr),
    group = as.factor(group)
  ) %>%
  # 过滤掉转换过程中可能产生的任何 NA 空值（安全保护）
  filter(!is.na(expr)) %>%
  group_by(gene) %>%
  t_test(expr ~ group, comparisons = my_comparisons) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance(
    p.col = "p.adj",
    symbols = c("****", "***", "**", "*", "ns")
  ) %>%
  # 这里的 group 必须与上面 mutate 的列以及 x 轴变量完全对齐
  add_xy_position(x = "group", dodge = 0.8)
# ── 5. 单基因绘图函数（修正美化版） ────────────────────────────────────────>
plot_one_gene <- function(gene_name, df, stat, colors) {
  
  df_g   <- df   %>% filter(gene == gene_name)
  stat_g <- stat %>% filter(gene == gene_name)
  
  if(nrow(df_g) == 0) return(ggplot() + theme_void()) # 空数据防崩保护
  
  y_max <- max(df_g$expr, na.rm = TRUE)
  y_min <- min(df_g$expr, na.rm = TRUE)
  y_rng <- y_max - y_min
  if (y_rng == 0) y_rng <- 1 # 防止数值完全一致时除以 0
  
  # ── 智能计算每条显著性比较线的高度，防止重叠 ──
  if (nrow(stat_g) > 0) {
    group_levels <- levels(factor(df_g$group))
    if (is.null(group_levels)) group_levels <- unique(df_g$group)
    
    stat_g <- stat_g %>%
      mutate(
        x1_pos = match(group1, group_levels),
        x2_pos = match(group2, group_levels),
        span   = abs(x1_pos - x2_pos)
      ) %>%
      arrange(span, x1_pos) %>% # 先按跨度排序，再按起点排序
      mutate(
        step_height = y_rng * 0.08,
        y.position  = y_max + (y_rng * 0.10) + (row_number() - 1) * step_height
      )
  }
  
  # 动态计算 y 轴上限，防止最高的显著性线和星号被顶部切掉
  y_upper <- if(nrow(stat_g) > 0) max(stat_g$y.position) + (y_rng * 0.15) else y_max + (y_rng * 0.2)
  
  ggplot(df_g, aes(x = group, y = expr, fill = group)) +
    geom_boxplot(
      color          = "black",
      linewidth      = 0.4,
      outlier.shape  = NA,
      outlier.size   = 1.0,
      outlier.stroke = 0.3,
      outlier.color  = "black",
      outlier.fill   = "white",
      width          = 0.45,
      alpha          = 0.85
    ) +
    
    # ── 显著性连接线标记 ──
    stat_pvalue_manual(
      stat_g,
      label           = "p.adj.signif",
      tip.length      = 0.01,         # 竖向短线下垂长度
      hide.ns         = TRUE,         # 隐藏不显著的组 (ns)
      size            = 2.5,
      color           = "grey10",
      linewidth       = 0.35
    ) +
    
    scale_fill_manual(values = colors) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.05)),
      limits = c(NA, y_upper)
    ) +
    labs(
      title = stringr::str_wrap(gene_name, width = 25), # 每行约25字符自动换行
      x = NULL, 
      y = "GSVA Enrichment"
    ) +
    theme_classic(base_size = 6) +
    theme(
      plot.title          = element_text(face = "bold", size = 8, hjust = 0.5, margin = margin(b = 3)),
      
      # 💡 核心修改：隐藏横坐标的分组文字，但保留并显示 X 轴刻度线
      axis.text.x         = element_blank(),                                 # 隐藏 X 轴文本
      axis.ticks.x        = element_line(linewidth = 0.35, color = "black"), # 显示 X 轴刻度线
      axis.ticks.length.x = unit(1, "mm"),                                  # 刻度线长度
      
      # Y轴设置
      axis.text.y         = element_text(size = 6.5, color = "black"),
      axis.title.y        = element_text(size = 6.5, color = "black"),
      axis.line           = element_line(linewidth = 0.35, color = "black"),
      axis.ticks.y        = element_line(linewidth = 0.35, color = "black"),
      
      legend.position     = "none",
      plot.background     = element_blank(),
      panel.background    = element_blank(),
      plot.margin         = margin(4, 4, 4, 4)
    )
}

# ── 2. 批量生成子图列表 ──────────────────────────────────────────────────
target_genes <- row.names(ssgsea_es)
plots <- lapply(target_genes, function(g) {
  plot_one_gene(g, df_long, stat_df, group_colors)
})

# ── 3. 处理 Y 轴标题（仅每行第 1 个图保留说明） ─────────────────────────────
num_cols <- 3

for (i in seq_along(plots)) {
  is_first_in_row <- ((i - 1) %% num_cols == 0)
  
  if (is_first_in_row) {
    plots[[i]] <- plots[[i]] + 
      labs(y = "GSVA Enrichment") +
      theme(
        axis.title.y = element_text(
          size = 8, 
          color = "black", 
          face = "plain", 
          angle = 90, 
          hjust = 0.5, 
          vjust = 0.5, 
          margin = margin(r = 4)
        )
      )
  } else {
    plots[[i]] <- plots[[i]] + 
      labs(y = NULL) + 
      theme(
        axis.title.y = element_blank()
      )
  }
}

# ── 4. 共享图例构建 (右侧垂直布局) ─────────────────────────────────────────
legend_plot <- ggplot(
  data.frame(group = factor(group_order, levels = group_order), x = 1, y = 1),
  aes(x = x, y = y, fill = group)
) +
  geom_point(shape = 22, size = 3.5, color = "black", stroke = 0.3) +
  scale_fill_manual(
    values = group_colors,
    labels = c("MFCB" ,"MFC" ,"MB" , 
               "FCB" ,"YB")
  ) +
  theme_void() +
  theme(
    legend.title     = element_blank(),
    legend.text      = element_text(size = 7.5, color = "grey20"),
    legend.key.size  = unit(0.35, "cm"),
    legend.direction = "vertical",
    legend.position  = "right"
  )

shared_legend <- cowplot::get_legend(legend_plot)

# ── 5. 使用 patchwork 进行左右拼接（网格矩阵 + 右侧图例） ───────────────────
grid_plots <- patchwork::wrap_plots(plots, ncol = num_cols)

# 左右组合：主图矩阵在左，图例在右
combined <- grid_plots | shared_legend

combined <- combined + 
  patchwork::plot_layout(widths = c(10, 2.2)) + # 可微调右侧图例占用的宽度比例
  patchwork::plot_annotation(
    title   = NULL,
    caption = NULL,
    theme   = theme(
      plot.title   = element_text(size = 11, face = "bold", hjust = 0.5, margin = margin(b = 6, t = 4)),
      plot.caption = element_text(size = 6.5, color = "grey40", hjust = 0, margin = margin(t = 4))
    )
  )
# ── 8. 输出保存 ──────────────────────────────────────────────────────────────

ggsave(
  filename = "./results/260901/260901-GSVA-PMFCT和对照15samples-boxplot_target_自定义3geneset.pdf",
  plot     = combined,
  width    = 6,  
  height   = 3,
  device   = cairo_pdf
)


#箱线图
#6基因表达情况
#特征基因，IHC marker####
library(tidyr)
library(dplyr)
library(tidyverse)
exp.mouse<- fread("../HBdataset/data/260305-论文/XHX跑/最新数据备用/292samples/328/260723-最全340samples-public+XHX-mouse-salmon_Symbol_tpm.txt",
                  header = T,data.table = F)
row.names(exp.mouse) <- exp.mouse[, 1]
exp.mouse<-exp.mouse[,-1]


p = identical(colnames(exp.mouse),rownames(clinical3));p
if(!p) {
  s = intersect(rownames(clinical3),colnames(exp.mouse))
  exp.mouse = exp.mouse[,s]
  clinical3 = clinical3[s,]
}

gene.expr.mouse <-log2(exp.mouse+1)%>%as.matrix

target_genes<-c( "Nkd1","Tnfrsf19", 
                 "Dlk1","Sall4",
                 "Slc4a1" ,"Trim58")#wnt-干性，红系
plot_matrix <- gene.expr.mouse[target_genes, ]
wnt_genes      <- c("Nkd1","Tnfrsf19")
stemness_genes <- c("Dlk1","Sall4")
erythroid_genes<-c("Slc4a1" ,"Trim58")

# 3. 构建行注释数据框 (row_annotation)
# 核心逻辑：先创建一个空的 Data Frame，行名设为基因名
row_annotation <- data.frame(
  Gene_Class = rep(NA, length(target_genes)),
  row.names = target_genes
)
row_annotation[wnt_genes, "Gene_Class"]      <- "WNT pathway"
row_annotation[stemness_genes, "Gene_Class"] <- "Fetal liver"
row_annotation[erythroid_genes, "Gene_Class"] <- "Erythropoiesis"
row_annotation$Gene_Class <- factor(
  row_annotation$Gene_Class,
  levels = c("Fetal liver", "WNT pathway", "Erythropoiesis"))
ha.r <- rowAnnotation(
  Gene_Class = row_annotation$Gene_Class, # 行注释的数据源
  col = list(
    Gene_Class = c(
      "Fetal liver"  = "#7876B1",   
      "WNT pathway"  = "#20854E",   
      "Erythropoiesis"  = "#BC3C29"   
    )
  ),
  show_annotation_name = FALSE,show_legend = FALSE,
  simple_anno_size = unit(1.5, "mm")# 是否在行注释上方显示“Pathway”这个标题
)

plot_matrix_scaled <- t(scale(t(plot_matrix)))
quantile(plot_matrix_scaled, c(0.01, 0.05, 0.95, 0.99), na.rm = TRUE)

set.seed(123578)
pdf("./results/260822/260822新-PMFCT和它的对照-ICR-21samples-6特征基因新.pdf", width = 6,height =4.5 )
ht <-Heatmap(plot_matrix_scaled,
             col = colorRamp2(c(-2, 0, 3), c("blue", "white", "red")),
             column_order = order(col_factor), 
             column_split = sampleAnno$treatment5,column_title = NULL,
             cluster_rows = FALSE,     column_gap = unit(1.5, "mm"),
             show_row_dend = FALSE,
             height = unit(nrow(gsva_es1) * 6.5, "mm"),
             width = unit(ncol(gsva_es1) * 3, "mm"),
             row_split = row_annotation$Gene_Class,
             row_title_side = "left",                 # 文字显示在最左侧（色块旁边）
             row_title_rot = 90,                      # 💡 旋转 90 度（如需不旋转横向显示，可改为 0）
             row_title_gp = gpar(
               fontsize = 9, 
               fontface = "bold", 
               col = "black"),
             row_names_gp = gpar(fontsize = 9),
             show_column_names = FALSE,
             row_dend_reorder = TRUE,
             row_gap = unit(1.5, "mm"),
             border = TRUE,
             name = "Z_Score",heatmap_legend_param = list(
               title_gp = gpar(fontsize = 9, fontface = "bold"), # 标题字号（如 9 pt）
               labels_gp = gpar(fontsize = 9)                    # 刻度数字字号（如 8 pt）
             ),
             top_annotation = ha.b,
             left_annotation = ha.r
)

draw(ht, 
     heatmap_legend_side = "right",      # 主图例在右侧
     annotation_legend_side = "right",  # 注释图例全部移到底部
     merge_legend = TRUE
)
dev.off()


#IHC marker,ICR+C57####
#增加例数每组4个样####
exp.mouse<- fread("../HBdataset/data/260305-论文/XHX跑/最新数据备用/292samples/328/260723-最全340samples-public+XHX-mouse-salmon_Symbol_tpm.txt",
                  header = T,data.table = F)
row.names(exp.mouse) <- exp.mouse[, 1]
exp.mouse<-exp.mouse[,-1]
#用clinical3
clinical<- read.csv("../HBdataset/data/260305-论文/XHX跑/最新数据备用/clinical-newGSVA3-最新newsubtype2-4rd-260626.csv",header = T)
#XHX构的10个模型和对照，加上PMFCT14-D28，FCT1008,YB1008，删除MYC-B-YAP-28
clinical3<-clinical[c(39:40,42:43,23:25,6:8,27:28,32,34:35,
                      78,81:83,276:278,102:107),]#ICR+两质粒+YB,XHX,172换176,加上赵斌的三组

clinical3<-clinical[c(39:40,42:43,23:25,6:8,
                      78,81:83,102:107),]#四组ICR
clinical3<-clinical[c(325:327,170:173,176,165:167
                      ,168:169,174),]#C57
clinical3<-clinical[c(325:327,170,172,173,165:167
                      ,168:169,174),]

rownames(clinical3) <- clinical3$sampleName
#用clinical3
p = identical(colnames(exp.mouse),rownames(clinical3));p
if(!p) {
  s = intersect(rownames(clinical3),colnames(exp.mouse))
  exp.mouse = exp.mouse[,s]
  clinical3 = clinical3[s,]
}

gene.expr.mouse <-log2(exp.mouse+1)%>%as.matrix



IHC_genes<-c( "Mki67","Afp", 
                      "Vim",
                      "Cyp2e1","Sall4","Glul","Gypa")
IHC_genes<-c( "Afp",
              "Vim",
              "Cyp2e1","Gata1")
IHC_genes<-c( "Mki67","Sall4",
              "Krt19",
              "Glul","Tfrc")#wnt-干性，红系
plot_matrix <- gene.expr.mouse[IHC_genes, ]

#箱线图及气泡图
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)      # stat_compare_means
library(rstatix)     # wilcox_test / add_xy_position
library(patchwork)   # 拼图
color2 <- c(
  "MFCB"= "#EA4335", 
  "FCB"="#34A853","Control"="#4285F4",
  "MFC"               = "#c4a7e7",
  "MB"            = "#FFB5A7",  
  "YB"="#FBBC05"
) 
color2 <- c(
  "MFCB"= "#EA4335", 
  "FCB"="#34A853","Control"="#4285F4",
  "YB"="#FBBC05"
) 
group_colors <- color2

# 分组顺序（影响x轴排列）
group_order <- c("MFCB" ,"MFC" ,"MB" ,"FCB" , "YB",
                 "Control"  )
group_order <- c("MFCB" , "YB","FCB" ,
                 "Control"  )

# 统计比较对（根据你的生物学问题按需增减）
my_comparisons <- list(
  c("MFCB", "MFC" ),
  c("MFCB", "MB"),
  c("MFCB", "FCB"),c("MFCB", "YB"),
c("FCB", "YB")
)
my_comparisons <- list(
  c("MFCB", "YB" ),
  c("MFCB", "FCB"),
  c("MFCB", "Control"), c("YB", "FCB"),c("Control", "YB"),
  c("Control", "FCB")
)
#加上组别
# 1. 转置并合并分组信息
#最全注释####
sampleAnno <- setDF(clinical3)
sampleAnno$Type[sampleAnno$Pathology == ""] <- "Control"
sampleAnno$MYC[sampleAnno$MYC == ""] <- "No"
sampleAnno$CTNNB1.mut[sampleAnno$CTNNB1.mut == ""] <- "No"
sampleAnno$FGF19.CCND1[sampleAnno$FGF19.CCND1 == ""] <- "No"
sampleAnno$YAP1.mut[sampleAnno$YAP1.mut == ""] <- "No"
sampleAnno$NFE2L2.mut[sampleAnno$NFE2L2.mut == ""] <- "No"

sample_info<-subset(clinical3[,c(1:2,6)])
df_long <- as.data.frame(t(plot_matrix)) %>%
  rownames_to_column("sampleName") %>%
  left_join(sample_info, by = "sampleName") %>%
  pivot_longer(cols = -c(sampleName, GSE,treatment5), names_to = "gene", values_to = "expr")
df_long$treatment5 <- factor(df_long$treatment5, levels = c("MFCB" ,"YB","FCB" , 
                                                            "Control"))#"MFC" ,"MB" 

colnames(df_long)[3]<-"group"

# ── 4. 统计检验（Wilcoxon，FDR校正）────────────────────────
stat_df <- df_long %>%
  mutate(
    expr = as.numeric(expr),
    group = as.factor(group)
  ) %>%
  # 过滤掉转换过程中可能产生的任何 NA 空值（安全保护）
  filter(!is.na(expr)) %>%
  group_by(gene) %>%
  t_test(expr ~ group, comparisons = my_comparisons) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance(
    p.col = "p.adj",
    symbols = c("****", "***", "**", "*", "ns")
  ) %>%
  # 这里的 group 必须与上面 mutate 的列以及 x 轴变量完全对齐
  add_xy_position(x = "group", dodge = 0.8)
# ── 5. 单基因绘图函数（修正美化版） ────────────────────────────────────────>
plot_one_gene <- function(gene_name, df, stat, colors) {
  
  df_g   <- df   %>% filter(gene == gene_name)
  stat_g <- stat %>% filter(gene == gene_name)
  
  if(nrow(df_g) == 0) return(ggplot() + theme_void()) # 空数据防崩保护
  
  y_max <- max(df_g$expr, na.rm = TRUE)
  y_min <- min(df_g$expr, na.rm = TRUE)
  y_rng <- y_max - y_min
  if (y_rng == 0) y_rng <- 1 # 防止数值完全一致时除以 0
  
  # ── 智能计算每条显著性比较线的高度，防止重叠 ──
  if (nrow(stat_g) > 0) {
    group_levels <- levels(factor(df_g$group))
    if (is.null(group_levels)) group_levels <- unique(df_g$group)
    
    stat_g <- stat_g %>%
      mutate(
        x1_pos = match(group1, group_levels),
        x2_pos = match(group2, group_levels),
        span   = abs(x1_pos - x2_pos)
      ) %>%
      arrange(span, x1_pos) %>% # 先按跨度排序，再按起点排序
      mutate(
        step_height = y_rng * 0.08,
        y.position  = y_max + (y_rng * 0.10) + (row_number() - 1) * step_height
      )
  }
  
  # 动态计算 y 轴上限，防止最高的显著性线和星号被顶部切掉
  y_upper <- if(nrow(stat_g) > 0) max(stat_g$y.position) + (y_rng * 0.15) else y_max + (y_rng * 0.2)
  
  ggplot(df_g, aes(x = group, y = expr, fill = group)) +
    geom_boxplot(
      color          = "black",
      linewidth      = 0.4,
      outlier.shape  = NA,
      outlier.size   = 1.0,
      outlier.stroke = 0.3,
      outlier.color  = "black",
      outlier.fill   = "white",
      width          = 0.45,
      alpha          = 0.85
    ) +
    
    # ── 显著性连接线标记 ──
    stat_pvalue_manual(
      stat_g,
      label           = "p.adj.signif",
      tip.length      = 0.01,         # 竖向短线下垂长度
      hide.ns         = TRUE,         # 隐藏不显著的组 (ns)
      size            = 3,
      color           = "grey10",
      linewidth       = 0.35
    ) +
    
    scale_fill_manual(values = colors) +
    scale_y_continuous(
      limits = c(0, y_upper),
      expand = expansion(mult = c(0, 0.05))
    )  +
    labs(
      title = stringr::str_wrap(gene_name, width = 25), # 每行约25字符自动换行
      x = NULL, 
      y = expression(log[2](TPM + 1))) +
    theme_classic(base_size = 6) +
    theme(
      plot.title          = element_text(face = "italic", size = 10, hjust = 0.5, margin = margin(b = 3)),
      
      # 💡 核心修改：隐藏横坐标的分组文字，但保留并显示 X 轴刻度线
      axis.text.x         = element_text(size = 7.5, color = "black", angle = 45, hjust = 1, vjust = 1),                              # 隐藏 X 轴文本
      axis.ticks.x        = element_line(linewidth = 0.35, color = "black"), # 显示 X 轴刻度线
      axis.ticks.length.x = unit(1, "mm"),                                  # 刻度线长度
      
      # Y轴设置
      axis.text.y         = element_text(size = 10, color = "black"),
      axis.title.y        = element_text(size = 10, color = "black"),
      axis.line           = element_line(linewidth = 0.35, color = "black"),
      axis.ticks.y        = element_line(linewidth = 0.35, color = "black"),
      
      legend.position     = "none",
      plot.background     = element_blank(),
      panel.background    = element_blank(),
      plot.margin         = margin(4, 4, 4, 4)
    )
}

# ── 2. 批量生成子图列表 ──────────────────────────────────────────────────
target_genes <- row.names(plot_matrix)
plots <- lapply(target_genes, function(g) {
  plot_one_gene(g, df_long, stat_df, group_colors)
})

# ── 3. 处理 Y 轴标题（仅每行第 1 个图保留说明） ─────────────────────────────
num_cols <- 7

for (i in seq_along(plots)) {
  is_first_in_row <- ((i - 1) %% num_cols == 0)
  
  if (is_first_in_row) {
    plots[[i]] <- plots[[i]] + 
      labs(y = expression(log[2](TPM + 1))) +
      theme(
        axis.title.y = element_text(
          size = 10, 
          color = "black", 
          face = "plain", 
          angle = 90, 
          hjust = 0.5, 
          vjust = 0.5, 
          margin = margin(r = 4)
        )
      )
  } else {
    plots[[i]] <- plots[[i]] + 
      labs(y = NULL) + 
      theme(
        axis.title.y = element_blank()
      )
  }
}
combined <- patchwork::wrap_plots(plots, ncol = num_cols) +
  patchwork::plot_annotation(
    title   = NULL,
    caption = NULL,
    theme   = theme(
      plot.title   = element_text(size = 10, face = "bold", hjust = 0.5, margin = margin(b = 6, t = 4)),
      plot.caption = element_text(size = 6.5, color = "grey40", hjust = 0, margin = margin(t = 4))
    )
  )

print(combined)



# ── 8. 输出保存 ──────────────────────────────────────────────────────────────
ggsave(
  filename = "./results/261012/261108新-PMFCT和它的对照-21samples-7IHCmarker-ICR高.pdf",
  plot     = combined,
  width    = 7.2,  
  height   = 4.2,
  device   = cairo_pdf
)#6组用8.2
#不显示右边注释，每组显示分组，codex####
pair_key <- function(a, b) {
  paste(pmin(as.character(a), as.character(b)),
        pmax(as.character(a), as.character(b)),
        sep = "|")
}

# 指定每组比较所在的层
comparison_layers <- c(
  setNames(1, pair_key("MFCB", "YB")),
  setNames(1, pair_key("Control", "FCB")),
  setNames(2, pair_key("YB", "FCB")),
  setNames(3, pair_key("MFCB", "FCB")),
  setNames(4, pair_key("Control", "YB")),
  setNames(5, pair_key("MFCB", "Control"))
)

plot_one_gene <- function(gene_name, df, stat, colors) {
  
  df_g <- df %>%
    filter(gene == gene_name) %>%
    mutate(group = factor(as.character(group), levels = group_order))
  
  stat_g <- stat %>%
    filter(gene == gene_name) %>%
    mutate(
      comparison = pair_key(group1, group2),
      layer = unname(comparison_layers[comparison])
    ) %>%
    filter(!is.na(layer), !is.na(p.adj.signif),
           p.adj.signif != "ns")
  
  if (nrow(df_g) == 0 || all(is.na(df_g$expr))) {
    return(ggplot() + theme_void())
  }
  
  y_max <- max(df_g$expr, na.rm = TRUE)
  y_min <- min(df_g$expr, na.rm = TRUE)
  
  # 控制括号之间的距离；可按需要调整 0.14
  y_span <- max(y_max, y_max - y_min, 0.5)
  bracket_base <- y_max + 0.10 * y_span
  bracket_step <- 0.12 * y_span
  
  if (nrow(stat_g) > 0) {
    stat_g <- stat_g %>%
      mutate(
        y.position = bracket_base + (layer - 1) * bracket_step
      )
    
    y_upper <- max(stat_g$y.position) + 0.20 * y_span
  } else {
    y_upper <- y_max + 0.20 * y_span
  }
  
  p <- ggplot(df_g, aes(x = group, y = expr, fill = group)) +
    geom_boxplot(
      color = "black",
      linewidth = 0.4,
      outlier.shape = NA,
      width = 0.45,
      alpha = 0.85
    ) +
    scale_x_discrete(limits = group_order, drop = FALSE) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(
      limits = c(0, y_upper),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      title = stringr::str_wrap(gene_name, width = 25),
      x = NULL,
      y = expression(log[2](TPM + 1))
    ) +
    theme_classic(base_size = 6) +
    theme(
      plot.title = element_text(
        face = "italic", size = 10, hjust = 0.5,
        margin = margin(b = 3)
      ),
      axis.text.x = element_text(
        size = 7.5, color = "black",
        angle = 45, hjust = 1, vjust = 1
      ),
      axis.ticks.x = element_line(
        linewidth = 0.35, color = "black"
      ),
      axis.ticks.length.x = grid::unit(1, "mm"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title.y = element_text(size = 10, color = "black"),
      axis.line = element_line(
        linewidth = 0.35, color = "black"
      ),
      axis.ticks.y = element_line(
        linewidth = 0.35, color = "black"
      ),
      legend.position = "none",
      plot.background = element_blank(),
      panel.background = element_blank(),
      plot.margin = margin(4, 4, 4, 4)
    )
  
  if (nrow(stat_g) > 0) {
    p <- p + ggpubr::stat_pvalue_manual(
      stat_g,
      label = "p.adj.signif",
      tip.length = 0.01,
      bracket.shorten = 0.06,
      step.increase = 0,
      size = 3.5,
      color = "grey10",
      linewidth = 0.35
    )
  }
  
  p
}
target_genes <- row.names(plot_matrix)
plots <- lapply(target_genes, function(g) {
  plot_one_gene(g, df_long, stat_df, group_colors)
})

# ── 3. 处理 Y 轴标题（仅每行第 1 个图保留说明） ─────────────────────────────
num_cols <- 7

for (i in seq_along(plots)) {
  is_first_in_row <- ((i - 1) %% num_cols == 0)
  
  if (is_first_in_row) {
    plots[[i]] <- plots[[i]] + 
      labs(y = expression(log[2](TPM + 1))) +
      theme(
        axis.title.y = element_text(
          size = 10, 
          color = "black", 
          face = "plain", 
          angle = 90, 
          hjust = 0.5, 
          vjust = 0.5, 
          margin = margin(r = 4)
        )
      )
  } else {
    plots[[i]] <- plots[[i]] + 
      labs(y = NULL) + 
      theme(
        axis.title.y = element_blank()
      )
  }
}
combined <- patchwork::wrap_plots(plots, ncol = num_cols) +
  patchwork::plot_annotation(
    title   = NULL,
    caption = NULL,
    theme   = theme(
      plot.title   = element_text(size = 10, face = "bold", hjust = 0.5, margin = margin(b = 6, t = 4)),
      plot.caption = element_text(size = 6.5, color = "grey40", hjust = 0, margin = margin(t = 4))
    )
  )

print(combined)



# ── 8. 输出保存 ──────────────────────────────────────────────────────────────
ggsave(
  filename = "./results/261012/261108新-PMFCT和它的对照-7IHCmarker-C57-12sample.pdf",
  plot     = combined,
  width    = 7.5,  
  height   = 2,
  device   = cairo_pdf
)#6组用8.2




#bubble####
library(ggh4x)
library(dplyr)
library(tibble)
sample_info<-subset(clinical3[,c(1:3)])
df_long <- as.data.frame(t(plot_matrix)) %>%
  rownames_to_column("sampleName") %>%
  left_join(sample_info, by = "sampleName") %>%
  pivot_longer(cols = -c(sampleName, GSE,treatment5), names_to = "gene", values_to = "expr")
df_long$treatment5 <- factor(df_long$treatment5, levels = c("MFCB" ,"MFC" ,"MB" , 
                                                            "FCB" ,"YB",
                                                            "Control" ))

colnames(df_long)[3]<-"group"
row_anno_df <- row_annotation %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene")
#对df_long的expr进行标准化
df_long1 <- df_long %>%
  group_by(gene) %>%                          # 按基因分组
  mutate(expr_scaled = as.vector(scale(expr))) %>% # 计算Z-score标准化后的值
  ungroup()
bubble_data <- df_long1 %>%
  group_by(group, gene) %>%
  summarise(
    mean_expr = mean(expr_scaled, na.rm = TRUE), # 用于颜色
    pct_expressed = mean(expr > 0) * 100 ,
    .groups = "drop"
  )

bubble_data <- bubble_data %>%
  left_join(row_anno_df, by = "gene") %>% 
  select(everything(), gene_class = Gene_Class)
# 2. 根据基因分类对基因进行排序，锁定绘图顺序
gene_order <- bubble_data %>%
  distinct(gene, gene_class) %>%
  arrange(gene_class, gene) %>% # 按分类排序，同分类内按字母排序
  pull(gene)

bubble_data$gene <- factor(bubble_data$gene, levels = gene_order)


# 3. 锁定 5 个组的横坐标顺序
bubble_data$group <- factor(bubble_data$group, 
                            levels = c("MFCB" ,"MFC" ,"MB" , 
                                       "FCB" ,"YB",
                                       "Control"))
gene_order <- bubble_data %>%
  distinct(gene, gene_class) %>%
  arrange(gene_class, gene) %>% # 按分类排序，同分类内按字母排序
  pull(gene)

bubble_data$gene <- factor(bubble_data$gene, levels = gene_order)
Gene_Class_Colors <- c(
  "Erythropoiesis"  = "#BC3C29"  ,"Fetal liver"  = "#7876B1",   
  "WNT pathway"  = "#20854E"
)
bubble_data$gene_class <- factor(
  bubble_data$gene_class, 
  levels = c("Erythropoiesis", "Fetal liver", "WNT pathway")
)

p_bubble <- ggplot(bubble_data, aes(x = group, y = gene)) +
  # 绘制气泡
  geom_point(aes(size = pct_expressed, fill = mean_expr), shape = 21, color = "grey20", stroke = 0.4) +
  
  # 颜色渐变
  scale_fill_gradient2(
    low = "#313695", mid = "#F7F7F7", high = "#A50026", 
    midpoint = mean(bubble_data$mean_expr, na.rm = TRUE),
    name = "Z_scores"
  ) +
  
  # 气泡大小范围
  scale_size_continuous(range = c(1, 6), name = "Expressed\nSamples (%)") +
  
  # 横坐标标签
  scale_x_discrete(labels = c("MFCB" ,"MFC" ,"MB" , 
                              "FCB" ,"YB",
                              "Control"
  )) + 
  
  # ── 核心：高级分面设置（含柱状背景框与文字旋转） ──────────────────────
  facet_grid2(
    gene_class ~ ., 
    scales = "free_y", 
    space = "free_y", 
    switch = "y", # 将分面标签放左侧
    strip = strip_themed(
      # background_y: 设置柱状背景框（"33"代表20%透明度的浅色填充）
      background_y = list(
        element_rect(fill = paste0(Gene_Class_Colors["Erythropoiesis"], "33"), color = "NA", linewidth = 0.4),
        element_rect(fill = paste0(Gene_Class_Colors["Fetal liver"], "33"), color = "NA", linewidth = 0.4),
        element_rect(fill = paste0(Gene_Class_Colors["WNT pathway"], "33"), color = "NA", linewidth = 0.4)
      ),
      # text_y: 设置文字样式（颜色对应、旋转-90度、居中对齐）
      text_y = list(
        element_text(color = Gene_Class_Colors["Erythropoiesis"], face = "bold", size = 8.5, angle = 90, hjust = 0.5, vjust = 0.5),
        element_text(color = Gene_Class_Colors["Fetal liver"], face = "bold", size = 8.5, angle = 90, hjust = 0.5, vjust = 0.5),
        element_text(color = Gene_Class_Colors["WNT pathway"], face = "bold", size = 8.5, angle = 90, hjust = 0.5, vjust = 0.5)
      )
    )
  )+
  
  # ── 标题与主题修饰 ──────────────────────────────────────────────────
  labs(
    title = "Expression Pattern of Target Genes Across Mouse Models",
    x = NULL, 
    y = NULL
  ) +
  theme_bw(base_size = 8) +
  theme(
    plot.title         = element_text(face = "bold", size = 11, hjust = 0.5, margin = margin(b = 10)),
    axis.text.x        = element_text(size = 8, color = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y        = element_text(size = 8, color = "black", face = "italic"),
    
    # 分面位置与间距
    strip.placement    = "outside",           # 将分面标签置于坐标轴外侧
    strip.background   = element_blank(),    # 必须保留此项或删除原本的 element_blank()，使 strip_themed 生效
    panel.spacing      = unit(0.4, "lines"),  
    panel.border       = element_rect(color = "grey80", fill = NA, linewidth = 0.4), 
    
    # 图例与网格线
    legend.title       = element_text(size = 7.5, face = "bold"),
    legend.text        = element_text(size = 7),
    legend.key.size    = unit(0.35, "cm"),
    panel.grid.major   = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    plot.margin        = margin(10, 10, 10, 15)
  )

# 查看图形
print(p_bubble)

# PDF 保存
ggsave(
  filename = "./results/260822/260822新-PMFCT和它的对照_21sample-bubble_target_6genes.pdf",
  plot     = p_bubble,
  width    = 6,   # 气泡图比较紧凑，不需要之前13.5那么宽了
  height   = 4.5,
  device   = cairo_pdf
)
#
##横纵坐标调换
p_bubble <- ggplot(bubble_data, aes(x = group, y = gene)) +
  # 1. 绘制气泡
  geom_point(aes(size = pct_expressed, fill = mean_expr), shape = 21, color = "grey20", stroke = 0.4) +
  
  # 2. 颜色渐变
  scale_fill_gradient2(
    low = "#313695", mid = "#F7F7F7", high = "#A50026", 
    midpoint = mean(bubble_data$mean_expr, na.rm = TRUE),
    name = "Z_scores"
  ) +
  
  # 3. 气泡大小范围
  scale_size_continuous(range = c(1, 3), name = "Expressed\nSamples (%)") +
  
  # 4. 原 group 轴标签
  scale_x_discrete(labels = c(
    "MFCB" ,"MFC" ,"MB" , 
    "FCB" ,"YB",
    "Control"
  )) + 
  
  # 5. 按 gene_class 分面
  facet_grid(. ~ gene_class, scales = "free_x", space = "free_x") +
  
  # 6. 交换 X 轴与 Y 轴
  coord_flip() +
  
  # 7. 标题与标签
  labs(
    title = NULL,
    x = NULL, 
    y = NULL
  ) +
  
  # 💡 8. 更改为经典半框主题，并进行微调
  theme_classic(base_size = 8) +
  theme(
    plot.title       = element_text(face = "bold", size = 11, hjust = 0.5, margin = margin(b = 10)),
    
    # 转置后文本对齐
    axis.text.x      = element_text(size = 8, color = "black", face = "italic", angle = 45, hjust = 1, vjust = 1),
    axis.text.y      = element_text(size = 8, color = "black", face = "plain"),
    
    # 💡 保留半框的坐标轴线
    axis.line        = element_line(color = "black", linewidth = 0.4),
    
    # 💡 移除默认的全外框，只保留分面背景设置
    panel.border     = element_blank(),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 8.5),
    
    # 图例与网格线
    legend.title     = element_text(size = 7.5, face = "bold"),
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(0.35, "cm"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(10, 10, 10, 15)
  )

print(p_bubble)
ggsave(
  filename = "./results/260822/260822新-PMFCT和它的对照_21sample-bubble_target_6genes-倒置.pdf",
  plot     = p_bubble,
  width    = 4.5,   # 气泡图比较紧凑，不需要之前13.5那么宽了
  height   = 3,
  device   = cairo_pdf
)  


#人源基因和内源基因表达情况,删除YAP1,NFE2L2####
exp.mouse<- fread("../HBdataset/data/260305-论文/XHX跑/最新数据备用/292samples/260723-最全295samples-public+XHX-mouse-salmon_Symbol_tpm.txt",
                  header = T,data.table = F)
row.names(exp.mouse) <- exp.mouse[, 1]
exp.mouse<-exp.mouse[,-1]
#用clinical3
p = identical(colnames(exp.mouse),rownames(clinical3));p
if(!p) {
  s = intersect(colnames(exp.mouse),rownames(clinical3))
  exp.mouse = exp.mouse[,s]
  clinical3 = clinical3[s,]
}

gene.expr.mouse <-log2(exp.mouse+1)%>%as.matrix
target_genes<-c( "Myc","Fgf15","Ccnd1","Ctnnb1","Yap1")#wnt-干性，红系
plot_matrix <- gene.expr.mouse[target_genes, ]

##读入人源基因的tpm再转化为log2
exp1<-fread('../HBdataset/data/260305-论文/XHX跑/最新数据备用/292samples/260723-最全295samples-public+XHX-HUMANgene-salmon_Symbol_tpm.txt',
            header = T,data.table = F)
row.names(exp1)<-exp1$Symbol

exp1<-exp1[,-1]
p = identical(colnames(exp1),rownames(clinical3));p
if(!p) {
  s = intersect(colnames(exp1),rownames(clinical3))
  exp1 = exp1[,s]
  clinical3 = clinical3[s,]
}
gene.expr.human <-log2(exp1+1)%>%as.matrix
gene.expr.human<-gene.expr.human[c("MYC","FGF19","CTNNB1","CCND1","YAP1"),]


gene.expr.hm<-rbind(plot_matrix,gene.expr.human)
#画图
Mouse_genes      <- c( "Myc","Fgf15","Ccnd1","Ctnnb1","Yap1")
Human_genes <- c("MYC","FGF19","CTNNB1","CCND1","YAP1")

# 2. 合并所有你想要展示在热图中的基因（作为最终的行名）
all_target_genes <- c(Mouse_genes,Human_genes)
# 如果有重复的基因，用 unique() 自动去重
all_target_genes <- unique(all_target_genes)

# 3. 干净利落地构建全新的 row_annotation
row_annotation <- data.frame(row.names = all_target_genes)

# 4. 批量赋值分类标签
row_annotation$Gene_Class <- case_when(
  rownames(row_annotation) %in% Mouse_genes     ~ "Mouse genes",
  rownames(row_annotation) %in% Human_genes       ~ "Human genes",
  TRUE                                          ~ "Others"
)

# 5. 转为因子并指定你在热图里从上到下的分块顺序
desired_row_order <- c("Mouse genes", "Human genes")
row_annotation$Gene_Class <- factor(row_annotation$Gene_Class, levels = desired_row_order)

# 行注释保持原样
ha.r <- rowAnnotation(
  Gene_Class = row_annotation$Gene_Class, # 行注释的数据源
  col = list(
    Gene_Class = c(
      "Mouse genes"  = "#CD5C5C",   
      "Human genes"  = "#20854E"
    )
  ),
  show_annotation_name = FALSE,show_legend = FALSE,
  simple_anno_size = unit(1.5, "mm")# 是否在行注释上方显示“Pathway”这个标题
)
plot_matrix_scaled <- t(scale(t(gene.expr.hm)))
quantile(plot_matrix_scaled, c(0.01, 0.05, 0.95, 0.99), na.rm = TRUE)
plot_matrix_scaled[is.na(plot_matrix_scaled)] <- 0

set.seed(123578)
pdf("./results/260822/260822新-PMFCT和它的对照-ICR-21samples-内外源特征基因新.pdf", width = 6,height =4.5 )
ht <-Heatmap(plot_matrix_scaled,
             col = colorRamp2(c(-2, 0, 3), c("blue", "white", "red")),
             column_order = order(col_factor), 
             column_split = sampleAnno$treatment5,column_title = NULL,
             cluster_rows = FALSE,     column_gap = unit(1.5, "mm"),
             show_row_dend = FALSE,
             height = unit(nrow(plot_matrix_scaled) * 6.5, "mm"),
             width = unit(ncol(plot_matrix_scaled) * 3, "mm"),
             row_split = row_annotation$Gene_Class,
             row_title_side = "left",                 # 文字显示在最左侧（色块旁边）
             row_title_rot = 90,                      # 💡 旋转 90 度（如需不旋转横向显示，可改为 0）
             row_title_gp = gpar(
               fontsize = 9, 
               fontface = "bold", 
               col = "black"),
             row_names_gp = gpar(fontsize = 9),
             show_column_names = FALSE,
             row_dend_reorder = TRUE,
             row_gap = unit(1.5, "mm"),
             border = TRUE,
             name = "Z_Score",heatmap_legend_param = list(
               title_gp = gpar(fontsize = 9, fontface = "bold"), # 标题字号（如 9 pt）
               labels_gp = gpar(fontsize = 9)                    # 刻度数字字号（如 8 pt）
             ),
             top_annotation = ha.b,
             left_annotation = ha.r
)

draw(ht, 
     heatmap_legend_side = "right",      # 主图例在右侧
     annotation_legend_side = "right",  # 注释图例全部移到底部
     merge_legend = TRUE
)
dev.off()

#气泡图,箱线图不合适

group_colors <- color2
# 分组顺序（影响x轴排列）
group_order <- c("MFCB" ,"MFC" ,"MB" , 
                 "FCB" ,"YB",
                 "Control")



#加上组别
# 1. 转置并合并分组信息
sample_info<-subset(clinical3[,c(1:3)])
df_long <- as.data.frame(t(gene.expr.hm)) %>%
  rownames_to_column("sampleName") %>%
  left_join(sample_info, by = "sampleName") %>%
  pivot_longer(cols = -c(sampleName, GSE,treatment5), names_to = "gene", values_to = "expr")
df_long$treatment5 <- factor(df_long$treatment5, levels = c("MFCB" ,"MFC" ,"MB" , 
                                                            "FCB" ,"YB",
                                                            "Control"))

colnames(df_long)[3]<-"group"

#bubble####
library(ggh4x)
library(dplyr)
library(tibble)
row_anno_df <- row_annotation %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene")
#对df_long的expr进行标准化
df_long1 <- df_long %>%
  group_by(gene) %>%                          # 按基因分组
  mutate(expr_scaled = as.vector(scale(expr))) %>% # 计算Z-score标准化后的值
  ungroup()
bubble_data <- df_long1 %>%
  group_by(group, gene) %>%
  summarise(
    mean_expr = mean(expr_scaled, na.rm = TRUE), # 用于颜色
    pct_expressed = mean(expr > 0) * 100 ,
    .groups = "drop"
  )

bubble_data <- bubble_data %>%
  left_join(row_anno_df, by = "gene") %>% 
  select(everything(), gene_class = Gene_Class)
# 2. 根据基因分类对基因进行排序，锁定绘图顺序
gene_order <- bubble_data %>%
  distinct(gene, gene_class) %>%
  arrange(gene_class, gene) %>% # 按分类排序，同分类内按字母排序
  pull(gene)

bubble_data$gene <- factor(bubble_data$gene, levels = gene_order)


# 3. 锁定 5 个组的横坐标顺序
bubble_data$group <- factor(bubble_data$group, 
                            levels = c("MFCB" ,"MFC" ,"MB" , 
                                       "FCB" ,"YB",
                                       "Control" ))

Gene_Class_Colors <- c(
  "Mouse genes"  = "#CD5C5C",   
  "Human genes"  = "#20854E" 
)

bubble_data$gene_class <- factor(
  bubble_data$gene_class, 
  levels = c( "Mouse genes", "Human genes")
)
p_bubble <- ggplot(bubble_data, aes(x = group, y = gene)) +
  # 绘制气泡
  geom_point(aes(size = pct_expressed, fill = mean_expr), shape = 21, color = "grey20", stroke = 0.4) +
  
  # 颜色渐变
  scale_fill_gradient2(
    low = "#313695", mid = "#F7F7F7", high = "#A50026", 
    midpoint = mean(bubble_data$mean_expr, na.rm = TRUE),
    name = "Z_scores"
  ) +
  
  # 气泡大小范围
  scale_size_continuous(range = c(1, 6), name = "Expressed\nSamples (%)") +
  
  # 横坐标标签
  scale_x_discrete(labels = c("MFCB" ,"MFC" ,"MB" , 
                              "FCB" ,"YB",
                              "Control"
  )) + 
  
  # ── 核心：高级分面设置（含柱状背景框与文字旋转） ──────────────────────
  facet_grid2(
    gene_class ~ ., 
    scales = "free_y", 
    space = "free_y", 
    switch = "y", # 将分面标签放左侧
    strip = strip_themed(
      # background_y: 设置柱状背景框（"33"代表20%透明度的浅色填充）
      background_y = list(
        element_rect(fill = paste0(Gene_Class_Colors["Mouse genes"], "33"), color = "NA", linewidth = 0.4),
        element_rect(fill = paste0(Gene_Class_Colors["Human genes"], "33"), color = "NA", linewidth = 0.4)
      ),
      # text_y: 设置文字样式（颜色对应、旋转-90度、居中对齐）
      text_y = list(
        element_text(color = Gene_Class_Colors["Mouse genes"], face = "bold", size = 8.5, angle = 90, hjust = 0.5, vjust = 0.5),
        element_text(color = Gene_Class_Colors["Human genes"], face = "bold", size = 8.5, angle = 90, hjust = 0.5, vjust = 0.5)
      )
    )
  ) +
  
  # ── 标题与主题修饰 ──────────────────────────────────────────────────
  labs(
    title = "Expression Pattern of Target Genes Across Mouse Models",
    x = NULL, 
    y = NULL
  ) +
  theme_bw(base_size = 8) +
  theme(
    plot.title         = element_text(face = "bold", size = 11, hjust = 0.5, margin = margin(b = 10)),
    axis.text.x        = element_text(size = 8, color = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y        = element_text(size = 8, color = "black", face = "italic"),
    
    # 分面位置与间距
    strip.placement    = "outside",           # 将分面标签置于坐标轴外侧
    strip.background   = element_blank(),    # 必须保留此项或删除原本的 element_blank()，使 strip_themed 生效
    panel.spacing      = unit(0.4, "lines"),  
    panel.border       = element_rect(color = "grey80", fill = NA, linewidth = 0.4), 
    
    # 图例与网格线
    legend.title       = element_text(size = 7.5, face = "bold"),
    legend.text        = element_text(size = 7),
    legend.key.size    = unit(0.35, "cm"),
    panel.grid.major   = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    plot.margin        = margin(10, 10, 10, 15)
  )

# 查看图形
print(p_bubble)

# PDF 保存
ggsave(
  filename = "./results/260822/260822新-PMFCT和它的对照-21samples_bubble_外源内源基因genes-少NRF2.pdf",
  plot     = p_bubble,
  width    = 4.5,   # 气泡图比较紧凑，不需要之前13.5那么宽了
  height   = 4,
  device   = cairo_pdf
)
# 💡 再次运行绘图代码,调换：####
p_bubble <- ggplot(bubble_data, aes(x = group, y = gene)) +
  # 1. 绘制气泡
  geom_point(aes(size = pct_expressed, fill = mean_expr), shape = 21, color = "grey20", stroke = 0.4) +
  
  # 2. 颜色渐变
  scale_fill_gradient2(
    low = "#313695", mid = "#F7F7F7", high = "#A50026", 
    midpoint = mean(bubble_data$mean_expr, na.rm = TRUE),
    name = "Z_scores"
  ) +
  
  # 3. 气泡大小范围
  scale_size_continuous(range = c(1, 3), name = "Expressed\nSamples (%)") +
  
  # 4. 原 group 轴标签
  scale_x_discrete(labels = c(
    "MFCB" ,"MFC" ,"MB" , 
    "FCB" ,"YB",
    "Control"
  )) + 
  
  # 5. 按 gene_class 分面
  facet_grid(. ~ gene_class, scales = "free_x", space = "free_x") +
  
  # 6. 交换 X 轴与 Y 轴
  coord_flip() +
  
  # 7. 标题与标签
  labs(
    title = NULL,
    x = NULL, 
    y = NULL
  ) +
  
  # 💡 8. 更改为经典半框主题，并进行微调
  theme_classic(base_size = 8) +
  theme(
    plot.title       = element_text(face = "bold", size = 11, hjust = 0.5, margin = margin(b = 10)),
    
    # 转置后文本对齐
    axis.text.x      = element_text(size = 8, color = "black", face = "italic", angle = 45, hjust = 1, vjust = 1),
    axis.text.y      = element_text(size = 8, color = "black", face = "plain"),
    
    # 💡 保留半框的坐标轴线
    axis.line        = element_line(color = "black", linewidth = 0.4),
    
    # 💡 移除默认的全外框，只保留分面背景设置
    panel.border     = element_blank(),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 8.5),
    
    # 图例与网格线
    legend.title     = element_text(size = 7.5, face = "bold"),
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(0.35, "cm"),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(10, 10, 10, 15)
  )

print(p_bubble)
ggsave(
  filename = "./results/260822/260822新-PMFCT和它的对照-21samples_bubble_外源内源基因genes-少NRF2-倒置.pdf",
  plot     = p_bubble,
  width    = 6,   # 气泡图比较紧凑，不需要之前13.5那么宽了
  height   = 3,
  device   = cairo_pdf
)  

