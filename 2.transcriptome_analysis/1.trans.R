rm(list=ls())
##
library("RColorBrewer")
library("tidyverse")
library("DESeq2")
library("ggcorrplot")
library("amap")
library("pheatmap")
library("gplots")
library(ggnewscale)

dataset <- read.table("bsta.tpm.txt",header = TRUE, row.names = 1)
pearson_cor <- as.matrix(cor(dataset, method="spearman"))
hc <- hcluster(t(dataset), method="spearman")
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
t <- function(){
  hclust(method = "average")
}
heatmap.2(pearson_cor, 
          Rowv=as.dendrogram(hc), 
          symm=T, 
          trace="none",
          dendrogram = c("row"),
          distfun=dist(method = "canberra"),
          col=hmcol, 
          #colsep=c(seq(0,15,1)),
          #rowsep=c(seq(0,15,1)),
          #sepcolor="gray90",
          #sepwidth=c(0.001,0.001),
          #margins=c(11,11),
          main="The spearman correlation of each sample")
rm(list=ls())
?heatmap.2

##load data
rm(list=ls())
countData <- as.matrix(read.csv("gene_count_matrix.csv", row.names="gene_id"))
colData <- read.csv("group_info.txt", sep="\t", row.names=1)
##check
all(rownames(colData) %in% colnames(countData))
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))
##construct martix and get res
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData , design = ~ Condition)
dds <- DESeq(dds)
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
rld <- vst(dds, blind=FALSE)
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
pearson_cor <- as.matrix(cor(rlogMat, method="pearson"))
hc <- hcluster(t(rlogMat), method="pearson")
heatmap.2(pearson_cor, 
          Rowv=as.dendrogram(hc), 
          symm=T, 
          trace="none",
          col=hmcol, 
          colsep=c(6,12,19),
          rowsep=c(6,13,19),
          sepcolor="gray90",
          sepwidth=c(0.001,0.001),
          #margins=c(11,11),
          main="The pearson correlation of each sample")
#res <- results(dds)
pcadata <- plotPCA(rld,intgroup = c("Condition"),returnData = TRUE)
pcadata$Popluation <- colData$Group
pcadata$Tissue <- colData$Tissue
percentVar <- round(100 * attr(pcadata, "percentVar"))
ggplot(pcadata, aes(PC1, PC2, color=Condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

res_deseq2 <- function(x,y,dds2 = dds,p = 0.01,lfc = 1,condition = "Condition"){
    out_name = paste(x,y,sep = "-")
  out_csv <- paste("deg",out_name,"csv",sep = ".")
  out_list <- paste("deg",out_name,"genelst.tsv",sep = ".")
  out_csv2 <- paste("all",out_name,"csv",sep = ".")
  n_lfc = 0-lfc
  res <- results(dds2, contrast=c(condition,x,y))
  res_v <- res
  res <- res[order(res$padj < p),]
  diff_gene_deseq2 <-subset(res,padj < p)
  resdata <- merge(na.omit(as.data.frame(diff_gene_deseq2)),
                   as.data.frame(counts(dds2,normalize=TRUE)),
                   by="row.names",sort = T)
  resdata2 <- merge(na.omit(as.data.frame(res_v)),
                   as.data.frame(counts(dds2,normalize=TRUE)),
                   by="row.names",sort = T)
  resdata <- resdata[resdata$log2FoldChange < n_lfc|resdata$log2FoldChange > lfc,]
  resdata <- resdata[order(resdata$log2FoldChange),]
  write_csv(resdata,out_csv)
  write_csv(resdata2[,1:7],out_csv2)
  gene_list <- as.data.frame(resdata[,1])
  colnames(gene_list) <- "gene"
  write_tsv(gene_list,out_list)
}
compare_list <- read.table("comparelist.txt")
for(i in 1:nrow(compare_list)){
  print (paste(compare_list[i,1],compare_list[i,2]))
  res_deseq2(as.character(compare_list[i,1]),as.character(compare_list[i,2],dds2=dds))  
}

rm(list=ls())
library("clusterProfiler")
GOannotation <- read.delim("merge.go.clustprofile.out.txt", stringsAsFactors=FALSE)
GOannotation = split(GOannotation, with(GOannotation, level))
KOannotation <- read.delim("KOannotation.tsv", stringsAsFactors=FALSE)
GOinfo <- read.delim("go.tb", stringsAsFactors=FALSE)
go_enrich <- function(x,y,p = 0.05){
  out_name = paste(x,y,sep = "-")
  out_csv <- paste("deg",out_name,"clustprofilego.csv",sep = ".")
  out_csv2 <- paste("deg",out_name,"clustprofilekegg.csv",sep = ".")
  g_List <- paste("deg",out_name,"genelst.tsv",sep = ".")
  gene_list<- read.csv(g_List, sep="\t",header = F)# 你的gene lis
  egocc <- enricher( gene_list[,1],
                     TERM2GENE=GOannotation[['cellular_component']][c(2,1)],
                     TERM2NAME=GOinfo[1:2],qvalueCutoff = 1,pvalueCutoff = p)
  egobp <- enricher( gene_list[,1],
                     TERM2GENE=GOannotation[['biological_process']][c(2,1)],
                     TERM2NAME=GOinfo[1:2],qvalueCutoff = 1,pvalueCutoff = p)
  egomf <- enricher( gene_list[,1],
                     TERM2GENE=GOannotation[['molecular_function']][c(2,1)],
                     TERM2NAME=GOinfo[1:2],qvalueCutoff = 1,pvalueCutoff = p)
  mf <- egomf@result
  bp <- egobp@result
  cc <- egocc@result
  
  bp$level <- c("BP")
  cc$level <- c("CC")
  mf$level <- c("MF")
  
  res <- rbind(bp,mf,cc)
  jud <- res[res$pvalue<p,]
  kegg <- enricher(gene_list[,1],
          TERM2GENE=KOannotation[c(3,1)],
          TERM2NAME=KOannotation[c(3,4)]
          )
  dir.create(out_name)
  out_f <- file.path(out_name,out_csv)
  out_f2 <- file.path(out_name,out_csv2)
  kegg_res <- kegg@result
  write_csv(jud,out_f)
  write_csv(kegg_res,out_f2)
}
plot_va <- function(x,y){
  out_name = paste(x,y,sep = "-")
  out_plot <- paste("all",out_name,"va.pdf",sep = ".")
  data_n <- paste("all",out_name,"csv",sep = ".")
  data <- read.csv(data_n)
  data$color <- "a"
  data[which(data$padj < 0.01 & data$log2FoldChange > 1),8] = "b"
  data[which(data$padj < 0.01 & data$log2FoldChange < -1),8] = "c"
  p <- ggplot(data)+
    geom_point(aes(x=log2FoldChange,y=-log10(padj),color=color),alpha=0.8)+
    theme_test()+
    geom_vline(aes(xintercept = 1),linetype="dashed",colour="grey80")+
    geom_vline(aes(xintercept = -1),linetype="dashed",colour="grey80")+
    geom_hline(aes(yintercept = 2),linetype="dashed",colour="grey80")+
    scale_color_manual(
      values = c("grey70","#FEA3CA","#77CCCC")
    )+
    ylim(0,25)+
    xlim(-15,15)+
    theme(
      legend.position = "none"
    )
  ggsave(out_plot,p,width = 4,height = 6)
}
compare_list <- read.table("comparelist.txt")
for(i in 1:nrow(compare_list)){
  print (paste(compare_list[i,1],compare_list[i,2]))
  go_enrich(as.character(compare_list[i,1]),as.character(compare_list[i,2]))
  #plot_va(as.character(compare_list[i,1]),as.character(compare_list[i,2]))
}

########
dataA <- as.data.frame(read.csv("Cma_Go_plot.csv",header = T))
dataA <- dataA[order(dataA$level,dataA$pvalue),]
y<-factor(dataA$Description,levels = rev(dataA$Description))
x=-(log10(dataA$pvalue))
ggplot(dataA,aes(x,y))+
  geom_point(aes(size=Count,color=-log10(dataA$pvalue),shape=level))+
  scale_color_gradient(low = "blue", high = "red")+ 
  labs(color=expression(-Log.P.value.),x="-LogP",y="",title="fix it manually")+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"), 
        axis.text = element_text(color = "black",size = 14),
        legend.text = element_text(size = 14),
        legend.title=element_text(size=14),
        axis.title.x = element_text(size = 14))+
  scale_size_continuous(range=c(4,8))

go_p <- function(x){
  title = paste(x," Go enrichment")
  name = paste(x,"_Go_plot.csv",sep = "")
  dataA <- as.data.frame(read.csv(name,header = T))
  dataA <- dataA[order(dataA$level,dataA$pvalue),]
  y<-factor(dataA$Description,levels = rev(dataA$Description))
  x=-(log10(dataA$pvalue))
  ggplot(dataA,aes(x,y))+
    geom_point(aes(size=Count,color=-log10(dataA$pvalue),shape=level))+
    scale_color_gradient(low = "blue", high = "red")+ 
    labs(color=expression(-Log(P.value)),x="-Log(P.value)",y="",title="fix it manually")+
    theme_bw()+
    theme(axis.line = element_line(colour = "black"), 
          axis.text = element_text(color = "black",size = 14),
          legend.text = element_text(size = 14),
          legend.title=element_text(size=14),
          axis.title.x = element_text(size = 14))+
    scale_size_continuous(range=c(4,8))
}
go_p("Cma")
library(ggrastr)
dat<-data.frame(x=factor(c("ACR-ECR","ACR-ATR","ECR-ETR","ATR-ETR"),
                         levels =c("ACR-ECR","ACR-ATR","ECR-ETR","ATR-ETR")),
                y=0,
                label=c("ACR-ECR","ACR-ATR","ECR-ETR","ATR-ETR"))
data <- read_delim(file = "R.format.txt",delim="\t")
data$col2 <- factor(data$col2,levels = c("ACR-ECR","ACR-ATR","ECR-ETR","ATR-ETR"))
data <- filter(data,abs(LFC) < 15)
ggplot()+
  #geom_col(data=datbar,aes(x=x,y=y),fill="grey",alpha=0.5)+
  #geom_col(data=datbar,aes(x=x,y=-y),fill="grey",alpha=0.5)+
  geom_jitter_rast(data=data %>% filter(col1 == 'more'),
              aes(x=col2,y=LFC,color=col1),size=0.5,height = 0,alpha=0.5)+
  geom_jitter_rast(data=data %>% filter(col1 == 'less'),
              aes(x=col2,y=LFC,color=col1),size=0.5,height = 0,alpha=0.5)+
  scale_color_manual(name=NULL,
                     values = c("#f1434a","grey60")
                     )+
  scale_y_continuous(
    breaks = (seq(-15,15,3))
  )+
  theme_minimal()+
  labs(x="Cluster",y="average logFC")+
  geom_hline(yintercept = 1,linetype = "dashed")+
  geom_hline(yintercept = -1,linetype = "dashed")+
  geom_tile(data=dat,
            aes(x=x,y=y,fill=x),
            height=2,
            alpha=0.8,show.legend = F)+
  theme(axis.line.y = element_line(),
        axis.ticks.y = element_line(),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = c(1,0),
        legend.direction = "vertical",
        #axis.text.x = element_blank()
  )+
  scale_fill_manual(values = c("#ADE2D0","#4DBBD5","#00A087","#FAE48B"))+
  geom_text(data=dat,aes(x=x,y=y,label=label))


dat<-data.frame(x=factor(c("ACL-ECL","ACL-ATL","ECL-ETL","ATL-ETL"),
                         levels =c("ACL-ECL","ACL-ATL","ECL-ETL","ATL-ETL")),
                y=0,
                label=c("ACL-ECL","ACL-ATL","ECL-ETL","ATL-ETL"))
data <- read_delim(file = "L.format.txt",delim="\t")
data$col2 <- factor(data$col2,levels = c("ACL-ECL","ACL-ATL","ECL-ETL","ATL-ETL"))
data <- filter(data,abs(LFC) < 15)
ggplot()+
  #geom_col(data=datbar,aes(x=x,y=y),fill="grey",alpha=0.5)+
  #geom_col(data=datbar,aes(x=x,y=-y),fill="grey",alpha=0.5)+
  geom_jitter_rast(data=data %>% filter(col1 == 'more'),
                   aes(x=col2,y=LFC,color=col1),size=0.5,height = 0,alpha=0.5)+
  geom_jitter_rast(data=data %>% filter(col1 == 'less'),
                   aes(x=col2,y=LFC,color=col1),size=0.5,height = 0,alpha=0.5)+
  scale_color_manual(name=NULL,
                     values = c("#f1434a","grey60")
  )+
  scale_y_continuous(
    breaks = (seq(-15,15,3))
  )+
  theme_minimal()+
  labs(x="Cluster",y="average logFC")+
  geom_hline(yintercept = 1,linetype = "dashed")+
  geom_hline(yintercept = -1,linetype = "dashed")+
  geom_tile(data=dat,
            aes(x=x,y=y,fill=x),
            height=2,
            alpha=0.8,show.legend = F)+
  theme(axis.line.y = element_line(),
        axis.ticks.y = element_line(),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = c(1,0),
        legend.direction = "vertical",
        #axis.text.x = element_blank()
  )+
  scale_fill_manual(values = c("#ADE2D0","#4DBBD5","#00A087","#FAE48B"))+
  geom_text(data=dat,aes(x=x,y=y,label=label))

dat<-data.frame(x=factor(c("ACR-ECR","ACR-ATR","ECR-ETR","ATR-ETR"),
                         levels =c("ACR-ECR","ACR-ATR","ECR-ETR","ATR-ETR")),
                y=0,
                label=c("ACR-ECR","ACR-ATR","ECR-ETR","ATR-ETR"))
data <- read_delim(file = "R.format.txt",delim="\t")
data$col2 <- factor(data$col2,levels = c("ACR-ECR","ACR-ATR","ECR-ETR","ATR-ETR"))
data <- filter(data,abs(LFC) < 15)
ggplot()+
  #geom_col(data=datbar,aes(x=x,y=y),fill="grey",alpha=0.5)+
  #geom_col(data=datbar,aes(x=x,y=-y),fill="grey",alpha=0.5)+
  geom_jitter_rast(data=data %>% filter(col1 == 'more'),
                   aes(x=col2,y=LFC,color=col1),size=0.5,height = 0,alpha=0.5)+
  geom_jitter_rast(data=data %>% filter(col1 == 'less'),
                   aes(x=col2,y=LFC,color=col1),size=0.5,height = 0,alpha=0.5)+
  scale_color_manual(name=NULL,
                     values = c("#f1434a","grey60")
  )+
  scale_y_continuous(
    breaks = (seq(-15,15,3))
  )+
  theme_minimal()+
  labs(x="Cluster",y="average logFC")+
  geom_hline(yintercept = 1,linetype = "dashed")+
  geom_hline(yintercept = -1,linetype = "dashed")+
  geom_tile(data=dat,
            aes(x=x,y=y,fill=x),
            height=2,
            alpha=0.8,show.legend = F)+
  theme(axis.line.y = element_line(),
        axis.ticks.y = element_line(),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = c(1,0),
        legend.direction = "vertical",
        #axis.text.x = element_blank()
  )+
  scale_fill_manual(values = c("#ADE2D0","#4DBBD5","#00A087","#FAE48B"))+
  geom_text(data=dat,aes(x=x,y=y,label=label))


dat<-data.frame(x=factor(c("ACR-ECR","ACR-ATR","ECR-ETR","ATR-ETR"),
                         ReveRs =c("ACR-ECR","ACR-ATR","ECR-ETR","ATR-ETR")),
                y=0,
                RabeR=c("ACR-ECR","ACR-ATR","ECR-ETR","ATR-ETR"))
data <- read_deRim(file = "R.format.txt",delim="\t")
data$col2 <- factor(data$col2,levels = c("ACR-ECR","ACR-ATR","ECR-ETR","ATR-ETR"))
data <- filter(data,abs(LFC) < 15)
ggplot()+
  #geom_col(data=datbar,aes(x=x,y=y),fill="grey",alpha=0.5)+
  #geom_col(data=datbar,aes(x=x,y=-y),fill="grey",alpha=0.5)+
  geom_jitter_rast(data=data %>% filter(col1 == 'more'),
                   aes(x=col2,y=LFC,color=col1),size=0.5,height = 0,alpha=0.5)+
  geom_jitter_rast(data=data %>% filter(col1 == 'less'),
                   aes(x=col2,y=LFC,color=col1),size=0.5,height = 0,alpha=0.5)+
  scale_color_manual(name=NULL,
                     values = c("#f1434a","grey60")
  )+
  scale_y_continuous(
    breaks = (seq(-15,15,3))
  )+
  theme_minimal()+
  labs(x="Cluster",y="average logFC")+
  geom_hline(yintercept = 1,linetype = "dashed")+
  geom_hline(yintercept = -1,linetype = "dashed")+
  geom_tile(data=dat,
            aes(x=x,y=y,fill=x),
            height=2,
            alpha=0.8,show.legend = F)+
  theme(axis.line.y = element_line(),
        axis.ticks.y = element_line(),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = c(1,0),
        legend.direction = "vertical",
        #axis.text.x = element_blank()
  )+
  scale_fill_manual(values = c("#ADE2D0","#4DBBD5","#00A087","#FAE48B"))+
  geom_text(data=dat,aes(x=x,y=y,label=label))

library(tidyverse)
data <- read_delim("deg_stat.txt",delim = "\t")
data$Group = factor(data$Group,levels=(c("ASC-ESC","ASC-AST","ESC-EST","AST-EST")))
data <- filter(data,Tissue=="R")
data$Number <- log2(data$Number)
data[c(2,4,6,8),3] = -data[c(2,4,6,8),3]
ggplot(data)+
  geom_bar(aes(x=Group,y=Number,fill=Class),
           stat = 'identity',
           position = 'identity',
           alpha=0.3,width=0.9)+
  theme_minimal()+
  theme(axis.line.y = element_line(),
        axis.ticks.y = element_line(),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = c(1,0),
        legend.direction = "vertical",
        #axis.text.x = element_blank()
  )+
  scale_fill_manual(
    values = c("#A1D8B1","#F88AAF")
    #values = c("grey95","grey95")
  )+
  scale_y_continuous(breaks=seq(-12,12,3),
                     position = "right",
                     limits = c(-12,12)
                     )
  #geom_text(aes(x = Group ,
  #              label=Number, 
   #             y=Number + 200), 
  #          position=position_dodge(0.5), vjust=0.5)
  #facet_grid(Tissue~.)
