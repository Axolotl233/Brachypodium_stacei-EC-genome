rm(list=ls())
library(tidyverse)
library(ggrepel)
data <- read.csv("z.filter_TE.pca.txt",row.names=1)
data2 <- data[,colMeans(data[,-ncol(data)]) > 0]
data2 <- data2[ , which(apply(data2, 2, var) != 0)]
data2[,ncol(data2)] <- data[,ncol(data)]
pca1 <- prcomp(data2[,-ncol(data2)],center = TRUE,scale. = TRUE)
df1 <- pca1$x
df1 <- as.data.frame(df1)
summ1 <- summary(pca1)
xlab1 <- paste0("PC1(",round(summ1$importance[2,1]*100,2),"%)")
ylab1 <- paste0("PC2(",round(summ1$importance[2,2]*100,2),"%)")
ggplot(data = df1,aes(x = PC1,y = PC2,color =data2[,ncol(data)])+
  stat_ellipse(aes(fill = data2[,ncol(data)]),
               type = "norm",geom = "polygon",alpha = 0.25,color = NA)+ # 添加置信椭圆
  geom_point(size = 3.5)+
  labs(x = xlab1,y = ylab1,color = "Condition",title = "PCA Scores Plot")+
  guides(fill = "none")+
  theme_bw()+
  scale_fill_manual(values = c("#EB977D","#818DB8"))+
  scale_colour_manual(values = c("#EB977D","#818DB8"))+
  theme(
    panel.grid=element_blank()
  )

rm(list=ls())
library(ggpubr)
library(ggsignif)
AS_co <- c("ACL","ATL","ACR","ATR")
ES_co <- c("ECL","ETL","ECR","ETR")
L <- c("ACL","ECL","ATL","ETL")
R <- c("ACR","ECR","ATR","ETR")

data <- read_table("insertions.fre.stat.gene.diff.tpm.txt")
data_l <- select(data,gene,L,highfre)
data_l$mean <-rowMeans(data_l[,L])
data_l$max <-apply(data_l[,L],1,max)
data_r <- select(data,gene,R,highfre)
data_r$mean <-rowMeans(data_r[,R])
data_r$max <-apply(data_r[,R],1,max)

data_l <- filter(data_l,mean>0.5) %>% select (-mean,-max)
data_r <- filter(data_r,mean>0.5) %>% select (-mean,-max)

data_t <- gather(data=data_l,key="class",value="value",-gene,-highfre)
data_t$class2 <- rep(c("Control","Treatment"),each=12)
data_t$class3 <- rep(c("AS","ES"),2,each = 6 )

ggplot(data_t,aes(class2,value))+
  geom_bar(aes(fill=class3),width = 0.7,stat = "identity",position="dodge")+
  facet_wrap(~gene,scales = "free")+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid = element_blank()
  )+
  labs(x="Leaf",y="FPKM")+
  scale_fill_manual(
    values = c("#E9967C","#808BB7")
  )

data_t <- gather(data=data_r,key="class",value="value",-gene,-highfre)
data_t$class2 <- rep(c("Control","Treatment"),each=16)
data_t$class3 <- rep(c("AS","ES"),2,each = 8 )

ggplot(data_t,aes(class2,value))+
  geom_bar(aes(fill=class3),width = 0.7,stat = "identity",position="dodge")+
  facet_wrap(~gene,scales = "free")+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid = element_blank()
  )+
  labs(x="Root",y="FPKM")+
  scale_fill_manual(
    values = c("#E9967C","#808BB7")
  )

data_t <- data_l[,2:3]
data_t <- as.data.frame(t(scale(t(data_t),center = F)))
data_t <- dplyr::bind_cols(data_t,data_l[,c(1,6)])
data_t <- tidyr::gather(data=data_t,key="class",value="value",-gene,-highfre)

ggplot(data_t,(aes(class,value)))+
  geom_line(aes(group=gene),linetype="dashed",color="grey80")+
  geom_point(aes(color=class,fill=class,group=gene,shape=highfre),
             position = position_dodge(0.05),size=2)+
  labs(x="class",y="normalized TPM",title = "comparsion")+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid = element_blank()
  )