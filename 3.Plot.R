library(tidyverse)
library(cowplot)
library(ggrastr)
library(RColorBrewer)
library(ggpubr)

data <- read_table("Pop.stat.summary.txt")
data2 <- gather(data, key = "class",value = "count",-sample)
data2$class <- factor(data2$class,levels=rev(c("REF_homo","heter","ALT_homo","miss")))
c <-  colorRampPalette(brewer.pal(9, "GnBu"))(8)
ggplot(data2,aes(sample,count,fill= class))+
  geom_bar(stat='identity', position='stack')+
  theme_bw() +labs(x = 'Sample',y = 'Count')+
  theme(axis.title =element_text(size = 12),
        axis.text =element_text(size = 12, color = 'black'))+
  scale_fill_manual(
    values = c(c[8],c[6],c[4],c[2])
  )+
  theme(
    axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
    panel.grid = element_blank()
  )

data <- read.table("pixy_fst.phase.txt")
cut_line = 10
#colnames(data) <- c("chr","start","end","var_num","w_fst","m_fst")
colnames(data) <- c("chr","start","end","w_fst","var_num","fix_fst","color","pos")

for(i in 1:nrow(data)){
  if(data[i,5]<cut_line){
    data[i,6] = 0
  }
  if(data[i,6] < 0){
    data[i,6] = 0
  }
}

res <- data.frame(i = seq(1,max(data$var_num),1),s = 0,r=0)
a_s <- sum(data$var_num)
data$cum_group <- findInterval(data$var_num,seq(1,1000,1))
for(i in seq(1,max(data$var_num),1)){
  tmp_d <- data[data$cum_group == i,]
  tmp_s <- sum(tmp_d$var_num)
  tmp_r <- tmp_s/a_s
  res[i,2] = tmp_s
  res[i,3] = tmp_r
}
res$c_r <- cumsum(res$r)

ggplot(res)+
  geom_line(aes(x=i,y=c_r))+
  xlab("var_num_in_window")+
  ylab("cumulative_percentage")+
  scale_x_continuous(breaks = seq(0,max(data$var_num),20))+
  scale_y_continuous(breaks = seq(0,1,0.1))+
  theme_bw()
ggplot(data,aes(x=var_num))+
  geom_histogram(binwidth = 1,colour= "black",position = "identity",boundary=0)+
  scale_x_continuous(breaks = seq(0,30,5),limits = c(0,30))+
  theme_bw()
ggplot(data,aes(x=w_fst))+
  geom_histogram(binwidth = 0.01,colour= "black",position = "identity",boundary=0)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  scale_x_continuous(
    breaks = c(
      seq(0,1,0.5),
      round(as.numeric(quantile(data$w_fst,probs = c(0.95))),2),
      round(as.numeric(quantile(data$w_fst,probs = c(0.5))),2)
    )
  )+
  geom_vline(aes(xintercept = quantile(data$w_fst,probs = c(0.95))),colour="firebrick3",size=0.4,alpha = 2/3,linetype="dashed")+
  geom_vline(aes(xintercept = quantile(data$w_fst,probs = c(0.5))),colour="navy",size=0.4,alpha = 2/3,linetype="dashed")

line <- as.numeric(quantile(data$fix_fst,probs = c(0.95)))
datatop5 <- data[data$fix_fst > line,1:8]
#write_tsv("pixy_fst.phase.top5.txt",x =datatop5)

p1 <- ggplot(data)+
  geom_col(aes(x = pos,y=fix_fst ,fill = color),size = 1)+
  #geom_point_rast(aes(x = pos,y=fix_fst ,color = color))+
  theme_cowplot()+
  scale_fill_manual(  values = c("#BC6367","#ACCECD","#B7AE9D"))+
  theme(
    legend.position = "none"
  ) +
  labs(x = NULL,y="FST")+
  scale_x_continuous(breaks = c(1490,4306,6934,9433,11714,13909,16078,18213,20174,22011),
                     labels = c("Bs01","Bs02","Bs03","Bs04","Bs05",
                                "Bs06","Bs07","Bs08","Bs09","Bs10"))+
  geom_hline(aes(yintercept = quantile(w_fst,probs = c(0.95)),linetype = "c"),colour="gray20",size=0.7,alpha = 2/3)+
  scale_linetype_manual(values=c("dashed"))
p1
data <- read.table("pixy_dxy.phase.txt")
cut_line = 10
colnames(data) <- c("chr","start","end","dxy","var_num","fix_dxy","color","pos")
data$fix_dxy <- (data$dxy*data$var_num)/10000
for(i in 1:nrow(data)){
  if(data[i,5]<cut_line){
    data[i,6]=0
  }
}
line <- as.numeric(quantile(data$fix_dxy,probs = c(0.95)))
datatop5 <- data[data$fix_dxy > line,1:8]
#write_tsv("pixy_dxy.phase.top5.txt",x =datatop5)

data$dxy = 0-data$fix_dxy
p2 <- ggplot(data)+
  geom_col(aes(x = pos,y=dxy ,fill = color),size = 1)+
  theme_cowplot()+
  scale_fill_manual(  values = c("#BC6367","#ACCECD","#B7AE9D"))+
  theme(
    legend.position = "none"
  ) +
  labs(x = NULL,y="DXY")+
  scale_x_continuous(breaks = c(1490,4306,6934,9433,11714,13909,16078,18213,20174,22011),
                     labels = c("Bs01","Bs02","Bs03","Bs04","Bs05",
                                "Bs06","Bs07","Bs08","Bs09","Bs10"))+
  geom_hline(aes(yintercept = quantile(dxy,probs = c(0.05)),linetype = "c"),colour="gray20",size=0.7,alpha = 2/3)+
  scale_linetype_manual(values=c("dashed"))

plot_grid(p1,p2,nrow = 2,align = "hv")

data <- read_table("permute_window_snp_fst.0.05.region.merge",col_names = F)
colnames(data) <- c("chr","start","end")
data <- data %>%
  mutate(len=end-start)
data <- data %>% mutate(data,
         class = case_when(
           len <= 100000 ~ "short",
           TRUE ~ "long"
         )
  )

data %>% 
  group_by(chr) %>%
  summarise(count = sum(len/10000))
data %>% group_by(class) %>%
  summarise(count2 = sum(len),count = n())
data %>% group_by(class,chr) %>%
  summarise(tmp = n()) %>%
  spread(key=class,value=tmp)

t <- matrix(c(0,5,13,271),nrow = 2,byrow = T)
fisher.test(t)

ggplot(data= data)+
  geom_bar(aes(x = len))
155 + 58 + 190
data1 <- read_table("island.dxy.txt")
data1 %>% filter(class == "Island")
ggplot(data1,aes(x = class,y = avg_dxy,fill=class)) + 
  #geom_violin(trim = T)+
  geom_boxplot(width = 0.7)+
  stat_compare_means(comparisons =  list(c("Background","Island")) ,
                     method = "wilcox.test", 
                     label = "p.format") +
  labs(x="",y="value",title = "Dxy")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(
    values = c("#92c5de","#b2abd2")
  )

data2 <- read_table("island.pi.txt")
data2$fix_pi = (data2$avg_pi * data2$no_sites) / 10000
data2 <- data2 %>% filter(fix_pi < 0.015)
ggplot(data2,
       aes(x = class,y = fix_pi,fill=class)
       ) + 
  #geom_violin(trim = T)+
  geom_boxplot(width = 0.7)+
  stat_compare_means(comparisons =  list(c("Background","Island")) ,
                     method = "wilcox.test", 
                     label = "p.format") +
  labs(x="",y="value",title = "Pi")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(
    values = c("#92c5de","#b2abd2")
  )+
  scale_y_continuous(
    breaks = seq(0,0.015,0.002),limits = c(0,0.016)
  )+
  facet_wrap(pop~.)

data3 <- read_table("island.rho.txt",col_names = F)
colnames(data3) <- c("Chr","Start","End","Rho","Pop","Class")

ggplot(data3,
       aes(x = Class,y = Rho,fill=Class)
  )+ 
  #geom_violin(trim = T)+
  geom_boxplot()+
  stat_compare_means(comparisons =  list(c("Background","Island")) ,
                     method = "wilcox.test", 
                     label = "p.format") +
  labs(x="",y="value",title = "Pi")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(
    values = c("#92c5de","#b2abd2")
  )+
  scale_y_continuous(
  #  breaks = seq(0,0.015,0.002),limits = c(0,0.016)
  )+
  facet_wrap(Pop~.)

ggplot(data3,
       aes(x = Pop,y = Rho,fill=Pop)
)+ 
  #geom_violin(trim = T)+
  geom_boxplot()+
  stat_compare_means(comparisons =  list(c("AS","ES")) ,
                     method = "wilcox.test", 
                     label = "p.format") +
  labs(x="",y="value",title = "Pi")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(
    values = c("#92c5de","#b2abd2")
  )+
  scale_y_continuous(
    #  breaks = seq(0,0.015,0.002),limits = c(0,0.016)
  )

ggplot(data3,aes(x = Rho,fill =Pop,color = Pop))+ 
  #geom_violin(trim = T)+
  geom_histogram(position = "identity",boundary=0,binwidth = 1)+
  labs(x="Rho",y="Count")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(
    values = c("#EB977D","#818DB8")
  )+
  scale_color_manual(
    values = c("#EB977D","#818DB8")
  )+
  scale_x_continuous(
    limits = c(0,100),
    breaks = seq(0,100,10)
  )+
  facet_wrap(.~Pop,nrow = 2)

data <- read_table("2.Pop.PCA.evec.ggplot2.fix")
a <- ggplot(data,aes(x=PC1,y=PC2,color = species)) +
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon",
               aes(fill=species),
               alpha=0.1)+
  scale_fill_manual(values = c("#EB977D","#818DB8"))+
  scale_color_manual(values = c("#EB977D","#818DB8"))+
  theme_bw()+
  theme(panel.grid = element_blank())
b <- ggplot(data,aes(x=PC1,y=PC3,color = species)) +
  geom_point(size = 2)+
  stat_ellipse(geom = "polygon",
               aes(fill=species),
               alpha=0.1)+
  scale_fill_manual(values = c("#EB977D","#818DB8"))+
  scale_color_manual(values = c("#EB977D","#818DB8"))+
  theme_bw()+
  theme(panel.grid = element_blank())
plot_grid(a,b,nrow = 1)

data1 <- read_table("Pop.extract.2.result.ggplot")
data_tmp <- read_table("sort.txt",col_names = F)
data1$id <- factor(data1$id,levels = rev(data_tmp$X1))
a <- ggplot(data1)+
  geom_col(aes(x=id,y=percent,fill=k))+
  scale_fill_manual(values = c("#818DB8","#EB977D"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))

data2 <- read_table("Pop.extract.3.result.ggplot")
data_tmp <- read_table("sort.txt",col_names = F)
data2$id <- factor(data2$id,levels = rev(data_tmp$X1))
b <- ggplot(data2)+
  geom_col(aes(x=id,y=percent,fill=k))+
  scale_fill_manual(values = c("#818DB8","#EB6264","#EB977D"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))

data3 <- read_table("Pop.extract.4.result.ggplot")
data_tmp <- read_table("sort.txt",col_names = F)
data3$id <- factor(data3$id,levels = rev(data_tmp$X1))
c <- ggplot(data3)+
  geom_col(aes(x=id,y=percent,fill=k))+
  scale_fill_manual(values = c("#EB6264","#6ABFE0","#818DB8","#EB977D"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))


plot_grid(a,b,c,nrow = 3)
