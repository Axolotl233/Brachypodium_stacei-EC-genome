args <- commandArgs(TRUE)
data <- read.table(args[1],header = F)
cut_line <- args[2]
TEST1 <- function(list_){
  m1 <- matrix(as.numeric(list_[c(2:5)]),ncol = 2)
  test <- chisq.test(m1,correct = T)
  p1 <- test$p.value
  return(p1)
}
TEST2 <- function(list_){
  m1 <- matrix(as.numeric(list_[c(2:5)]),ncol = 2)
  test <- fisher.test(m1)
  p1 <- test$p.value
  return(p1)
}
data$p_chisq <- apply(data,1,TEST1)
data$p_fisher <- apply(data,1,TEST2)
data$fdr_chisq <- p.adjust(data$p_chisq, method = "fdr")
data$fdr_fisher <- p.adjust(data$p_fisher, method = "fdr")
colnames(data) <- c("#gene","fixed_count","polymorphic_count","all_fixed","all_polymorphic","p_chisq","p_fisher","fdr_chisq","fdr_fisher")
write.table(data,file="select_gene_HKA.txt",row.names = F,col.names = T,quote = F,sep="\t")
data_sub <- data[data$fdr_fisher < cut_line,]
data_sub <- data_sub[data_sub$fixed_count > 0,]
write.table(data_sub,file="select_gene_HKA_sig.txt",row.names = F,col.names = T,quote = F,sep="\t")