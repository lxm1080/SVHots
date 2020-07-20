

library('argparse')
parser <- ArgumentParser(description="plot distriubution.jpg")
parser$add_argument('--path', type="character")
parser$add_argument('--file', type="character")
parser$add_argument('--superenhancertitle', type="character")
parser$add_argument('--enhancertitle', type="character")
parser$add_argument('--drivergenetitle', type="character")
parser$add_argument('--gwastitle', type="character")
parser$add_argument('--eQTLtitle', type="character")
args <- parser$parse_args()

path<-args$path
file<-args$file
superenhancertitle<-args$superenhancertitle
enhancertitle<-args$enhancertitle
drivergenetitle<-args$drivergenetitle
gwastitle<-args$gwastitle
eQTLtitle<-args$eQTLtitle



#用来画单个signature/type的格子图
#setwd('C:/Users/Lv/Desktop/annotation')


#path<-'Results/signature3/'
#file<-'Results/signature3/annotation_signature3.csv'
#superenhancertitle<-'gastric related super enhancer'
#enhancertitle<-'gastric related enhancer'
#drivergenetitle<-'STAD driver gene'
#eQTLtitle<-'gastric related eQTL'
#gwastitle<-'gastric cancer snp'

#annotation<-read.csv('Results/rs3_2_annotation.csv')
annotation<-read.csv(file)

name<-c() #colname
anno_heatmap<-cbind(annotation$num.gene,annotation$num.cancer.census.gene)
name<-c(name,'gene','cancer census gene')

#driver gene
anno_heatmap<-cbind(anno_heatmap,annotation$num.all.cancer.driver.gene)
name<-c(name,'all cancer driver gene')

if( gsub(" ",".",drivergenetitle) %in% colnames(annotation)){
  
  anno_heatmap<-cbind(anno_heatmap,annotation[,paste('num.',gsub(" ",".",drivergenetitle),sep='')])
  name<-c(name,drivergenetitle)
}

#super enhancer
anno_heatmap<-cbind(anno_heatmap,annotation$num.all.super.enhancer)
name<-c(name,'all super enhancer')

if( gsub(" ",".",superenhancertitle) %in% colnames(annotation)){
  
  anno_heatmap<-cbind(anno_heatmap,annotation[,paste('num.',gsub(" ",".",superenhancertitle),sep='')])
  name<-c(name,superenhancertitle)
}

#enhancer
anno_heatmap<-cbind(anno_heatmap,annotation$num.all.enhancer)
name<-c(name,'all enhancer')

if( gsub(" ",".",enhancertitle) %in% colnames(annotation)){
  
  anno_heatmap<-cbind(anno_heatmap,annotation[,paste('num.',gsub(" ",".",enhancertitle),sep='')])
  name<-c(name,enhancertitle)
}



#gwas
anno_heatmap<-cbind(anno_heatmap,annotation$num.all.gwas.snp)
name<-c(name,'all gwas snp')

if( gsub(" ",".",gwastitle) %in% colnames(annotation)){
  
  anno_heatmap<-cbind(anno_heatmap,annotation[,paste('num.',gsub(" ",".",gwastitle),sep='')])
  name<-c(name,gwastitle)
}

#eQTL

if( gsub(" ",".",eQTLtitle) %in% colnames(annotation)){
  
  anno_heatmap<-cbind(anno_heatmap,annotation[,paste('num.',gsub(" ",".",eQTLtitle),sep='')])
  name<-c(name,eQTLtitle)
}


#anno_heatmap<-cbind(annotation$num.allSEA, annotation$num.SEA_stomach.smooth.muscle,annotation$num.gene, annotation$num.drivergene, annotation$num.stomachdrivergene, annotation$num.cancer.census.gene, annotation$num.gastric.cancer.snp, annotation$num.eQTL)
#colnames(anno_heatmap)<-c("all super enhancer","related super enhancer","gene","all cancer driver gene","related driver gene","cancer census gene","GWAS SNP","eQTL")
colnames(anno_heatmap)<-name
sum=colSums(anno_heatmap)
#把每一列除以这一列的和
for( i in 1:ncol(anno_heatmap))
  if(sum[i]==0) {
    next
  }else{
    anno_heatmap[,i]<-anno_heatmap[,i]/sum[i]
  }
#rownames(anno_heatmap)<-annotation$cancer.census.gene #行名变成cancercensusgene,但是有的地方没有
rownames(anno_heatmap)<-annotation$hotspot.id #行名变成hotspot.id,但是有的地方没有
anno_heatmap <- t(anno_heatmap)#转置,转置以后class竟然变成了matrix
anno_heatmap<-as.data.frame(anno_heatmap)#转成数据框

#接下来按照那个热图教程画图
library(reshape2)
library(ggplot2)
pdf(paste(path,"distribution.pdf"), width=9, height=5)
# 转换前，先增加一列ID列，保存行名字
anno_heatmap$ID <- rownames(anno_heatmap)
# melt：把正常矩阵转换为长表格模式的函数。工作原理是把全部的非id列的数值列转为1列，命名为value；所有字符列转为variable列。
# id.vars 列用于指定哪些列为id列；这些列不会被merge，会保留为完整一列。
anno_heatmap <- melt(anno_heatmap, id.vars=c("ID"))
colnames(anno_heatmap)<-c('element','hotspot.id','proportion')
#ggplot(anno_heatmap, aes(x=hotspot.id,y=element))+ geom_tile(aes(fill=proportion)) + theme(axis.text.x=element_text(size=8,angle=90, hjust=1, vjust=1),axis.text.y=element_text(size=12))+ scale_fill_gradient(low = "white", high = 'springgreen4')+labs(title="distribution of elements in hotspot intervals")
ggplot(anno_heatmap, aes(x=hotspot.id,y=element))+ geom_tile(aes(fill=proportion)) + theme(axis.text.x=element_text(size=8,angle=90, hjust=1, vjust=1),axis.text.y=element_text(size=12))+ scale_fill_gradient(low = "white", high = 'springgreen4')+labs(title="distribution of elements in hotspot intervals")+ggsave( file = paste(path,"distribution.jpg"), width = 9, height = 5, type = "cairo", dpi = 600)
#ggplot(anno_heatmap, aes(x=variable,y=ID))+ geom_tile(aes(fill=value)) + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=1))+ scale_fill_gradient(low = "white", high = "springgreen4")+labs(title="distribution of elements in hotspot intervals")
dev.off()



# # data_m: 是前面费了九牛二虎之力得到的数据表
# # aes: aesthetic的缩写，一般指定整体的X轴、Y轴、颜色、形状、大小等
# # 在最开始读入数据时，一般只指定x和y，其它后续指定
# p <- ggplot(anno_heatmap, aes(x=variable,y=ID)) 
# 
# # 热图就是一堆方块根据其值赋予不同的颜色，所以这里使用fill=value, 用数值做填充色。
# p <- p + geom_tile(aes(fill=value))
# 
# # theme: 是处理图美观的一个函数，可以调整横纵轴label的选择、图例的位置等
# # 这里选择X轴标签45度。
# # hjust和vjust调整标签的相对位置，具体见 
# # 简单说，hjust是水平的对齐方式，0为左，1为右，0.5居中，0-1之间可以取任意值。vjust是垂直对齐方式，0底对齐，1为顶对齐，0.5居中，0-1之间可以取任意值
# p <- p + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=1))
# p <- p + scale_fill_gradient(low = "white", high = "black")
# dev.off()
