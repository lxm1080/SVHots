source('code/findhotspot/preprocess.R')
library('argparse')
parser <- ArgumentParser(description="Calculation of the posterior pairing probability")
parser$add_argument('-f', '--find', type="character", help="type/category/signature need to find hotspot")
parser$add_argument('-k', '--kmin', type="integer", default=8,help="smallest sample number in every hotspot region")
parser$add_argument('-i', '--input', type="character", help="input file")
parser$add_argument('-c','--cluster', type="character", default='',help="1:clustered 0:nonclustered")

args <- parser$parse_args()

find<-args$find
kmin<-args$kmin
sv.vcf<-args$input



#kmin.my <- 4
#ÊäÈësv.vcf
sv.vcf<-read.csv(sv.vcf,stringsAsFactors=FALSE)


#clusterÊôÐÔ
if(args$cluster!=''){
  if(args$cluster=="0"){
    sv.vcf<-subset(sv.vcf, Clustered ==0 )#nonclustered
    
  }
  if(args$cluster=="1"){
    sv.vcf<-subset(sv.vcf, Clustered ==1 )
  }
}
#gatric_rs1<-subset(sv.vcf,as.logical(Signature.1.prob>0.5&Clustered ==0) )#Signature.1.prob>0.5&Clustered ==0




# where data is stored and output written
fp <- paste('Results/',find,'/',sep='')
# parameters of the PCF algorithm
gamma <- 8
imd.factor <- 2
name<-find

process(sv.vcf, fp, gamma, kmin, imd.factor,name)