
preprocess<-function(gatric_rs1){

source('code/findhotspot/drawRandomRearrs.R')
source('code/findhotspot/runLightPcf.R')
source('code/findhotspot/drawSegmentation.R')
source('code/findhotspot/fastPCF.R')
source('code/findhotspot/extract.kat.regions2.R')
source('code/findhotspot/annotateBgModel.R')
source('code/findhotspot/calcIntermutDist.R')

gatric_rs1 <- gatric_rs1[c("Sample", "Type", "CHROM_1", "POS_1", "CHROM_2","POS_2")]
names(gatric_rs1)<-c("sample","Type","Chromosome.1","pos.1","Chromosome.2","pos.2")
gatric_rs1$Chromosome.1<-gsub("chr","",gatric_rs1$Chromosome.1)
gatric_rs1$Chromosome.2<-gsub("chr","",gatric_rs1$Chromosome.1)

x <- gatric_rs1$Type
x[which(x=='DEL')] <-2
x[which(x=='DUP')] <-4
x[which(x=='INV')] <-8
x[which(x=='BND')] <-32
gatric_rs1$pf<-x
#gatric_rs1$pf<-as.integer(gatric_rs1$pf)

#gatric_rs1$pf<-rep(8,nrow(gatric_rs1))

gatric_rs1$Type<-as.factor(gatric_rs1$Type)
gatric_rs1$sample<-as.factor(gatric_rs1$sample)
gatric_rs1$Chromosome.1<-as.factor(gatric_rs1$Chromosome.1)
gatric_rs1$Chromosome.2<-as.factor(gatric_rs1$Chromosome.2)


return(gatric_rs1)
}

#输入参数gatric_rs1, fp, gamma, kmin.my, imd.factor,name
process<-function(gatric_rs1, fp, gamma, kmin.my, imd.factor,name){
  gatric_rs1<-preprocess(gatric_rs1)
  bg.rate=nrow(gatric_rs1)/1400000000
  exp.num=bg.rate*0.5*1000000
  allBins<-read.csv('code/findhotspot/data/allBins.csv')
  x<-rep(exp.num,nrow(allBins))
  allBins$exp.col<-x   #在allBins里面添加一列，allBins还有待完善
  
  
  hotspots.gatric.rs1<- runLightPcf(rearr.df=gatric_rs1, # RS1 rearrangements 
                                    kmin=kmin.my,# kmin parameter for PCF
                                    gamma=gamma, # gamma parameter for PCF
                                    bg.rate=bg.rate, # genome-wide rate of RS1 breakpoints
                                    logScale=TRUE, # log-scale transformation of inter-mutation disances
                                    doMerging=TRUE, # merge adjacent hotspots
                                    obs.exp.thresh=imd.factor, # paramter for PCF
                                    allBins = allBins, # expected breakpoints/bin
                                    exp.col='exp.col', # use the RS1 column from allBins
                                    plot.path='../data/' # where to save a plot
  )
  hotspots.gatric.rs1$hotspot.id <- paste0('chr',hotspots.gatric.rs1$chr,'_', round(hotspots.gatric.rs1$start.bp/1e6,1), 'Mb' ) # give IDs to each hotspots

  drawSegmentation(gatric_rs1,  hotspots.gatric.rs1 , fn=paste0(fp, name,'.pdf'), name) # draw the rainfall plots
  #drawSegmentation(gatric_rs1,  hotspots.gatric.rs1 , fn=paste0(fp, 'obs.rs1.segmentation.i=',imd.factor , '-g=',gamma,'.pdf'), 'observed rs1') # draw the rainfall plots
 
  write.csv(hotspots.gatric.rs1, file=paste0(fp, name,'_hotspot', '.csv'), row.names = FALSE) # save table with hotspots
  #write.csv(hotspots.gatric.rs1, file=paste0(fp, 'obs.rs1.hotspots.i=',imd.factor, '-g=',gamma,'.csv'), row.names = FALSE) # save table with hotspots
}
