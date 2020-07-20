palimpsest_main <- function (sv_data)
{
source("code/classify/data.R")
source("code/classify/palimpsest_clonality.R")
source("code/classify/palimpsest_deconvolution.R")
source("code/classify/palimpsest_plotting.R")
source("code/classify/palimpsest_preprocess.R")
source("code/classify/palimpsest_timings.R")
source("code/classify/palimpsest_utils.R")
library(gtools)
library('ggplot2')
library('reshape2')
library(RColorBrewer)


#BiocManager::install("Biobase")
library(NMF)

datadir <- "/" # Path to directory containing data files
resdir <- "Results/classify";if(!file.exists(resdir))	dir.create(resdir) # Path to directory where to export results
#load(file.path(datadir,"sv_data.RData"))
#library(BSgenome.Hsapiens.UCSC.hg19) # Reference genome of choice
#ref_genome <- BSgenome.Hsapiens.UCSC.hg19
#sv_data<-read.csv('sv_data.csv',colClasses=c('character','character','character','numeric','character','numeric'))
load("code/classify/data/ensgene_hg19.RData") # Ensembl genes table (ensgene)
load("code/classify/data/cytoband_hg19.RData") # cytoband table (cyto)

resdir. <- file.path(resdir,"SV_signatures");if(!file.exists(resdir.)){dir.create(resdir.)}# Defining the results directory

library(bedr)
library(RCircos) # Loading dependencies necessary for SV annotation and CIRCOS plots
# To gain the functionality of bedr package you will need to have the BEDTools program installed and in your default PATH
data("UCSC.HG19.Human.CytoBandIdeogram")
# Preprocessing SV inputs and annotating for further analysis:
sv.vcf <- preprocessInput_sv(input_data =  sv_data,ensgene = ensgene,resdir = resdir.)
propSVsByCat <- palimpsestInput(vcf = sv.vcf,sample.col = "Sample",type="SV",mutcat.col = "Category",proportion = FALSE)
denovo_signatures <- deconvolution_nmf(input_data = propSVsByCat,type = "SV",range_of_sigs = 2:12,nrun =20,method = "brunet",resdir = resdir.)

# define list of colors for visualizing mutational signatures. Selecting default colors
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
mycol.sv <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
mycol.sv <- mycol.sv[sample.int(length(mycol.sv),nrow(denovo_signatures))];names(mycol.sv) <- rownames(denovo_signatures)

# Calculating contributions(exposures) of signatures in each sample:
SVsignatures_exp <- deconvolution_fit(vcf=sv.vcf,type = "SV",input_data = propSVsByCat,threshold = 6,input_signatures = denovo_signatures,sig_cols = mycol.sv,plot = T,resdir = resdir.)
#这里会报warnings()

# Estimate the probability of each event being due to each process
sv.vcf <- palimpsestOrigin(vcf=sv.vcf, type = "SV", sample.col = "Sample", mutcat.col = "Category", signature_contribution=SVsignatures_exp$sig_nums, input_signatures=denovo_signatures)
write.table(sv.vcf,file='Results/signature.csv',sep=',',row.names=F)

  
# Plotting the exposures of signatures across the series:
pdf(file.path(resdir.,"signature_content_plot.pdf"),width=10,height=7)
signature_content_plot <- deconvolution_exposure(SVsignatures_exp$sig_nums,SVsignatures_exp$sig_props,sig_cols = mycol.sv)#这一句会报错:Error in theme_bw(base_size = 14) : could not find function "theme_bw"
print(signature_content_plot)
dev.off()

return(sv.vcf)
}