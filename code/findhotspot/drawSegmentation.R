drawSegmentation <- function  (rearr.df,  hotspots.df,  fn, main.text) {

                                        # Draw the rainfall plots for rearrangement breakpoints
                                        #   rearr.df: data frame with rearrangements, one rearrangement per row
                                        #   hotspots.df: hotspots identified
                                        #   fn: file for saving PDF
                                        #   main.text: tile for PDF

    # make a data frame: one breakpoint per line
    rearr.bps <- rbind(data.frame(sample=rearr.df$sample, chr=rearr.df$Chromosome.1, position=rearr.df$pos.1, pf=rearr.df$pf),
                       data.frame(sample=rearr.df$sample, chr=rearr.df$Chromosome.2, position=rearr.df$pos.2, pf=rearr.df$pf)
                       )
    rearr.bps <- rearr.bps[order(rearr.bps$chr, rearr.bps$position),]

    intermut.dist <- calcIntermutDist(rearr.bps)

    pdf(fn, width=20, height=5)

    for (loop.chr in as.character(1:22)) {
        intermut.dist.chr <- subset(intermut.dist,chr==loop.chr)

        interDist = log10(intermut.dist.chr$distPrev);

        plot(intermut.dist.chr$position, interDist, main=paste(main.text, 'chromosome', loop.chr), pch=19, col='aquamarine3', cex=0.5)

        hotspots.chr <- subset(hotspots.df, chr==loop.chr )

        if (nrow(hotspots.chr)>0) {
            segments(x0=hotspots.chr$start.bp, y0=hotspots.chr$avgDist.bp,x1=hotspots.chr$end.bp, y1=hotspots.chr$avgDist.bp, lwd=2 )
        }

    }


    dev.off()

}

# 
# fn=paste0(fp, 'obs.rs1.segmentation.i=',imd.factor , '-g=',gamma,'.pdf')
# > rearr.df<-gatric_rs1#全部
# > hotspots.df<-hotspots.gatric.rs1
# > main.text<-'observed rs1'
# 
# rearr.bps <- rbind(data.frame(sample=rearr.df$sample, chr=rearr.df$Chromosome.1, position=rearr.df$pos.1, pf=rearr.df$Type),
#                    data.frame(sample=rearr.df$sample, chr=rearr.df$Chromosome.2, position=rearr.df$pos.2, pf=rearr.df$Type)
# )
# rearr.bps <- rearr.bps[order(rearr.bps$chr, rearr.bps$position),]
# 
# intermut.dist <- calcIntermutDist(rearr.bps)
# pdf("Results/tets3.pdf", width=10, height=5)
# 
# loop.chr="12"
# intermut.dist.chr <- subset(intermut.dist,chr==loop.chr)
# 
# intermut.dist.chr.del <- subset(intermut.dist.chr,pf==2)
# interDist.del = log10(intermut.dist.chr.del$distPrev)
# 
# intermut.dist.chr.dup <- subset(intermut.dist.chr,pf==4)
# interDist.dup = log10(intermut.dist.chr.dup$distPrev)
# 
# intermut.dist.chr.inv <- subset(intermut.dist.chr,pf==1|pf==8)
# interDist.inv = log10(intermut.dist.chr.inv$distPrev)
# 
# intermut.dist.chr.trn <- subset(intermut.dist.chr,pf==32)
# interDist.trn = log10(intermut.dist.chr.trn$distPrev)
# 
# points(intermut.dist.chr.del$position, interDist.del, xlim=c(0,1.5e+08),ylim=c(0,7),main=paste('nonclustered SV', 'chromosome', loop.chr), pch=19, col='blue', cex=0.5)
# plot(intermut.dist.chr.dup$position, interDist.dup, xlim=c(6.0e+07,1.1e+08),ylim=c(0,7),main=paste('nonclustered dup ', 'chromosome', loop.chr), pch=19, col='red', cex=0.5)
# plot(intermut.dist.chr.inv$position, interDist.inv, xlim=c(6.0e+07,1.1e+08),ylim=c(0,7),main=paste('nonclustered SV', 'chromosome', loop.chr), pch=19, col='yellow', cex=0.5)
# plot(intermut.dist.chr.trn$position, interDist.trn, xlim=c(6.0e+07,1.1e+08),ylim=c(0,7),main=paste('nonclustered trn', 'chromosome', loop.chr), pch=19, col='pink', cex=0.5)
# plot(intermut.dist.chr.del$position, interDist.del, xlim=c(6.0e+07,1.1e+08),ylim=c(0,7),main=paste('nonclusered del', 'chromosome', loop.chr), pch=19, col='blue', cex=0.5)
# plot(intermut.dist.chr$position, interDist, xlim=c(6.0e+07,1.1e+08),ylim=c(0,7),main=paste(main.text, 'chromosome', loop.chr), pch=19, col='blue', cex=0.5)
# dev.off()
