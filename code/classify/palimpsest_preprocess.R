#' preprocessInput_sv
#'
#' Annotating the mutation data with necessary fields for further analysis
#' @param input_data Table describing somatic structural rearrangements.
#' @param ensgene Gene table for annotations. A table of Ensembl genes is provided with the package.
#' @param resdir Results directory where graphical outputs should be exported.
#'
#' @return
#' @export
#'
#' @examples
preprocessInput_sv <- function(input_data = NULL,ensgene = ensgene, resdir = resdir)
{
  Sample.col = "Sample"; CHROM_1.col = "CHROM_1"; CHROM_2.col = "CHROM_2"; POS_1.col = "POS_1"; POS_2.col = "POS_2"; type.col = "Type";
  chroms <- unique(input_data[, CHROM_1.col])
  if (1 %in% chroms == TRUE) {
    input_data[, CHROM_1.col] <- paste("chr", input_data[,CHROM_1.col], sep = "")
    input_data[, CHROM_2.col] <- paste("chr", input_data[,CHROM_2.col], sep = "")
  }
  input_data <- input_data[mixedorder(input_data[, CHROM_1.col]),
                           ]
  input_data$strand.mut <- "+"
  input_data$strand.gene <- NA
  print("Annotating mutation data:")
  vcf <- palimpsest_dfPosXSegm(input_data, dfPos.chrom.col = CHROM_1.col,
                               dfPos.pos.col = POS_1.col, ensgene, dfSegm.chrom.col = "Chromosome.Name",
                               dfSegm.start.col = "Gene.Start..bp.", dfSegm.end.col = "Gene.End..bp.",
                               colsToAdd = "Associated.Gene.Name", namesColsToAdd = "Associated.Gene.Name")
  print("Adding mutation categories:")
  vcf$strand.gene <- c("-", NA, "+")[ensgene[match(vcf$Associated.Gene.Name,
                                                   ensgene$Associated.Gene.Name), "Strand"] + 2]
  vcf <- palimpsest_addSVcategoriesToVcf(vcf, type.col = type.col,
                                             sample.col = Sample.col, CHROM_1.col = CHROM_1.col, CHROM_2.col = CHROM_2.col,
                                             POS_1.col = POS_1.col, POS_2.col = POS_2.col, resdir = resdir)
  return(vcf)
}
