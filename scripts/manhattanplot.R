# title: 'Generate Manhattan plots and a table for GWAS summary data'
# author: 'Mei-Yu Lai'
# date: '08/08/2022'

library(getopt)
library(calibrate)
library(biomaRt)
library(vautils) # search nearest genes
library(tidyverse)

spec <- matrix(c('help', 'h', 0, 'logical', 'Help documentation',
                 'input.path', 'i', 1, 'character', 'Path of input files',
                 'plots.output.path', 'o', 1, 'character', 'Path for saving manhattan plots',
                 'table.output', 't', 1, 'character', 'Output file for a table'),
               byrow = T, ncol = 5)
opt <- getopt(spec)
if (!is.null(opt$help)){
  cat(getopt(spec, usage = T))
  q(status = 1)
}
if (!is.null(opt$input.path)) input.path <- opt$input.path
if (!is.null(opt$plots.output.path)) plots.output.path <- opt$plots.output.path
if (!is.null(opt$table.output)) table.output <- opt$table.output

## Use the default ENSEMBL Variation Mart & Human dataset
snpMart = useEnsembl(biomart = "snps", 
                     dataset = "hsapiens_snp",
                     version = 'GRCh37')

setwd('/Volumes/Benitez Lab/Meiyu/gwas_4785')

# files <- list.files(path = "rank8C/", pattern = "\\.linear$")
# proteins <- sub('.assoc.linear', '', files)
# files <- "ENTPD6.assoc.linear"
# proteins <- "ENTPD6"

manhattan.plot <- function(input.path, plots.output.path, table.output){
  files <- list.files(path = input.path, pattern = "\\.linear$")
  proteins <- sub('.assoc.linear', '', files)
  
  hits_table <- NULL
  for (p in 1:length(files)){
    x <- read.table(paste0(input.path, files[p]), header = T)
    d = data.frame(CHR = x$CHR, BP = x$BP, BETA = x$BETA, P = x$P, 
                   pos = NA, index = NA, SNP = x$SNP, stringsAsFactors = FALSE)
    d <- d[order(d$CHR, d$BP), ]
    d$logp <- -log10(d$P)
    d$index = rep.int(seq_along(unique(d$CHR)), times = tapply(d$SNP, d$CHR, length))
    nchr = length(unique(d$CHR))
    
    lastbase = 0
    ticks = NULL
    for (i in unique(d$index)) {
      if (i == 1) {
        d[d$index == i, ]$pos = d[d$index == i, ]$BP
      }
      else {
        lastbase = lastbase + max(d[d$index == (i - 1), "BP"])
        d[d$index == i, "BP"] = d[d$index == i, "BP"] - min(d[d$index == i, "BP"]) + 1
        d[d$index == i, "pos"] = d[d$index == i, "BP"] + lastbase
      }
    }
    ticks <- tapply(d$pos, d$index, quantile, probs = 0.5)
    xlabel = "Chromosome"
    labs <- unique(d$CHR)  
    
    xmax = ceiling(max(d$pos) * 1.03)
    xmin = floor(max(d$pos) * -0.03)
    def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", las = 1, pch = 20, 
                     xlim = c(xmin, xmax), ylim = c(0, max(8, ceiling(max(d$logp[is.finite(d$logp)])))), 
                     xlab = xlabel, ylab = expression(-log[10](italic(p))),
                     main = proteins[p])
    dotargs <- list()
    
    tiff(paste0(plots.output.path, proteins[p], '.tiff'), res = 300, width = 2500, height = 2000)
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
    axis(1, at = ticks, labels = labs)
    col = c("blue", "steelblue")
    col = rep_len(col, max(d$index))
    icol = 1
    for (i in unique(d$index)) {
      points(d[d$index == i, "pos"], d[d$index == i, "logp"], 
             col = col[icol], pch = 20)
      icol = icol + 1
    }
    
    # genomewideline
    abline(h = -log10(5e-08), col = "red", lty = 2)
    annotatePval <- 5e-08
    topHits = subset(d, P <= annotatePval)
    
    par(xpd = TRUE)
    topHits <- topHits[order(topHits$P), ]
    topSNPs <- NULL
    for (i in unique(topHits$CHR)) {
      chrSNPs <- topHits[topHits$CHR == i, ]
      topSNPs <- rbind(topSNPs, chrSNPs[1, ])
    }
    
    rs_table <- NULL
    if(!is.null(topSNPs)) {
      ## Submit the query
      for (j in 1:nrow(topSNPs)){
        snp <- unlist(strsplit(topSNPs[j, 'SNP'], ":"))
        CHR <- snp[1]; POS <- snp[2]; A1 <- snp[3]; A2 <- snp[4]
        BETA <- topSNPs[j, 'BETA']
        P <- topSNPs[j, 'P']
        res <- getBM(attributes = c('refsnp_id', 'chr_name', 'chrom_start', 'chrom_end', 'allele', 'associated_gene'),
                     filters = c('chr_name','start','end'), 
                     values = list(CHR, POS, POS), 
                     mart = snpMart)
        snp_1 <- res[which(res$chrom_start == POS & res$chrom_end == POS),]  # may have multiple SNPs at the same position
        near_genes <- find_nearest_gene(data = snp_1, flanking = 1000, build = "hg19",
                                        collapse = TRUE, snp = "refsnp_id", chr = "chr_name",
                                        bp = "chrom_start")
        nearest_genes <- ifelse(length(near_genes) != 0, near_genes$GENES, 'NA')
        rs_table <- rbind(rs_table, data.frame(Symbol = proteins[p], SNP = topSNPs[j, 'SNP'], snp_1[1,], BETA, P, nearest_genes))
      }
      
      symbol <- gsub("\\..*", "", rs_table$Symbol)
      df <- data.frame(symbol, rs_table$nearest_genes)
      rs_table$cis <- apply(df, 1, function(x){
        ifelse(grepl(x[1], x[2]), 1, 0)
      })
      
      hits_table <- rbind(hits_table, rs_table)
      
      textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625, 
             labs = rs_table[, 'refsnp_id'], cex = 0.5, srt = 90)
    }
    dev.off()
  }
  write.table(hits_table, table.output, sep = '\t', row.names = F, quote = F)
}
manhattan.plot(input.path, plots.output.path, table.output)




## manhattan and QQ plots
# for (i in 1:9){
#   ptm <- proc.time()
#   assc_results <- read.table(paste0("wgs_8_lrrk2_status/", proteins[i], "_lrrk2.assoc.linear"), head=TRUE)
#   
#   tiff(paste0("wgs_8_lrrk2_status/", proteins[i], "_lrrk2.tiff"), res = 300, height = 2000, width = 2500)
#   manhattan(assc_results, chr="CHR", bp="BP", p="P", snp="SNP", main = proteins[i],
#             col = c("blue4", "orange3"), ylim = c(0, max(10, max(-log10(assc_results$P)), na.rm = T)), annotatePval = 0.01)
#   dev.off()  
#   
#   # jpeg(paste0("wgs_8_results/QQplot_", i, ".jpeg"))
#   # qq(assc_results$P, main = "Q-Q plot of GWAS p-values : log")
#   # dev.off()
#   print(proc.time() - ptm)
# }
