#!/usr/bin/env Rscript
# @Date    : 2017-12-12 12:49:13
# @Author  : JT Guo
# @Email   : guojt-4451@163.com

# Args <- commandArgs(T)

# input = Args[1]  # with header: sample chr pos ref alt
# mc = Args[2]    # number of process
# outname = Args[3]   # outname

library("optparse")

option_list = list(
  make_option(c("-i", "--input"), type = "character", default = NULL, help = "the input file with header: sample chr pos ref alt",
    metavar = "file", dest = "input"), 
  make_option(c("-o", "--outname"), type = "character", default = NULL, help = "the output name",
    metavar = "str", dest = "outname"),
  make_option(c("-p", "--processNumber"), type = "character", default = 1, help = "the number of process",
    metavar = "int", dest = "mc"),
  make_option(c("-r", "--refSignature"), type = "character", default = NULL, help = "the file of reference signature",
    metavar = "file", dest = "refSignature"),
  make_option(c("-c", "--cutoff"), type = "character", default = 0.06, help = "Discard any signature contributions with a weight less than this amount[0.06]",
    metavar = "float", dest = "cutoff"))

opt_parser = OptionParser(usage = "usage: deconstructSigs [options]", 
    option_list = option_list, 
    add_help_option = TRUE, prog = NULL, description = "", epilogue = "")
opt = parse_args(opt_parser, print_help_and_exit=TRUE)

if (length(opt)==1){
  print_help(opt_parser)
  q()
}

library(deconstructSigs)
library(parallel)
library(stringr)

input = opt$input
outname = opt$outname
mc = opt$mc
if(is.null(opt$refSignature)){
    refSignture = signatures.cosmic
}else {
    refSignature = read.table(opt$refSignature, header=T, check.names=F, sep=",", stringsAsFactors=F)
}
cutoff = opt$cutoff
sigs.input <- mut.to.sigs.input(mut.ref = input, 
                                sample.id = "sample", 
                                chr = "chr", 
                                pos = "pos", 
                                ref = "ref", 
                                alt = "alt")

sigs <- mclapply(rownames(sigs.input), function(x){
    s = whichSignatures(tumor.ref = sigs.input, 
                        signatures.ref = refSignature, 
                        sample.id = x, 
                        contexts.needed = TRUE,
                        signature.cutoff = cutoff,
                        tri.counts.method = 'default')
}, mc.cores = mc)

pdf(paste0(outname,".pdf"), width = 6, height = 6)
for(i in 1:length(sigs)){
    plotSignatures(sigs[[i]], sub = outname)
    makePie(sigs[[i]], sub = outname)
}
dev.off()

sigs.weights <- sapply(sigs, function(x){
    unlist(x$weights)
})

rownames(sigs.weights) <- gsub("\\.","", rownames(sigs.weights))
colnames(sigs.weights) <- rownames(sigs.input)

save.image(paste0(outname, ".RData"))

