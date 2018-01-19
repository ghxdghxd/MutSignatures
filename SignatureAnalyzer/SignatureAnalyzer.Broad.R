#!/usr/bin/env Rscript
# @Date    : 2018-01-17 17:00:02
# @Author  : JTGuo
# @Email   : guojt-4451@163.com

############################################################################################
############################################################################################
#### Copyright (c) 2017, Broad Institute
#### Redistribution and use in source and binary forms, with or without
#### modification, are permitted provided that the following conditions are
#### met:

####     Redistributions of source code must retain the above copyright
####     notice, this list of conditions and the following disclaimer.

####     Redistributions in binary form must reproduce the above copyright
####     notice, this list of conditions and the following disclaimer in
####     the documentation and/or other materials provided with the
####     distribution.

####     Neither the name of the Broad Institute nor the names of its
####     contributors may be used to endorse or promote products derived
####     from this software without specific prior written permission.

#### THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#### "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#### LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#### A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#### HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#### SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#### LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#### DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#### THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#### (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#### OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
############################################################################################
############################################################################################

######################################################################################################
####### SignatureAnalyzer - Mutation Signature Profiling using Bayesian NMF algorithms 
######################################################################################################
####### For details on the implementation 
####### see J Kim, Mouw K, P Polak et al, Somatic ERCC2 mutations are associated with a distinct genomic signature in urothelial tumors 
####### Nat. Genet. (2016) DOI: 10.1038/ng.3557
####### For details on the original algorithms 
####### see Tan, V.Y. & Févotte, C. Automatic relevance determination in nonnegative matrix factorization with the beta-divergence.
####### IEEE Trans. Pattern Anal. Mach. Intell. 35, 1592–1605 (2013).

######################################################################################################
####### The official name of this program is "SignatureAnalyzer" so please cite our ERCC2 paper above with this official name 
####### when you need to cite this program. 
######################################################################################################

library("optparse")

option_list = list(
  make_option(c("-i", "--input"), type = "character", default = NULL, help = "the maf file",
  metavar = "file", dest = "maf"), 
  make_option(c("-n", "--name"), type = "character", default = NULL, help = "the disease name",
  metavar = "str", dest = "name"),
  make_option(c("-o", "--outdir"), type = "character", default = NULL, help = "the OUTPUT dir",
  metavar = "dir", dest = "outdir"))

opt_parser = OptionParser(usage = "usage: %prog [options]", 
    option_list = option_list, 
    add_help_option = TRUE, prog = NULL, description = "", epilogue = "")
opt = parse_args(opt_parser, print_help_and_exit=TRUE)

if (length(opt)==1){
  print_help(opt_parser)
  q()
}

getScriptLocation = function() {
  cmd.args = commandArgs(trailingOnly = FALSE)
  cmd.args = cmd.args[grep("--file=", cmd.args)]
  return(dirname(gsub("--file=", "", cmd.args)))
}

locateDir = paste0(getScriptLocation(), "/")

source(paste(locateDir, "Functions.R", sep = ""))
source(paste(locateDir, "get.context96_annotated.from.maf.R", sep = ""))

if(is.null(opt$outdir)){
  CURRENT <- paste(getwd(), "/", sep = "")
  opt$outdir <- paste(CURRENT, "OUTPUT_lego96/", sep = "")
}



system(paste("mkdir", opt$outdir, sep = " "))

library(gplots)
library(ggplot2)
library(reshape2)
############## read COSMIC signatures ######################################
sanger <- read.delim(paste0(locateDir, "signature.new.COSMIC.txt"), header = T, sep = "\t", as.is = T)
base.change <- paste(substring(sanger[, 1], 1, 1), substring(sanger[, 1], 3, 3), sep = "")
context <- paste(substring(sanger[, 2], 1, 1), substring(sanger[, 2], 3, 3), sep = "")
rownames(sanger) <- paste(base.change, context, sep = "")
sanger <- sanger[, c(4:33)]
sanger <- sanger[match(context96.label, rownames(sanger), nomatch = 0), ]

##########################################################

############## INPUT ######################################
## maf file (tab-delimited) should have the following column headers - "sample", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele2", "ref_context"
## "sample" - tumor_sample_name, "Variant_Type" == SNP , "ref_context" = reference sequence at mutated base (at least tri-nucleotide sequence)
## lego96 - mutation counts matrix (96 by # of samples) 
## This matrix should contain mutation counts along 96 tri-nucleotide mutation contexts (rows) across samples (columns). 
## Rownames of the lego matrix should be 4-letters ex) CGTA (C to G mutation at 5'-TCA-3'contexts) (see the acoompanied example lego matrix).
###########################################################

maf <- read.delim(opt$maf, header = T, sep = "\t", as.is = T, comment = "#")
x <- get.spectrum96.from.maf(maf)  ### getting input lego matrix (96 by # of samples)
maf <- x[[1]]
lego96 <- x[[2]]

#####################################################
############## BayesNMF parameters
############## n.iter = number of independent simulations
############## Kcol = number of initial signatures
############## tol = tolerance for convergence
############## a0 = hyper-parameter
n.iter <- 10
Kcol <- 96
tol <- 1e-07
a0 <- 10
tumor.type <- opt$name  ### please specify your cohort name here
##################################

##################################
############### Choose pirors for W and H
############### Default = L1KL (expoential priors); L2KL (half-normal priors)
##################################
prior <- "L1KL"
if (prior == "L1KL") {
  method <- paste("L1KL.lego96", tumor.type, sep = ".")
} else {
  method <- paste("L2KL.lego96", tumor.type, sep = ".")
}
##################################

##################################
############## Default = FALSE ; TRUE - to reduce the effect of hyper-mutant samples in the signature discovery
##################################
hyper <- FALSE
if (hyper) {
  lego96 <- get.lego96.hyper(lego96)
  method <- paste(method, "hyper", sep = ".")
}
##################################

##########################################################
###################### Running the algorithms ############
##########################################################
for (i in 1:n.iter) {
  if (prior == "L1KL") {
    res <- BayesNMF.L1KL(as.matrix(lego96), 1e+05, a0, tol, Kcol, Kcol, 1)
  } else {
    res <- BayesNMF.L2KL(as.matrix(lego96), 1e+05, a0, tol, Kcol, Kcol, 1)
  }
  save(res, file = paste(opt$outdir, paste(method, a0, i, "RData", sep = "."), sep = ""))
}

##########################################################
###################### Analysis ##########################
##########################################################
res.WES <- get.stats.simulation(tumor.type, n.iter, opt$outdir)

############## frequency figure - this plot shows a distribution of # of signatures..
pdf(file = paste(opt$outdir, paste(method, a0, "signature.freq.pdf", sep = "."), sep = ""), width = 4, height = 5)
s1 <- 1.5
s2 <- 2
par(mfrow = c(1, 1))
par(mar = c(5, 5, 2, 1))
barplot(table(unlist(res.WES[[4]])), cex = s1, cex.axis = s1, cex.main = s1, cex.names = s1, cex.lab = s1,
  xlab = "# of signatures", ylab = "Freq.", main = paste(tumor.type, sep = "."))
dev.off()

##########################################################
############## select the best solution (maximum posteria solution) for given K
##########################################################
tmpK <- unlist(res.WES[[4]])
unique.K <- sort(unique(tmpK))
n.K <- length(unique.K)
MAP <- list()
for (i in 1:n.K) {
  tmpX <- res.WES[[5]]
  tmpX[tmpK != unique.K[i]] <- 1e+06
  MAP[[i]] <- which.min(tmpX)
}

tmpK <- unlist(res.WES[[4]])
unique.K <- sort(unique(tmpK))
n.K <- length(unique.K)
MAP <- list()
for (i in 1:n.K) {
  tmpX <- res.WES[[5]]
  tmpX[tmpK != unique.K[i]] <- 1e+06
  MAP[[i]] <- which.min(tmpX)
}
MAP <- unlist(MAP)
names(MAP) <- unique.K
MAP.nontrivial <- MAP[names(MAP) != 1]
##########################################################
##########################################################

n.K <- length(MAP.nontrivial)
if (n.K > 0) { 
for (j in 1:n.K) {
  load(paste(opt$outdir,paste(method,a0,MAP.nontrivial[[j]],"RData",sep="."),sep=""))
  W <- res[[1]]
  H <- res[[2]]
  W1 <- W
  H1 <- H
  W.norm <- apply(W,2,function(x) x/sum(x))
  for (i in 1:ncol(W)) {
          H1[i,] <- H1[i,]*colSums(W)[i]
          W1[,i] <- W1[,i]*rowSums(H)[i]
  }
  lambda <- res[[5]]
  df <- get.df.solution(W1,H,lambda,tumor.type)
  K <- length(table(df$class))

  ############# Signature plot
  width <- 16
  height <- ifelse(K==1,3,K*2)
  pdf(file = paste0(opt$outdir, paste(method, a0, paste0("MAP", K), "signature.pdf", sep = ".")), width = width,
    height = height)
  p <- plot.lego.observed.barplot(df, paste("Mutation Signatures in ", tumor.type, sep = ""))
  plot(p)
  dev.off()

  index.nonzero <- colSums(W) != 0
  lambda <- unlist(res[[5]][length(res[[5]])])
  lambda <- lambda/min(lambda)
  lambda <- lambda[index.nonzero]
  names(lambda) <- paste("W", seq(1:length(lambda)), sep = "")
  K <- sum(index.nonzero)
  if (K != 1) {
    W0.tumor <- W[, index.nonzero]
    H0.tumor <- H[index.nonzero, ]
    x <- W0.tumor
    y <- H0.tumor
    x.norm <- apply(x, 2, function(x) x/sum(x))
    W1.tumor <- x.norm
    for (j in 1:K) y[j, ] <- y[j, ] * colSums(x)[j]
    H1.tumor <- y
    y.norm <- apply(y, 2, function(x) x/sum(x))
    H2.tumor <- y.norm
  } else {
    stop("No non-trivial solutations; All simulations converged to a trivial solution with one signature")
  }

  ############# Reconstructing the activity of signatures
  W.mid <- W1.tumor
  H.mid <- H1.tumor
  H.norm <- H2.tumor
  if (length(grep("__", colnames(H.mid))) != 0) {
    hyper <- colnames(H.mid)[grep("__", colnames(H.mid))]
    H.hyper <- H.mid[, colnames(H.mid) %in% hyper]
    H.nonhyper <- H.mid[, !(colnames(H.mid) %in% hyper)]
    sample.hyper <- colnames(H.hyper)
    sample.hyper <- sapply(sample.hyper, function(x) strsplit(x, "__")[[1]][[1]])
    unique.hyper <- unique(sample.hyper)
    n.hyper <- length(unique.hyper)
    x.hyper <- array(0, dim = c(nrow(H.hyper), n.hyper))
    for (i in 1:n.hyper) {
      x.hyper[, i] <- rowSums(H.hyper[, sample.hyper %in% unique.hyper[i]])
    }
    colnames(x.hyper) <- unique.hyper
    rownames(x.hyper) <- rownames(H.mid)
    H.mid <- (cbind(H.nonhyper, x.hyper))
    H.norm <- apply(H.mid, 2, function(x) x/sum(x))
  }
  W.norm <- apply(W.mid, 2, function(x) x/sum(x))

  ##########################################################
  ############# W.norm = extracted signatures normalized to one
  ############# H.mid = activity of signatures across samples (expected mutations associated with signatures)
  ############# H.norm = normalized signature activity 
  ##########################################################
  WH <- list(W.norm, H.mid, H.norm)
  save(WH, file = paste(opt$outdir, paste(method, a0, paste("MAP", K, sep = ""), "WH.RData", sep = "."),
    sep = ""))

  ############# Activity plot
  p1 <- plot.activity.barplot(H.mid,H.norm,1.0,tumor.type)
  pdf(file = paste0(opt$outdir, paste(method, a0, "activity.barplot1", K, "pdf", sep = ".")), width = 15,
    height = 12)
  plot(p1)
  dev.off()

  main <- paste("Cosine similarity;",tumor.type,sep="")
  pdf(file = paste0(opt$outdir, paste(method, a0, "signature.comprison.sanger", K, "pdf", sep = ".")),
    width = 4, height = 6)
  s1 <- 1.5
  s2 <- 2
  par(mfrow = c(1, 1))
  par(mar = c(0, 2, 1, 1))
  x <- W.mid[1:96, ]
  colnames(x) <- paste("W", c(1:K), sep = "")
  plot.heatmap.2(plot.W.correlation(sanger, x), T, T, "ward.D", main)
  dev.off()
  }
}
