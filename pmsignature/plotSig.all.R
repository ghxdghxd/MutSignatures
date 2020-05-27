Args <- commandArgs(T)

sig_RData = Args[1]       # {Cancer_Type}/{Cancer_Type}.sig.{Signature_num}.RData

out_pdf = gsub("RData", "pdf", sig_RData)

# sig_RData = "LUAD/LUAD.sig.5.RData"

library(GenomicRanges, BSgenome.Hsapiens.UCSC.hg19)
library(pmsignature)
library(ggplot2)
library(cowplot)
library(reshape)
library(lsa)
library(dplyr)
library(egg)
library(RColorBrewer)

load(sig_RData)
fontsize = 20


getMutSignatures <- function(Param, sinInd) {
  numBases = Param@flankingBasesNum
  centerBase <- (1 + numBases) / 2
  vF <- Param@signatureFeatureDistribution[sinInd, , ]
  v1 <- vF[1, 1:6]
  v2 <- vF[2:(numBases), 1:4]
  A <- matrix(0, numBases, 4)
  B <- matrix(0, 4, 4)
  
  for (l in 1:numBases) {
    if (l < centerBase) {
      A[l, ] <- v2[l, ]
    } else if (l > centerBase) {
      A[l, ] <- v2[l - 1, ]
    }
  }
  A[centerBase, 2] <- sum(v1[1:3])
  A[centerBase, 4] <- sum(v1[4:6])
  
  B[2, c(1, 3, 4)] <- v1[1:3] / sum(v1[1:3])
  B[4, c(1, 2, 3)] <- v1[4:6] / sum(v1[4:6])
  
  A <- as.data.frame(A[c(2, 3, 4), ])
  B <- as.data.frame(B)
  colnames(A) = colnames(B) = rownames(B) = c("A", "C", "G", "T")
  AC <- c()
  AT <- c()
  for(i in 1:4){
    for(j in 1:4){
      AC <- c(AC, A[1, i]*A[3, j]*A[2, 2])
      AT <- c(AT, A[1, i]*A[3, j]*A[2, 4])
    }
  }
  BC <- c()
  BT <- c()
  for(k in c(1, 3, 4)){
    BC <- c(BC, B[2, k]*AC)
  }
  for(k in c(1, 2, 3)){
    BT <- c(BT, B[4, k]*AT)
  }
  D <- c(BC, BT)
  D <- as.data.frame(D)
  colnames(D) <- paste0("S", sinInd)
  rownames(D) <- c("CA A.A", "CA A.C", "CA A.G", "CA A.T", "CA C.A", "CA C.C", "CA C.G", "CA C.T", "CA G.A", "CA G.C", "CA G.G", 
                   "CA G.T", "CA T.A", "CA T.C", "CA T.G", "CA T.T", "CG A.A", "CG A.C", "CG A.G", "CG A.T", "CG C.A", "CG C.C", 
                   "CG C.G", "CG C.T", "CG G.A", "CG G.C", "CG G.G", "CG G.T", "CG T.A", "CG T.C", "CG T.G", "CG T.T", "CT A.A", 
                   "CT A.C", "CT A.G", "CT A.T", "CT C.A", "CT C.C", "CT C.G", "CT C.T", "CT G.A", "CT G.C", "CT G.G", "CT G.T", 
                   "CT T.A", "CT T.C", "CT T.G", "CT T.T", "TA A.A", "TA A.C", "TA A.G", "TA A.T", "TA C.A", "TA C.C", "TA C.G",
                   "TA C.T", "TA G.A", "TA G.C", "TA G.G", "TA G.T", "TA T.A", "TA T.C", "TA T.G", "TA T.T", "TC A.A", "TC A.C",
                   "TC A.G", "TC A.T", "TC C.A", "TC C.C", "TC C.G", "TC C.T", "TC G.A", "TC G.C", "TC G.G", "TC G.T", "TC T.A", 
                   "TC T.C", "TC T.G", "TC T.T", "TG A.A", "TG A.C", "TG A.G", "TG A.T", "TG C.A", "TG C.C", "TG C.G", "TG C.T", 
                   "TG G.A", "TG G.C", "TG G.G", "TG G.T", "TG T.A", "TG T.C", "TG T.G", "TG T.T")
  return(as.matrix(D))
}

getMatSignatures <- function(Param, BG = NULL, reorderSig = NULL, Signame = NULL){
    signatureNum <- Param@signatureNum
  if (Param@isBackGround){
    BG_signatureNum = signatureNum - 1
  }else{
    BG_signatureNum = signatureNum
  }
  
  if (is.null(reorderSig)) {
    reorderSig <- 1:BG_signatureNum
  }
  sigOrder <- reorderSig
  
  sig_mm <- getMutSignatures(Param, sigOrder[1])
  if(BG_signatureNum > 1){
      for(i in 2:BG_signatureNum){
        sig_mm <- cbind(sig_mm, getMutSignatures(Param, sigOrder[i]))
      }
  }
  if(!is.null(BG_prob)){
    bg_name = unique(substr(names(BG_prob), 1, 9))
    bg_name_2 = sort(c(paste(bg_name, "1", sep=","), paste(bg_name, "2", sep=",")))
    BG_new = BG_prob[match(bg_name_2, names(BG_prob))]
    names(BG_new) = bg_name_2
    BG <- matrix(BG_new, ncol=2, byrow=T)
    rownames(BG) <- bg_name
    BG[is.na(BG)] <- 0
    BG <- BG/sum(BG) # normalize to 1
    BG <- as.data.frame(apply(BG, 1, sum))
    BG$name <- sapply(rownames(BG),function(x){
      name = as.numeric(unlist(strsplit(x, ",")))
      a = c("CA", "CG", "CT", "TA", "TC", "TG")
      b = c("A", "C", "G", "T")
      c = c("+", "-")
      return(paste(a[name[1]], paste(b[name[median(1:length(name))]], b[name[median(1:length(name))+1]], sep = "."), sep = " "))
    })
    BG <- sapply(unique(BG$name), function(x){
      return(sum(BG[which(BG$name==x),1]))
    })
    sig_mm <- cbind(sig_mm, Sbg=BG[match(rownames(sig_mm), names(BG))])
  }
  
  if (!is.null(Signame)){
    colnames(sig_mm) <- Signame
  }
  return(sig_mm)
}

sig_mat <- getMatSignatures(Param, BG = BG_prob)

plotSignatures <- function(sig_mat, charSize = 10, scales = "fixed", scale_y_sqrt = FALSE){
  level = colnames(sig_mat)
  sig_mat <- cbind(sig_mat, Total_SNVs = rowSums(sig_mat))
  w_df = melt(sig_mat, varnames = c("motif", "signature"))
  w_df$alteration = sub("([ACGTN])([ACGTN]) .+", "\\1>\\2", w_df$motif)
  # w_df$context = sub("[ACGTN][ACGTN] (.+)", "\\1", w_df$motif)
  w_df$context = sub("([ACGTN])([ACGTN]) ([ACGTN]).([ACGTN])", "\\3\\1\\4", w_df$motif)
  w_df$signature = factor(w_df$signature, levels = c("Total_SNVs", level))
  
  p = ggplot(w_df)
  p = p + geom_bar(aes_string(x = "context", y = "value", fill = "alteration"),
                   stat = "identity", position = "identity")
  # stat = "identity", position = "identity",  colour="black",size=0.1)
  p = p + facet_grid(signature ~ alteration, scales = scales, switch = "y")
  p = p + theme_bw()
  p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = charSize*0.6),
                axis.text.y = element_text(hjust = 0.5, size = charSize),
                axis.title=element_text(size = charSize*1.2),
                panel.background = element_blank(),
                strip.text = element_text(size = charSize),
                axis.ticks.x = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                legend.position = "none",
                panel.spacing.y = unit(1, "lines"),
                plot.margin=unit(c(0.5,0.5,0.5,0.5), "lines"))
  p = p + scale_fill_manual(values = alpha(c("Sky Blue","black", "firebrick3", 
                                             "gray","lightgreen", "lightpink")))
  p = p + xlab("Tri-nucleotide Sequence Motifs") + ylab("Activity")
  if(scale_y_sqrt == TRUE){
    p = p + scale_y_sqrt() 
  }
  return(p)
}

p_sigs = plotSignatures(sig_mat, charSize = fontsize, scales = "free", scale_y_sqrt = FALSE)


rownames(sig_mat)<- sub("([ACGTN])([ACGTN]) ([ATCG]).([ATCG])", "\\3[\\1>\\2]\\4", rownames(sig_mat))

cosmicSig <- read.table("/home/jintao/Projects/TCGA/SignatureAnalyzer/signature.new.COSMIC.txt", header = T, sep = "\t",
    stringsAsFactors = F)
rownames(cosmicSig) <- cosmicSig$Somatic.Mutation.Type

sig_mat <-cbind(sig_mat, cosmicSig[match(rownames(sig_mat), rownames(cosmicSig)),])

corr.sig <- sapply(colnames(sig_mat)[1:Signature_num], function(x){
  return(sapply(paste0("Signature.", 1:30), function(y){
    return(cosine(sig_mat[, x], sig_mat[, y]))
  }))
})

m <- melt(corr.sig)
m$X1 <- gsub("Signature\\.", "", m$X1)
m$X1 <- factor(m$X1, levels = 1:30)
m$X2 <- gsub("S", "", m$X2)
# hm.palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')), space='Lab')
actual <- c('#2F2C62', '#42399B', '#4A52A7', '#59AFEA', '#7BCEB8', '#A7DA64',
            '#EFF121', '#F5952D', '#E93131', '#D70131', '#D70131', '#000000')
hm.palette <- colorRampPalette(actual, space='Lab')

mm <- m$value
mm[which(mm < 0.65)] <- NA

for(i in unique(m$X2)){
  if(sum(mm[which(m$X2 == i)], na.rm = T) == 0){
    MAX = which.max(m$value[which(m$X2 == i)])
    mm[which(m$X2 == i)][MAX] = m$value[which(m$X2 == i)][MAX]
  }
}

mm <-as.numeric(format(mm,digits = 2))
# m$X1 = factor(paste(m$X1," "), levels = paste(m$X1," "))
p_heat <- ggplot(m, aes(x = X2, y = X1, fill = value)) + geom_tile(color="white", size=0.1) +
  coord_fixed(ratio = 3/6) + 
  # geom_text(aes(label=mm), size = fontsize * 5/14) +
  geom_text(aes(label=mm), size = fontsize * 8/14) +
  scale_fill_gradientn(colours = hm.palette(10), limit = c(0,1), space = "Lab", name=NULL) +
  ylab("COSMIC Mutational Signatures") + xlab("") +
  scale_x_discrete(breaks = levels(factor(m$X2)),
                   labels = paste0("S", levels(factor(m$X2)))) +
  scale_y_discrete(position = "right") +
  theme(axis.text.x = element_text(hjust = 0.5, size = fontsize),
        axis.text.y = element_text(hjust = 0.5, size = fontsize),
        axis.title = element_text(size = fontsize),
        axis.line = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title=element_blank(),
        legend.position="top",
        legend.text=element_text(size = fontsize),
        legend.key.height=unit(25,"pt"),
        legend.key.width=unit(50,"pt"),
        plot.margin=unit(c(0,0,0,0),"lines"))

pdf(out_pdf, width = 25, height = 30) 
p_sigs
p_heat
dev.off()

# Sigcolor = c("#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02","#666666")[1:Param@signatureNum]
# Signame = c(paste0("S",1:(Param@signatureNum-1)), "BG")

# getMutSpetra <- function(maf, numbase=1){
#   m_raw = as.data.frame.matrix(table(maf[, c("Tumor_Sample_Barcode", "Variant_Type")]))
#   m <- maf[which(maf$Variant_Type=="SNP"), c("Tumor_Sample_Barcode","ref_context", "Tumor_Seq_Allele2")]
#   m$ref_context <- substr(m$ref_context, 10, 12)
#   m$ref_context <- toupper(paste(m$ref_context, m$Tumor_Seq_Allele2, sep = ""))
#   changeBase <- function(x){
#     if(x == "A"){return("T")}else if(x == "C"){return("G")}else if(x == "G"){return("C")}else{return("A")}
#   }
#   m$ref_context <- sapply(m$ref_context, function(x){
#     name = unlist(strsplit(x, ""))
#     if(name[2] %in% c("C", "T")){
#       return(paste(paste(name[2], name[4], sep = ""), paste(name[1], name[3],sep = ".")))
#     }else{
#       name = sapply(name, changeBase)
#       return(paste(paste(name[2], name[4], sep = ""), paste(name[3], name[1],sep = ".")))
#     }
#   })
#   m$alteration = sub("([ACGTN])([ACGTN]) .+", "\\1>\\2", m$ref_context)
#   m$context = sub("[ACGTN][ACGTN] (.+)", "\\1", m$ref_context)
#   m_sub <- as.data.frame.matrix(table(m[, c("Tumor_Sample_Barcode", "alteration")]))
#   return(cbind(sampleID = rownames(m_raw), m_raw, m_sub[match(rownames(m_raw), rownames(m_sub)),]))
# }

# mem <- getMutSpetra(maf,1)

# mem[, c("C>A_fq", "C>G_fq", "C>T_fq", 
#   "T>A_fq", "T>C_fq", "T>G_fq")] <- mem[, c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")]/mem$SNP

# mem_sig = getMembershipValue(Param)
# colnames(mem_sig) = Signame

# mem = cbind(mem, mem_sig[match(rownames(mem), rownames(mem_sig)), ])

# maf$Variant1 <- sapply(maf$Variant_Classification, function(x){
#   if(x == "Silent"){
#     return("Syn.")
#   }else if(x %in% c("3'UTR", "5'Flank", "5'UTR", "RNA", "Intron", "IGR", "lincRNA")){
#     return(NA)
#   }else{
#     return("Non_syn.")
#   }
# })

# Mat0 = as.data.frame.matrix(table(maf[, c("Tumor_Sample_Barcode", "Variant1")]))
# Mat0$sample = rownames(Mat0)
# Mat0$sum = Mat0$Non_syn. + Mat0$Syn.
# Mat0 = arrange(Mat0, desc(sum), desc(`Non_syn.`), desc(`Syn.`))
# Mat0$sample = factor(Mat0$sample, levels = Mat0$sample)
# mem$sampleID = factor(mem$sampleID, levels = Mat0$sample)

# Mat0 = melt(Mat0, id.vars = c("sample","sum"))
# Mat0$variable <- factor(Mat0$variable, levels = c('Syn.', 'Non_syn.'))

# p_snv <- ggplot(Mat0) + geom_bar(aes(x = sample, y = value, fill = variable), 
#   stat = "identity", colour="black", size = 0.1) +
#   ylab("Mutation Count") + 
#   scale_fill_manual(values = brewer.pal(12,"Paired")[c(3,2)]) +
#   theme(strip.background = element_blank(),
#         strip.text = element_blank(),
#         axis.line = element_blank(),
#         # axis.text.x = element_text(size=10, vjust = 0.5, angle = 90),
#         axis.text.x = element_blank(),
#         axis.text.y = element_text(size=fontsize),
#         axis.ticks.x = element_blank(),
#         axis.title.x = element_blank(),
#         axis.title.y = element_text(size=fontsize),
#         legend.position = c(0.5, 0.5),
#         legend.direction = "horizontal",
#         legend.title = element_blank(),
#         legend.text = element_text(size=fontsize),
#         legend.key.height=unit(10,'mm'),
#         legend.key.width=unit(4,'mm'),
#         plot.margin = unit(c(0, 0, 0, 0), "lines"))


# plotSamSig <- function(mem, Sigcolor = NULL, charSize = 15, rm.BG = FALSE, legend = TRUE,
#                        legendposition = "bottom") {
#   mm <- mem[, c("sampleID", Signame)]
#   if(rm.BG == TRUE){
#     mm[, Signame] <- mm[, Signame]/(1-mm$BG)
#     mm$BG <- NULL
#   }
#   mm <- mm %>% gather(signature, intensity, -sampleID)
#   mm$signature = factor(mm$signature, levels = Signame)
#   p2 = ggplot(mm)
#   p2 = p2 + geom_bar(aes(sampleID, intensity, fill = signature),
#                      stat = "identity")
#   p2 = p2 + theme_bw()
#   p2 = p2 + theme( axis.text.x = element_blank(),
#                    # axis.text.x = element_text(size=10, vjust = 0,5, angle = 90),
#                    axis.text.y = element_text(hjust = 0.5, size = charSize),
#                    axis.title.y=element_text(size=charSize*1.2),
#                    axis.title.x=element_blank(),
#                    axis.ticks.x=element_blank(),
#                    panel.background = element_blank(),
#                    panel.grid = element_blank(),
#                    panel.border = element_blank(),
#                    strip.text = element_text(size=charSize))
#   if(legend){
#     p2 = p2 + theme(legend.position = legendposition,
#                     legend.text = element_text(size=charSize),
#                     legend.title = element_blank())
#   }else{
#     p2 = p2 + theme(legend.position = "none")
#   }
#   if(!is.null(Sigcolor)){
#     p2 = p2 + scale_fill_manual(values = Sigcolor)
#   }
#   p2 = p2 + scale_y_continuous(expand=c(0, 0))
#   p2 = p2 + xlab("") + ylab("Activity")
#   return(p2)
# }


# p_sample = plotSamSig(mem, charSize = fontsize, legend = T, rm.BG = F, 
#     Sigcolor = Sigcolor, legendposition = "right") + 
# theme(legend.key.height=unit(10,'mm'), legend.key.width=unit(4,'mm'),
#     # panel.spacing.x=unit(-0.1, "lines"),
#     plot.margin = unit(c(0, 0, 0.25, 0), "lines"))

# pdf(gsub("RData", "vis.pdf", sig_RData), width = 25, height = 12) 
# ggarrange(p_snv, p_sample, ncol = 1, heights = c(2,2))
# dev.off()

# save.image(gsub("RData", "vis.RData", sig_RData))