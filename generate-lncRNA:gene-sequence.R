#this procedure is to get the sequence for XH lncRNA and gene promoter 
#for human
load("~/Documents/divergent-lncRNA-gene/human/human-XH-OPENNESS-open-data.RData")
rm(list=setdiff(ls(),c("XH_lncrna_info_open1","XH_nearbyGene_info_open1")))
#load lncpromoter.txt and promoter.txt
XH_lncrna_promoter_info<-lncpromoter[match(XH_lncrna_info_open1$XH_lncrna_noncode_openness,lncpromoter$X6),c(2,3,4,5,6)]
colnames(XH_lncrna_promoter_info)<-lncpromoter[2,c(2,3,4,5,6)]
XH_nearbyGene_promoter_info<-promoter[match(XH_nearbyGene_info_open1$gene_symbol,promoter$X6),c(2,3,4,5,6)]
colnames(XH_nearbyGene_promoter_info)<-promoter[2,c(2,3,4,5,6)]

library(Biostrings)
s = readDNAStringSet("~/Documents/Homo_sapiens_UCSC_hg19/Sequence/WholeGenomeFasta/genome.fa",format="fasta")
chrom<-names(s)
XH_genepromoter_seq<-lapply(1:nrow(XH_nearbyGene_promoter_info),function(x) substring(as.character(s[XH_nearbyGene_promoter_info$chrom[x]]),as.numeric(XH_nearbyGene_promoter_info$start[x]),as.numeric(XH_nearbyGene_promoter_info$end[x])))
XH_lncpromoter_seq<-lapply(1:nrow(XH_lncrna_promoter_info),function(x) substring(as.character(s[XH_lncrna_promoter_info$chrom[x]]),as.numeric(XH_lncrna_promoter_info$start[x]),as.numeric(XH_lncrna_promoter_info$end[x])))
#write sequence in fasta format
library(seqinr)
write.fasta(sequences = XH_lncpromoter_seq, names = XH_lncrna_info_open1$gene_symbol, nbchar = 80, file.out = "~/Documents/divergent-lncRNA-gene/human/1MBrandom-lc/human-lncrnapromoter-seq.fasta")
write.fasta(sequences = XH_genepromoter_seq, names = XH_nearbyGene_promoter_info$gene, nbchar = 80, file.out = "~/Documents/divergent-lncRNA-gene/human/1MBrandom-lc/human-genepromoter-seq.fasta")

#for mouse
load("~/Documents/divergent-lncRNA-gene/mouse/XH-lncRNA-CG-OPENNESS-open-classification.RData")
rm(list=setdiff(ls(),c("XH_lncRNA_info_open1","XH_nearbyGene_info_open1")))
#load pormoter.txt and lncpromoter.txt
XH_lncrna_promoter_info<-lncpromoter[match(XH_lncRNA_info_open2$XH_lncrna_noncode.X2,lncpromoter$X6),c(2,3,4,5,6)]
colnames(XH_lncrna_promoter_info)<-lncpromoter[2,c(2,3,4,5,6)]

XH_nearbyGene_promoter_info<-promoter[match(XH_nearbyGene_info_open1$gene_symbol,promoter$X6),c(2,3,4,5,6)]
colnames(XH_nearbyGene_promoter_info)<-promoter[2,c(2,3,4,5,6)]

library("BSgenome.Mmusculus.UCSC.mm9") #load mm9 UCSC whole genome sequence
genome<-BSgenome.Mmusculus.UCSC.mm9
chrom<-names(genome)
XH_genepromoter_seq<-lapply(1:nrow(XH_nearbyGene_promoter_info),function(x) substring(as.character(genome[[XH_nearbyGene_promoter_info$chrom[x]]]),as.numeric(XH_nearbyGene_promoter_info$start[x]),as.numeric(XH_nearbyGene_promoter_info$end[x])))
XH_lncpromoter_seq<-lapply(1:nrow(XH_lncrna_promoter_info),function(x) substring(as.character(genome[[XH_lncrna_promoter_info$chrom[x]]]),as.numeric(XH_lncrna_promoter_info$start[x]),as.numeric(XH_lncrna_promoter_info$end[x])))

library(seqinr)
write.fasta(sequences = XH_genepromoter_seq, names = XH_nearbyGene_promoter_info$gene, nbchar = 80, file.out = "~/Documents/divergent-lncRNA-gene/mouse/1MBrandom-lc/mouse-genepromoter-seq.fasta")
write.fasta(sequences = XH_lncpromoter_seq, names = XH_lncRNA_info_open2$gene_symbol, nbchar = 80, file.out = "~/Documents/divergent-lncRNA-gene/mouse/1MBrandom-lc/mouse-lncpromoter-seq.fasta")

# for lncRNA/gene equences, just replace XH_lncRNA_promoter_info with XH_lncRNA_info for lncRNA, replace XH_nearbyGene_promoter_info with XH_nearbyGene_info for gene.