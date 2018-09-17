# this procedure is to get the lncRNA list for a given gene round by at least one divergent lncRNA
library("randomForest")
#take mouse data as an example
# load("~/Documents/divergent-lncRNA-gene/in Hongkong/mouse/XH-divergent-lncrna-gene-openness-open-regression.RData")

# generate the nearby lncRNA within 1MB for a given gene
XH_gene_nearbyL<-lapply(1:nrow(XH_gene_info_open),function(x) which(lncrna_chrom == XH_gene_chrom[x] & abs(lncrna_start-XH_gene_start[x])<=1e+6))

lncRNA_list<-list()
for (i in 1:nrow(XH_gene_open)) {
  Y<-XH_gene_open[i,]
  X<-t(lncrna_open[XH_gene_nearbyL[[i]],])
  model<-randomForest(X,Y,xtest=X,ytest=Y)
  importance<-model$importance
  sr<-sort(importance,decreasing = T,index.return=T)
  lncRNA_list[[i]]<-XH_gene_nearbyL[[i]][sr$ix]
}

# get location of divergent lncRNA in gene's nearby lncRNAs list
match_divergent<-match(XH_lncRNA_info_open$NONCODE.G,rownames(lncrna_open))

# get the location of divergent lncRNA in lncRNA rank list
divergene_randk_loc<-lapply(1:length(lncRNA_list),function(x) which(lncRNA_list[[x]]==match_divergent[[x]])/length(lncRNA_list[[x]]))