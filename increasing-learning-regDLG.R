# this procedure is for generating regulatory network for divergent lncRNA
# using increasing learning framework

load("~/Documents/divergent-lncRNA-gene/Increasing-learning-mouse/initial-data.RData")
# initial positive is those eight experimental validated regulatory interactiond
# using correlation under each genomic feature as the input feature
# each input is a pair of lncrna-gene
# the negative is generated from cloasely positioned random lncRNA-gene pair
source('~/Documents/CCLE/cell line sensitivity via kernel fusion/crossvalSVM.R')
source('~/Documents/CCLE/cell line sensitivity via kernel fusion/getperf.R')

star_sample_lab<-match(c("Foxd3","Evx1","Ccnyl1","Sox3","Nr2f1","Rab11b","Zfp687","Gata3"),rownames(gene_A))
C<-which(Adj==1,arr.ind = T)
C1<-C[-star_sample_lab,]# all candidate regulatory lncRNA-gene associations

dist_loc<-matrix(0,nrow(gene_A),nrow(lncrna_A))
for (i in 1:nrow(XH_nearbyGene_info_open_exp_seq)) {for (j in 1:nrow(XH_lncRNA_info_open_exp_seq)) {a<-XH_nearbyGene_info_open_exp_seq$chromosome[i];b<-XH_lncRNA_info_open_exp_seq$chromosome[j]; if (a==b) {dist_loc[i,j]<-abs(as.numeric(XH_nearbyGene_info_open_exp_seq$start[i])-as.numeric(XH_lncRNA_info_open_exp_seq$start[j]))} else {dist_loc[i,j]<-Inf}}}
rownames(dist_loc)<-rownames(gene_A)
colnames(dist_loc)<-rownames(lncrna_A)

XH_dist<-diag(dist_loc)
Xst<-cbind(cor_seq,cor_exp,cor_open,XH_dist/1e+6)[-star_sample_lab,]
colnames(Xst)<-c("cor_seq","cor_exp","cor_open","dist")
Xst[is.na(Xst)]<-0


gene_seq<-gene_A[,1:84]
lncrna_seq<-lncrna_A[,1:84]

gene_exp<-gene_A[,(84+1):(84+146)]
lncrna_exp<-lncrna_A[,(84+1):(84+146)]
 
gene_open<-gene_A[,(84+146+1):(84+146+55)]
lncrna_open<-lncrna_A[,(84+146+1):(84+146+55)]

# generating the candicates for negatives
A<-dist_loc
diag(A)<-0
M<-which(A<=1e+7 & A!=0,arr.ind = T)

# the first round
Xp<-data.matrix(star_sample1_info[,c(1,2,3,4)])
colnames(Xp)<-c("cor_seq","cor_exp","cor_open","dist")
Xp[,4]<-abs(Xp[,4])/1e+6

library(pracma)
p1<-randperm(nrow(M),5*nrow(Xp))

cs<-unlist(lapply(1:length(p1),function(x) cor(gene_seq[M[p1[x],1],],lncrna_seq[M[p1[x],2],])))
ce<-unlist(lapply(1:length(p1),function(x) cor(gene_exp[M[p1[x],1],],lncrna_exp[M[p1[x],2],])))
co<-unlist(lapply(1:length(p1),function(x) cor(gene_open[M[p1[x],1],],lncrna_open[M[p1[x],2],],method = "spearman")))
dd<-unlist(lapply(1:length(p1),function(x) dist_loc[M[p1[x],1],M[p1[x],2]]))

Xn<-cbind(cs,ce,co,dd/1e+6)
colnames(Xn)<-c("cor_seq","cor_exp","cor_open","dist")
X<-rbind(Xp,Xn)
X[is.na(X)]<-0
Ylab<-c(rep("a",nrow(Xp)),rep("b",nrow(Xn)))
library(e1071)
crossV<-crossvalSVM(X,Ylab,5,10,0.1)
perf<-getperf(crossV[[1]],crossV[[2]])

model<-svm(X,Ylab,type="C-classification",kernel="radial",cost=10,gamma=0.1,class.weights=c(a=1, b=1),probability=TRUE)
predictedY<-attr(predict(model, Xst,probability=TRUE),"probabilities")[,1]

# the second round
Xp1<-rbind(Xp,Xst[which(predictedY>0.96),]) # adding 55 postives
Xst1<-Xst[-which(predictedY>0.96),]

p2<-randperm(nrow(M),5*nrow(Xp1))
cs1<-unlist(lapply(1:length(p2),function(x) cor(gene_seq[M[p2[x],1],],lncrna_seq[M[p2[x],2],])))
ce1<-unlist(lapply(1:length(p2),function(x) cor(gene_exp[M[p2[x],1],],lncrna_exp[M[p2[x],2],])))
co1<-unlist(lapply(1:length(p2),function(x) cor(gene_open[M[p2[x],1],],lncrna_open[M[p2[x],2],],method = "spearman")))
dd1<-unlist(lapply(1:length(p2),function(x) dist_loc[M[p2[x],1],M[p2[x],2]]))

Xn1<-cbind(cs1,ce1,co1,dd1/1e+6)
colnames(Xn1)<-c("cor_seq","cor_exp","cor_open","dist")
X<-rbind(Xp1,Xn1)
X[is.na(X)]<-0
Ylab<-c(rep("a",nrow(Xp1)),rep("b",nrow(Xn1)))
crossV<-crossvalSVM(X,Ylab,5,10,0.01)
perf<-getperf(crossV[[1]],crossV[[2]])

model1<-svm(X,Ylab,type="C-classification",kernel="radial",cost=10,gamma=0.01,class.weights=c(a=1, b=1),probability=TRUE)
predictedY1<-attr(predict(model1, Xst1,probability=TRUE),"probabilities")[,1]

# the thrid round
Xp2<-rbind(Xp1,Xst1[which(predictedY1>0.9),]) # adding 71 postives
Xst2<-Xst1[-which(predictedY1>0.9),]

p3<-randperm(nrow(M),5*nrow(Xp2))
cs2<-unlist(lapply(1:length(p3),function(x) cor(gene_seq[M[p3[x],1],],lncrna_seq[M[p3[x],2],])))
ce2<-unlist(lapply(1:length(p3),function(x) cor(gene_exp[M[p3[x],1],],lncrna_exp[M[p3[x],2],])))
co2<-unlist(lapply(1:length(p3),function(x) cor(gene_open[M[p3[x],1],],lncrna_open[M[p3[x],2],],method = "spearman")))
dd2<-unlist(lapply(1:length(p3),function(x) dist_loc[M[p3[x],1],M[p3[x],2]]))

Xn2<-cbind(cs2,ce2,co2,dd2/1e+6)
colnames(Xn2)<-c("cor_seq","cor_exp","cor_open","dist")
X<-rbind(Xp2,Xn2)
X[is.na(X)]<-0
Ylab<-c(rep("a",nrow(Xp2)),rep("b",nrow(Xn2)))
crossV<-crossvalSVM(X,Ylab,5,10,0.1)
perf<-getperf(crossV[[1]],crossV[[2]])

model2<-svm(X,Ylab,type="C-classification",kernel="radial",cost=10,gamma=0.1,class.weights=c(a=1, b=1),probability=TRUE)
predictedY2<-attr(predict(model2, Xst2,probability=TRUE),"probabilities")[,1]

# the forth round
Xp3<-rbind(Xp2,Xst2[which(predictedY2>0.8),]) # adding 81 postives
Xst3<-Xst2[-which(predictedY2>0.8),]

p4<-randperm(nrow(M),5*nrow(Xp3))
cs3<-unlist(lapply(1:length(p4),function(x) cor(gene_seq[M[p4[x],1],],lncrna_seq[M[p4[x],2],])))
ce3<-unlist(lapply(1:length(p4),function(x) cor(gene_exp[M[p4[x],1],],lncrna_exp[M[p4[x],2],])))
co3<-unlist(lapply(1:length(p4),function(x) cor(gene_open[M[p4[x],1],],lncrna_open[M[p4[x],2],],method = "spearman")))
dd3<-unlist(lapply(1:length(p4),function(x) dist_loc[M[p4[x],1],M[p4[x],2]]))

Xn3<-cbind(cs3,ce3,co3,dd3/1e+6)
colnames(Xn3)<-c("cor_seq","cor_exp","cor_open","dist")
X<-rbind(Xp3,Xn3)
X[is.na(X)]<-0
Ylab<-c(rep("a",nrow(Xp3)),rep("b",nrow(Xn3)))
crossV<-crossvalSVM(X,Ylab,5,10,0.1)
perf<-getperf(crossV[[1]],crossV[[2]])

model3<-svm(X,Ylab,type="C-classification",kernel="radial",cost=10,gamma=0.1,class.weights=c(a=1, b=1),probability=TRUE)
predictedY3<-attr(predict(model3, Xst3,probability=TRUE),"probabilities")[,1]

# the fifth round
Xp4<-rbind(Xp3,Xst3[which(predictedY3>0.7),]) # adding 83 postives
Xst4<-Xst3[-which(predictedY3>0.7),]

p5<-randperm(nrow(M),5*nrow(Xp4))
cs4<-unlist(lapply(1:length(p5),function(x) cor(gene_seq[M[p5[x],1],],lncrna_seq[M[p5[x],2],])))
ce4<-unlist(lapply(1:length(p5),function(x) cor(gene_exp[M[p5[x],1],],lncrna_exp[M[p5[x],2],])))
co4<-unlist(lapply(1:length(p5),function(x) cor(gene_open[M[p5[x],1],],lncrna_open[M[p5[x],2],],method = "spearman")))
dd4<-unlist(lapply(1:length(p5),function(x) dist_loc[M[p5[x],1],M[p5[x],2]]))

Xn4<-cbind(cs4,ce4,co4,dd4/1e+6)
colnames(Xn4)<-c("cor_seq","cor_exp","cor_open","dist")
X<-rbind(Xp4,Xn4)
X[is.na(X)]<-0
Ylab<-c(rep("a",nrow(Xp4)),rep("b",nrow(Xn4)))
crossV<-crossvalSVM(X,Ylab,5,10,0.1)
perf<-getperf(crossV[[1]],crossV[[2]])

model4<-svm(X,Ylab,type="C-classification",kernel="radial",cost=10,gamma=0.1,class.weights=c(a=1, b=1),probability=TRUE)
predictedY4<-attr(predict(model4, Xst4,probability=TRUE),"probabilities")[,1]

# the sixth round
Xp5<-rbind(Xp4,Xst4[which(predictedY4>0.7),]) # adding 72 postives
Xst5<-Xst4[-which(predictedY4>0.7),]

p6<-randperm(nrow(M),5*nrow(Xp5))
cs5<-unlist(lapply(1:length(p6),function(x) cor(gene_seq[M[p6[x],1],],lncrna_seq[M[p6[x],2],])))
ce5<-unlist(lapply(1:length(p6),function(x) cor(gene_exp[M[p6[x],1],],lncrna_exp[M[p6[x],2],])))
co5<-unlist(lapply(1:length(p6),function(x) cor(gene_open[M[p6[x],1],],lncrna_open[M[p6[x],2],],method = "spearman")))
dd5<-unlist(lapply(1:length(p6),function(x) dist_loc[M[p6[x],1],M[p6[x],2]]))

Xn5<-cbind(cs5,ce5,co5,dd5/1e+6)
colnames(Xn5)<-c("cor_seq","cor_exp","cor_open","dist")
X<-rbind(Xp5,Xn5)
X[is.na(X)]<-0
Ylab<-c(rep("a",nrow(Xp5)),rep("b",nrow(Xn5)))
crossV<-crossvalSVM(X,Ylab,5,10,0.1)
perf<-getperf(crossV[[1]],crossV[[2]])

model5<-svm(X,Ylab,type="C-classification",kernel="radial",cost=10,gamma=0.1,class.weights=c(a=1, b=1),probability=TRUE)
predictedY5<-attr(predict(model5, Xst5,probability=TRUE),"probabilities")[,1]

# the seventh round
Xp6<-rbind(Xp5,Xst5[which(predictedY5>0.7),]) # adding 26 postives
Xst6<-Xst5[-which(predictedY5>0.7),]

p7<-randperm(nrow(M),5*nrow(Xp6))
cs6<-unlist(lapply(1:length(p7),function(x) cor(gene_seq[M[p7[x],1],],lncrna_seq[M[p7[x],2],])))
ce6<-unlist(lapply(1:length(p7),function(x) cor(gene_exp[M[p7[x],1],],lncrna_exp[M[p7[x],2],])))
co6<-unlist(lapply(1:length(p7),function(x) cor(gene_open[M[p7[x],1],],lncrna_open[M[p7[x],2],],method = "spearman")))
dd6<-unlist(lapply(1:length(p7),function(x) dist_loc[M[p7[x],1],M[p7[x],2]]))

Xn6<-cbind(cs6,ce6,co6,dd6/1e+6)
colnames(Xn6)<-c("cor_seq","cor_exp","cor_open","dist")
X<-rbind(Xp6,Xn6)
X[is.na(X)]<-0
Ylab<-c(rep("a",nrow(Xp6)),rep("b",nrow(Xn6)))
crossV<-crossvalSVM(X,Ylab,5,10,0.1)
perf<-getperf(crossV[[1]],crossV[[2]])

model6<-svm(X,Ylab,type="C-classification",kernel="radial",cost=10,gamma=0.1,class.weights=c(a=1, b=1),probability=TRUE)
predictedY6<-attr(predict(model6, Xst6,probability=TRUE),"probabilities")[,1]

# the eighth round
Xp7<-rbind(Xp6,Xst6[which(predictedY6>0.7),]) # adding 25 postives
Xst7<-Xst6[-which(predictedY6>0.7),]

p8<-randperm(nrow(M),5*nrow(Xp7))
cs7<-unlist(lapply(1:length(p8),function(x) cor(gene_seq[M[p8[x],1],],lncrna_seq[M[p8[x],2],])))
ce7<-unlist(lapply(1:length(p8),function(x) cor(gene_exp[M[p8[x],1],],lncrna_exp[M[p8[x],2],])))
co7<-unlist(lapply(1:length(p8),function(x) cor(gene_open[M[p8[x],1],],lncrna_open[M[p8[x],2],],method = "spearman")))
dd7<-unlist(lapply(1:length(p8),function(x) dist_loc[M[p8[x],1],M[p8[x],2]]))

Xn7<-cbind(cs7,ce7,co7,dd7/1e+6)
colnames(Xn7)<-c("cor_seq","cor_exp","cor_open","dist")
X<-rbind(Xp7,Xn7)
X[is.na(X)]<-0
Ylab<-c(rep("a",nrow(Xp7)),rep("b",nrow(Xn7)))
crossV<-crossvalSVM(X,Ylab,5,10,0.1)
perf<-getperf(crossV[[1]],crossV[[2]])

model7<-svm(X,Ylab,type="C-classification",kernel="radial",cost=10,gamma=0.1,class.weights=c(a=1, b=1),probability=TRUE)
predictedY7<-attr(predict(model7, Xst7,probability=TRUE),"probabilities")[,1]

# the ninth round
Xp8<-rbind(Xp7,Xst7[which(predictedY7>0.7),]) # adding 12 postives
Xst8<-Xst7[-which(predictedY7>0.7),]

p9<-randperm(nrow(M),5*nrow(Xp8))
cs8<-unlist(lapply(1:length(p9),function(x) cor(gene_seq[M[p9[x],1],],lncrna_seq[M[p9[x],2],])))
ce8<-unlist(lapply(1:length(p9),function(x) cor(gene_exp[M[p9[x],1],],lncrna_exp[M[p9[x],2],])))
co8<-unlist(lapply(1:length(p9),function(x) cor(gene_open[M[p9[x],1],],lncrna_open[M[p9[x],2],],method = "spearman")))
dd8<-unlist(lapply(1:length(p9),function(x) dist_loc[M[p9[x],1],M[p9[x],2]]))

Xn8<-cbind(cs8,ce8,co8,dd8/1e+6)
colnames(Xn8)<-c("cor_seq","cor_exp","cor_open","dist")
X<-rbind(Xp8,Xn8)
X[is.na(X)]<-0
Ylab<-c(rep("a",nrow(Xp8)),rep("b",nrow(Xn8)))
crossV<-crossvalSVM(X,Ylab,5,10,0.1)
perf<-getperf(crossV[[1]],crossV[[2]])

model8<-svm(X,Ylab,type="C-classification",kernel="radial",cost=10,gamma=0.1,class.weights=c(a=1, b=1),probability=TRUE)
predictedY8<-attr(predict(model8, Xst8,probability=TRUE),"probabilities")[,1]

# the last round
Xp9<-rbind(Xp8,Xst8[which(predictedY8>0.7),]) # adding 12 postives
Xst9<-Xst8[-which(predictedY8>0.7),]

p10<-randperm(nrow(M),5*nrow(Xp9))
cs9<-unlist(lapply(1:length(p10),function(x) cor(gene_seq[M[p10[x],1],],lncrna_seq[M[p10[x],2],])))
ce9<-unlist(lapply(1:length(p10),function(x) cor(gene_exp[M[p10[x],1],],lncrna_exp[M[p10[x],2],])))
co9<-unlist(lapply(1:length(p10),function(x) cor(gene_open[M[p10[x],1],],lncrna_open[M[p10[x],2],],method = "spearman")))
dd9<-unlist(lapply(1:length(p10),function(x) dist_loc[M[p10[x],1],M[p10[x],2]]))

Xn9<-cbind(cs9,ce9,co9,dd9/1e+6)
colnames(Xn9)<-c("cor_seq","cor_exp","cor_open","dist")
X<-rbind(Xp9,Xn9)
X[is.na(X)]<-0
Ylab<-c(rep("a",nrow(Xp9)),rep("b",nrow(Xn9)))
crossV<-crossvalSVM(X,Ylab,5,10,0.1)
perf<-getperf(crossV[[1]],crossV[[2]])

model9<-svm(X,Ylab,type="C-classification",kernel="radial",cost=10,gamma=0.1,class.weights=c(a=1, b=1),probability=TRUE)
predictedY9<-attr(predict(model9, Xst9,probability=TRUE),"probabilities")[,1]

# the final round
Xp_f<-rbind(Xp9,Xst9[which(predictedY9>0.7),]) # adding 1 postive
Xst_f<-Xst9[-which(predictedY9>0.7),]

pf<-randperm(nrow(M),5*nrow(Xp_f))
cs_f<-unlist(lapply(1:length(pf),function(x) cor(gene_seq[M[pf[x],1],],lncrna_seq[M[pf[x],2],])))
ce_f<-unlist(lapply(1:length(pf),function(x) cor(gene_exp[M[pf[x],1],],lncrna_exp[M[pf[x],2],])))
co_f<-unlist(lapply(1:length(pf),function(x) cor(gene_open[M[pf[x],1],],lncrna_open[M[pf[x],2],],method = "spearman")))
dd_f<-unlist(lapply(1:length(pf),function(x) dist_loc[M[pf[x],1],M[pf[x],2]]))

Xn_f<-cbind(cs_f,ce_f,co_f,dd_f/1e+6)
colnames(Xn_f)<-c("cor_seq","cor_exp","cor_open","dist")
X<-rbind(Xp_f,Xn_f)
X[is.na(X)]<-0
Ylab<-c(rep("a",nrow(Xp_f)),rep("b",nrow(Xn_f)))
crossV<-crossvalSVM(X,Ylab,5,10,0.01)
perf<-getperf(crossV[[1]],crossV[[2]])

model_f<-svm(X,Ylab,type="C-classification",kernel="radial",cost=10,gamma=0.01,class.weights=c(a=1, b=1),probability=TRUE)
predictedY_f<-attr(predict(model_f, Xst_f,probability=TRUE),"probabilities")[,1]

Xp_f1<-rbind(Xp_f,Xst_f[which(predictedY_f>0.7),])# adding another two with low co-exp and high co-open


#enhancer region checking
mg<-unlist(lapply(1:length(strsplit(rownames(Xp_f)," ")), function(x) strsplit(rownames(Xp_f)," ")[[x]][1]))
ml<-unlist(lapply(1:length(strsplit(rownames(Xp_f)," ")), function(x) strsplit(rownames(Xp_f)," ")[[x]][2]))
Xpf_gene_info<-subset(XH_nearbyGene_info_open_exp_seq,gene_symbol %in% mg)
Xpf_lncrna_info<-XH_lncRNA_info_open_exp_seq[match(ml,XH_lncRNA_info_open_exp_seq$gene_symbol),]
Xpf_lncRNA_en<-lapply(1:nrow(Xpf_lncrna_info),function(x) which(mouse_enhancer_info$chrom==Xpf_lncrna_info$chromosome[x] & as.numeric(mouse_enhancer_info$start)<as.numeric(Xpf_lncrna_info$start[x]) & as.numeric(mouse_enhancer_info$end)>as.numeric(Xpf_lncrna_info$start[x])))
Xpf_lncRNA_en1<-lapply(1:nrow(Xpf_lncrna_info),function(x) which(mouse_enhancer_info$chrom==Xpf_lncrna_info$chromosome[x] & as.numeric(mouse_enhancer_info$start)<as.numeric(Xpf_lncrna_info$end[x]) & as.numeric(mouse_enhancer_info$end)>as.numeric(Xpf_lncrna_info$end[x])))


#tissue specific generating
# for mouse
Xpf.gene<-unlist(lapply(1:length(strsplit(rownames(Xpf_r1)," ")),function(x) strsplit(rownames(Xpf_r1)," ")[[x]][1]))
Xpf.lncrna<-unlist(lapply(1:length(strsplit(rownames(Xpf_r1)," ")),function(x) strsplit(rownames(Xpf_r1)," ")[[x]][2]))

Xpf.tissue<-list()
for (i in 1:length(Xpf.gene)) {
  a<-colnames(gene_exp)[which(gene_exp[which(rownames(gene_exp) %in% Xpf.gene[i]),]>6)]
  t<-mouse_XH_LandG_exp_sample$tissue[match(a,mouse_XH_LandG_exp_sample$sampleID)]
  b<-colnames(lncrna_exp)[which(lncrna_exp[which(rownames(lncrna_exp) %in% Xpf.lncrna[i]),]>4)]
  t1<-mouse_XH_LandG_exp_sample$tissue[match(b,mouse_XH_LandG_exp_sample$sampleID)]
  l<-intersect(t,t1)
  s<-table(as.factor(c(t,t1)))[l]
  Xpf.tissue[[i]]<-names(s)[which(s==max(s))]
  rm(a,b,t,t1,l,s)
}

#for human
regX.gene<-unlist(lapply(1:length(strsplit(rownames(regX_r)," ")),function(x) strsplit(rownames(regX_r)," ")[[x]][1]))
regX.lncrna<-unlist(lapply(1:length(strsplit(rownames(regX_r)," ")),function(x) strsplit(rownames(regX_r)," ")[[x]][2]))


regX.tissue<-list()
for (i in 1:length(regX.gene)) {
  a<-colnames(gene_exp)[which(gene_exp[which(rownames(gene_exp) %in% regX.gene[i]),]>6)]
  t<-sample_TT[a]
  b<-colnames(lncrna_exp)[which(lncrna_exp[which(rownames(lncrna_exp) %in% regX.lncrna[i]),]>4)]
  t1<-sample_TT[b]
  l<-intersect(t,t1)
  l[is.na(l)]<-"NA"
  l<-l[-which(l %in% "NA")]
  s<-table(as.factor(c(t,t1)))[l]
  regX.tissue[[i]]<-names(s)[which(s==max(s))]
  rm(a,b,t,t1,l,s)
}






  
  