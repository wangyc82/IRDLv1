IRDLv1 is procedure for identify the regulaoty elements fro divergent lncRNAs.
It includes the MATLAB code for generation of the 3mer sequence feature from lncRNA/gene nucleic acid sequences (codonComposition.m),
the R code for random forest crossvall validation (crossvalRF.R),
the R code for generation the nucleic acid sequences for lncRNA/gene and their promoters (generate-lncRNA/gene-sequence).
To get the divergent lncRNA that regulate a given gene, using RF regression model on lncRNA/gene features, 
and rank those feature in descending order based on feature importance generated from RF regression model: R code is given by generate-lncRNAlist.

The corresponding data used is required by contacting ycwang@nwipb.cas.cn
