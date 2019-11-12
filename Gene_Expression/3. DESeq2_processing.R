#########################GENE EXPRESSION PROCESSING########################################################################

#ssh -X sec

#cd /nfs/secure/scratch6/nekrut_gtex/rxv923/Transcript_Assembly/GreatApe_Summary
#source activate renv



library(tximport)
##Refseq id to gene annotation
tar2gen<-read.table("/nfs/brubeck.bx.psu.edu/scratch6/rahul/Apes/data/Amplicon/ReadCountSummary/hg38_refFlat_10252016_GenetoID_appendCDS.txt", header=F, sep="\t", stringsAsFactors=F)
names(tar2gen)=c("target_id","gene")


#Salmon processed output
base_dir <- "/nfs/brubeck.bx.psu.edu/scratch6/rahul/Apes/analysis/Bonobo/salmon_expression/quants/"
samples=c("160817_M02286","160906_M02286","161104_M02286","BonHL")
sal_dirsBo <- sapply(samples, function(id) file.path(base_dir, id,"quant.sf"))

base_dir <- "/nfs/brubeck.bx.psu.edu/scratch6/rahul/Apes/analysis/Chimpanzee/salmon_expression/quants/"
samples=c("SRR2040590","SRR2040591","ChimpHL")
sal_dirsCh <- sapply(samples, function(id) file.path(base_dir, id,"quant.sf"))

base_dir <- "/nfs/secure/project/nekrut_gtex/rxv923/Transcript_Assembly/Expression_African/salmon_expression/quants/"
samples=c("GTEX-QLQW","GTEX-P4QS","GTEX-UPJH")
sal_dirsHS <- sapply(samples, function(id) file.path(base_dir, id,"quant.sf"))

base_dir <- "/nfs/brubeck.bx.psu.edu/scratch6/rahul/Apes/analysis/Gorilla/salmon_expression/quants/"
samples=c("6-6-2013","1-13-2014","GorHL")
sal_dirsGr <- sapply(samples, function(id) file.path(base_dir, id,"quant.sf"))

base_dir <- "/nfs/brubeck.bx.psu.edu/scratch6/rahul/Apes/analysis/BorneanOrangutan/salmon_expression/quants/"
samples=c("3405R1","3405R2","3405R3","BOrHLR1","BOrHLR2")
sal_dirsBOr <- sapply(samples, function(id) file.path(base_dir, id,"quant.sf"))

base_dir <- "/nfs/brubeck.bx.psu.edu/scratch6/rahul/Apes/analysis/SumatranOrangutan/salmon_expression/quants/"
samples=c("1991R1","1991R2","1991R3")
sal_dirSOr <- sapply(samples, function(id) file.path(base_dir, id,"quant.sf"))

###Generate a summary matrix
list_dir <- c(sal_dirsBo,sal_dirsCh,sal_dirsHS,sal_dirsGr,sal_dirsBOr,sal_dirSOr)
gene_level=tximport(list_dir, type = "salmon",countsFromAbundance = "no", tx2gene = tar2gen, geneIdCol="gene",txIdCol="target_id" )

ApesData=gene_level$counts
ApesData<-as.matrix(ApesData)
ApesData=apply(ApesData, 1:2, round)
head(ApesData)


#From the summary table pick one replicate per sample and omit all the other samples.
ad<-ApesData[,c(1,4,5,6,7,8,9,10,11,13,14,17,19)]

##Data description for normalization.
species <-factor(c(rep("Bonobo",2),rep("Chimpanzee",3),rep("Human",3),rep("Gorilla",2),rep("Borang",2),rep("Sorang",1)))
type <- factor(c("paired_end","single_end","paired_end","paired_end","single_end","paired_end","paired_end","paired_end","paired_end","single_end","paired_end","paired_end","paired_end"))
sample<-c("B1","B2","C1","C2","C3","H1","H2","H3","G1","G2","BO1","BO2","S01")
run<- paste0("run",1:12)

cd=data.frame(cbind(species,type,sample),stringsAsFactors = FALSE)
row.names(cd)=colnames(ad)
head(cd)

#Normalization
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = ad, colData = cd, design = ~ type+species)
dds <- DESeq(dds)

#Final read counts used in downstream analysis
d<-counts(dds, normalized=TRUE)


#PCA specific normalization
rld <- rlog(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)

r<-assay(rld)
v<-assay(vsd)
head(r)
head(v)

#save(dds,d,v,r, file = "DESeq_greatApe_normalized_Salmon.RData")
#write.table(d,file="greatApe_Counts_Testis_2-3samples.tab", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
#write.table(r,file="greatApe_RLDCounts_Testis_2-3samples.tab", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
#write.table(v,file="greatApe_VSTCounts_Testis_2-3samples.tab", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
#save(dds,rld,vsd, file = "DESeq2_greatApe_normalized_dds_rld_vsd_Salmon.RData")



##List of ampliconic genes, the great ape assembled ampliconic transcripts are labelled as BPY2_CDS, CDY_CDS,.... .
##If for a species, we did not finde a transcript in our assembly, we used a human transcript for the same gene family as a dummy .
ampliconList=c("BPY2B","BPY2","BPY2C","BPY2_CDS","CDY1","CDY1B","CDY2A","CDY2B","CDY_CDS","DAZ1","DAZ2","DAZ3","DAZ4","DAZ_CDS","HSFY1","HSFY2","HSFY_CDS","PRY","PRY2","PRY_CDS","RBMY1A1","RBMY1J","RBMY1F","RBMY1E","RBMY1D","RBMY1B","RBMY_CDS","TSPY10","TSPY1","TSPY3","TSPY8","TSPY4","TSPY2","TSPY_CDS","VCY","VCY1B","VCY1_CDS","XKRY","XKRY2")


#Parsing ampliconic gene expression counts.
ampD=NULL
for (gene in ampliconList){
  print(gene)
  ampD=rbind(ampD,d[which(row.names(d)==gene),])
}
amp_normD=cbind(ampliconList,ampD)


Amplicon_summaryD=rbind(
	  BPY=apply(ampD[1:4,],2,sum),
	  CDY=apply(ampD[5:9,],2,sum),
	  DAZ=apply(ampD[10:14,],2,sum),
	  HSFY=apply(ampD[15:17,],2,sum),
	  PRY=apply(ampD[18:20,],2,sum),
	  RBMY=apply(ampD[21:27,],2,sum),
	  TSPY=apply(ampD[28:34,],2,sum),
	  VCY=apply(ampD[35:37,],2,sum),
	  XKRY=apply(ampD[38:39,],2,sum)
)


#Generate output which can be used for EVEmodel analysis and summary of expression figure and table.
write.table(Amplicon_summaryD,file="Apes_normalizedDDS_DESeq2_Amp_Salmon.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
fh <- file("Apes_WithoutNA_normalized_DESeq2_Amp_Salmon.dat", open="wt")
writeLines("5", fh)
write.table(Amplicon_summaryD[c(1,2,3,6,7),], fh ,sep=" ", quote=FALSE, row.names=TRUE,col.names=FALSE)
close(fh)