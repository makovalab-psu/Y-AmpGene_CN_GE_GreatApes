#cd /nfs/brubeck.bx.psu.edu/scratch6/rahul/Apes/analysis/ape_compare/
#source activate renv

library(ggplot2)
library(gridExtra)
library(reshape2)
library(ggsignif)
library(plyr)
library(ggpubr)
library(grid) 
library("pheatmap")
library("RColorBrewer")
library(ape)
library(dplyr)


greatA<-read.table(file="New_Normalized_all_updated",header=T,sep="\t",stringsAsFactors=F, fill=TRUE)
greatA_CN<-cbind(greatA[,c(1,2)],apply(greatA[,-c(1,2)],1,function(x){mean(x,na.rm=T)}))
colnames(greatA_CN)[3]="CopyNumber"

annotation<-read.table(file="New_great_ape_info",header=T,sep="\t",stringsAsFactors=F)

greatAPE_summary<-cbind.data.frame(greatA_CN,annotation[match(greatA_CN$IID,annotation$IID),2],stringsAsFactors=F)
colnames(greatAPE_summary)[4]="Species"
greatAPE_summary$Gene[which(greatAPE_summary$Gene=="BPY")]="BPY2"
greatAPE_summary$Species[which(greatAPE_summary$Species=="Chimp")]="Chimpanzee"
greatAPE_summary$Species[which(greatAPE_summary$Species=="Borangutan")]="B.orangutan"
greatAPE_summary$Species[which(greatAPE_summary$Species=="Sorangutan")]="S.orangutan"

GA_table<-dcast(greatAPE_summary,Gene~IID,value.var="CopyNumber")

species <- split(greatAPE_summary  , f = greatAPE_summary$Species)
geneFamily<-split(greatAPE_summary , f = greatAPE_summary$Gene )


#Summary of Copy Number
#Median, Variance and RANGE

##CAFE output using median values 
Species_level_median<-function(data){
	test<-dcast(data,Gene~IID,value.var="CopyNumber")
	avg_val<-apply(test[,-1],1,function(x){median(x,na.rm=T)})
	out<-as.data.frame(avg_val)
	rownames(out)<-test[,1]
	colnames(out)<-data$Species[1]
	return(out)
}

Species_median_copynumber=llply(species,Species_level_median)

Species_level_range<-function(data){
	test<-dcast(data,Gene~IID,value.var="CopyNumber")
	var_val<-apply(test[,-1],1,function(x){range(x,na.rm=T)})
	out<-as.data.frame(var_val)
	rownames(out)<-c("Min","Max")
	colnames(out)<-test[,1]
	return(out)
}

Species_range_copynumber=llply(species,Species_level_range)


Species_level_variance<-function(data){
	test<-dcast(data,Gene~IID,value.var="CopyNumber")
	var_val<-apply(test[,-1],1,function(x){var(x,na.rm=T)})
	out<-as.data.frame(var_val)
	rownames(out)<-test[,1]
	colnames(out)<-data$Species[1]
	return(out)
}

Species_variance_copynumber=llply(species,Species_level_variance)

tempBCv<-merge(Species_variance_copynumber$Bonobo,Species_variance_copynumber$Chimpanzee,by='row.names', all=TRUE)
rownames(tempBCv) <- tempBCv$Row.names; tempBCv$Row.names <- NULL

tempBCHv<-merge(tempBCv,Species_variance_copynumber$Human,by='row.names', all=TRUE)
rownames(tempBCHv) <- tempBCHv$Row.names; tempBCHv$Row.names <- NULL

tempBCHGv<-merge(tempBCHv,Species_variance_copynumber$Gorilla,by='row.names', all=TRUE)
rownames(tempBCHGv) <- tempBCHGv$Row.names; tempBCHGv$Row.names <- NULL

tempBCHGBov<-merge(tempBCHGv,Species_variance_copynumber$B.orangutan,by='row.names', all=TRUE)
rownames(tempBCHGBov) <- tempBCHGBov$Row.names; tempBCHGBov$Row.names <- NULL

GreatAPE_var<-merge(tempBCHGBov,Species_variance_copynumber$S.orangutan,by='row.names', all=TRUE)
rownames(GreatAPE_var) <- GreatAPE_var$Row.names; GreatAPE_var$Row.names <- NULL


GreatAPE_var[is.na(GreatAPE_var)] <- 0
GreatAPE_var=round(GreatAPE_var,digits = 2)


tempBCr<-merge(t(Species_range_copynumber$Bonobo),t(Species_range_copynumber$Chimpanzee),by='row.names', all=TRUE)
rownames(tempBCr) <- tempBCr$Row.names; tempBCr$Row.names <- NULL

tempBCHr<-merge(tempBCr,t(Species_range_copynumber$Human),by='row.names', all=TRUE)
rownames(tempBCHr) <- tempBCHr$Row.names; tempBCHr$Row.names <- NULL

tempBCHGr<-merge(tempBCHr,t(Species_range_copynumber$Gorilla),by='row.names', all=TRUE)
rownames(tempBCHGr) <- tempBCHGr$Row.names; tempBCHGr$Row.names <- NULL

tempBCHGBor<-merge(tempBCHGr,t(Species_range_copynumber$B.orangutan),by='row.names', all=TRUE)
rownames(tempBCHGBor) <- tempBCHGBor$Row.names; tempBCHGBor$Row.names <- NULL

GreatAPE_range<-merge(tempBCHGBor,t(Species_range_copynumber$S.orangutan),by='row.names', all=TRUE)
rownames(GreatAPE_range) <- GreatAPE_range$Row.names; GreatAPE_range$Row.names <- NULL

GreatAPE_range[is.na(GreatAPE_range)] <- 0
GreatAPE_range=round(GreatAPE_range,digits = 2)



tempBC<-merge(Species_median_copynumber$Bonobo,Species_median_copynumber$Chimpanzee,by='row.names', all=TRUE)
rownames(tempBC) <- tempBC$Row.names; tempBC$Row.names <- NULL

tempBCH<-merge(tempBC,Species_median_copynumber$Human,by='row.names', all=TRUE)
rownames(tempBCH) <- tempBCH$Row.names; tempBCH$Row.names <- NULL

tempBCHG<-merge(tempBCH,Species_median_copynumber$Gorilla,by='row.names', all=TRUE)
rownames(tempBCHG) <- tempBCHG$Row.names; tempBCHG$Row.names <- NULL

tempBCHGBo<-merge(tempBCHG,Species_median_copynumber$B.orangutan,by='row.names', all=TRUE)
rownames(tempBCHGBo) <- tempBCHGBo$Row.names; tempBCHGBo$Row.names <- NULL

GreatAPE_CAFE<-merge(tempBCHGBo,Species_median_copynumber$S.orangutan,by='row.names', all=TRUE)
rownames(GreatAPE_CAFE) <- GreatAPE_CAFE$Row.names; GreatAPE_CAFE$Row.names <- NULL

GreatAPE_CAFE[is.na(GreatAPE_CAFE)] <- 0
GreatAPE_CAFE=round(GreatAPE_CAFE,digits = 2)

#write.table(round(GreatAPE_CAFE),file="GreatAPE_CN_median_CAFE.txt", sep="\t", quote=F, col.names=T, row.names=T)

#write.table(round(tempBCHGBo),file="GreatAPEminusSorang_CN_median_CAFE.txt", sep="\t", quote=F, col.names=T, row.names=T)


#Species level information
Species_median_var<-function(data){
	test<-dcast(data,Gene~IID,value.var="CopyNumber")
	sum_val<-apply(test[,-1],2,sum)
	out<-c(median(sum_val),var(sum_val),range(sum_val))
	return(round(out,digits = 2))
}

Species_overall=llply(species,Species_median_var)
Species_overall

### Table S1 :Summary of ampliconic gene copy number across great apes. The median, variance, and range of each individual gene family and  all families together in each each species studied.

#Populate the table using the below variables.
GreatAPE_CAFE
GreatAPE_var
GreatAPE_range
Species_overall


##############PLOT

###COLOR CODE
ColorsSP.name <-  c(Bonobo="#E69F00", Chimpanzee="#F0E442" , Human="#D55E00", Gorilla="#009E73", B.orangutan="#0078D7", S.orangutan="#56B4E9")
ColorsDT.name <-  c(BPY2="#E69F00", CDY="#0078D7", DAZ="#009E73", HSFY="#D55E00", PRY="#F0E442", RBMY="#56B4E9", TSPY="#CC79A7", VCY="#9999CC", XKRY="#63a884")
ColorsSP5.name <-  c(Bonobo="#E69F00", Chimpanzee="#F0E442" , Human="#D55E00", Gorilla="#009E73", Orangutan="#0078D7")



###FIG 1B. Plot of the first two principal components (PCs) of Y ampliconic gene copy numbers across great ape species
#Read the great ape ampliconic copy number data
GA_dat<-read.table('../Data_files/GA_copy_number_ddPCR_replicates',header=T,sep="\t")

#Calculate the mean copy number across all the ddPCR reaplicates for each great ape male individual
GA_dat$mean<-apply(GA_dat[,c(3:7)],1,mean,na.rm=T)

#Read the table with great ape species information for each individual ID
great_ape<-read.table('../Data_files/great_ape_species_info',sep="\t",header=T)

#Merge the copy number information with great ape species information into one data frame
GA_dat2<-merge(GA_dat,great_ape,by="IID",sort=F)

#Prepare a data frame that displays IIDs areas rows and genes as columns and average copy number across ddPCR replicates
GA_dat3<-dcast(GA_dat2,IID~Gene,value.var="mean")
GA_dat3$IID<-as.character(GA_dat3$IID)
GA_dat3<-join(GA_dat3,great_ape,by="IID")

#Run Principal Component Analysis on copy numbers of ampliconic genes shared across great apes
pca<-prcomp(GA_dat3[,c(2:4, 7:8)],center=T,scale=T)

#Load eigenvalues and eigenvectors for all ampliconic gene copy numbers
eigval.cnv<-data.frame(V1=pca$sdev)
eigvec.cnv<-data.frame(pca$x)
eigvec.cnv$IID<-as.character(pca$IID)
eigvec.cnv$Species<-as.character(GA_dat3$Species)

#Plot PC1 versus PC2 for all ampliconic gene copy numbers
cnv.p1vp2<-ggplot(eigvec.cnv,aes(PC1,PC2,color=Species, size = 3))+geom_point(shape=21,color="black",aes(fill=Species))+theme_bw()+scale_fill_manual(values = c("#009E73", "#D55E00", "#56B4E9","#E69F00","#0078D7","#F0E442")) + labs (x = "PC1 (68.7% explained var.)", y = "PC2 (22.8% explained var.)")
ggsave('Figure_1B.pdf',cnv.p1vp2,height=7,width=7)


###Figure S1. Relationship between median and variance of copy number across all great apes. 
#Plot mean and variance for all ampliconic gene copy numbers for all great ape species
GA_dat4<-data.frame(Gene=colnames(GA_dat3)[c(2:10)],Median=apply(GA_dat3[,c(2:10)],2,median,na.rm=T),Variance=apply(GA_dat3[,c(2:10)],2,var,na.rm=T))
fig_S1<-ggplot(GA_dat4,aes(log(Median),log(Variance),color=Gene))+geom_point()+stat_smooth(method="lm",se=F,color="grey")+theme_bw()+geom_text_repel(aes(log(Median),log(Variance),label=Gene)) 
ggsave('Figure_S1.pdf',fig_S1,height=7,width=7)



###Figure S2. Principal components analysis of the overall copy number of Y ampliconic gene families across great apes.
#Plot proportion of variance explained by each PC
eigval.cnv$prop<-eigval.cnv$V1/sum(eigval.cnv$V1)
eigval.cnv$PC<-seq(1,nrow(eigval.cnv),1)

#Limit no. of PCs on the plot to 5
screecnv<-ggplot(eigval.cnv,aes(PC,prop))+geom_point()+geom_line()+theme_bw()+labs(y="Prop. of variance explained")+scale_x_continuous(breaks=seq(1,9,1),limits=c(1,5)) 
ggsave('Figure_S2.pdf',screecnv,height=7,width=7)



###Figure 2. Variation in copy number of Y ampliconic gene families in great apes.

greatAPE_summary$Species<- factor(greatAPE_summary$Species, levels = c("Bonobo","Chimpanzee","Human","Gorilla","B.orangutan","S.orangutan"))
boxP<- ggplot(greatAPE_summary, aes( x=Species, y=CopyNumber),size=0.5) +
  geom_boxplot( aes(fill = Species,color=Species), size=0.25,outlier.size=0.25 ) +
  stat_summary(fun.y=median, geom="point", size=0.2, shape=23 )+
  scale_y_continuous(breaks = round(seq(0,60, by = 5),1))+
  theme_bw() + 
  theme(                              
    axis.title.x = element_text(color="black",size=9),  axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    axis.title.y = element_text(color="black",size=9),
    plot.title = element_text(face="bold", color = "black",size=9)) +
  labs(title= "")+
  labs(x = "Ampliconic gene family")+
  labs(y = "Copy number")+labs(fill = "Species")+
  scale_fill_manual(values=ColorsSP.name)+ scale_colour_manual(values=ColorsSP.name)+
  facet_wrap( ~ Gene, nrow = 1)


grid.newpage()
pdf('Fig_SummaryCN_Boxplot.pdf', width=6.5, height=4)
boxP+theme(legend.title=element_text(size=8.5),legend.text=element_text(size=8), legend.key.size = unit(0.7,"line"),legend.position="bottom" )+ theme(strip.text = element_text(face = "italic"))+guides(fill=guide_legend(nrow=1,byrow=TRUE))
dev.off()
#NOTE: Manually added the sample size using Adobe Illustrator.




######################################################

#Figure 3. Larger Y ampliconic gene families are more variable across great apes.

####var vs copy number

CopyNumber_Median_Var_bySpecies<-function (species_data)
{
cn_table<-dcast(species_data,Gene~IID,value.var="CopyNumber")
sp_name<-species_data$Species[1]
rownames(cn_table)<-cn_table$Gene
cn_table<-cn_table[,-1] #Remove Gene column 
medianAmpliconCN=apply(cn_table,1,median)
varAmpliconCN=apply(cn_table,1,var)
geneAmpliconCN=rownames(cn_table)
summaryAmpliconCN=cbind(geneAmpliconCN,medianAmpliconCN,varAmpliconCN)
colnames(summaryAmpliconCN)=c("Gene","Median","Variance")
ACN_M_V=as.data.frame(summaryAmpliconCN,stringsAsFactors = FALSE)
ACN_M_V$Median=log(as.numeric(ACN_M_V$Median))
ACN_M_V$Variance=log(as.numeric(ACN_M_V$Variance))
id=which(ACN_M_V$Median!=-Inf)
r=cor.test(ACN_M_V$Median[id],ACN_M_V$Variance[id],method = "spearman")
print(c(round(r$estimate,digits = 2),round(r$p.value,digits = 4)))
pCNV<- ggplot(ACN_M_V, aes(y=Variance, x=Median, color=Gene)) +
  geom_smooth(data=ACN_M_V, mapping=aes(x = Median, y = Variance), col="black",method=lm,se=FALSE, size=0.5) +  
  geom_point(size=2,aes(color=Gene),shape = 18)+
  geom_text(size=2,aes(label=Gene,color=Gene),hjust=0, vjust=0, fontface="italic") +
  xlim(0,5)+ 
  ylim(-7,5)+
  theme_bw() + 
  theme(                              
    axis.title.x = element_text(color="black",size=9),
    axis.title.y = element_text(color="black",size=9),
    plot.title = element_text(face="bold", color = "black",size=10)) +     ##axis.text.x = element_text(face = "bold")
  labs(title= sp_name) +
  labs(x = "log(median copy number)") +
  labs(y = "log(variance in copy number)")+ 
  annotate("text",x=0,y=-Inf,size=2.5,vjust=-1, hjust=-0.1,label=paste("rho =",round(r$estimate,digits = 2),"(p =",round(r$p.value,digits = 4),")",by=""))+
  scale_fill_manual(values=ColorsDT.name)+ scale_colour_manual(values=ColorsDT.name)+
  theme(legend.position="none")
return(pCNV)
}
plots_species_cn_medianVSvar<-llply(species,CopyNumber_Median_Var_bySpecies)


g1=ggplot_gtable(ggplot_build(plots_species_cn_medianVSvar$Bonobo))
g2=ggplot_gtable(ggplot_build(plots_species_cn_medianVSvar$Chimpanzee))
g3=ggplot_gtable(ggplot_build(plots_species_cn_medianVSvar$Human))
g4=ggplot_gtable(ggplot_build(plots_species_cn_medianVSvar$Gorilla))
g5=ggplot_gtable(ggplot_build(plots_species_cn_medianVSvar$B.orangutan))
g6=ggplot_gtable(ggplot_build(plots_species_cn_medianVSvar$S.orangutan))


maxWidth = grid::unit.pmax(g1$widths[2:3], g2$widths[2:3], g3$widths[2:3], g4$widths[2:3], g5$widths[2:3], g6$widths[2:3])
g1$widths[2:3] <- as.list(maxWidth)
g2$widths[2:3] <- as.list(maxWidth)
g3$widths[2:3] <- as.list(maxWidth)
g4$widths[2:3] <- as.list(maxWidth)
g5$widths[2:3] <- as.list(maxWidth)
g6$widths[2:3] <- as.list(maxWidth)

grid.newpage()
pdf('Fig_logCN_Median_Var.pdf' ,width=6, height=5)
grid.arrange(arrangeGrob(g1,g2,g3,g4,g5,g6,ncol=3))
dev.off()

####NOTE: REMOVE the internal axis and arrange the gene names using Adobe Illustrator. Also, updated the p as italic. Moved labels of dots away from the line to make them visible. 



####Figure S3. Across species, gene families with higher copy number have higher variance. 

CopyNumber_Median_Var_byFamily<-function (geneFamily_data)
{
gf_table<-split(geneFamily_data, geneFamily_data$Species)
gene_name<-geneFamily_data$Gene[1]
medianAmpliconCN=lapply(gf_table,function(x){median(x$CopyNumber)})
medianAmpliconCN=ldply (medianAmpliconCN, data.frame)
varAmpliconCN=lapply(gf_table,function(x){var(x$CopyNumber)})
varAmpliconCN=ldply (varAmpliconCN, data.frame)
summaryAmpliconCN=cbind(medianAmpliconCN,varAmpliconCN[,2])
colnames(summaryAmpliconCN)=c("Species","Median","Variance")
ACN_M_V=as.data.frame(summaryAmpliconCN,stringsAsFactors = FALSE)
ACN_M_V$Median=log(as.numeric(ACN_M_V$Median))
ACN_M_V$Variance=log(as.numeric(ACN_M_V$Variance))
id=which(ACN_M_V$Median!=-Inf)
r=cor.test(ACN_M_V$Median[id],ACN_M_V$Variance[id],method = "spearman")
print(c(round(r$estimate,digits = 2),round(r$p.value,digits = 4)))
pCNV<- ggplot(ACN_M_V, aes(y=Variance, x=Median, color=Species)) +
  geom_smooth(data=ACN_M_V, mapping=aes(x = Median, y = Variance), col="black",method=lm,se=FALSE, size=0.5) +  
  geom_point(size=1,aes(color=Species),shape = 18)+
  geom_text(size=2,aes(label=Species,color=Species),hjust=0, vjust=0, fontface="italic") +
  xlim(0,5)+ 
  ylim(-7,5)+
  theme_bw() + 
  theme(                              
    axis.title.x = element_text(color="black",size=9),
    axis.title.y = element_text(color="black",size=9),
    plot.title = element_text(face="bold.italic", color = "black",size=10)) +
  labs(title= gene_name) +
  labs(x = "log(median copy number)") +
  labs(y = "log(variance in copy number)")+ 
  annotate("text",x=0,y=-Inf,size=2.5,vjust=-1, hjust=-0.1,label=paste("rho =",round(r$estimate,digits = 2),"(p =",round(r$p.value,digits = 4),")",by=""))+
  scale_fill_manual(values=ColorsSP.name)+ scale_colour_manual(values=ColorsSP.name)+
  theme(legend.position="none")
return(pCNV)
}
plots_geneFamily_cn_medianVSvar<-llply(geneFamily,CopyNumber_Median_Var_byFamily)

pdf('Fig_GF_logCN_Median_Var.pdf',width=6, height=7)
grid.arrange(plots_geneFamily_cn_medianVSvar$BPY,plots_geneFamily_cn_medianVSvar$CDY,plots_geneFamily_cn_medianVSvar$DAZ,plots_geneFamily_cn_medianVSvar$HSFY,plots_geneFamily_cn_medianVSvar$PRY, plots_geneFamily_cn_medianVSvar$RBMY, plots_geneFamily_cn_medianVSvar$TSPY,plots_geneFamily_cn_medianVSvar$XKRY, ncol=3)
dev.off()
####NOTE: REMOVE the internal axis and arrange the gene names using Adobe Illustrator. Also, updated the p as italic. Moved labels of dots away from the line to make them visible. 



#TABLE 1 ANOVA values
GeneLevel_SpeciesvsCN_ANOVA<-function(genelevel){
  anova_S <- aov(CopyNumber ~ Species, data = genelevel)
  out<-c(summary(anova_S),TukeyHSD(anova_S),pairwise.t.test(genelevel$CopyNumber,genelevel$Species, p.adjust.method = "bonferroni"))
  return(out)
}

familylevel_test=llply(geneFamily,GeneLevel_SpeciesvsCN_ANOVA)

CNpvalueANOVA<-as.data.frame(llply(familylevel_test, function(x){c(x[[1]][["F value"]][[1]],x[[1]][["Pr(>F)"]][[1]])}))






###########PERMUTATION TEST

#Table S2. P-values from permutation tests for copy number differences between Sumatran and Bornean orangutans. 
#Results Pvalue

oneMean.test <- function(data_GF,nA,nB) {
datashuffeled<-sample(data_GF)
sAmean<-apply(datashuffeled[,1:nA],1,mean)
sBmean<-apply(datashuffeled[,(nA+1):ncol(datashuffeled)],1,mean)
abs(sAmean-sBmean)
}
GreatApe_SpeciesvsCopyNumberPermutationTest<-function(speciesA,speciesB){
	set.seed(9)
	sA<-dcast(speciesA,Gene~IID,value.var="CopyNumber")
	sB<-dcast(speciesB,Gene~IID,value.var="CopyNumber")
	id=match(sA$Gene,sB$Gene)
	data_GF<-cbind(sA,sB[id,-1])
	sAmean<-apply(sA[,-1],1,mean)
	sBmean<-apply(sB[,-1],1,mean)
	trueM<-abs(sAmean-sBmean[id])
	nA<-ncol(sA)-1
	nB<-ncol(sB)-1
	permutation <- replicate(1000000, oneMean.test(data_GF[,-1], nA, nB))
	pvalue=matrix(,nrow(data_GF))
	p=matrix(,nrow(data_GF))
	for(i in seq(nrow(data_GF)))
	{
		pvalue[i]<-mean(abs(trueM[i]) < abs(permutation[i,]))
		# p[i]<-ggplot(as.data.frame(permutation[i,]), aes(x=permutation[i,])) +  geom_histogram(color="black", fill="white")+geom_vline(aes(xintercept=trueM[i]))+
			# annotate("text",x=-Inf,y=Inf,vjust=1, hjust=0,label=paste("p= ",pvalue,by=""))+
			# theme_bw() + labs(title= data_GF$Gene[i], hjust = 0) +labs(x = "Diff. mean CN")
	
	}
	out<-as.data.frame(pvalue)
	rownames(out)<-data_GF[,1]
	return(out)
	}

GreatApe_SpeciesvsCopyNumberPermutationTest(species$Bonobo,species$Chimp)
           V1
BPY  0.001668
CDY  0.000000
DAZ  0.000000
RBMY 0.000000
TSPY 0.000000

pvalue_susbspecies<-GreatApe_SpeciesvsCopyNumberPermutationTest(species$Borangutan,species$Sorangutan)
BPY  0.863196
CDY  0.747264
DAZ  0.020307
HSFY 0.194625
PRY  0.041594
RBMY 0.094442
TSPY 0.027756
XKRY 0.001324




#########CAFE analysis

#####Submit CAFE Job

##APE_CN_median_CAFE.tab is obtained from GreatAPE_CAFE rounded to nearest number.

#!cafe
#version
#date
load -i data_APE/APE_CN_median_CAFE.tab -t 10 -l reports/log_run_6APE_KYA.txt -p 0.05
tree ((((Bonobo:2432,Chimp:2432):5641,Human:8073):4633,Gorilla:12706):15346.38,(Borangutan:578,Sorangutan:578):27474.38);
lambda -s -t ((((1,1)1,1)1,1)1,(1,1)1)
report reports/report_6APE_KYA
########################################

##Table S3. The branch-level p-values showing the presence of significant shift in copy number when compared to its immediate ancestor in the great ape phylogenetic tree. 
##Table 1. CAFE P-value


data<-read.table(file="cafe_output/report_6APE_KYA.cafe", sep="\t", skip=11,stringsAsFactors=FALSE)
gf<-data$V3
CAFE_pvalue<-data$V3
x=unlist(strsplit(data$V4,"[,]"))
x=gsub("[(]","",x)
x=gsub("[)]","",x)
y<-matrix(x,ncol=10,byrow=T)
out<-apply(y,1:2,function(x){as.double(x)})
out[is.na(out)]=1
out_binary<-apply(out,1:2,function(x){if(x<=0.005){x=1}else{x=0}})
for(k in seq(length(gf))){if(as.double(gf[k])<=0.05){gf[k]=1}else{gf[k]=0}}

CAFE_pvalue<-data$V3
tableS3=out[c(2,6,7,9),]


#Generating the trees for each gene family with the size of internal nodes.

fileConn<-file("CDY_Apetree.txt")
writeLines(paste(data$V2[2],";",sep=""), fileConn)
close(fileConn)

fileConn<-file("RBMY_Apetree.txt")
writeLines(paste(data$V2[6],";",sep=""), fileConn)
close(fileConn)

fileConn<-file("TSPY_Apetree.txt")
writeLines(paste(data$V2[7],";",sep=""), fileConn)
close(fileConn)

fileConn<-file("XKRY_Apetree.txt")
writeLines(paste(data$V2[9],";",sep=""), fileConn)
close(fileConn)



####Update the species name before generating figure. You open each file and update the names.
((((Bonobo_3:1759,Chimpanzee_5:1759)_4:4079,Human_4:5838)_5:3350,Gorilla_8:9188)_7:11099.9,(Bornean_orangutan_35:418,Sumatran_orangutan_36:418)_35:19869.9)_15;

###Figure 4. Results of CAFE analysis identifying Y ampliconic gene families with significant shifts in gene copy number when compared to their ancestors.
pdf('Fig_CAFE_TREE .pdf',width=6, height=6)
library(ape)
par(mfrow=c(2,2))
CDYapes<-read.tree(file="CDY_Apetree.txt")
plot(CDYapes, show.tip.label=TRUE, show.node.label=TRUE, edge.color = c("black","black","black","black","black","black","black","red","black","black"),edge.width=c(1,1,1,1,1,1,1,3,1,1), main=substitute(paste(italic("CDY"))))

RBMYapes<-read.tree(file="RBMY_Apetree.txt")
plot(RBMYapes, show.tip.label=TRUE, show.node.label=TRUE, edge.color = c("black","black","black","red","blue","black","black","black","black","black"),edge.width=c(1,1,1,3,3,1,1,1,1,1), main=substitute(paste(italic("RBMY"))))

TSPYapes<-read.tree(file="TSPY_Apetree.txt")
plot(TSPYapes, show.tip.label=TRUE, show.node.label=TRUE, edge.color = c("black","black","black","red","blue","black","blue","black","red","blue"),edge.width=c(1,1,1,3,3,1,3,1,3,3), main=substitute(paste(italic("TSPY"))))

XKRYapes<-read.tree(file="XKRY_Apetree.txt")
plot(XKRYapes, show.tip.label=TRUE, show.node.label=TRUE, edge.color = c("black","black","black","black","black","black","black","red","blue","red"),edge.width=c(1,1,1,1,1,1,1,3,3,3), main=substitute(paste(italic("XKRY"))))

dev.off()

###NOTE: Changed the font italic for gene names and non italic for species name and rearranged the figures to look compact.


####Loading Gene Expression data

####################Gene Expression
Amp_GE<-read.table("GreatAPE_GE_Salmon_allSamples.txt",sep="\t", header=TRUE, stringsAsFactor=FALSE)
Amp_GE[8,c(2,3)]=NA
Amp_GE[4,c(2,3)]=NA
Amp_GE[5,c(2,3)]=NA
Amp_GE[9,c(2,3)]=NA
Amp_GE[4,c(4,5,6)]=NA
Amp_GE[5,c(4,5,6)]=NA
Amp_GE[9,c(4,5,6)]=NA
Amp_GE[8,c(10,11)]=NA
Amp_GE[8,c(12,13,14)]=NA

data_ApeGE<-melt(Amp_GE)


sampleID<-colnames(Amp_GE)[-1]
species <-factor(c(rep("Bonobo",2),rep("Chimpanzee",3),rep("Human",3),rep("Gorilla",2),rep("B.orangutan",2),rep("S.orangutan",1)))


data_ApeGE<-cbind(data_ApeGE,species[match(data_ApeGE$variable,sampleID)])
colnames(data_ApeGE)=c("Gene","SampleID","Expression", "Species")


speciesGE <- split(data_ApeGE  , f = data_ApeGE$Species)
geneFamilyGE<-split(data_ApeGE , f = data_ApeGE$Gene )


##Figure 5. Summary of gene expression levels across great apes.

logdata_ApeGE<-data_ApeGE
logdata_ApeGE$Expression<-log(logdata_ApeGE$Expression)
pGE<- ggplot(logdata_ApeGE, aes(y=Expression, x=Gene, color=Species)) + 
	geom_point(size=2,aes(color=Species), shape=18) +
	#ylim(0,6000)+
	theme_bw() + 
	theme(                              
	axis.title.x = element_text( color="black", size=12),
	axis.title.y = element_text( color="black", size=12),
	plot.title = element_text(face="bold", color = "black", size=12),
	axis.text.x = element_text(face = "italic")) +
	#labs(title= "Gene expression of ampliconic gene families") +
	labs(x = "Gene family") +
	labs(y = "log(gene expression)")+
	scale_fill_manual(values=ColorsSP.name)+ scale_colour_manual(values=ColorsSP.name)


pdf('Fig5_GE_summary.pdf',width=5.5, height=5)
pGE
dev.off()
#NOTE:Add sample count next to species in legend

##########################################Gene Expression and Copy Number
Amplicon_summaryD<-read.table("GreatAPE_GE_Salmon_allSamples.txt", sep="\t", stringsAsFactors=F, header=T, row.names=1)


Summary_Expression<-cbind(
Bonobo=apply(Amplicon_summaryD[,1:2],1,median),
Chimp=apply(Amplicon_summaryD[,3:5],1,median),
Human=apply(Amplicon_summaryD[,6:8],1,median),
Gorilla=apply(Amplicon_summaryD[,9:10],1,median),
Borangutan=apply(Amplicon_summaryD[,11:12],1,median)
)

Summary_Expression["VCY","Bonobo"]=0
Summary_Expression["HSFY","Bonobo"]=0
Summary_Expression["PRY","Bonobo"]=0
Summary_Expression["XKRY","Bonobo"]=0
Summary_Expression["HSFY","Chimp"]=0
Summary_Expression["PRY","Chimp"]=0
Summary_Expression["XKRY","Chimp"]=0
Summary_Expression["VCY","Gorilla"]=0
Summary_Expression["VCY","Borangutan"]=0

GreatAPE_GE=round(Summary_Expression,digits = 2)

#########################################################################################Compare both
temp_cn=melt(as.matrix(GreatAPE_CAFE[,-6]))
temp_ge=melt(as.matrix(GreatAPE_GE))

temp_ge$Var2 = as.character(temp_ge$Var2)
temp_ge$Var2[which(temp_ge$Var2=="Chimp")]="Chimpanzee"
temp_ge$Var2[which(temp_ge$Var2=="Borangutan")]="B.orangutan"


temp_ge[,1:2]==temp_cn[,1:2]

greatAPE_CN_GE<-cbind(temp_cn,temp_ge[,3])
colnames(greatAPE_CN_GE)=c("Gene", "Species","CopyNumber","Expression")
greatAPE_CN_GE_log<-greatAPE_CN_GE
greatAPE_CN_GE_log$CopyNumber=log(greatAPE_CN_GE_log$CopyNumber)
greatAPE_CN_GE_log$Expression=log(greatAPE_CN_GE_log$Expression)

GAgeneFamilies <- split( greatAPE_CN_GE , f = greatAPE_CN_GE$Gene )
GAspecies <- split( greatAPE_CN_GE , f = greatAPE_CN_GE$Species )


###PLOT CN vs GE

##Figure 6. Relationship between copy number and gene expression of Y ampliconic gene families in great ape species.

logCopyNumber_GeneExpression_bySpecies<-function (species_data)
{
species_data$CopyNumber=log(species_data$CopyNumber)
species_data$Expression=log(species_data$Expression)
id=which(species_data$CopyNumber!=-Inf)
r=cor.test(species_data$CopyNumber[id],species_data$Expression[id],method = "spearman")
sp_name=species_data$Species[1]
idg=which(species_data$Expression>0) #filter genes with less than zero expression level.
species_data_sub<-species_data[idg,]
print(c(round(r$estimate,digits = 2),round(r$p.value,digits = 4)))
pCNV<- ggplot(species_data_sub, aes(y=Expression, x=CopyNumber, color=Gene)) +
  geom_smooth(data=species_data_sub, mapping=aes(x = CopyNumber, y = Expression), col="black",method=lm,se=FALSE, size=0.5) +  
  geom_point(size=2,aes(color=Gene), shape=18)+
  geom_text(size=2,aes(label=Gene,color=Gene),hjust=0, vjust=0, fontface="italic") +
  xlim(0,5)+ 
  ylim(0,8)+
  theme_bw() + 
  theme(                              
    axis.title.x = element_text(color="black",size=9),
    axis.title.y = element_text(color="black",size=9),
    plot.title = element_text(face="bold", color = "black",size=10)) +
  labs(title= sp_name) +
  labs(x = "log(copy number)") +
  labs(y = "log(gene expression)")+ 
  annotate("text",x=1,y=-Inf,size=2.5,vjust=-1, hjust=-0.1,label=paste("rho =",round(r$estimate,digits = 2),"(p =",round(r$p.value,digits = 4),")",by=""))+
  scale_fill_manual(values=ColorsDT.name)+ scale_colour_manual(values=ColorsDT.name)  
return(pCNV)
}


plots_species_cn_ge<-llply(GAspecies,logCopyNumber_GeneExpression_bySpecies)

pdf('Fig_logCN_logGE_Median.pdf' ,width=6, height=5)
grid.arrange(plots_species_cn_ge$Bonobo+ theme(legend.position="none"), plots_species_cn_ge$Chimpanzee+ theme(legend.position="none"), plots_species_cn_ge$Human+ theme(legend.position="none"), plots_species_cn_ge$Gorilla+ theme(legend.position="none"), plots_species_cn_ge$B.orangutan+ theme(legend.position="none"), ncol=3)
dev.off()
####NOTE: REMOVE the internal axis and arrange the gene names using Adobe Illustrator. Also, updated the p as italic. Moved labels of dots away from the line to make them visible. 




##Figure 7. Relationship between copy number and gene expression across species.
logCopyNumber_GeneExpression_byGeneFamily<-function (genefam_data)
{
genefam_data$CopyNumber=log(genefam_data$CopyNumber)
genefam_data$Expression=log(genefam_data$Expression)
id=which(genefam_data$CopyNumber!=-Inf)
r=cor.test(genefam_data$CopyNumber[id],genefam_data$Expression[id],method = "spearman")
#r=cor.test(genefam_data$CopyNumber,genefam_data$Expression,method = "spearman")
gf_name=genefam_data$Gene[1]
print(c(round(r$estimate,digits = 2),round(r$p.value,digits = 4)))
pCNV<- ggplot(genefam_data, aes(y=Expression, x=CopyNumber, color=Species)) +
  geom_smooth(data=genefam_data, mapping=aes(x = CopyNumber, y = Expression), col="black",method=lm,se=FALSE, size=0.5) +  
  geom_point(size=2,aes(color=Species), shape=18)+
  geom_text(size=2,aes(label=Species,color=Species),hjust=0, vjust=0, fontface="italic") +
  xlim(0,5)+ 
  ylim(0,8)+
  theme_bw() + 
  theme(                              
    axis.title.x = element_text(color="black",size=9),
    axis.title.y = element_text(color="black",size=9),
    plot.title = element_text(face="bold", color = "black",size=10)) +
  labs(title= gf_name) +
  labs(x = "log(copy number)") +
  labs(y = "log(gene expression)")+ 
  annotate("text",x=1,y=-Inf,size=2.5,vjust=-1, hjust=-0.1,label=paste("rho =",round(r$estimate,digits = 2),"(p =",round(r$p.value,digits = 4),")",by=""))+
  scale_fill_manual(values=ColorsSP.name)+ scale_colour_manual(values=ColorsSP.name)  +  theme(legend.position="none")
return(pCNV)
}
plots_genefam_cn_ge<-llply(GAgeneFamilies,logCopyNumber_GeneExpression_byGeneFamily)

pdf('Fig_byGF_logCN_logGE_Median.pdf' ,width=6, height=5)
grid.arrange(plots_genefam_cn_ge$BPY2, plots_genefam_cn_ge$CDY, plots_genefam_cn_ge$DAZ, plots_genefam_cn_ge$RBMY, plots_genefam_cn_ge$TSPY, ncol=3)
dev.off()

####NOTE: REMOVE the internal axis and arrange the gene names using Adobe Illustrator. Also, updated the p as italic. Moved labels of dots away from the line to make them visible. 


####Supplemental Figure SA. The great ape phylogenetic tree with five individuals sampled per species (star phylogeny) used in CAFE analysis.
library(ape)
test<-read.tree("tree_star5_phylogeny.nwk")

pdf('Tree_5star_phylogeny.pdf',width=4, height=6)
plot(test,cex=0.8)
dev.off()
####NOTE: Add numbers on edges and remove white space.




#########################GENE EXPRESSION PROCESSING########################################################################

#ssh -X sec

#cd /nfs/secure/scratch6/nekrut_gtex/rxv923/Transcript_Assembly/GreatApe_Summary
#source activate renv



library(tximport)
tar2gen<-read.table("/nfs/brubeck.bx.psu.edu/scratch6/rahul/Apes/data/Amplicon/ReadCountSummary/hg38_refFlat_10252016_GenetoID_appendCDS.txt", header=F, sep="\t", stringsAsFactors=F)
names(tar2gen)=c("target_id","gene")

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



#list_dir <- c(sal_dirsBo,sal_dirsCh,sal_dirsHS,sal_dirsGr,sal_dirsBOr)
list_dir <- c(sal_dirsBo,sal_dirsCh,sal_dirsHS,sal_dirsGr,sal_dirsBOr,sal_dirSOr)
gene_level=tximport(list_dir, type = "salmon",countsFromAbundance = "no", tx2gene = tar2gen, geneIdCol="gene",txIdCol="target_id" )

ApesData=gene_level$counts
ApesData<-as.matrix(ApesData)
ApesData=apply(ApesData, 1:2, round)
head(ApesData)

ad<-ApesData[,c(1,4,5,6,7,8,9,10,11,13,14,17,19)]
#species <-factor(c(rep("Bonobo",2),rep("Chimpanzee",3),rep("Human",3),rep("Gorilla",2),rep("Borang",2)))
species <-factor(c(rep("Bonobo",2),rep("Chimpanzee",3),rep("Human",3),rep("Gorilla",2),rep("Borang",2),rep("Sorang",1)))
#type <- factor(c("paired_end","single_end","paired_end","paired_end","single_end","paired_end","paired_end","paired_end","paired_end","single_end","paired_end","paired_end"))
type <- factor(c("paired_end","single_end","paired_end","paired_end","single_end","paired_end","paired_end","paired_end","paired_end","single_end","paired_end","paired_end","paired_end"))

#sample<-c("B1","B2","C1","C2","C3","H1","H2","H3","G1","G2","O1","O2")
sample<-c("B1","B2","C1","C2","C3","H1","H2","H3","G1","G2","BO1","BO2","S01")
run<- paste0("run",1:12)

cd=data.frame(cbind(species,type,sample),stringsAsFactors = FALSE)
row.names(cd)=colnames(ad)
head(cd)

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = ad, colData = cd, design = ~ type+species)
dds <- DESeq(dds)

rld <- rlog(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)

#save(dds,rld,vsd, file = "DESeq2_greatApe_normalized_dds_rld_vsd_Salmon.RData")

d<-counts(dds, normalized=TRUE)
r<-assay(rld)
v<-assay(vsd)
head(r)
head(v)

#save(dds,d,v,r, file = "DESeq_greatApe_normalized_Salmon.RData")
#write.table(d,file="greatApe_Counts_Testis_2-3samples.tab", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
#write.table(r,file="greatApe_RLDCounts_Testis_2-3samples.tab", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
#write.table(v,file="greatApe_VSTCounts_Testis_2-3samples.tab", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)


ampliconList=c("BPY2B","BPY2","BPY2C","BPY2_CDS","CDY1","CDY1B","CDY2A","CDY2B","CDY_CDS","DAZ1","DAZ2","DAZ3","DAZ4","DAZ_CDS","HSFY1","HSFY2","HSFY_CDS","PRY","PRY2","PRY_CDS","RBMY1A1","RBMY1J","RBMY1F","RBMY1E","RBMY1D","RBMY1B","RBMY_CDS","TSPY10","TSPY1","TSPY3","TSPY8","TSPY4","TSPY2","TSPY_CDS","VCY","VCY1B","VCY1_CDS","XKRY","XKRY2")

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

write.table(Amplicon_summaryD,file="Apes_normalizedDDS_DESeq2_Amp_Salmon.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
fh <- file("Apes_WithoutNA_normalized_DESeq2_Amp_Salmon.dat", open="wt")
#fh <- file("Apes_WithoutNA_normalized_DESeq2_Amp_k25.dat", open="wt")
writeLines("5", fh)
write.table(Amplicon_summaryD[c(1,2,3,6,7),], fh ,sep=" ", quote=FALSE, row.names=TRUE,col.names=FALSE)
close(fh)
		 

###To run EVE

##
##cat data_Apes/exampleNindivs_2-3each.nindiv
##2 3 3 2 2

#EVE run
./EVEmodel -S -n 5 -t data_Apes/Great_apest_Ychr_PPPY.awk -i data_Apes/exampleNindivs_2-3each.nindiv -d data_Apes/Apes_WithoutNA_normalized_DESeq2_Amp_Salmon.dat -f _DESeq2RealDataGreatApes_Salmon -v 10

R
x=read.table("results/betaTestLRTs_DESeq2RealDataGreatApes_Salmon.res",sep=" ", header=FALSE, stringsAsFactors=FALSE)
x
pchisq(as.double(x), df=1, lower.tail=FALSE)
b=read.table("results/indivBetaMLparams_DESeq2RealDataGreatApes.res",sep=" ", header=FALSE, stringsAsFactors=FALSE)
log(b$V4)



################################################JUST to plot

load("DESeq2_greatApe_normalized_dds_rld_vsd_Salmon.RData")



library("DESeq2")
library("pheatmap")
library("RColorBrewer")

pdf('SupFig_GE_rld_PCAplot.pdf',width=6, height=6)
plotPCA(rld, intgroup=c("species", "sample"))
dev.off()


pdf('SupFig_GE_vsd_PCAplot.pdf',width=6, height=6)
plotPCA(vsd, intgroup=c("species", "sample"))
dev.off()

sampleDists <- dist(t(assay(rld)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$sample, rld$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf('SupFig_GE_rld_heatmap.pdf',width=6, height=6)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()


sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$sample, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf('SupFig_GE_vsd_heatmap.pdf',width=6, height=6)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()



####################Phenotype PLOTS

########Figure S4. Phenotypes related to sperm competition in great apes.
phenotype_APE<-read.table("Phenotype_sperm.txt",sep="\t", stringsAsFactors=F, header=T)

gMV<-ggplot(data=phenotype_APE, aes(x=ID, y=Midpiece_volume,color=ID)) +
  theme_bw() +
  geom_point()+
  labs(y = "MidpieceVol.")+
  labs(x = "Species")+
  scale_x_discrete(limits=c("Bonobo","Chimpanzee","Human","Orangutan","Gorilla"))+theme(legend.position="none") +
  scale_fill_manual(values=ColorsSP5.name)+ scale_colour_manual(values=ColorsSP5.name) 

gTW<-ggplot(data=phenotype_APE, aes(x=ID, y=Residual_testes_weight,color=ID)) +
  theme_bw() +
  geom_point()+
  labs(y = "ResTW")+
  labs(x = "Species")+
  scale_x_discrete(limits=c("Bonobo","Chimpanzee","Human","Orangutan","Gorilla"))+theme(legend.position="none")  +
  scale_fill_manual(values=ColorsSP5.name)+ scale_colour_manual(values=ColorsSP5.name) 

gSI<-ggplot(data=phenotype_APE, aes(x=ID, y=Spermatogenic_Index,color=ID)) +
  theme_bw() +
  geom_point()+
  labs(y = "Sperm. index")+
  labs(x = "Species")+
  scale_x_discrete(limits=c("Bonobo","Chimpanzee","Human","Orangutan","Gorilla"))+theme(legend.position="none") +
  scale_fill_manual(values=ColorsSP5.name)+ scale_colour_manual(values=ColorsSP5.name) 
  
gSC<-ggplot(data=phenotype_APE, aes(x=ID, y=Sperm_conc,color=ID)) +
  theme_bw() +
  geom_point()+
  labs(y = "Sperm. conc")+
  labs(x = "Species")+
  scale_x_discrete(limits=c("Bonobo","Chimpanzee","Human","Orangutan","Gorilla"))+theme(legend.position="none") +
  scale_fill_manual(values=ColorsSP5.name)+ scale_colour_manual(values=ColorsSP5.name) 

gSM<-ggplot(data=phenotype_APE, aes(x=ID, y=Sperm_Motility,color=ID)) +
  theme_bw() +
  geom_point()+
  labs(y = "Sperm. motility")+
  labs(x = "Species")+
  scale_x_discrete(limits=c("Bonobo","Chimpanzee","Human","Orangutan","Gorilla"))+theme(legend.position="none") +
  scale_fill_manual(values=ColorsSP5.name)+ scale_colour_manual(values=ColorsSP5.name) 



g1=ggplot_gtable(ggplot_build(gMV))
g2=ggplot_gtable(ggplot_build(gTW))
g3=ggplot_gtable(ggplot_build(gSC))
g4=ggplot_gtable(ggplot_build(gSM))

g = rbind(g1, g2, g3, g4, size = "last")
grid.newpage()
pdf('FigSup_phenotype_four.pdf',width=6, height=6)
grid.draw(g)
dev.off()

####Figure S5. Copy number variation in CDY, RBMY, TSPY, and XKRY.
CopyNumber_Var_figure_genefamily<-function (genefam_data)
{
copynumber_GeneFamily<-genefam_data
copynumber_GeneFamily$Species = as.character(copynumber_GeneFamily$Species)
copynumber_GeneFamily$Species[which(copynumber_GeneFamily$Species=="B.orangutan")]="Orangutan"
copynumber_GeneFamily<-copynumber_GeneFamily[ ! copynumber_GeneFamily$Species %in% c("S.orangutan"), ]
x<-copynumber_GeneFamily %>% group_by(Species) %>% summarize(Variance = var(CopyNumber, na.rm = TRUE))
x<-as.data.frame(x)
gCN<-ggplot(data=copynumber_GeneFamily, aes(x=Species, y=CopyNumber,color=Species)) +
  theme_bw() +
  geom_boxplot()+
  geom_point(data=x,aes(x=Species, y=Variance),color="black", shape=18,size=3)+
  labs(y = "Copy number")+
  labs(x = "Species")+
  scale_x_discrete(limits=c("Bonobo","Chimpanzee","Human","Orangutan","Gorilla"))+theme(legend.position="none") +
  scale_fill_manual(values=ColorsSP5.name)+ scale_colour_manual(values=ColorsSP5.name)
return(gCN)  
}

plots_genefam_cnVar_dotplot<-llply(geneFamily,CopyNumber_Var_figure_genefamily)

gC3=ggplot_gtable(ggplot_build(plots_genefam_cnVar_dotplot$CDY))
gR3=ggplot_gtable(ggplot_build(plots_genefam_cnVar_dotplot$RBMY))
gT3=ggplot_gtable(ggplot_build(plots_genefam_cnVar_dotplot$TSPY))
gX3=ggplot_gtable(ggplot_build(plots_genefam_cnVar_dotplot$XKRY))

grid.newpage()
pdf('Fig_CN_andVAR_phenotype.pdf',width=6, height=6)
gVar = rbind(gC3,gR3,gT3,gX3, size = "last")
grid.draw(gVar)
dev.off()

# var_data<-melt(as.matrix(GreatAPE_var))
# VARgeneFamilies <- split( var_data , f = var_data$Var1 )
# Var_figure_genefamily<-function (genefam_data)
# {
# copynumber_GeneFamily<-genefam_data[-6,]
# copynumber_GeneFamily$Var2 = as.character(copynumber_GeneFamily$Var2)
# copynumber_GeneFamily$Var2[which(copynumber_GeneFamily$Var2=="B.orangutan")]="Orangutan"
# gCN<-ggplot(data=copynumber_GeneFamily, aes(x=Var2, y=value,color=Var2)) +
  # theme_bw() +
  # geom_point()+
  # labs(y = "var(CN)")+
  # labs(x = "Species")+
  # scale_x_discrete(limits=c("Bonobo","Chimpanzee","Human","Orangutan","Gorilla"))+theme(legend.position="none") +
  # scale_fill_manual(values=ColorsSP5.name)+ scale_colour_manual(values=ColorsSP5.name)
# return(gCN)  
# }

# plots_genefam_cn_var_dotplot<-llply(VARgeneFamilies,Var_figure_genefamily)

# GE_figure_genefamily<-function (genefam_data) {
# geneexp_GeneFamily<-genefam_data
# geneexp_GeneFamily$Species = as.character(geneexp_GeneFamily$Species)
# geneexp_GeneFamily$Species[which(geneexp_GeneFamily$Species=="B.orangutan")]="Orangutan"
# gGE<-ggplot(data=geneexp_GeneFamily, aes(x=Species, y=Expression,color=Species)) +
  # theme_bw() +
  # geom_boxplot()+
  # labs(y = "Gene expression")+
  # labs(x = "Species")+
  # scale_x_discrete(limits=c("Bonobo","Chimpanzee","Human","Orangutan","Gorilla"))+theme(legend.position="none")  +
  # scale_fill_manual(values=ColorsSP5.name)+ scale_colour_manual(values=ColorsSP5.name)  
 # return(gGE)
# }  


# plots_genefam_ge_dotplot<-llply(geneFamilyGE,GE_figure_genefamily)

CopyNumber_figure_genefamily<-function (genefam_data)
{
copynumber_GeneFamily<-genefam_data
copynumber_GeneFamily$Species = as.character(copynumber_GeneFamily$Species)
copynumber_GeneFamily$Species[which(copynumber_GeneFamily$Species=="B.orangutan")]="Orangutan"
copynumber_GeneFamily<-copynumber_GeneFamily[ ! copynumber_GeneFamily$Species %in% c("S.orangutan"), ]
gCN<-ggplot(data=copynumber_GeneFamily, aes(x=Species, y=CopyNumber,color=Species)) +
  theme_bw() +
  geom_boxplot()+
  labs(y = "Copy number")+
  labs(x = "Species")+
  scale_x_discrete(limits=c("Bonobo","Chimpanzee","Human","Orangutan","Gorilla"))+theme(legend.position="none") +
  scale_fill_manual(values=ColorsSP5.name)+ scale_colour_manual(values=ColorsSP5.name)
return(gCN)  
}
plots_genefam_cn_dotplot<-llply(geneFamily,CopyNumber_figure_genefamily)

gR3=ggplot_gtable(ggplot_build(plots_genefam_cn_dotplot$RBMY))
gT4=ggplot_gtable(ggplot_build(plots_genefam_cn_dotplot$TSPY))
g2 = rbind(g3, g4, gR3, gT4, size = "last")
grid.newpage()
pdf('FigSup_phenotype_RBMYTSPY.pdf',width=6, height=6)
grid.draw(g2)
dev.off()


# gT4=ggplot_gtable(ggplot_build(plots_genefam_cn_dotplot$TSPY))
# gT5=ggplot_gtable(ggplot_build(plots_genefam_ge_dotplot$TSPY))
# gT = rbind(g1, g2, g3,g4, g5, size = "last")
# grid.newpage()
# grid.draw(gT)

# grid.newpage()
# pdf('Fig_BPY2_phenotype.pdf',width=6, height=6)
# gB3=ggplot_gtable(ggplot_build(plots_genefam_cn_dotplot$BPY2))
# gB4=ggplot_gtable(ggplot_build(plots_genefam_cn_var_dotplot$BPY2))
# gB5=ggplot_gtable(ggplot_build(plots_genefam_ge_dotplot$BPY2))
# gB = rbind(gB3,gB4, gB5, size = "last")
# grid.draw(gB)
# dev.off()

# grid.newpage()
# pdf('Fig_CDY_phenotype.pdf',width=6, height=6)
# gC3=ggplot_gtable(ggplot_build(plots_genefam_cn_dotplot$CDY))
# gC4=ggplot_gtable(ggplot_build(plots_genefam_cn_var_dotplot$CDY))
# gC5=ggplot_gtable(ggplot_build(plots_genefam_ge_dotplot$CDY))
# gC = rbind(gC3,gC4, gC5, size = "last")
# grid.draw(gC)
# dev.off()

# grid.newpage()
# pdf('Fig_DAZ_phenotype.pdf',width=6, height=6)
# gD4=ggplot_gtable(ggplot_build(plots_genefam_cn_dotplot$DAZ))
# gD5=ggplot_gtable(ggplot_build(plots_genefam_ge_dotplot$DAZ))
# gD = rbind(g1, g2, g3,gD4, gD5, size = "last")
# grid.draw(gD)
# dev.off()

# grid.newpage()
# pdf('Fig_PRY_phenotype.pdf',width=6, height=6)
# gP4=ggplot_gtable(ggplot_build(plots_genefam_cn_dotplot$PRY))
# gP5=ggplot_gtable(ggplot_build(plots_genefam_ge_dotplot$PRY))
# gP = rbind(g1, g2, g3,gP4, gP5, size = "last")
# grid.draw(gP)
# dev.off()

# grid.newpage()
# pdf('Fig_HSFY_phenotype.pdf',width=6, height=6)
# gH4=ggplot_gtable(ggplot_build(plots_genefam_cn_dotplot$HSFY))
# gH5=ggplot_gtable(ggplot_build(plots_genefam_ge_dotplot$HSFY))
# gH = rbind(g1, g2, g3,gH4, gH5, size = "last")
# grid.draw(gH)
# dev.off()

# grid.newpage()
# pdf('Fig_RBMY_phenotype.pdf',width=6, height=6)
# gR3=ggplot_gtable(ggplot_build(plots_genefam_cn_dotplot$RBMY))
# gR4=ggplot_gtable(ggplot_build(plots_genefam_cn_var_dotplot$RBMY))
# gR5=ggplot_gtable(ggplot_build(plots_genefam_ge_dotplot$RBMY))
# gR = rbind(gR3,gR4, gR5, size = "last")
# grid.draw(gR)
# dev.off()

# grid.newpage()
# pdf('Fig_TSPY_phenotype.pdf',width=6, height=6)
# gT3=ggplot_gtable(ggplot_build(plots_genefam_cn_dotplot$TSPY))
# gT4=ggplot_gtable(ggplot_build(plots_genefam_cn_var_dotplot$TSPY))
# gT5=ggplot_gtable(ggplot_build(plots_genefam_ge_dotplot$TSPY))
# gT = rbind(gT3,gT4, gT5, size = "last")
# grid.draw(gT)
# dev.off()


# grid.newpage()
# pdf('Fig_VCY_phenotype.pdf',width=6, height=6)
# gV4=ggplot_gtable(ggplot_build(plots_genefam_cn_dotplot$VCY))
# gV5=ggplot_gtable(ggplot_build(plots_genefam_ge_dotplot$VCY))
# gV = rbind(g1, g2, g3,gV4, gV5, size = "last")
# grid.draw(gV)
# dev.off()

# grid.newpage()
# pdf('Fig_XKRY_phenotype.pdf',width=6, height=6)
# gX3=ggplot_gtable(ggplot_build(plots_genefam_cn_dotplot$XKRY))
# gX4=ggplot_gtable(ggplot_build(plots_genefam_cn_var_dotplot$XKRY))
# gX5=ggplot_gtable(ggplot_build(plots_genefam_ge_dotplot$XKRY))
# gX = rbind(gX3,gX4, gX5, size = "last")
# grid.draw(gX)
# dev.off()


# grid.newpage()
# pdf('Fig_VARCN_phenotype.pdf',width=6, height=6)
# gVar = rbind(gC4,gR4,gT4,gX4, size = "last")
# grid.draw(gVar)
# dev.off()

# grid.newpage()
# pdf('Fig_ALLGENECN_phenotype.pdf',width=6, height=6)
# gall = rbind(gC3,gC4,gR3,gR4, gR5, gC5,gT3,gT4, gT5,gX3,gX4, gX5, size = "last")
# grid.draw(gall)
# dev.off()



