############################################################################################################
############################################################################################################
############################################################################################################
###in below folder the slurm scripts were run
#/nfs/brubeck.bx.psu.edu/scratch6/rahul/Apes/simulationRori_EVE/EVE_release/slurm_jobs_simulation/Jobs_5genes/
#SAMPLE CODE
#for i in `seq 1 50`
#do
#  ./EVEmodel -H -n 5 -t data_Apes/Great_apest.awk -i data_Apes/exampleNindivs_2-3each.nindiv -h 100 -g 5 -a 3 -j 24 -f _5genes_SimRoripapertau40_simulationReplicate${i} -v 10
#  echo "output: $i"
#done
#Once we have the output run the below code to parse the output files and generate plots.
############################################################################################################

##The sample size for each species used
#cat exampleNindivs_2-3each.nindiv
#2 3 3 2 2

##The phylogenetic tree used (mean from Fig 1 -> Locke, Devin P., et al. "Comparative and demographic analysis of orang-utan genomes." Nature 469.7331 (2011): 529.)
#cat Great_apest.awk
#5
#((((1:100,2:100):425,3:525):175,4:700):800,5:1500);


#Folder name might change keep the path updated
path="/nfs/brubeck.bx.psu.edu/scratch6/rahul/Apes/simulationRori_EVE/EVE_release_simulation/results_simulation2-3each5species5gene"  #results_simulation2-3each_withoutNAgenes
setwd("/nfs/brubeck.bx.psu.edu/scratch6/rahul/Apes/simulationRori_EVE/EVE_release_simulation/results_simulation2-3each5species5gene")

logRatios_sig0.5 <- list.files(path, pattern="betaTestLRTs5genes_theta100.00_sig20.50_alp3.00_beta60.00_5genes_SimRoripapersigma0.5_simulationReplicate")
logRatios_sig1 <- list.files(path, pattern="betaTestLRTs5genes_theta100.00_sig21.00_alp3.00_beta30.00_5genes_SimRoripapersigma1_simulationReplicate")
logRatios_sig3 <- list.files(path, pattern="betaTestLRTs5genes_theta100.00_sig23.00_alp3.00_beta10.00_5genes_SimRoripapersigma3_simulationReplicate")
logRatios_sig5 <- list.files(path, pattern="betaTestLRTs5genes_theta100.00_sig25.00_alp3.00_beta6.00_5genes_SimRoripapersigma5_simulationReplicate")
logRatios_sig10 <- list.files(path, pattern="betaTestLRTs5genes_theta100.00_sig210.00_alp3.00_beta3.00_5genes_SimRoripapersigma10_simulationReplicate")
logRatios_sig15 <- list.files(path, pattern="betaTestLRTs5genes_theta100.00_sig215.00_alp3.00_beta2.00_5genes_SimRoripapersigma15_simulationReplicate")
logRatios_sig20 <- list.files(path, pattern="betaTestLRTs5genes_theta100.00_sig220.00_alp3.00_beta1.33_5genes_SimRoripapersigma20_simulationReplicate")
logRatios_sig40 <- list.files(path, pattern="betaTestLRTs5genes_theta100.00_sig240.00_alp3.00_beta0.75_5genes_SimRoripapersigma40_simulationReplicate")
logRatios_sig60 <- list.files(path, pattern="betaTestLRTs5genes_theta100.00_sig260.00_alp3.00_beta0.50_5genes_SimRoripapersigma60_simulationReplicate")

logRatios_alpha0.5 <- list.files(path, pattern="betaTestLRTs5genes_theta100.00_sig25.00_alp0.50_beta1.00_5genes_SimRoripaperalpha0.5_simulationReplicate")
logRatios_alpha1 <- list.files(path, pattern="betaTestLRTs5genes_theta100.00_sig25.00_alp1.00_beta2.00_5genes_SimRoripaperalpha1_simulationReplicate")
logRatios_alpha3 <- list.files(path, pattern="betaTestLRTs5genes_theta100.00_sig25.00_alp3.00_beta6.00_5genes_SimRoripaperalpha3_simulationReplicate")
logRatios_alpha5 <- list.files(path, pattern="betaTestLRTs5genes_theta100.00_sig25.00_alp5.00_beta10.00_5genes_SimRoripaperalpha5_simulationReplicate")
logRatios_alpha10 <- list.files(path, pattern="betaTestLRTs5genes_theta100.00_sig25.00_alp10.00_beta20.00_5genes_SimRoripaperalpha10_simulationReplicate")
logRatios_alpha15 <- list.files(path, pattern="betaTestLRTs5genes_theta100.00_sig25.00_alp15.00_beta30.00_5genes_SimRoripaperalpha15_simulationReplicate")
logRatios_alpha20 <- list.files(path, pattern="betaTestLRTs5genes_theta100.00_sig25.00_alp20.00_beta40.00_5genes_SimRoripaperalpha20_simulationReplicate")
logRatios_alpha40 <- list.files(path, pattern="betaTestLRTs5genes_theta100.00_sig25.00_alp40.00_beta80.00_5genes_SimRoripaperalpha40_simulationReplicate")
logRatios_alpha60 <- list.files(path, pattern="betaTestLRTs5genes_theta100.00_sig25.00_alp60.00_beta120.00_5genes_SimRoripaperalpha60_simulationReplicate")


logRatios_tau0.5 <- list.files(path, pattern="betaTestLRTs5genes_theta100.00_sig25.00_alp3.00_beta0.30_5genes_SimRoripapertau0.5_simulationReplicate")
logRatios_tau1 <- list.files(path, pattern="betaTestLRTs5genes_theta100.00_sig25.00_alp3.00_beta0.60_5genes_SimRoripapertau1_simulationReplicate")
logRatios_tau3 <- list.files(path, pattern="betaTestLRTs5genes_theta100.00_sig25.00_alp3.00_beta1.80_5genes_SimRoripapertau3_simulationReplicate")
logRatios_tau5 <- list.files(path, pattern="betaTestLRTs5genes_theta100.00_sig25.00_alp3.00_beta3.00_5genes_SimRoripapertau5_simulationReplicate")
logRatios_tau10 <- list.files(path, pattern="betaTestLRTs5genes_theta100.00_sig25.00_alp3.00_beta6.00_5genes_SimRoripapertau10_simulationReplicate")
logRatios_tau15 <- list.files(path, pattern="betaTestLRTs5genes_theta100.00_sig25.00_alp3.00_beta9.00_5genes_SimRoripapertau15_simulationReplicate")
logRatios_tau20 <- list.files(path, pattern="betaTestLRTs5genes_theta100.00_sig25.00_alp3.00_beta12.00_5genes_SimRoripapertau20_simulationReplicate")
logRatios_tau40 <- list.files(path, pattern="betaTestLRTs5genes_theta100.00_sig25.00_alp3.00_beta24.00_5genes_SimRoripapertau40_simulationReplicate")
logRatios_tau60 <- list.files(path, pattern="betaTestLRTs5genes_theta100.00_sig25.00_alp3.00_beta36.00_5genes_SimRoripapertau60_simulationReplicate")



readLogRatio_SummaryCounts <- function(logRatios, threshold){
	listOflogRatios <- lapply(logRatios, function(x) read.table(x,header=FALSE, sep=" ", stringsAsFactors=FALSE)) 
	df_logRation <- data.frame(matrix(unlist(listOflogRatios),nrow=6,byrow=F))
	df_logRation_final <- df_logRation[complete.cases(df_logRation), ]
	logRation.pvalue=apply(df_logRation_final,1:2,function(x){pchisq(x, df=1, lower.tail=FALSE)})
	threshold_binary=apply(logRation.pvalue,1:2,function(x){ifelse(x<threshold,0,1) })
	apply(threshold_binary,1,function(x){sum(x)/length(x)})
}

thres<-0.05

counts_sigma_pvalue0.5<-cbind(readLogRatio_SummaryCounts(logRatios_sig0.5,thres),
readLogRatio_SummaryCounts(logRatios_sig1,thres),
readLogRatio_SummaryCounts(logRatios_sig3,thres),
readLogRatio_SummaryCounts(logRatios_sig5,thres),
readLogRatio_SummaryCounts(logRatios_sig10,thres),
readLogRatio_SummaryCounts(logRatios_sig15,thres),
readLogRatio_SummaryCounts(logRatios_sig20,thres),
readLogRatio_SummaryCounts(logRatios_sig40,thres),
readLogRatio_SummaryCounts(logRatios_sig60,thres))
colnames(counts_sigma_pvalue0.5)=c("0.5","1","3","5","10","15","20","40","60")


counts_alpha_pvalue0.5<-cbind(
readLogRatio_SummaryCounts(logRatios_alpha0.5,thres),
readLogRatio_SummaryCounts(logRatios_alpha1,thres),
readLogRatio_SummaryCounts(logRatios_alpha3,thres),
readLogRatio_SummaryCounts(logRatios_alpha5,thres),
readLogRatio_SummaryCounts(logRatios_alpha10,thres),
readLogRatio_SummaryCounts(logRatios_alpha15,thres),
readLogRatio_SummaryCounts(logRatios_alpha20,thres),
readLogRatio_SummaryCounts(logRatios_alpha40,thres),
readLogRatio_SummaryCounts(logRatios_alpha60,thres))
colnames(counts_alpha_pvalue0.5)=c("0.5","1","3","5","10","15","20","40","60")


counts_tau_pvalue0.5<-cbind(
readLogRatio_SummaryCounts(logRatios_tau0.5,thres),
readLogRatio_SummaryCounts(logRatios_tau1,thres),
readLogRatio_SummaryCounts(logRatios_tau3,thres),
readLogRatio_SummaryCounts(logRatios_tau5,thres),
readLogRatio_SummaryCounts(logRatios_tau10,thres),
readLogRatio_SummaryCounts(logRatios_tau15,thres),
readLogRatio_SummaryCounts(logRatios_tau20,thres),
readLogRatio_SummaryCounts(logRatios_tau40,thres),
readLogRatio_SummaryCounts(logRatios_tau60,thres))
colnames(counts_tau_pvalue0.5)=c("0.5","1","3","5","10","15","20","40","60")

breaks <- c(0.0,0.02,0.05,0.1,0.2,0.4,0.6,0.8,0.9,0.95,0.98,1.00)
col = c("deepskyblue","dodgerblue3","blue","darkblue","darkgreen","green2","greenyellow","yellow","orange","red","red4")
leg.txt <- c("0-0.02","0.02-0.05","0.05-0.10","0.10-0.20", "0.20-0.40", "0.40-0.60", "0.60-0.80","0.80-0.90","0.90-0.95","0.95-0.98","0.98-1.00")

pdf("SupFig_P-value_heatmap_LEGEND.pdf", width=6, height=6)
matplot(c(1, 8), c(0, 4.5), type = "n")
legend(6, 3, leg.txt, fill=col,cex=0.9,title="%Pval > Thresh")
dev.off()
#The legend from above plot was added manually to the figures below.
				  
pdf("SupFig_P-value_heatmap_Sigmaparameter.pdf", width=6, height=6)
heatmap(counts_sigma_pvalue0.5, Colv = NA, Rowv = NA, scale="none",main="\n Sigma parameter simulation",breaks=breaks, col=col,xlab="Sigma \n", ylab="Genes")
#legend(7.3, 4.5, leg.txt, fill=col,cex=0.9,title="%Pval > Thresh")
dev.off()

pdf("SupFig_P-value_heatmap_Alphaparameter.pdf", width=6, height=6)
heatmap(counts_alpha_pvalue0.5, Colv = NA, Rowv = NA, scale="none",main="\n Alpha parameter simulation",breaks=breaks, col=col,xlab="Alpha \n", ylab="Genes")
#legend(7.3, 4.5, leg.txt, fill=col,cex=0.9,title="%Pval > Thresh")
dev.off()

pdf("SupFig_P-value_heatmap_Tauparameter.pdf", width=6, height=6)
heatmap(counts_tau_pvalue0.5, Colv = NA, Rowv = NA, scale="none",main="\n P-value above 0.05 for Tau parameter variation",breaks=breaks, col=col,xlab="Tau \n", ylab="Genes")
#legend(7.3, 4.5, leg.txt, fill=col,cex=0.9,title="%Pval > Thresh")
dev.off()
				  
