####Generate input files
#cd /nfs/brubeck.bx.psu.edu/scratch6/rahul/Apes/analysis/ddPCR_CN/CAFE_simulation_star5indiv
R

library(reshape2)
library(plyr)

greatA<-read.table(file="../New_Normalized_all_updated",header=T,sep="\t",stringsAsFactors=F, fill=TRUE)
greatA_CN<-cbind(greatA[,c(1,2)],apply(greatA[,-c(1,2)],1,function(x){mean(x,na.rm=T)}))
colnames(greatA_CN)[3]="CopyNumber"

annotation<-read.table(file="../New_great_ape_info",header=T,sep="\t",stringsAsFactors=F)

greatAPE_summary<-cbind(greatA_CN,annotation[match(greatA_CN$IID,annotation$IID),2])
colnames(greatAPE_summary)[4]="Species"

####generate input files for CAFE to test for different combinations.
set.seed(9)
test<-dcast(greatAPE_summary,Gene~IID,value.var="CopyNumber")
id<-match(annotation$IID,colnames(test))
greatAPE_df<-test[,id]
rownames(greatAPE_df)<-test$Gene
greatAPE_df[is.na(greatAPE_df)] <- 0


#Sample summary and order in the dataframe
#Gorilla	Chimp	Borangutan	Bonobo	Sorangutan	Human
# 14		9		7          	7       5			10


###Randomly pick 5 individuals per species
for(i in 1:100){
genefamily=c("BPY2","CDY","DAZ","HSFY","PRY","RBMY","TSPY","VCY","XKRY")
each5_randGA=greatAPE_df[,c(sample(31:37,5),sample(15:23,5),sample(43:52,5),sample(1:14,5),sample(24:30,5),sample(38:42,5))]
output=round(each5_randGA,digits = 2)
output=cbind(genefamily,genefamily,output)
colnames(output)=c("Description","ID","B1","B2","B3","B4","B5","C1","C2","C3","C4","C5","H1","H2","H3","H4","H5","G1","G2","G3","G4","G5","BO1","BO2","BO3","BO4","BO5","SO1","SO2","SO3","SO4","SO5")
write.table(output,file=paste("Random_5each_great_apes_100/GreatApes_CopyNumber_5each_",i,".tab",sep=""), quote=FALSE , sep ="\t",col.names = TRUE, row.names = FALSE)
}


###WE use PPY for final analysis
# for(i in 1:100){
# genefamily=c("BPY2","CDY","DAZ","HSFY","PRY","RBMY","TSPY","VCY","XKRY")
# each5_randGA=greatAPE_df[,c(sample(31:37,4),sample(15:23,4),sample(43:52,4),sample(1:14,4),sample(24:30,4),sample(38:42,4))]
# output=round(each5_randGA,digits = 2)
# output=cbind(genefamily,genefamily,output)
# colnames(output)=c("Description","ID","B1","B2","B3","B4","C1","C2","C3","C4","H1","H2","H3","H4","G1","G2","G3","G4","BO1","BO2","BO3","BO4","SO1","SO2","SO3","SO4")
# write.table(output,file=paste("Random_4each_great_apes_100_PPG/GreatApes_CopyNumber_4each_",i,".tab",sep=""), quote=FALSE , sep ="\t",col.names = TRUE, row.names = FALSE)
# }


####Generate CAFÉ scripts

#R test the tree
library(ape)
test<-read.tree("tree_star5_phylogeny.nwk")
plot(test)


#Generate scripts
python

import os

work_dir="/nfs/brubeck.bx.psu.edu/scratch6/rahul/Apes/analysis/ddPCR_CN/CAFE_simulation_star5indiv/Random_5each_great_apes_100/"
for I in range(1,102):
	sl=open(work_dir+"script_CAFE_each5_specieslevelnode_KYA_1k_"+str(I)+".sh","w")
	sl.write("#!cafe\n")
	sl.write("#version\n")
	sl.write("#date\n")
	sl.write("load -i GreatApes_CopyNumber_5each_"+str(I)+".tab  -t 10 -l reports/log_run_each5_specieslevel_KYA_1K_inter_"+str(I)+".txt -p 0.05\n")
	sl.write("#1k for all internal node\n")
	#sl.write("tree ((((((((B1:288,B2:288):1,B3:289):1,B4:290):1,B5:291):1469,((((C1:548,C2:548):1,C3:549):1,C4:550):1,C5:551):1209):4079,((((H1:130,H2:130):1,H3:131):1,H4:132):1,H5:133):5706):3350,((((G1:148,G2:148):1,G3:149):1,G4:150):1,G5:150):9038):11099.87,(((((BO1:74,BO2:74):1,BO3:75):1,BO4:76):1,BO5:77):342,((((SO1:11,SO2:11):1,SO3:12):1,SO4:13):1,SO5:14):405):19869.87)\n")  ##PPG
	sl.write("tree ((((((((B1:208,B2:208):1,B3:209):1,B4:210):1,B5:211):2222,((((C1:396,C2:396):1,C3:397):1,C4:398):1,C5:399):2034):5641,((((H1:93,H2:93):1,H3:94):1,H4:95):1,H5:96):7978):4633,((((G1:108,G2:108):1,G3:109):1,G4:110):1,G5:111):12596):15346.38,(((((BO1:54,BO2:54):1,BO3:55):1,BO4:56):1,BO5:57):522,((((SO1:8,SO2:8):1,SO3:9):1,SO4:10):1,SO5:11):568):27474.38)\n")
	sl.write("lambda -s -t ((((((((1,1)1,1)1,1)1,1)1,((((1,1)1,1)1,1)1,1)1)1,((((1,1)1,1)1,1)1,1)1)1,((((1,1)1,1)1,1)1,1)1)1,(((((1,1)1,1)1,1)1,1)1,((((1,1)1,1)1,1)1,1)1)1)\n")
	sl.write("report reports/report_each5APE_specieslevel_KYA_1K_inter_"+str(I)+"\n")


####Run CAFÉ

cat run_CAFE_100_job.sh

#!/bin/bash
for i in `seq 1 100`
do
        I=$i
        echo" processing : cafe script_CAFE_each5_specieslevelnode_KYA_1k_${I}.sh"
        cafe script_CAFE_each5_specieslevelnode_KYA_1k_${I}.sh
done



#tail -n9 reports/report_each4APE_specieslevel_KYA_1K_inter_1.cafe |awk '{split($4,a,"\),\("); print $1,a[4],a[8],a[12],a[16]}'


nohup sh run_CAFE_100_job.sh &



# IDs of nodes:((((((((B1<0>,B2<2>)<1>,B3<4>)<3>,B4<6>)<5>,B5<8>)<7>,((((C1<10>,C2<12>)<11>,C3<14>)<13>,C4<16>)<15>,C5<18>)<17>)<9>,((((H1<20>,H2<22>)<21>,H3<24>)<23>,H4<26>)<25>,H5<28>)<27>)<19>,((((G1<30>,G2<32>)<31>,G3<34>)<33>,G4<36>)<35>,G5<38>)<37>)<29>,(((((BO1<40>,BO2<42>)<41>,BO3<44>)<43>,BO4<46>)<45>,BO5<48>)<47>,((((SO1<50>,SO2<52>)<51>,SO3<54>)<53>,SO4<56>)<55>,SO5<58>)<57>)<49>)<39>
# Output format for: ' Average Expansion', 'Expansions', 'No Change', 'Contractions', and 'Branch-specific P-values' = (node ID, node ID): (0,2) (1,4) (3,6) (5,8) (7,17) (10,12) (11,14) (13,16) (15,18) (9,27) (20,22) (21,24) (23,26) (25,28) (19,37) (30,32) (31,34) (33,36) (35,38) (29,49) (40,42) (41,44) (43,46) (45,48) (47,57) (50,52) (51,54) (53,56) (55,58)

#Parse p-values
R

sumGF<-matrix(0,nrow=9,ncol=10)
sumGL<-rep(0,9)
for( i in 1:100){
data<-read.table(file=paste("Random_5each_great_apes_100/reports/report_each5APE_specieslevel_KYA_1K_inter_",i,".cafe",sep=""), sep="\t", skip=11,stringsAsFactors=FALSE)
gf<-data$V3
x=unlist(strsplit(data$V4,"[,]"))
x=gsub("[(]","",x)
x=gsub("[)]","",x)
y<-matrix(x,ncol=58,byrow=T)
y_subset<-as.data.frame(y[,c(9,10,19,20,29,30,39,40,49,50)])
out<-apply(y_subset,1:2,function(x){as.double(x)})
out[is.na(out)]=1
out_binary<-apply(out,1:2,function(x){if(x<=0.005){x=1}else{x=0}})
for(k in seq(length(gf))){if(as.double(gf[k])<=0.05){gf[k]=1}else{gf[k]=0}}
sumGF=sumGF+out_binary
sumGL=sumGL+gf
}


colnames(sumGF)<-seq(1,10)
rownames(sumGF)<-c("BPY2","CDY","DAZ","HSFY","PRY","RBMY","TSPY","VCY","XKRY")

#Supplemetary table
sumGF





