## Script to accompany "Significance Analysis of Prognostic Signatures" ##
## This script is run to obtain a null distribution of SAPS component p values by running the SAPS procedure on "random" gene sets
## This is typically run after saps.R and before computeSAPS.Permute.PValue.R
## Questions: abeck2@bidmc.harvard.edu
rm(list=ls())
#set home dir
home<-"D:\\Dropbox\\SAPS_Revisions\\RevisedFilesToSubmitToPLOSCB\\December.4.2012_Revised.Final\\R-Script-Workspace-Data_rev\\"
setwd(home)
library(impute)
library(survcomp)
library(genefu)


### Are you doing ovary or breast analysis?
### Run the correct line and comment out the other
anType<-"Br"
anType<-"Ov"

if(anType=="Br"){
   load("Breast.RData")
}

if(anType=="Ov"){
   load("Ovary.RData")
}


### RandSurv Function returns a matrix of log rank p values
### for survival differences, based on clustering data with random gene sets
### of the size specified in "geneSeq"
randSurv<-function(dat,numClusts=2,iter=25,time,event,geneSeq){
  # Matrix to hold log-rank p values
  lrp<-matrix(NA,nrow=iter,ncol=length(geneSeq),dimnames=list(1:iter,geneSeq))
  for(iii in 1:length(geneSeq)){
    for (jjj in 1:iter){
      geneSamp<-sample(colnames(dat),geneSeq[iii])
      dat1<-scale(dat[,geneSamp])
        lab<-kmeans(dat1[,geneSamp],numClusts)$cluster
        survtest<-survdiff(Surv(time,event)~lab)
        lrp[jjj,iii]<-1 - pchisq(survtest$chisq, 1)
     }
       cat("loop", iii, " of ", length(geneSeq), "\n")
       flush.console()
}
     return (lrp)
}

# #############################################
# ######## To compute survival statistics for random gene sets#########
# #############################################
maxGeneSetSize<-250
iter=10000 # # of samplings of random gene sets
geneSeq<-c(5,10,25,50,100,150,200,250) # Smaller sequence to cover molSigDB
#For Ovary
if(anType=="Ov"){
   randMat_global<-randSurv(dat,numClusts=2,iter=iter,time,event,geneSeq)
   randMat_ang<-randSurv(dat[st=="Angiogenic",],numClusts=2,iter=iter,time[st=="Angiogenic"],event[st=="Angiogenic"],geneSeq)
   randMat_nonAng<-randSurv(dat[st=="Non-angiogenic",],numClusts=2,iter=iter,time[st=="Non-angiogenic"],event[st=="Non-angiogenic"],geneSeq)
}
 
 # For Breast
 if(anType=="Br"){
    randMat_global<-randSurv(dat.st,numClusts=2,iter=iter,time,event,geneSeq)
    randMat_er.h<-randSurv(dat.st[st=="ER+/HER2- High Prolif",],numClusts=2,iter=iter,time[st=="ER+/HER2- High Prolif"],event[st=="ER+/HER2- High Prolif"],geneSeq)
    randMat_er.l<-randSurv(dat.st[st=="ER+/HER2- Low Prolif",],numClusts=2,iter=iter,time[st=="ER+/HER2- Low Prolif"],event[st=="ER+/HER2- Low Prolif"],geneSeq)
    randMat_h2<-randSurv(dat.st[st=="HER2+",],numClusts=2,iter=iter,time[st=="HER2+"],event[st=="HER2+"],geneSeq)
    randMat_tn<-randSurv(dat.st[st=="ER-/HER2-",],numClusts=2,iter=iter,time[st=="ER-/HER2-"],event[st=="ER-/HER2-"],geneSeq)
 }


########################
#### P_pure prognosis ############
########################
genes<-colnames(dat.st)
randSets<-vector("list",length=length(geneSeq))
for(iii in 1:length(geneSeq)){
   for (jjj in 1:iter){
      randSets[[iii]][[jjj]]<-sample(colnames(dat),geneSeq[iii])
   }
}

# For breast
if(anType=="Br"){
   p_pure<- array(NA,dim=c(length(geneSeq),iter,6),dimnames=list(geneSeq,1:iter,c("Size","Global","ER_H","ER_L","H2","TN")))
   dat.x<-dat.st # For Breast
}

# For ovary
if(anType=="Ov"){
   p_pure<- array(NA,dim=c(length(geneSeq),iter,4),dimnames=list(geneSeq,1:iter,c("Size","Global","Angiogenic","Non-angiogenic")))
   dat.x<-dat # For ovary
}


for(i in 1:dim(p_pure)[1]){
   for(j in 1:dim(p_pure)[2]){
   comGenes<-intersect(genes,unique(as.character(randSets[[i]][[j]])))
   if(length(comGenes)){
      p_pure[i,j,"Size"]<-length(comGenes)
         
      # Global
      dats<-scale(dat.x[,is.element(genes,comGenes)])
      lab<-kmeans(dats,2)$cluster
      survtest<-survdiff(Surv(time,event)~lab)
      p_pure[i,j,"Global"]<-   1 - pchisq(survtest$chisq, 1)
      
      
      # For ovary
      if(anType=="Ov"){
         dats<-scale(dat.x[st=="Angiogenic",is.element(genes,comGenes)])
         lab<-kmeans(dats,2)$cluster
         survtest<-survdiff(Surv(time[st=="Angiogenic"],event[st=="Angiogenic"])~lab)
         p_pure[i,j,"Angiogenic"]<-   1 - pchisq(survtest$chisq, 1)
         
         dats<-scale(dat.x[st=="Non-angiogenic",is.element(genes,comGenes)])
         lab<-kmeans(dats,2)$cluster
         survtest<-survdiff(Surv(time[st=="Non-angiogenic"],event[st=="Non-angiogenic"])~lab)
         p_pure[i,j,"Non-angiogenic"]<-   1 - pchisq(survtest$chisq, 1)
      }
      
      # For Breast
      if(anType=="Br"){
         dats<-scale(dat.x[st=="ER+/HER2- High Prolif",is.element(genes,comGenes)])
         lab<-kmeans(dats,2)$cluster
         survtest<-survdiff(Surv(time[st=="ER+/HER2- High Prolif"],event[st=="ER+/HER2- High Prolif"])~lab)
         p_pure[i,j,"ER_H"]<-   1 - pchisq(survtest$chisq, 1)
         
         dats<-scale(dat.x[st=="ER+/HER2- Low Prolif",is.element(genes,comGenes)])
         lab<-kmeans(dats,2)$cluster
         survtest<-survdiff(Surv(time[st=="ER+/HER2- Low Prolif"],event[st=="ER+/HER2- Low Prolif"])~lab)
         p_pure[i,j,"ER_L"]<-   1 - pchisq(survtest$chisq, 1)
         
         dats<-scale(dat.x[st=="HER2+",is.element(genes,comGenes)])
         lab<-kmeans(dats,2)$cluster
         survtest<-survdiff(Surv(time[st=="HER2+"],event[st=="HER2+"])~lab)
         p_pure[i,j,"H2"]<-   1 - pchisq(survtest$chisq, 1)
         
         dats<-scale(dat.x[st=="ER-/HER2-",is.element(genes,comGenes)])
         lab<-kmeans(dats,2)$cluster
         survtest<-survdiff(Surv(time[st=="ER-/HER2-"],event[st=="ER-/HER2-"])~lab)
         p_pure[i,j,"TN"]<-   1 - pchisq(survtest$chisq, 1)
      
      }
   }
      cat("Inside loop", j, " of ", dim(p_pure)[2], "\n")
      flush.console()
   }
   
   cat("Outside loop", i, " of ", dim(p_pure)[1], "\n")
   flush.console()
}

########################
#### P_random prognosis 
########################
### Now, we want to compute P random prognosis
p_rand<- array(NA,dim=dim(p_pure),dimnames=dimnames(p_pure))
p_rand[,,"Size"]<-p_pure[,,"Size"]

for(i in 1:dim(p_rand)[1]){
   for(j in 1:dim(p_rand)[2]){
      if(!is.na(p_rand[i,j,"Size"])){
         if(p_rand[i,j,"Size"] <= max(geneSeq)){
         col.temp<-as.character(min(geneSeq[geneSeq>=p_rand[i,j,"Size"]]))
      if(!is.na(col.temp)){
         # Compute non-parametric raw p value for each gene set
         if(anType == "Ov"){
            p_rand[i,j,"Global"]<-sum(randMat_global[,col.temp]<=p_pure[i,j,"Global"])/nrow(randMat_global)
            p_rand[i,j,"Angiogenic" ]<-sum(randMat_ang[,col.temp]<=p_pure[i,j,"Angiogenic"])/nrow(randMat_ang)
            p_rand[i,j,"Non-angiogenic"]<-sum(randMat_nonAng[,col.temp]<=p_pure[i,j,"Non-angiogenic"])/nrow(randMat_nonAng)
         }
         if(anType== "Br"){
            p_rand[i,j,"Global"]<-sum(randMat_global[,col.temp]<=p_pure[i,j,"Global"])/nrow(randMat_global)
            p_rand[i,j,"ER_H"]<-sum(randMat_er.h[,col.temp]<=p_pure[i,j,"ER_H"])/nrow(randMat_er.h)
            p_rand[i,j,"ER_L"]<-sum(randMat_er.l[,col.temp]<=p_pure[i,j,"ER_L"])/nrow(randMat_er.l)
            p_rand[i,j,"H2"]<-sum(randMat_h2[,col.temp]<=p_pure[i,j,"H2"])/nrow(randMat_h2)
            p_rand[i,j,"TN"]<-sum(randMat_tn[,col.temp]<=p_pure[i,j,"TN"])/nrow(randMat_tn)
               }
            }
         }
      }
   cat("Inside loop", j, " of ", dim(p_rand)[2], "\n")
   flush.console()
   }
   cat("Outside loop", i, " of ", dim(p_rand)[1], "\n")
   flush.console()
}
p_rand<-replace(p_rand,p_rand==0,1/iter)

################################
################################
### Compute P_enrichment########
################################
################################

## 1. Generate pre-ranked lists, based on concordance index
if(anType=="Ov"){
   data.s<-dat
   rr.1<- t(apply(data.s[, ,drop=F],2, function(x, y, z) {
      tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE)
      return(c("cindex"=tt$c.index, "z"=(tt$c.index - 0.5)/tt$se))
   }, y=time, z=event))

   data.s<-dat[st=="Angiogenic",]
   rr.2<- t(apply(data.s[, ,drop=F],2, function(x, y, z) {
      tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE)
      return(c("cindex"=tt$c.index, "z"=(tt$c.index - 0.5)/tt$se))
   }, y=time[st=="Angiogenic"], z=event[st=="Angiogenic"]))


   data.s<-dat[st=="Non-angiogenic",]
   rr.3<- t(apply(data.s[, ,drop=F],2, function(x, y, z) {
      tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE)
      return(c("cindex"=tt$c.index, "z"=(tt$c.index - 0.5)/tt$se))
   }, y=time[st=="Non-angiogenic"], z=event[st=="Non-angiogenic"]))
}
}


if(anType=="Br"){
   data.s=dat.st
   rr.1<- t(apply(data.s[, ,drop=F],2, function(x, y, z) {
      tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE)
      return(c("cindex"=tt$c.index, "z"=(tt$c.index - 0.5)/tt$se))
   }, y=time, z=event))
   
   data.s<-dat.st[st=="ER+/HER2- High Prolif",]
   rr.2<- t(apply(data.s[, ,drop=F],2, function(x, y, z) {
      tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE)
      return(c("cindex"=tt$c.index, "z"=(tt$c.index - 0.5)/tt$se))
   }, y=time[st=="ER+/HER2- High Prolif"], z=event[st=="ER+/HER2- High Prolif"]))


   data.s<-dat.st[st=="ER+/HER2- Low Prolif",]
   rr.3<- t(apply(data.s[, ,drop=F],2, function(x, y, z) {
      tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE)
      return(c("cindex"=tt$c.index, "z"=(tt$c.index - 0.5)/tt$se))
   }, y=time[st=="ER+/HER2- Low Prolif"], z=event[st=="ER+/HER2- Low Prolif"]))

   data.s<-dat.st[st=="HER2+",]
   rr.4<- t(apply(data.s[, ,drop=F],2, function(x, y, z) {
      tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE)
      return(c("cindex"=tt$c.index, "z"=(tt$c.index - 0.5)/tt$se))
   }, y=time[st=="HER2+"], z=event[st=="HER2+"]))

   data.s<-dat.st[st=="ER-/HER2-",]
   rr.5<- t(apply(data.s[, ,drop=F],2, function(x, y, z) {
      tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE)
      return(c("cindex"=tt$c.index, "z"=(tt$c.index - 0.5)/tt$se))
   }, y=time[st=="ER-/HER2-"], z=event[st=="ER-/HER2-"]))
}

## 2. Save .RNK files for GSEA
rnk.genes<-as.character(rownames(rr.1))
rnk.dir<-"XXXX" ## Set the directory to save files
if(!file.exists(rnk.dir)) {dir.create(rnk.dir)}
setwd(rnk.dir)

if(anType=="Ov"){
   write.table(cbind(rnk.genes,rr.1[,2]), file="ovary_nss_global.rnk", col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
   write.table(cbind(rnk.genes,rr.2[,2]), file="ovary_nss_angio.rnk", col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
   write.table(cbind(rnk.genes,rr.3[,2]), file="ovary_nss_non.angio.rnk", col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
}

if(anType=="Br"){
   write.table(cbind(rnk.genes,rr.1[,2]), file="breast_global.rnk", col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
   write.table(cbind(rnk.genes,rr.2[,2]), file="breast_erh.rnk", col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
   write.table(cbind(rnk.genes,rr.3[,2]), file="breast_erl.rnk", col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
   write.table(cbind(rnk.genes,rr.4[,2]), file="breast_h2.rnk", col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
   write.table(cbind(rnk.genes,rr.5[,2]), file="breast_tn.rnk", col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
}


## Now, use the rnk files in a pre-ranked GSEA using the JAVA gsea application
## This application can be called from within R
# Example (change file locations to match your system):

## We will have to generate a huge file with all our random gene sets
## Our random gene sets are in randSets, a list of length 8, with 100 gene sets each
if(1){

randForGSEA<-matrix(NA,nrow=max(geneSeq),ncol=length(geneSeq)*iter)
for(i in 1:(length(geneSeq))){
   for(j in 1:iter){
     randForGSEA[1:geneSeq[i],(i-1)*iter + j]<-randSets[[i]][[j]]
   }
}

cols.temp<-vector(length=dim(randForGSEA)[2])
for(i in 1:length(geneSeq)){
   for(j in 1:iter){
      cols.temp[(i-1)*iter + j]=paste(geneSeq[i],j,sep=".")
   }
}
colnames(randForGSEA)<-cols.temp
randForGSEA<-rbind(cols.temp,randForGSEA)
write.table(randForGSEA,"randForGSEA.10K.gmx",quote=F,sep="\t",row.names=F,na="")
}

# # Split up GMX Files into 8 subfiles so the GSEA doesn't crash
setwd(home)
tempGMX<-read.table("randForGSEA.10K.gmx",sep="\t",header=T)
dim(tempGMX)
tempGMX[1:10,1:10]
for(i in 1:length(geneSeq)){
   write.table(tempGMX[,(((i-1)*iter)+1):(i*iter)],file=paste("randForGSEA.",geneSeq[i],".gmx",sep=""),sep="\t",quote=F,row.names=F)
   cat(paste("iteration ", i))
}

## Now, we are ready to run GSEA in Java
# GLOBAL and 4 SUBTYPES on each of 8 .gmx files
## This is an example of how the analysis was run.
## To re-use the code,directories must be changed to point to the appropriate
## files on your directory system
##Breast Global
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.5.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_global.rnk -scoring_scheme weighted -rpt_label Breast_Global_5_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.10.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_global.rnk -scoring_scheme weighted -rpt_label Breast_Global_10_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.25.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_global.rnk -scoring_scheme weighted -rpt_label Breast_Global_25_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.50.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_global.rnk -scoring_scheme weighted -rpt_label Breast_Global_50_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.100.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_global.rnk -scoring_scheme weighted -rpt_label Breast_Global_100_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.150.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_global.rnk -scoring_scheme weighted -rpt_label Breast_Global_150_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.200.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_global.rnk -scoring_scheme weighted -rpt_label Breast_Global_200_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.250.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_global.rnk -scoring_scheme weighted -rpt_label Breast_Global_250_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)

# ER H
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.5.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_erh.rnk -scoring_scheme weighted -rpt_label Breast_ERH_5_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.10.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_erh.rnk -scoring_scheme weighted -rpt_label Breast_ERH_10_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.25.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_erh.rnk -scoring_scheme weighted -rpt_label Breast_ERH_25_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.50.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_erh.rnk -scoring_scheme weighted -rpt_label Breast_ERH_50_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.100.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_erh.rnk -scoring_scheme weighted -rpt_label Breast_ERH_100_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.150.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_erh.rnk -scoring_scheme weighted -rpt_label Breast_ERH_150_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.200.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_erh.rnk -scoring_scheme weighted -rpt_label Breast_ERH_200_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.250.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_erh.rnk -scoring_scheme weighted -rpt_label Breast_ERH_250_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)

# ER L
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.5.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_erl.rnk -scoring_scheme weighted -rpt_label Breast_ERL_5_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.10.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_erl.rnk -scoring_scheme weighted -rpt_label Breast_ERL_10_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.25.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_erl.rnk -scoring_scheme weighted -rpt_label Breast_ERL_25_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.50.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_erl.rnk -scoring_scheme weighted -rpt_label Breast_ERL_50_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.100.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_erl.rnk -scoring_scheme weighted -rpt_label Breast_ERL_100_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.150.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_erl.rnk -scoring_scheme weighted -rpt_label Breast_ERL_150_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.200.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_erl.rnk -scoring_scheme weighted -rpt_label Breast_ERL_200_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.250.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_erl.rnk -scoring_scheme weighted -rpt_label Breast_ERL_250_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)

# H2
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.5.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_h2.rnk -scoring_scheme weighted -rpt_label Breast_H2_5_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.10.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_h2.rnk -scoring_scheme weighted -rpt_label Breast_H2_10_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.25.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_h2.rnk -scoring_scheme weighted -rpt_label Breast_H2_25_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.50.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_h2.rnk -scoring_scheme weighted -rpt_label Breast_H2_50_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.100.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_h2.rnk -scoring_scheme weighted -rpt_label Breast_H2_100_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.150.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_h2.rnk -scoring_scheme weighted -rpt_label Breast_H2_150_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.200.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_h2.rnk -scoring_scheme weighted -rpt_label Breast_H2_200_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.250.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_h2.rnk -scoring_scheme weighted -rpt_label Breast_H2_250_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)

# TN
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.5.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_tn.rnk -scoring_scheme weighted -rpt_label Breast_TN_5_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.10.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_tn.rnk -scoring_scheme weighted -rpt_label Breast_TN_10_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.25.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_tn.rnk -scoring_scheme weighted -rpt_label Breast_TN_25_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.50.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_tn.rnk -scoring_scheme weighted -rpt_label Breast_TN_50_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.100.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_tn.rnk -scoring_scheme weighted -rpt_label Breast_TN_100_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.150.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_tn.rnk -scoring_scheme weighted -rpt_label Breast_TN_150_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.200.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_tn.rnk -scoring_scheme weighted -rpt_label Breast_TN_200_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.250.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/5.30.12.breast_subtype.specific_rank.rnkDir/breast_tn.rnk -scoring_scheme weighted -rpt_label Breast_TN_250_SS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)


##Global Ovary
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.5.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/6.3.12_Ovary_NonSS_RankDir/ovary_nss_global.rnk -scoring_scheme weighted -rpt_label Ovary_Global_5_NSS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.10.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/6.3.12_Ovary_NonSS_RankDir/ovary_nss_global.rnk -scoring_scheme weighted -rpt_label Ovary_Global_10_NSS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.25.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/6.3.12_Ovary_NonSS_RankDir/ovary_nss_global.rnk -scoring_scheme weighted -rpt_label Ovary_Global_25_NSS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.50.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/6.3.12_Ovary_NonSS_RankDir/ovary_nss_global.rnk -scoring_scheme weighted -rpt_label Ovary_Global_50_NSS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.100.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/6.3.12_Ovary_NonSS_RankDir/ovary_nss_global.rnk -scoring_scheme weighted -rpt_label Ovary_Global_100_NSS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.150.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/6.3.12_Ovary_NonSS_RankDir/ovary_nss_global.rnk -scoring_scheme weighted -rpt_label Ovary_Global_150_NSS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.200.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/6.3.12_Ovary_NonSS_RankDir/ovary_nss_global.rnk -scoring_scheme weighted -rpt_label Ovary_Global_200_NSS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.250.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/6.3.12_Ovary_NonSS_RankDir/ovary_nss_global.rnk -scoring_scheme weighted -rpt_label Ovary_Global_250_NSS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)

## Angio
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.5.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/6.3.12_Ovary_NonSS_RankDir/ovary_nss_angio.rnk -scoring_scheme weighted -rpt_label Ovary_Angio_5_NSS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.10.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/6.3.12_Ovary_NonSS_RankDir/ovary_nss_angio.rnk -scoring_scheme weighted -rpt_label Ovary_Angio_10_NSS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.25.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/6.3.12_Ovary_NonSS_RankDir/ovary_nss_angio.rnk -scoring_scheme weighted -rpt_label Ovary_Angio_25_NSS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.50.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/6.3.12_Ovary_NonSS_RankDir/ovary_nss_angio.rnk -scoring_scheme weighted -rpt_label Ovary_Angio_50_NSS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.100.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/6.3.12_Ovary_NonSS_RankDir/ovary_nss_angio.rnk -scoring_scheme weighted -rpt_label Ovary_Angio_100_NSS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.150.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/6.3.12_Ovary_NonSS_RankDir/ovary_nss_angio.rnk -scoring_scheme weighted -rpt_label Ovary_Angio_150_NSS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.200.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/6.3.12_Ovary_NonSS_RankDir/ovary_nss_angio.rnk -scoring_scheme weighted -rpt_label Ovary_Angio_200_NSS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.250.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/6.3.12_Ovary_NonSS_RankDir/ovary_nss_angio.rnk -scoring_scheme weighted -rpt_label Ovary_Angio_250_NSS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)


##Non.Angio
## Angio
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.5.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/6.3.12_Ovary_NonSS_RankDir/ovary_nss_non.angio.rnk -scoring_scheme weighted -rpt_label Ovary_Non.Angio_5_NSS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.10.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/6.3.12_Ovary_NonSS_RankDir/ovary_nss_non.angio.rnk -scoring_scheme weighted -rpt_label Ovary_Non.Angio_10_NSS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.25.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/6.3.12_Ovary_NonSS_RankDir/ovary_nss_non.angio.rnk -scoring_scheme weighted -rpt_label Ovary_Non.Angio_25_NSS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.50.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/6.3.12_Ovary_NonSS_RankDir/ovary_nss_non.angio.rnk -scoring_scheme weighted -rpt_label Ovary_Non.Angio_50_NSS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.100.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/6.3.12_Ovary_NonSS_RankDir/ovary_nss_non.angio.rnk -scoring_scheme weighted -rpt_label Ovary_Non.Angio_100_NSS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.150.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/6.3.12_Ovary_NonSS_RankDir/ovary_nss_non.angio.rnk -scoring_scheme weighted -rpt_label Ovary_Non.Angio_150_NSS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.200.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/6.3.12_Ovary_NonSS_RankDir/ovary_nss_non.angio.rnk -scoring_scheme weighted -rpt_label Ovary_Non.Angio_200_NSS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)
gsea1<-"java -Xmx8g -cp D:\\Dropbox\\random_pathway_bc\\gsea2-2.07.jar xtools.gsea.GseaPreranked -gmx D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\randForGSEA.250.gmx  -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk D:/Dropbox/random_pathway_bc/6.3.12_Ovary_NonSS_RankDir/ovary_nss_non.angio.rnk -scoring_scheme weighted -rpt_label Ovary_Non.Angio_250_NSS_Random -include_only_symbols true -make_sets true -plot_top_x 7000 -rnd_seed  -set_max 250 -set_min 1 -zip_report false -out D:\\Dropbox\\random_pathway_bc\\FilesToSubmitToPCB\\RScriptWorkspaceDataRevisions\\RandomGSEAv4 - Revisions -gui false"
system(gsea1)


############
##############
##########
##########
#####
### 


## Reading in GSEA results from these analyses
iter=10000 # overall iterations for SAPS
iter.gsea<-1000 # GSEA iterations


# In p_enrich, we'll put in raw p values from GSEA for all of the 10K gene sets
# of each size
if(anType=="Br"){
   p_enrich<-array(NA,dim=c(length(geneSeq),iter,6),dimnames=list(geneSeq,1:iter,c("Size","Global","ERH","ERL","H2","TN")))
}
if(anType=="Ov"){
   p_enrich<-array(NA,dim=dim(p_pure),dimnames=dimnames(p_pure))
   dimnames(p_enrich)[3]<- list(c("Size","Global","Angio","Non.Angio"))
}
p_dir=p_enrich
## bigDir is an output directory from GSEA containing all sub-folders for breast or ovary
if(anType=="Br"){
   bigDir<-"Breast GSEA output directory on random data" 
}
if(anType=="Ov"){
   # bigDir<-"Ovary GSEA output directory on random data"
}
setwd(bigDir)
dir1<-dir()

## This loop extracts the P_Enrich from the GSEA out put and puts it into the P_Enrich array
for(i in 1:length(dir1)){
   temp1<-strsplit(dir1,"_",fixed=T)   
   st.temp<-unlist(lapply(temp1,function(x)(x[2])))
   size.temp<-unlist(lapply(temp1,function(x)(x[3])))  
   setwd(paste(bigDir,dir1[i],sep="/"))
   upFile<-dir(pattern="gsea_report_for_na_pos_")[2]
   downFile<-dir(pattern="gsea_report_for_na_neg_")[2]
   upF<-read.table(upFile,sep="\t",header=T)
   downF<-read.table(downFile,sep="\t",header=T)
   tempTab<-rbind(upF[,c("NAME","NOM.p.val"   )],downF[,c("NAME","NOM.p.val"   )])
   tempTab[,2]<-replace(tempTab[,2],tempTab[,2]==0,1/iter.gsea)
   tempTab<-cbind(tempTab,c(rep(1,nrow(upF)),rep(-1,nrow(downF))))
   head(tempTab)
   gs.name<-as.numeric(unlist(lapply(strsplit(as.character(tempTab[,"NAME"]),"\\."),function(x)(x[2]))))
   tempTab<-tempTab[order(gs.name),]
   head(tempTab)
   p_enrich[size.temp[i],,st.temp[i]]<-tempTab[,"NOM.p.val"]
   p_dir[size.temp[i],,st.temp[i]]<-tempTab[,3]
   cat(paste("Iteration ",i, "/",length(dir1)))
}

### Save the 3 P Values for the computeSAPS.Permute.PValue.R Script
if(anType=="Ov"){
   save(p_pure,p_rand,p_enrich,file="Ovary.Ps.OnPermutedData.RData")
}
if(anType=="Br"){
   save(p_pure,p_rand,p_enrich,file="Breast.Ps.OnPermutedData.RData")
}

#######################
########################
#########################
##### Ready to use these p values to computeSAPS.Permute.PValue.R
######################
#############################
####################
###############
######
###
#