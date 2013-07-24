## Script to ccompay "Significance Analysis of Prognostic Signatures ##
## This script is run after running saps.R and runSAPSonPermutedData.R ##
## to obtain permutation-based p values for a set of SAPS scores ##
## abeck2@bidmc.harvard.edu

rm(list=ls())

# Select home directory
home<-"D:\\Dropbox\\SAPS_Revisions\\RevisedFilesToSubmitToPLOSCB\\December.4.2012_Revised.Final\\R-Script-Workspace-Data_rev\\"

setwd(home)

## Select analysis type
anType="Br"
anType="Ov"

if(anType=="Ov"){
load("Ovary.Ps.OnPermutedData.RData" )
load("OvaryOutput_TradScaled.RData")
}
if(anType=="Br"){
load("Breast.Ps.OnPermutedData.RData")
load("BreastOutput_SubScaled.RData")
}

saps.score.r<-array(NA,dim=dim(p_pure),dimnames=list(rownames(p_pure),colnames(p_pure),unlist(dimnames(p_pure)[3])))

for(i in 1:nrow(saps.score.r)){
   for(j in 1:ncol(saps.score.r)){
      
      if(anType=="Br"){
      saps.score.r[i,j,"Global"]<-max(p_enrich[i,j,"Global"],p_pure[i,j,"Global"],p_rand[i,j,"Global"],na.rm=T)
      saps.score.r[i,j,"ER_H"]<-max(p_enrich[i,j,"ERH"],p_pure[i,j,"ER_H"],p_rand[i,j,"ER_H"],na.rm=T)
      saps.score.r[i,j,"ER_L"]<-max(p_enrich[i,j,"ERL"],p_pure[i,j,"ER_L"],p_rand[i,j,"ER_L"],na.rm=T)
      saps.score.r[i,j,"H2"]<-max(p_enrich[i,j,"H2"],p_pure[i,j,"H2"],p_rand[i,j,"H2"],na.rm=T)
      saps.score.r[i,j,"TN"]<-max(p_enrich[i,j,"TN"],p_pure[i,j,"TN"],p_rand[i,j,"TN"],na.rm=T)
      }
   
      if(anType=="Ov"){
      saps.score.r[i,j,"Global"]<-max(p_enrich[i,j,"Global"],p_pure[i,j,"Global"],p_rand[i,j,"Global"],na.rm=T)
      saps.score.r[i,j,"Angiogenic"]<-max(p_enrich[i,j,"Angio"],p_pure[i,j,"Angiogenic"],p_rand[i,j,"Angiogenic"],na.rm=T)
      saps.score.r[i,j,"Non-angiogenic"]<-max(p_enrich[i,j,"Non.Angio"],p_pure[i,j,"Non-angiogenic"],p_rand[i,j,"Non-angiogenic"],na.rm=T)
      }
   }
   cat(paste("Iteration ", i, " of ", nrow(saps.score.r)),"\n")
  # cat("\\n")
}

if(anType=="Br"){
   saps.score<-matrix(NA,nrow=nrow(allPs),ncol=5,dimnames=list(rownames(allPs),c("Global","ER_H","ER_L","H2","TN")))
   saps.score[,"Global"]<-apply(allPs[,c('P_pure.Global','P_random.Global','P_gsea.Global')],1,max,na.rm=T)
   saps.score[,"ER_H"]<-apply(allPs[,c('P_pure.ER_H','P_random.ER_H','P_gsea.ER_H')],1,max,na.rm=T)
   saps.score[,"ER_L"]<-apply(allPs[,c('P_pure.ER_L','P_random.ER_L','P_gsea.ER_L')],1,max,na.rm=T)
   saps.score[,"TN"]<-apply(allPs[,c('P_pure.TN','P_random.TN','P_gsea.TN')],1,max,na.rm=T)
   saps.score[,"H2"]<-apply(allPs[,c('P_pure.H2','P_random.H2','P_gsea.H2')],1,max,na.rm=T)
   saps.p<-matrix(NA,nrow=nrow(allPs),ncol=5,dimnames=list(rownames(allPs),c("Global","ER_H","ER_L","H2","TN")))
}

if(anType=="Ov"){
   saps.score<-matrix(NA,nrow=nrow(allPs),ncol=3,dimnames=list(rownames(allPs),c("Global","Angiogenic","Non-angiogenic")))
   saps.score[,"Global"]<-apply(allPs[,c('P_pure.Global','P_random.Global','P_gsea.Global')],1,max,na.rm=T)
   saps.score[,"Angiogenic"]<-apply(allPs[,c('P_pure.Angio','P_random.Angio','P_gsea.Angio')],1,max,na.rm=T)
   saps.score[,"Non-angiogenic"]<-apply(allPs[,c('P_pure.Non-angio','P_random.Non-angio','P_gsea.Non-angio')],1,max,na.rm=T) 
   saps.p<-matrix(NA,nrow=nrow(allPs),ncol=3,dimnames=list(rownames(allPs),c("Global","Angiogenic","Non-angiogenic")))
}

geneSeq<-as.numeric(unlist(dimnames(saps.score.r)[1]))

size<-allPs[,"gsSizes"  ]

for(i in 1:nrow(saps.p)){
   for(j in 1:ncol(saps.p)){
      temp.ps<-saps.score.r[geneSeq==min(geneSeq[geneSeq>=size[i]]),,colnames(saps.p)[j]]
      saps.p[i,colnames(saps.p)[j]] = sum(temp.ps<=saps.score[i,colnames(saps.p)[j]])/length(temp.ps)
   }
   cat("Iteration ", i, " of ", nrow(saps.p),"\n")
}
min(saps.p)
saps.p<-replace(saps.p,saps.p==0,1/dim(saps.score.r)[2])
min(saps.p)

cor(saps.p,saps.score,method="sp")

saps.p.adj<-apply(saps.p,2,p.adjust,"BH")
# 
if(anType=="Ov"){
saps.score.adj<-matrix(NA,nrow=nrow(saps.p),ncol=ncol(saps.p),dimnames=dimnames(saps.p))
saps.score.adj[,"Global"]<- apply(allPs.adj[,c("P_pure.Global","P_random.Global"  ,"P_gsea.Global")],1,max,na.rm=T)
saps.score.adj[,"Angiogenic" ]<- apply(allPs.adj[,c("P_pure.Angio","P_random.Angio"  ,"P_gsea.Angio")],1,max,na.rm=T)
saps.score.adj[,"Non-angiogenic" ]<- apply(allPs.adj[,c("P_pure.Non-angio","P_random.Non-angio"  ,"P_gsea.Non-angio")],1,max,na.rm=T)
}
if(anType=="Br"){
saps.score.adj<-matrix(NA,nrow=nrow(saps.p),ncol=ncol(saps.p),dimnames=dimnames(saps.p))
saps.score.adj[,"Global"]<- apply(allPs.adj[,c("P_pure.Global","P_random.Global"  ,"P_gsea.Global")],1,max,na.rm=T)
saps.score.adj[,"ER_H" ]<- apply(allPs.adj[,c("P_pure.ER_H","P_random.ER_H"  ,"P_gsea.ER_H")],1,max,na.rm=T)
saps.score.adj[,"ER_L" ]<- apply(allPs.adj[,c("P_pure.ER_L","P_random.ER_L"  ,"P_gsea.ER_L")],1,max,na.rm=T)
saps.score.adj[,"H2" ]<- apply(allPs.adj[,c("P_pure.H2","P_random.H2"  ,"P_gsea.H2")],1,max,na.rm=T)
saps.score.adj[,"TN" ]<- apply(allPs.adj[,c("P_pure.TN","P_random.TN"  ,"P_gsea.TN")],1,max,na.rm=T)
}

if(anType=="Br"){
load("BreastDirs.RData")
save(dirs,saps.score.adj,allPs,allPs.adj,saps.p,saps.p.adj,saps.score,saps.score.r,size,file="FinalOutput_Breast.RData")
}

if(anType=="Ov"){
   load("OvaryDirs.RData")
   save(dirs,saps.score.adj,allPs,allPs.adj,saps.p,saps.p.adj,saps.score,saps.score.r,size,file="FinalOutput_Ovary.RData")
}

