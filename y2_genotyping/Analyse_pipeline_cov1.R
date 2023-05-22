#####
#Date : 28.01.2023
# Y2 helitron 
#####

#elimination des warning
options(warn=-1)

#recuperation des arguments
args <- commandArgs(trailingOnly = TRUE)

#TE family
te<-as.character(args[1])

inputFile<-try(system("ls *txt",intern=TRUE))
files <- unlist(strsplit(inputFile,split=" "))

length(files)
M<-data.frame()

all<-paste("all_insertion_",te,".names",sep="")
M<-read.table(all,sep="\n",h=F)

for (i in 1:length(files)){
 res<-read.table(files[i],sep="\n",h=F)
 tmp<-M$V1 %in% res$V1
 M<-cbind(M,tmp)
 }
colnames(M)<-c("insertion",files)

#pos cov1
pos<-paste("all_position_cov1_",te,".names",sep="")
POS<-read.table(pos,sep="\n",h=F)


#Matrice final
library(data.table)

M2<-M[which(M$insertion %in% POS$V1),]
#M3<-transpose(M2)
write.table(M2,file="matrice_final.csv",sep="\t",row.names=F,col.names=T,quote=F)
#write.table(M3,file="matrice_final.csv",sep="\t",row.names=F,col.names=T,quote=F)


