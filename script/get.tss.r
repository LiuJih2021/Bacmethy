Args<-commandArgs(TRUE)
Args<-as.numeric(Args)
para_pro<-Args[1]
para_cds<-Args[2]
positive<-read.table("positive.gtf",sep = "\t",quote = "")
TSS<-positive$V4
end<-positive$V5
promoter<-TSS-para_pro
end2<-c(0,end)
for(i in 1:length(promoter)){if(promoter[i]<end2[i]){ promoter[i]=end2[i]}}
promoter<-cbind(promoter,TSS+para_cds)
write.table(promoter,"tss.po.bed",sep = "\t", row.names =FALSE, col.names =FALSE,quote = FALSE)

negative<-read.table("negative.gtf",sep = "\t",quote = "")
 TSS<-negative$V4
 end<-negative$V5
 promoter<-end+para_pro
 for(i in 1:(length(promoter)-1)){if(promoter[i]>TSS[i+1])promoter[i]=TSS[i+1]}
promoter<-cbind(end-para_cds,promoter)
 write.table(promoter,"tss.ne.bed",sep = "\t", row.names =FALSE, col.names =FALSE,quote = FALSE)


