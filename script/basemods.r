library(ggplot2)
library(ggthemes)
Args<-commandArgs(TRUE)
file<-Args[1]
motif<-read.table(file,sep = "\t")
colnames(motif)<-c("coverage","fraction","identificationQV")

pdf("plot.pdf")
p<-ggplot(motif,aes(x=coverage,y=identificationQV,colour=fraction))+
  geom_point()+theme_hc()+theme(legend.position = "right")+geom_hline(aes(yintercept=40),linetype=5,col="red")+ scale_colour_gradient2(low = 'blue', high = 'red', midpoint = 0.75)
print(p)
dev.off()

motif$per<-as.factor(ifelse(motif$fraction>0.59, ifelse(motif$fraction>0.74,ifelse(motif$fraction>0.89,'90-100%','75-90%'),'60-75%'),'0-60%'))
data2<-as.data.frame(table(motif$per)) 
colnames(data2)<-c("fraction","Frequence")
pdf("fraction.plot.pdf")
p1=ggplot(data2,aes(x =fraction ,y = Frequence))+
  geom_bar(stat = 'identity',aes(fill = Frequence))+theme_hc()
print(p1)
dev.off()
