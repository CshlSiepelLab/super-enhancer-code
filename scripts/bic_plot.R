library(ggplot2)
library(data.table)
library(gridExtra)
library(plotrix)
setwd("C:/Users/Noah Dukler/OneDrive/Documents/Enhancer_letter")

dat=as.data.table(read.csv("model_bic.csv",nrows=12,header=TRUE))[,1:4,with=FALSE]

rbic2=dat[Error=="Log-normal",list(Model,Error,Rel.bic=BIC[1]-BIC),by="Dataset"]
rbic2=rbind(rbic2[Dataset=="Alpha-globin"],rbic2[Dataset=="Wap"])
rgap=ifelse(rbic2$Rel.bic < -5 , rbic2$Rel.bic +50, rbic2$Rel.bic )

pdf("plots/fancy_bic.pdf")
x<-barplot(height = rgap,width=1,  pch=16, xlab= 'Model', ylab='Relative BIC', 
        yaxt="n",col=c(4,4,4,3,3,3),cex.names=1.4,cex.lab=1.5)
text(cex=1.5, x=x+0.25, y=-16, rbic2$Model, xpd=TRUE, srt=45, pos=2)
xat <- pretty(rgap)
xat <- xat[xat!=-5]
xlab = xat 
xlab[xlab < -5]=xlab[xlab < -5]-50
axis(2,at=xat, labels=xlab,cex.axis=1.5)
axis.break(2,-5,style="slash") 
dev.off()

rbic=dat[,list(Model,Error,Rel.bic=BIC-BIC[1]),by="Dataset"]
rbic[,Model:=factor(Model,levels=c("Additive","Exponential","Logistic"))]

g=ggplot(rbic,aes(x=Model,y=-Rel.bic,fill=Error,color=Error))+
  facet_wrap(~Dataset,scales = "free_y")+
  geom_bar(stat = "identity",position="dodge")+
  theme_bw(base_size = 22)+
  ylab("Relative BIC")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "bottom")
  

pdf("plots/supp_bic_barplot.pdf",height=5,width=7)
plot(g)
dev.off()


dat2=as.data.table(read.csv("model_bic.csv",header=TRUE))[c(3,13:14),1:4,with=FALSE]
rbic2=dat2[,list(Model,Error,Rel.bic=BIC-min(BIC)),by="Dataset"]
rbic2[,Model:=factor(Model,levels=c("Logistic","Logistic-H12","Logistic-H23","Logistic-H12/13"))]


g2=ggplot(rbic2,aes(x=Model,y=-Rel.bic,fill=Error,color=Error))+
  facet_grid(~Dataset)+
  geom_bar(stat = "identity",position="dodge")+
  theme_bw(base_size = 22)+
  ylab("Relative BIC")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "bottom")

pdf("plots/bic_hierarchy_barplot.pdf",height=5,width=7)
plot(g2)
dev.off()
