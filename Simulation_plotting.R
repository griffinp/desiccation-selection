#! usr/bin/Rscript

###########################################################################################
# 
#                      R script written by Alex Fournier-Level 
#               supporting the results from Griffin et al. (submitted)
#                 freely available and modifiable under GNU licence
# requires -fields
# designed to reproduce the analysis of the result from the Simul_SelectionRegime.py script
###########################################################################################

library(fields)

#Replace by the name of the file were the simulation outputs were printed to:
Data=read.csv("OUT_NormCov.csv",header=F)
Data=Data[order(Data[,1],Data[,2]),]

#par(mfrow=c(1,2))
Data_sub=Data[Data[,6]>1,]
TP=Data_sub[,7]#/Data_sub[,6]
FP_St=Data_sub[,11]#/Data_sub[,12]
FP_Lg=Data_sub[,15]#/Data_sub[,14]

################################################################################################
#                Correlation btw True and False positives
################################################################################################

Table=matrix(NA,ncol=7,nrow=5)
Table[,1]=c("2L","2R","3L","3R","All")
for (X in Table[,1]) {
	if (X=="All") {
	index=c(1:nrow(Data_sub))
	} else index=c(Data_sub[,1]==X)
	Table[X==Table[,1],2]=mean(TP[index])
	Table[X==Table[,1],3]=sd(TP[index])
	Table[X==Table[,1],4]=mean(FP_St[index])
	Table[X==Table[,1],5]=sd(FP_St[index])
	Table[X==Table[,1],6]=mean(FP_Lg[index])
	Table[X==Table[,1],7]=sd(FP_Lg[index])
}

Table_St=xtabs(~TP+FP_St)
SUM_Table=sum(Table_St)
Table_St=log(1+Table_St)
Table_St=Table_St[1:15,1:15]
Table_Lg=xtabs(~TP+FP_Lg)
Table_Lg=log(1+Table_Lg)
Table_Lg=Table_Lg[1:15,1:15]

Table_all=unlist(list(Table_St,Table_Lg))
MAX_Table=max(Table_all)
MIN_Table=min(Table_all)
Yellow2Red=colorRampPalette(c("white","orange","red"))
pdf(file="HeatMap_TPvsFP.pdf",width=10,height=5)
par(mfrow=c(1,2),mar=c(5,4.5,4,6))

image(Table_St,col=Yellow2Red(10),xaxt='n',yaxt='n',xlab="Number of True Associations in Causal Gene",ylab="Number of False Positives in Linked Gene",main="Linked Gene ~10Kbp from Causal Gene",border=F)
axis(1,at=seq(0,1,length.out=8),label=round(seq(min(TP),15,length.out=8)))
axis(2,at=seq(0,1,length.out=8),label=round(seq(min(FP_St),15,length.out=8)))

image(Table_Lg,col=Yellow2Red(10),xaxt='n',yaxt='n',xlab="Number of True Associations in Causal Gene",ylab="Number of False Positives in Linked Gene",main="Linked Gene ~1Mbp from Causal Gene",border=F)
axis(1,at=seq(0,1,length.out=8),label=round(seq(min(TP),15,length.out=8)))
axis(2,at=seq(0,1,length.out=8),label=round(seq(min(FP_St),15,length.out=8)))
image.plot(exp(Table_Lg), col=Yellow2Red(10),zlim=c(exp(MIN_Table)/SUM_Table,exp(MAX_Table)/SUM_Table),lab.breaks=signif((exp(seq(MIN_Table,MAX_Table,length.out=11))-1)/SUM_Table,1),legend.only=T)
dev.off()

dev.new()
plot(TP~Data_sub[,4],xlab="Initial frequency of causal allele",ylab="True positive rate")


################################################################################################
#                FalseNeg/FalsePos over genomic windows
################################################################################################

slideFunct=function(data, coord, step){
	data=as.numeric(as.character(data))
	position=seq(from=1, to=(max(coord,na.rm=T)-step/2), by=step)
	result=c()
	for(i in position) {
		end=which(coord>(i+step/2))[1]
		if (end==1) {
		prop=NA
		result=c(result,prop)
		} else {
			chunk=data[1:(end-1)]
			size=length(chunk)
			prop=sum(chunk>=1)/size
			result=c(result,prop)
			total=length(coord)
			data=data[end:total]
			coord=coord[end:total]
		}
	}
	return(cbind(position,result))
}

slideCount=function(data, coord, step){
	data=as.numeric(as.character(data))
	position=seq(from=1, to=(max(coord,na.rm=T)-step/2), by=step)
	result=c()
	for(i in position) {
		end=which(coord>(i+step/2))[1]
		if (end==1) {
		prop=NA
		result=c(result,prop)
		} else {
			chunk=data[1:(end-1)]
			size=length(chunk)
			#prop=sum(chunk>=1)/size
			prop=mean(chunk)
			result=c(result,prop)
			total=length(coord)
			data=data[end:total]
			coord=coord[end:total]
		}
	}
	return(cbind(position,result))
}


#dev.new()
png("DesSelExp_TPvsFP.png",height=800,width=1200)
par(mfrow=c(2,2))
Data_compiled=matrix(NA,nrow=0,ncol=4)
for (Chr in c("2L","2R","3L","3R")) {
	Data_sub=Data[Data$V1==Chr,]
	TPos=as.numeric(Data_sub$V7>1)
	TPos_window=slideFunct(TPos, Data_sub[,2], 50000)
	FPos=as.numeric(Data_sub$V11>1|Data_sub$V15>1)
	FPos_window=slideFunct(FPos, Data_sub[,8], 50000)
	GeneDense_window=slideCount(Data_sub[,6], Data_sub[,2], 50000)
	Data_compiled=rbind(Data_compiled,cbind(rep(Chr,nrow(TPos_window)),TPos_window,FPos_window[,2]))
	print(Chr)
	print(summary(TPos_window))
	#plot(FPos_window[,2]~GeneDense_window[,2],main=Chr)
	plot(log(TPos_window[,2]/FPos_window[,2])~TPos_window[,1],ylim=c(-4,4),main=Chr,xlab="Position in Mbp",ylab="Rate",type="line",col=3,xaxt='n')
	abline(h=0,col=2,lty=2)
	#lines(FPos_window,main=Chr,col=2)
	axis(1,at=seq(min(TPos_window[,1]),max(TPos_window[,1]),by=5000000),label=round(seq(min(TPos_window[,1]),max(TPos_window[,1]),by=5000000)/1000000))
}
dev.off()

colnames(Data_compiled)=c("Chr","Pos","TP_rate","FP_rate")
write.csv(Data_compiled,file="DesSelExp_TPvsFP.csv",row.names=F)

################################################################################################
#                Association Detection as function of the freq
################################################################################################

TPos=as.numeric(Data$V7>1)
FPos_St=as.numeric(Data$V11>1)
FPos_Lg=as.numeric(Data$V15>1)
Freq=Data$V4

#dev.new()
png("Proba_TPvsFP_h01.png",height=500,width=1200)
par(mfrow=c(1,3),mar=c(5, 5, 4, 2))
model=glm(TPos~Freq,family=binomial(link='logit'))
summary(model)
plot(Freq,TPos+rnorm(length(TPos),0,0.01),xlab="Frequency",ylab="Probability of detecting an Association",main="True Association",col=grey(.6),cex.axis=2,cex.main=2,cex.lab=2) # plot with body size on x-axis 
curve(predict(model,data.frame(Freq=x),type="resp"),add=TRUE,col=2) # draws a curve based on prediction from logistic regression model

model=glm(FPos_St~Freq,family=binomial(link='logit'))
summary(model)
plot(Freq,FPos_St+rnorm(length(FPos_St),0,0.01),xlab="Frequency",ylab="Probability of detecting an Association",main="False Association 10Kbp apart",col=grey(.6),cex.axis=2,cex.main=2,cex.lab=2) # plot with body size on x-axis 
curve(predict(model,data.frame(Freq=x),type="resp"),add=TRUE,col=2) # draws a curve based on prediction from logistic regression model

model=glm(FPos_Lg~Freq,family=binomial(link='logit'))
summary(model)
plot(Freq,FPos_Lg+rnorm(length(FPos_Lg),0,0.01),xlab="Frequency",ylab="Probability of detecting an Association",main="False Association 1Mbp apart",col=grey(.6),cex.axis=2,cex.main=2,cex.lab=2) # plot with body size on x-axis 
curve(predict(model,data.frame(Freq=x),type="resp"),add=TRUE,col=2) # draws a curve based on prediction from logistic regression model

dev.off()



