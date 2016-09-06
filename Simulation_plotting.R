setwd("~/Dropbox/Desiccation\ Selection\ Paper/Dros\ Selection\ Revisions")

rates <- read.csv(file="DesSelExp_TPvsFP.csv", sep=",", quote="\"", header=TRUE)

for(i in levels(as.factor(rates$Chr))){
  templevel <- i
  tempname <- paste("sub", as.character(i), sep="")
  tempsubset <- subset(rates, rates$Chr==i)
  assign(tempname, tempsubset)
}


plot(log10(sub3R$TP_rate/sub3R$FP_rate)~sub3R$Pos, 
     type="l",
     ylab="log10(true positive/false positive)",
     xlab="Position (Mbp)",
     ylim=c(-1.5, 1.5))
abline(h=0, col="red")



D1_counts <- D1_pre_counts[,c("chr", "pos", "MB_counts", "C1_counts", "C2_counts", "C3_counts", "C4_counts", "C5_counts", "D1_counts", "D2_counts", "D3_counts", "D4_counts", "D5_counts")]
