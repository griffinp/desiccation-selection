# source("http://bioconductor.org/biocLite.R")
# biocLite("qvalue")

library(RSQLite)


Sample_code <- c("C1", "C2", "C3", "C4", "C5",
                 "D1", "D2", "D3", "D4", "D5")

#############
# FUNCTIONS #
#############

split_fst_field <- function(fst_field){
  fst <- as.numeric(strsplit(fst_field, split="=", fixed=TRUE)[[1]][2])
  return(fst)
}

########
# MAIN #
########

###############################
# Comparing each gen13 sample #
# to the gen0 MB sample       #
###############################


db <- dbConnect(SQLite(), dbname="fst.sqlite")

# import all Fst files

setwd("~/Documents/Drosophila Selection Experiment/fst_and_diversity_files")

#rm(all_C_data_frame)
for(i in Sample_code){
  input_file_name <- paste("MB_", i, "_100000bp_window_151216.fst", sep="")
  dbWriteTable(conn = db, name = i, value = input_file_name,
            sep="\t", row.names = FALSE, header = TRUE)
  dbListFields(db, i)
}

dbSendQuery(db, 'CREATE TABLE C1_C2 AS SELECT * FROM C1 LEFT OUTER JOIN C2 ON (C1.chr = C2.chr AND C1.midpoint_pos=C2.midpoint_pos)')
dbSendQuery(db, 'CREATE TABLE C3_C4 AS SELECT * FROM C3 LEFT OUTER JOIN C4 ON (C3.chr = C4.chr AND C3.midpoint_pos=C4.midpoint_pos)')
dbSendQuery(db, 'CREATE TABLE C5_D1 AS SELECT * FROM C5 LEFT OUTER JOIN D1 ON (C5.chr = D1.chr AND C5.midpoint_pos=D1.midpoint_pos)')
dbSendQuery(db, 'CREATE TABLE D2_D3 AS SELECT * FROM D2 LEFT OUTER JOIN D3 ON (D2.chr = D3.chr AND D2.midpoint_pos=D3.midpoint_pos)')
dbSendQuery(db, 'CREATE TABLE D4_D5 AS SELECT * FROM D4 LEFT OUTER JOIN D5 ON (D4.chr = D5.chr AND D4.midpoint_pos=D5.midpoint_pos)')

dbSendQuery(db, 'CREATE TABLE C1_C2_C3_C4 AS SELECT * FROM C1_C2 LEFT OUTER JOIN C3_C4 ON (C1_C2.chr = C3_C4.chr AND C1_C2.midpoint_pos=C3_C4.midpoint_pos)')
dbSendQuery(db, 'CREATE TABLE C5_D1_D2_D3 AS SELECT * FROM C5_D1 LEFT OUTER JOIN D2_D3 ON (C5_D1.chr = D2_D3.chr AND C5_D1.midpoint_pos=D2_D3.midpoint_pos)')

dbSendQuery(db, 'CREATE TABLE C5_to_D5 AS SELECT * FROM C5_D1_D2_D3 LEFT OUTER JOIN D4_D5 ON (C5_D1_D2_D3.chr = D4_D5.chr AND C5_D1_D2_D3.midpoint_pos=D4_D5.midpoint_pos)')

dbSendQuery(db, 'CREATE TABLE all_reps AS SELECT * FROM C1_C2_C3_C4 LEFT OUTER JOIN C5_to_D5 ON (C1_C2_C3_C4.chr = C5_to_D5.chr AND C1_C2_C3_C4.midpoint_pos=C5_to_D5.midpoint_pos)')

testquery <- dbGetQuery(db, "SELECT `chr`, `midpoint_pos`, `fst`, `fst:1`, `fst:2`, `fst:3`, `fst:4`, `fst:5`, `fst:6`, `fst:7`, `fst:8`, `fst:9` FROM all_reps")


testquery <- dbGetQuery(db, "SELECT `chr`, `midpoint_pos` FROM all_reps")

colnames(testquery) <- c("chr", "midpoint_pos", Sample_code)

fst_only <- apply(testquery[,3:12], MARGIN = c(1,2), FUN = split_fst_field)

new_data <- data.frame(testquery[,1:2], fst_only)

new_data$C_mean 

###################################
# Comparing gen13 sample pairwise #
# to each other                   #
###################################

#Reading pairwise Fst for the all-vs-all comparison file

#all_input <- read.csv("all_MB_gen3_gen13_window_coarse.fst", sep="\t", header=FALSE)
# all_input <- read.csv("all_MB_gen3_gen13_window_coarse_lax.fst", sep="\t", header=FALSE)
all_input <- read.csv("MB_gen3_gen13_all_100000bp_window_151217.fst", sep="\t", header=FALSE)

chrs <- c("X", "2L", "2R", "3L", "3R")
samples <- c("C1", "C2", "C3", "C4", "C5", "D1", "D2", "D3", "D4", "D5")
pairwise_samples <- c("MB", "C1.3", "C1.13", "C2.3", "C2.13", "C3.3", "C3.13", "C4.3", "C4.13", "C5.3", "C5.13", "D1.3", "D1.13", "D2.3", "D2.13", "D3.3", "D3.13", "D4.3", "D4.13", "D5.3", "D5.13")



justchrs <- subset(all_input, as.character(all_input[,1])%in%chrs)

fst_labels <- c()
fst_type_first <- c()
fst_type_second <- c()
fst_gen_first <- c()
fst_gen_second <- c()
for (i in 1:20){
  for (j in (i+1):21){
    tempfirst <- pairwise_samples[i]
    tempsecond <- pairwise_samples[j]
    typefirst <- substr(tempfirst, 1, 1)
    genfirst <- unlist(strsplit(tempfirst, split="[.]"))[2]
    typesecond <- substr(tempsecond, 1, 1)
    gensecond <- unlist(strsplit(tempsecond, split="[.]"))[2]
    templabel <- paste(tempfirst, tempsecond, sep=":")
    fst_labels <- c(fst_labels, templabel)
    fst_type_first <- c(fst_type_first, typefirst)
    fst_type_second <- c(fst_type_second, typesecond)
    fst_gen_first <- c(fst_gen_first, genfirst)
    fst_gen_second <- c(fst_gen_second, gensecond)
  }
}
fst_gen_first[1:20] <- "0"

colnames(justchrs) <- c("chr", "pos", "no_snps", "frac_covered", "av_min_coverage", fst_labels)
split_fst <- function(x){
  val <- as.numeric(unlist(strsplit(as.character(x), "="))[2])
  return(val)
}

justchrs[,6:215] <- as.data.frame(lapply(justchrs[,6:215],FUN = function(x) {sapply(x,FUN=split_fst)}))
justchrs$pos_in_mb <- justchrs$pos/1000000

#in theory should be able to simplify this and vectorize my original function instead
#using something like
#Vsplit_fst <- Vectorize(split_fst)
#but not sure that this was working

exclude_gen_3 <- c(rep("TRUE", times=5), fst_gen_first != "3" & fst_gen_second != "3")
among_D_13 <- c(rep("TRUE", times=5), fst_gen_first != "3" & fst_gen_second != "3" & fst_type_first == "D" & fst_type_second == "D")
among_C_13 <- c(rep("TRUE", times=5), fst_gen_first != "3" & fst_gen_second != "3" & fst_type_first == "C" & fst_type_second == "C")

no_gen_3 <- justchrs[,exclude_gen_3==TRUE]

D_13 <- justchrs[,among_D_13==TRUE]
D_mean <- rowMeans(D_13[,6:15], na.rm=TRUE)
D_comps <- colnames(D_13)[6:15]
C_13 <- justchrs[,among_C_13==TRUE]
C_mean <- rowMeans(C_13[,6:15], na.rm=TRUE)
C_comps <- colnames(C_13)[6:15]
both_means <- data.frame(justchrs[,1:5], D_13, D_mean, C_13, C_mean)

view_X <- subset(both_means, chr=="X")





try_l <- layout(matrix(c(1,2,3,4,5),5,1,byrow=TRUE),widths=lcm(19),heights=lcm(5))
layout.show(try_l)
filename <- "Fst_among_replicate_lines_less_stringent_coverage_threshold.jpg"

jpeg(file=filename, width=21, height=16, units="cm", res=300)
layout(matrix(c(1,2,3,4,5,6),3,2,byrow=TRUE),widths=lcm(10.5),heights=lcm(5))
par(mar=c(3,3,1,1)+0.1)
par(mgp=c(2,1,0))
par(oma=c(0,0,0,0))
#fill in the top left-hand corner with an empty plot
plot(0,type='n',axes=FALSE,ann=FALSE)
for (i in chrs[1:5]){
  tempchr <- subset(both_means, chr==i)
  #tempmain <- paste("Chr", i, sep=" ")
  # use the following commented-out lines if you want to make individual jpgs per chr
  #print(paste('Plotting for', i))
  #filename <- paste(i, "coarse mean pairwise Fst.jpg", sep=" ")
  #jpeg(file=filename, width=20, height=8, units="cm", res=300)
  if (i=="2L"| i=="3L"){
    temp_ylab <- expression('Pairwise F'['ST'])
  }
  else temp_ylab <- ""
  if (i=="3R"| i=="3L"){
    temp_xlab <- "Position (Mb)"
  }
  else temp_xlab <- ""
  plot(tempchr[,11]~tempchr$pos_in_mb, type="l", col="pink", 
       ylab=temp_ylab, xlab=temp_xlab, ylim=c(0, 0.25), 
       #main=tempmain, 
       bty="l", cex.axis=1.2, lwd=0.5)
  for (j in 12:20){
    points(tempchr[,j]~tempchr$pos_in_mb, type="l", col="pink", lwd=0.5)
  }
  for (k in 27:36){
    points(tempchr[,k]~tempchr$pos_in_mb, type="l", col="lightblue", lwd=0.5)
  }
  points(tempchr$D_mean~tempchr$pos_in_mb, type="l", col="red")
  points(tempchr$C_mean~tempchr$pos_in_mb, type="l", col="blue")
  text(x=max(tempchr$pos_in_mb)-1, y=0.22, label=i, cex=2)
  #dev.off()
}
dev.off()



##################################

for (i in chrs){
  tempchr <-subset(both_means, chr==i)
  print(paste('Plotting for', i))
  filename <- paste(i, "Mean pairwise Fst for Desiccation-selected vs Control Lines.jpg", sep = " ")
  jpeg(file=filename, width=20, height=20, units="cm", res=300)
  plot(tempchr$D_mean~tempchr$C_mean, cex=0.5, type="p", col=rgb(100, 100, 100, 100, maxColorValue=255), xlab="Mean Fst among Control Lines", ylab="Mean Fst among Selected Lines", main=i)
  abline(a = 0, b = 1, col="red", lty=2)
  dev.off()
}

####### WORKING ON RESHAPING DATA SO AS TO APPLY A MIXED LINEAR MODEL #######
both_means$chrpos <- paste(both_means$chr, both_means$pos, sep="_")
test_values <- both_means[,c(1,11:20,27:36,38)]


melted <- melt(test_values)
melted$type <- substr(melted$variable, 1, 1)

# measured <- c(colnames(D_13[6:15]), colnames(C_13[6:15]))
# collist <- paste("both_means$", c(colnames(D_13[6:15]), colnames(C_13[6:15])))
# 
# reshaped <- reshape(both_means, varying=collist, times=collist, ids=both_means$chrpos, direction="long")
# testcast <- recast(both_means, chrpos~variable, id.var=38, measure.var= c(11:20, 27:36))

library(lme4)

fst_model_0 <- lmer(value~type+(1|chrpos), data=melted)
fst_model <- lmer(value~chr+type+(1|chrpos), data=melted)

#convert position to distance from centromere and include as a fixed effect?



#Plotting pairwise Fst for each pair of samples from the .igv file
for (eachsample in samples){
  othersamples <- samples[samples != eachsample]
  tempfstname <- paste('MB_', eachsample, '_noN_window.fst.igv', sep="")
  
  temp_fst <- read.csv(tempfstname, sep="\t", header=TRUE)
  colnames(temp_fst) <- c("chr", "posstart", "pos", "type", "fst")
  temp_fst$chrpos <- paste(temp_fst$chr, temp_fst$pos, sep="_")
  
  for (other in othersamples){
    otherfstname <- paste('MB_', other, '_noN_window.fst.igv', sep="")
    
    other_fst <- read.csv(otherfstname, sep="\t", header=TRUE)
    colnames(other_fst) <- c("chr", "posstart", "pos", "type", "fst2")
    other_fst$chrpos <- paste(other_fst$chr, other_fst$pos, sep="_")
    
    temp_merged <- merge(temp_fst, other_fst)
    temp_merged <- subset(temp_merged, chr%in%chrs)
  
      for (i in chrs){
        tempchr <-subset(temp_merged, chr==i)
        print(paste('Plotting for', i, eachsample, "vs", other, sep=" "))
        filename <- paste(eachsample, other, i, "pairwise Fst.jpg", sep = " ")
        jpeg(file=filename, width=20, height=20, units="cm", res=300)
        plot(tempchr$fst~tempchr$fst2, cex=0.5, pch=19, col=rgb(100, 100, 100, 100, maxColorValue=255), xlab=other, ylab=eachsample, main=i)
        dev.off()
      }
   }
}

#Now want to make a 'database' with one entry for each genomic window
#and the columns containing chr, pos, C1 Fst, C2 Fst, ... D5 Fst

for (eachsample in samples){
  tempfstfilename <- paste('MB_', eachsample, '_noN_window.fst.igv', sep="")
  tempfstvarname <- paste(eachsample, 'fst', sep="_")
  
  temp_fst <- read.csv(tempfstfilename, sep="\t", header=TRUE)
  colnames(temp_fst) <- c("chr", "posstart", "pos", "type", tempfstvarname)
  temp_fst$chrpos <- paste(temp_fst$chr, temp_fst$pos, sep="_")
  
  if (eachsample == "C1"){
    alldata <- temp_fst
  }
  else {
    alldata <- merge(alldata, temp_fst, all = TRUE)
    alldata <- subset(alldata, chr %in% chrs)
#     differences <- setdiff(alldata$chrpos, temp_fst$chrpos)
#     augmented <- data.frame(chr)
  }
}

D_values <- subset(alldata, D1_fst > 0.15 & D2_fst > 0.15 & D3_fst > 0.15 & D4_fst > 0.15 & D5_fst > 0.15)
convert_to_tf <- data.frame(alldata[,1:5], alldata[6:15]>0.15)
counts_of_d <- apply(convert_to_tf[,c(11:15)], 1, FUN=sum, na.rm=TRUE)
counts_of_c <- apply(convert_to_tf[,c(6:10)], 1, FUN=sum, na.rm=TRUE)

alldata2 <- data.frame(alldata, counts_of_c, counts_of_d)

high_d <- subset(alldata2, counts_of_d >= 4)
high_c <- subset(alldata2, counts_of_c >= 2)
high_d_and_high_c <- subset(alldata2, counts_of_d >= 2 & counts_of_c >= 2)

for (ch in chrs){

d1 <- subset(alldata2, chr == ch & counts_of_d == 1)
d2 <- subset(alldata2, chr == ch & counts_of_d == 2)
d3 <- subset(alldata2, chr == ch & counts_of_d == 3)
d4 <- subset(alldata2, chr == ch & counts_of_d == 4)
d5 <- subset(alldata2, chr == ch & counts_of_d == 5)
c1 <- subset(alldata2, chr == ch & counts_of_c == 1)
c2 <- subset(alldata2, chr == ch & counts_of_c == 2)
c3 <- subset(alldata2, chr == ch & counts_of_c == 3)
c4 <- subset(alldata2, chr == ch & counts_of_c == 4)
c5 <- subset(alldata2, chr == ch & counts_of_c == 5)

filename <- paste(ch, 'Fst_above_0.15.jpeg', sep="_")
jpeg(file=filename, width=20, height=10, units="cm", res=300)
plot(apply(d1[,c(11:15)], 1, max, na.rm=TRUE)~d1$pos, pch=19, cex=0.5, col="grey", ylim=c(0, 0.6), xlab="Position (bp)", ylab="Fst (F13 vs MB)", main = ch)
points(apply(d2[,c(11:15)], 1, max, na.rm=TRUE)~d2$pos, pch=19, cex=0.5, col="orange")
points(apply(d3[,c(11:15)], 1, max, na.rm=TRUE)~d3$pos, pch=19, cex=0.5, col="red")
points(apply(d4[,c(11:15)], 1, max, na.rm=TRUE)~d4$pos, pch=19, cex=0.5, col="blue")
points(apply(d5[,c(11:15)], 1, max, na.rm=TRUE)~d5$pos, pch=19, cex=0.5, col="green")
points(apply(c1[,c(6:10)], 1, max, na.rm=TRUE)~c1$pos, pch=21, cex=0.35, col="black", bg="white")
points(apply(c2[,c(6:10)], 1, max, na.rm=TRUE)~c2$pos, pch=21, cex=0.35, col="black", bg="orange")
points(apply(c3[,c(6:10)], 1, max, na.rm=TRUE)~c3$pos, pch=21, cex=0.35, col="black", bg="red")
points(apply(c4[,c(6:10)], 1, max, na.rm=TRUE)~c4$pos, pch=21, cex=0.35, col="black", bg="blue")
points(apply(c5[,c(6:10)], 1, max, na.rm=TRUE)~c5$pos, pch=21, cex=0.35, col="black", bg="green")
dev.off()
}

#Now want to make a 'database' with one entry for each genomic window
#and the columns containing chr, pos, C1 theta, C2 theta, ... C5 theta

Csamples <- c("C1", "C2", "C3", "C4", "C5")
Dsamples <- c("D1", "D2", "D3", "D4", "D5")
MBsamples <- "MB"

alldata <- data.frame(chr=character(),
                      pos=integer(), 
                      no_snps=integer(), 
                      prop_covered=numeric(), temptheta=numeric()) 
for (eachsample in Csamples){
  tempthetafilename <- paste(eachsample, '.13_noN.theta', sep="")
  tempthetavarname <- paste(eachsample, 'theta', sep="_")
  
  temp_theta <- read.csv(tempthetafilename, sep="\t", header=FALSE)
  colnames(temp_theta) <- c("chr", "pos", "no_snps", "prop_covered", tempthetavarname)
  temp_theta$chrpos <- paste(temp_theta$chr, temp_theta$pos, sep="_")
  temp_theta <- temp_theta[,c("chr", "pos", "chrpos", tempthetavarname)]
  
  if (eachsample == "C1"){
    alldata <- temp_theta
  }
  else {
    alldata <- merge(alldata, temp_theta, all = TRUE)
    alldata <- subset(alldata, chr %in% chrs)
    #     differences <- setdiff(alldata$chrpos, temp_fst$chrpos)
    #     augmented <- data.frame(chr)
  }
}

head(alldata, 20)
write.table(alldata[,c(1,2,4:8)], file="Control_lines_theta.txt", sep="\t", row.names=FALSE, quote=FALSE)

setwd("~/Documents/Drosophila Selection Experiment/fst_and_diversity_files")
alldata <- read.table("Control_lines_theta.txt", sep="\t", header=TRUE)

alldata_2L <- subset(alldata, alldata$chr=="2L")
plot(alldata_2L$pos, as.numeric(as.character(alldata_2L$C1_theta)), type="l", cex=0.1, pch=19, col=rgb(1, 0.5, 0, 0.5))
points(alldata_2L$pos, as.numeric(as.character(alldata_2L$C2_theta)), type="l", cex=0.1, pch=19, col=rgb(0, 0.3, 0.8, 0.5))

loessC1 <- loess(as.numeric(as.character(alldata_2L$C1_theta))~alldata_2L$pos, surface="interpolate", span=1, trace.hat="approximate")
plot(loessC1, type="l")


####### OLD STUFF ########

# setwd("~/Documents/Drosophila Selection Experiment/fst_and_diversity_files")
# chrs <- c("X", "2L", "2R", "3L", "3R", "4")
# samples <- c("C1", "C2", "C3", "C4", "C5", "D1", "D2", "D3", "D4", "D5")
# pairwise_samples <- c("MB", "C1.3", "C1.13", "C2.3", "C2.13", "C3.3", "C3.13", "C4.3", "C4.13", "C5.3", "C5.13", "D1.3", "D1.13", "D2.3", "D2.13", "D3.3", "D3.13", "D4.3", "D4.13", "D5.3", "D5.13")
# 
# #Plotting Fst vs pi for each chromosome for each sample
# for (eachsample in samples){
#     temppiname <- paste(eachsample, '.13_noN.pi', sep="")
#     tempfstname <- paste('MB_', eachsample, '_noN_window.fst.igv', sep="")
# 
#     temp_pi <- read.csv(temppiname, sep="\t", header=FALSE)
#     colnames(temp_pi) <- c("chr", "pos", "no_snps", "frac", "pi")
#     temp_pi$chrpos <- paste(temp_pi$chr, temp_pi$pos, sep="_")
#     temp_pi$pi <- replace(temp_pi$pi, temp_pi$pi=="na", NA)
#     temp_pi$pi <- as.numeric(as.character(temp_pi$pi))
# 
#     temp_fst <- read.csv(tempfstname, sep="\t", header=TRUE)
#     colnames(temp_fst) <- c("chr", "posstart", "pos", "type", "fst")
#     temp_fst$chrpos <- paste(temp_fst$chr, temp_fst$pos, sep="_")
# 
#     temp_merged <- merge(temp_pi, temp_fst)
#     temp_merged <- subset(temp_merged, chr%in%chrs)
# 
#     for (i in chrs){
#         tempchr <-subset(temp_merged, chr==i)
#         print(paste('Plotting for', eachsample, i, sep=" "))
#         filename <- paste(eachsample, i, "plot of Fst vs pi.jpg", sep = " ")
#         jpeg(file=filename, width=20, height=20, units="cm", res=300)
#         plot(tempchr$fst~tempchr$pi, cex=0.5, pch=19, col=rgb(100, 100, 100, 100, maxColorValue=255), ylim=c(0, 0.55), xlim=c(0, 0.035))
#         dev.off()
#     }
# }
# 
# Mpi <- read.csv("MB_noN.pi", sep="\t", header=FALSE)
# for (i in chrs){
#   tempchr <- subset(Mpi, Mpi[1]==i)
#   avpi <- mean(as.numeric(as.character(tempchr[,5])), na.rm=TRUE)
#   varpi <- var(as.numeric(as.character(tempchr[,5])), na.rm=TRUE)
#   print(paste(i, "Av. value of pi =", avpi, "Variance =", varpi, sep=" "))
# }




#Reading pairwise Fst for the pairwise MB vs gen13 comparison files




