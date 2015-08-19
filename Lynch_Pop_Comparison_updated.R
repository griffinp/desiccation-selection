setwd("/vlsci/VR0267/pgriffin/hsm/sandra/output_wgs/variant_calls/")
test_file <- read.csv('MB_C2_noN_reduced.mpileup.sync', sep = "\t", header=FALSE)
output_file_prefix <- "MB_C2_drift_SNP_test_150224_"
plot1_prefix <- "Minuslog10_p_vs_starting_freq_MB_C2_"
plot2_prefix <- "End_freq_vs_starting_freq_MB_C2_"

library(polysat)

# setwd("/Users/pgriffin/Documents/Drosophila Selection Experiment/Pileup_Files")
# test_file <- read.csv('MB_C1_test.sync', sep = "\t", header=FALSE)
# output_file_name <- "MB_C1_test_drift_SNP_test_150209.txt"


#############################################
#                                           #
# Simulating allele frequencies over 13 gen #
# under drift                               #
#                                           #
#############################################

N=200
#for the nuclear genome
x=1:(2*N)
#for the mitochondrial genome
mit_x=1:N

GenNum=13
Niter=1000
pA_table_nuclear=list()
Freq_bin=paste(seq(.001,.501,by=.01))
pA_start=.5

for (pA_start in Freq_bin) {

	pA_table_nuclear[[paste(pA_start)]] = matrix(NA,nrow=GenNum+1)[,-1]

	for (i in 1:Niter) {
	pA=as.numeric(pA_start)
	pA_record=pA
  
		for (g in 1:GenNum) {
			Pop=sample(c("A","B"),2*N,prob=c(pA,1-pA),replace=T)
			pA=table(Pop)["A"]/(2*N)
			if (is.na(pA)) pA=0
			pA_record=c(pA_record,pA)
		}
  
	pA_table_nuclear[[pA_start]]=cbind(pA_table_nuclear[[pA_start]],pA_record)
	}
}

#now repeat for the mitochondrial genome
for (pA_start in Freq_bin) {
  
  pA_table_nuclear[[paste(pA_start)]] = matrix(NA,nrow=GenNum+1)[,-1]
  
  for (i in 1:Niter) {
    pA=as.numeric(pA_start)
    pA_record=pA
    
    for (g in 1:GenNum) {
      Pop=sample(c("A","B"),2*N,prob=c(pA,1-pA),replace=T)
      pA=table(Pop)["A"]/(2*N)
      if (is.na(pA)) pA=0
      pA_record=c(pA_record,pA)
    }
    
    pA_table_nuclear[[pA_start]]=cbind(pA_table_nuclear[[pA_start]],pA_record)
  }
}


#############################################
#                                           #
# Calculating the 95% CI around the max     #
# likelihood estimate of final allele       #
# frequency for each simulated starting     #
# allele frequency                          #
#                                           #
#############################################

#This is the log-likelihood function for a log normal distribution (to deal with the case pA=0, I had to build it with a log(x+1))
loglik <- function(par, y) {sum(dnorm(log(1+y), mean=par[1], sd=sqrt(par[2]), log = TRUE))}
y=seq(0,1,0.01)

MLest=list()
#dev.new()
#plot('n',xlim=c(0,max(unlist(pA_table))),ylim=c(0,25),main=paste("Frequency density for intital freq =",pA_start),xlab="Allele Frequency",ylab="Density")
for (pA_start in Freq_bin) {
	#lines(density(pA_table[[pA_start]][14,]))
	MLest[[pA_start]]<- optim(c(mean = .2, var = .1), fn = loglik,y = pA_table[[pA_start]][14,], control = list(fnscale = -1,reltol = 1e-16))$par
  #MLest[[pA_start]]<- optimize(f=loglik, y=pA_table[[pA_start]][14,], maximum=TRUE, interval=c(0,1))$par
  #lines(dnorm(log(1+y), mean=MLest[[pA_start]][1], sd=sqrt(MLest[[pA_start]][2]))~y,lty=2)
}

esti=unlist(MLest)
esti=cbind(esti[grep("mean",names(esti))],esti[grep("var",names(esti))])
rownames(esti)=gsub(".mean","",rownames(esti))
Out=matrix(NA,ncol=5,nrow=length(Freq_bin))
rownames(Out)=Freq_bin
colnames(Out)=c("val", "lower_est","upper_est","lower_sim","upper_sim")
for (pA_start in Freq_bin) {
  Out[pA_start,1]<-round(as.numeric(as.character(pA_start)), digits=2)
  Out[pA_start,2]=qlnorm(0.025,esti[pA_start,1],sqrt(esti[pA_start,2]))-1
  Out[pA_start,3]=qlnorm(0.975,esti[pA_start,1],sqrt(esti[pA_start,2]))-1
  Out[pA_start,4]=quantile(pA_table[[pA_start]][14,],.025)
  Out[pA_start,5]=quantile(pA_table[[pA_start]][14,],.975)
}

# dev.new()
# plot('n',xlim=c(0,.5),ylim=c(0,.5),)
# abline(0,1)
# lines(Out[,3]~as.numeric(rownames(Out)),col=2,lty=2)
# lines(Out[,4]~as.numeric(rownames(Out)),col=2,lty=2)

#############################################
#                                           #
# Functions for reading in data.            #
# Dealing with real data from a two-pop     #
# 'sync' file (created from a pileup file)  #
#                                           #
#############################################

# Considerations
# - have to make sure the maj and min alleles are the same in both pops
# (can use my GWAS pileup file processing code for this)
# - have to make zero- and one-bounded (as Alex did above)

library("bbmle")

convert_to_bases <- function(a_vector){
  #A function to convert column positions (1,2,3,4) into base calls (A,C,G,T)
  call <- c('test')
  for (j in a_vector){
    if (j == 1){
      call <- c(call, "A")
    }
    if (j == 2){
      call <- c(call, "C")
    }
    if (j == 3){
      call <- c(call, "G")
    }
    if (j == 4){
      call <- c(call, "T")
    }
    if (j %in% c(5,6)){
      call <- c(call, "N")
    }
  }
  return(call[-1])
}

convert_to_vector <- function(base_call_row){
  # A function to convert a row from a sync file into a vector
  return(as.numeric(as.character(unlist(base_call_row))))
}

is_triallelic <- function(base_calls_vector){
  #Test of whether a locus is triallelic site
  total <- sum(base_calls_vector)
  fivepc <- total/20
  no_above_fivepc <- length(subset(base_calls_vector, base_calls_vector >= fivepc))
  return(no_above_fivepc >= 3)
}

find_maj_and_min <- function(base_calls){
  # Function to find the major and minor allele counts 
  # and column location
  majcount <- apply(base_calls, 1, max)
  majlocation <- apply(base_calls, 1, which.max)
  minlocation <- 0
  mincount <- 0
  for (i in c(1:nrow(base_calls))){
    arow <- base_calls[i,]
    tempcount <- as.integer(majcount[i])
    templocation <- majlocation[i]
    newrow <- arow[-templocation]
    mincount <- c(mincount, max(newrow))
    if (length(which(arow==majcount[i]))==1){
      minlocation <- c(minlocation, match(max(newrow), arow))
    }
    else {
      minlocation <- c(minlocation, which(arow==majcount[i])[2])
    }
    #print(minlocation)
  }
  minlocation <- minlocation[-1]
  mincount <- mincount[-1]
  majbase <- convert_to_bases(majlocation)
  minbase <- convert_to_bases(minlocation)
  return(data.frame(majcount, majbase,mincount, minbase))
}

find_overall_maj_and_min_allele <- function(all_data){
  # Function to find the allele that is major overall in
  # a two-population comparison
  # NB This version of the script sets the allele minor at Gen0
  # to be minor overall
  maj <- 'test'
  min <- 'test'
  for (i in c(1:nrow(all_data))){
    temprow <- all_data[i,]
    conmajbase <- as.character(temprow$majbase)
    conminbase <- as.character(temprow$minbase)
    selmajbase <- as.character(temprow$majbase.1)
    selminbase <- as.character(temprow$minbase.1)
    conmajcount <- temprow$majcount
    conmincount <- temprow$mincount
    selmajcount <- temprow$majcount.1
    selmincount <- temprow$mincount.1
    
    if (conmincount <= 1){
      if (conmajbase == selmajbase){
        output_maj <- conmajbase
        output_min <- selminbase
      }
      else{
        output_maj <- conmajbase
        output_min <- selmajbase
      }
    }
    if (selmincount <= 1){
      if (conmajbase == selmajbase | conminbase == selmajbase){
        output_maj <- conmajbase
        output_min <- conminbase
      }
      else{
        output_maj <- conmajbase
        output_min <- selmajbase
      }
    }
    if (conmincount > 1 & selmincount > 1){
      if (conmajbase == selmajbase){
        if (conminbase == selminbase){
          output_maj <- conmajbase
          output_min <- conminbase
        }
        else{
          output_maj <- conmajbase
          output_min <- selminbase
        }
      }
      if (conmajbase == selminbase){
        if (conminbase == selmajbase){
          output_maj <- conmajbase
          output_min <- conminbase
        }
      }
      if (conmajbase != selmajbase & conmajbase != selminbase){
        if (conminbase == selminbase){
            output_maj <- selmajbase
          output_min <- conminbase
        }
      }
    }
    maj <- c(maj, output_maj)
    min <- c(min, output_min)
  }
  return(data.frame(maj[-1], min[-1]))
}

find_columns <- function(majandmin){
  # This function takes a dataframe 'overall maj and min'
  # and outputs the column numbers  of a 'basecalls' dataframe
  # in which to find the right counts
  calldf <- data.frame(matrix(ncol=2, nrow=nrow(majandmin)))
 for (i in c(1:nrow(majandmin))){
   temprow <- majandmin[i,]
   maj <- temprow[1]
   min <- temprow[2]
   call <- c()
   for (each in c(maj, min)){
     if (each == "A"){
       call <- c(call, 1)
     }
     if (each == "C"){
       call <- c(call, 2)
     }
     if (each == "G"){
       call <- c(call, 3)
     }
     if (each == "T"){
       call <- c(call, 4)
     }
     if (each == "N"){
       call <- c(call, 6)
     }
   }
   calldf[i,] <- call
 }
 return(calldf)
}

#############################################
#                                           #
# Functions for MLE estimation              #
# of allele frequencies, allele frequency   #
# CIs and for the Lynch et al. likelihood   #
# test of whether two pools are sig diff in #
# allele frequency.                         #
#                                           #
#############################################

startprobmaj <- function(y, N, epsi){
  # This function calculates the error term
  # from Lynch et al. 2014 
  # for the probability of observing the major allele
  outmaj <- (y/(2*N))*(1-(4*epsi/3))+(epsi/3)
  return(outmaj)
  print(N)
}

startprobmin <- function(y, N, epsi){
  # This function calculates the error term
  # from Lynch et al. 2014 
  # for the probability of observing the minor allele
  outmin <- (y/(2*N))*((4*epsi/3)-1)+(1-epsi)
  return(outmin)
  print(N)
}

loglik.binom <- function(p, y, N, epsi) {
  # Basic log likelihood function for a single pop Y
  epsi = epsi
  sum(dbinom(y, size=1, prob=p, log = TRUE)+log(startprobmaj(y, N, epsi))+log(startprobmin(y, N, epsi)))
  #sum(log(1+dbinom(y, size=1, prob=p, log = FALSE)*(startprobmaj(y, N, epsi))*(startprobmin(y, N, epsi))))
}

loglik.binom2 <- function(p, y, N, epsi) {
  # Same as above, but negative
  # (required for using the 'mle2' function)
  epsi=epsi
  -sum(dbinom(y, size=1, prob=p, log = TRUE)+log(startprobmaj(y, N, epsi))+log(startprobmin(y, N, epsi)))
  #-sum(log(1+dbinom(y, size=1, prob=p, log = FALSE)*(startprobmaj(y, N, epsi))*(startprobmin(y, N, epsi))))
}

#To avoid getting a lot of -Inf values,
# try factorising the log - and adding 1. (Alex)

#-sum(log(1+dbinom(XXX,log=F)*min*maj))


loglik.binom.homogeneous <- function(p, y, z, Ny, Nz, epsi) {
  # Log likelihood function for the 'homogeneous' case
  # from Lynch et al. 2014
  # i.e. 
  epsi=epsi
  #sum(dbinom((y), size=1, prob=p, log = TRUE)+log(startprobmaj(y, Ny, epsi))+log(startprobmin(y, Ny, epsi)))+
  #sum(dbinom((z), size=1, prob=p, log = TRUE)+log(startprobmaj(z, Nz, epsi))+log(startprobmin(z, Nz, epsi)))
  sum(log(1+dbinom((y), size=1, prob=p, log = FALSE)*(startprobmaj(y, Ny, epsi))*(startprobmin(y, Ny, epsi))))+
  sum(log(1+dbinom((z), size=1, prob=p, log = FALSE)*(startprobmaj(z, Nz, epsi))*(startprobmin(z, Nz, epsi))))
}


lynch_test <- function(Ny, Nz, epsi, input_data_row){
  # This function implements the method from the Lynch et al. 2014 paper (on p 1215)
  # to compare allele frequencies between two pools
  
  count_y <- sum(input_data_row[,c(1:2)])
  count_z <- sum(input_data_row[,c(3:4)])
  py <- input_data_row[,1]/count_y
  pz <- input_data_row[,3]/count_z
  Ny=Ny
  Nz=Nz
  epsi=epsi
  y <- c(rep(0,times=input_data_row[1]), rep(1, times=input_data_row[2]))
  z <- c(rep(0, times=input_data_row[3]), rep(1, times=input_data_row[4]))
  pall <- sum(input_data_row[,c(1,3)])/sum(input_data_row)
  
  L01L02 <- loglik.binom.homogeneous(p=pall, y=y, z=z, Ny=Ny, Nz=Nz, epsi=epsi)
  Lf1=loglik.binom(py, y, Ny, epsi)
  Lf2=loglik.binom(pz, z, Nz, epsi)
  chi_val <- 2*(Lf1+Lf2-L01L02)
  #print(chi_val)
  if (is.nan(chi_val)){
    sig_test <- 1
  }
  else if (chi_val < 0) {
    sig_test <- pchisq((-1)*chi_val, df=1, lower.tail=FALSE)
  }
  else {
    sig_test <- pchisq(chi_val, df=1, lower.tail=FALSE)
  }  
  return(c(chi_val, sig_test))
}

mle_estimation <- function(epsi, N, input_data_row){
  # This function estimates allele frequency for a single input row
  # (as opposed to the whole dataframe) which is useful for error catching
  epsi = epsi
  N = N
  majmin_counts <- input_data_row
  
  maf <- majmin_counts[2]/(majmin_counts[1]+majmin_counts[2])
  snp_serie <- c(rep(0,times=majmin_counts[1]), rep(1, times=majmin_counts[2]))
  MLest <- optimize(f=loglik.binom, maximum=TRUE, interval=c(0,1), y=snp_serie, N=N, epsi=epsi)$maximum
  testmle <- mle2(loglik.binom2, optimizer="optimize", lower=0, upper=1, start=list(p=MLest), data=list(y=snp_serie, N=N, epsi=epsi))
  mle_max <- round(coef(testmle)["p"], digits=2)
  CI_95_lower <- round(confint(testmle)[[1]], digits=2)
  CI_95_upper <- round(confint(testmle)[[2]], digits=2)
  obs <- round(majmin_counts[2]/(majmin_counts[1]+majmin_counts[2]), digits=2)
 
  output <- c(obs, mle_max, CI_95_lower, CI_95_upper)   
  return(output)
}

when_mle_wont_converge <- function(epsi, N, input_data_row){
  #This function can be applied when the other estimator fails
  # which happens with zero or very low minor allele counts
  # - won't give CIs but will give just MLE
  
  epsi = epsi
  N = N
  majmin_counts <- input_data_row
  snp_serie <- c(rep(0,times=majmin_counts[1]), rep(1, times=majmin_counts[2]))
  MLest <- optimize(f=loglik.binom, maximum=TRUE, interval=c(0,1), y=snp_serie, N=N, epsi=epsi)$maximum
  mle_max <- round(MLest, digits=2)
  obs <- round(majmin_counts[2]/(majmin_counts[1]+majmin_counts[2]), digits=2)
  
  output <- c(obs, mle_max, NA, NA)
  return(output)
}


#############################################
#                                           #
#     Now read in the real data             #
#                                           #
#############################################

# make a subset for temporary ease of handling
print("Subsetting input file by p-value")
prelowestp <- subset(test_file, test_file[,6] <= 0.00001)
print("Subsetting by p-value completed")
#lowestp <- test_file
rm(test_file)

chrs <- c("2L", "2R", "3L", "3R", "X")
print("starting loop")

for (tempchr in chrs){
  print(paste("About to subset for", tempchr))
  lowestp <- subset(prelowestp, prelowestp[,1]==tempchr)
  print(paste(nrow(lowestp), "rows in subset"))
  output_file_name <- paste(output_file_prefix, tempchr, '.txt', sep="")
  print(paste("Processing to produce", output_file_name))
  
  gen0_base_call_col <- lowestp[,4]
  gen0_base_calls <- data.frame(do.call('rbind', strsplit(as.character(gen0_base_call_col), ':', fixed=TRUE)))
  gen0_base_vector <- apply(gen0_base_calls, 2, convert_to_vector)
  gen0_triallelic_calls <- apply(gen0_base_vector, 1, is_triallelic)
  gen0_base_calls <- apply(gen0_base_calls, 2, as.numeric)
  gen0_exp <- find_maj_and_min(gen0_base_calls)
  
  gen13_base_call_col <- lowestp[,5]
  gen13_base_calls <- data.frame(do.call('rbind', strsplit(as.character(gen13_base_call_col), ':', fixed=TRUE)))
  gen13_base_vector <- apply(gen13_base_calls, 2, convert_to_vector)
  gen13_triallelic_calls <- apply(gen13_base_vector, 1, is_triallelic)
  gen13_base_calls <- apply(gen13_base_calls, 2, as.numeric)
  gen13_exp <- find_maj_and_min(gen13_base_calls)
  
  rm(alldata)
  alldata <- data.frame(lowestp, gen0_exp, gen13_exp)
  overall_majmin <- find_overall_maj_and_min_allele(alldata)
  cols_to_grab <- find_columns(overall_majmin)
  majmin_counts <- data.frame(matrix(nrow=nrow(cols_to_grab), ncol=4))
  colnames(majmin_counts) <- c("gen0_maj", "gen0_min", "gen13_maj", "gen13_min")
  
  for (eachrow in c(1:nrow(cols_to_grab))){
    tempcols <- as.numeric(cols_to_grab[eachrow,])
    majmin_counts[eachrow,1:2] <- gen0_base_calls[eachrow,tempcols]
    majmin_counts[eachrow,3:4] <- gen13_base_calls[eachrow,tempcols]
  }
  
  #majmin_counts <- majmin_counts[500:700,]
  
  #############################################
  #                                           #
  #  Perform the Lynch et al. population      #
  #  allele frequency comparison              #
  #                                           #
  #############################################
  
  #Lynch population comparison
  lynchoutput <- data.frame(matrix(ncol=2, nrow=nrow(majmin_counts)))
  for (eachrow in 1:nrow(majmin_counts)){
    temprow <- majmin_counts[eachrow,]
    templynch <- lynch_test(Ny=200, Nz=50, epsi=0.01, input_data_row=temprow)
    lynchoutput[eachrow,] <- templynch
  }
  
  #############################################
  #                                           #
  #  Calculate individual pop frequency       #
  #  MLE and ML 95% CI                        #
  #                                           #
  #############################################
  
  startinput <- majmin_counts[,1:2]
  endinput <- majmin_counts[,3:4]
  
  print(paste("Length of startinput = ", nrow(startinput)))
  
  startoutput <- data.frame(matrix(nrow=nrow(startinput), ncol=4))
  epsi=0.01
  N=200
  for (eachrow in 1:nrow(startinput)){
    tempstartinput <- startinput[eachrow,]
    
    row_result <- tryCatch({
      output_line <- mle_estimation(epsi=0.01, N=200, input_data_row=tempstartinput)
      #print(output_line)
    }, error = function(err) {
      print("MLE error detected; alternative script running")
      output_line <- when_mle_wont_converge(epsi=0.01, N=200, input_data_row=tempstartinput)
      return(output_line)
    })
    startoutput[eachrow,] <- row_result
    
  }
  
  endoutput <- data.frame(matrix(nrow=nrow(endinput), ncol=4))
  for (eachrow in 1:nrow(endinput)){
    tempendinput <- endinput[eachrow,]
    row_result <- tryCatch({
      mle_estimation(epsi=0.01, N=50, input_data_row=tempendinput)
    }, error = function(err) {
      print("MLE error detected; alternative script running")
      error_out_line <- when_mle_wont_converge(epsi=0.01, N=50, input_data_row=tempendinput)
      return(error_out_line)
    })
    endoutput[eachrow,] <- row_result
  }
  
  
  colnames(startoutput) <- c("start_obs", "start_mle", "start_lower", "start_upper")
  colnames(endoutput) <- c("end_obs", "end_mle", "end_lower", "end_upper")
  
  
  direction <- as.numeric(startoutput$start_mle) < as.numeric(endoutput$end_mle)
  
  start_limit <- c()
  end_limit <- c()
  for (i in 1:nrow(startoutput)){
    direc <- direction[i]
    if (direc=='TRUE') {
      start_limit_to_add <- as.numeric(as.character(startoutput$start_upper))[i]
      end_limit_to_add <- as.numeric(as.character(endoutput$end_lower))[i]
    } 
    if (direc=='FALSE') {
      start_limit_to_add <- as.numeric(as.character(startoutput$start_lower))[i]
      end_limit_to_add <- as.numeric(as.character(endoutput$end_upper))[i]
    }
    if (is.na(as.factor(startoutput$start_lower)[i])){
      start_limit_to_add <- as.numeric(as.character(startoutput$start_mle))[i]
    }
    if (is.na(as.factor(endoutput$end_lower)[i])){
      end_limit_to_add <- as.numeric(as.character(endoutput$end_mle))[i]
    }
    start_limit <- c(start_limit, start_limit_to_add)
    end_limit <- c(end_limit, end_limit_to_add)
  }
  
  final_table <- data.frame(startoutput, endoutput, start_limit, end_limit)
  
  sig_test <- c()
  # My homemade significance test of whether
  # the allele freq change as calculated by
  # starting at the upper bound of Pop 1 allele freq estimate
  # and ending the lower bound of Pop 2 allele freq CI
  # is in the 2.5% tail of the simulated distribution
  # (calculated way back at the start)
  for (i in 1:nrow(final_table)){
    temprow <- final_table[i,]
    tempstart <- temprow$start_limit
    tempend <- temprow$end_limit
    if (tempstart > 0.5){
      tempstart <- 1-tempstart
      tempend <- 1-tempend
    }
    temptest <- Out[as.character(Out[,1])==as.character(tempstart),]
    if (tempend < temptest[4] | tempend > temptest[5]){
      sig_test <- c(sig_test, "sig")
    }
    else {
      sig_test <- c(sig_test, "ns")
    }
  }
  
  out_tab <- data.frame(lowestp, final_table, sig_test, lynchoutput)
  
  #Finally, over-write the nonsense results for triallelic loci
  for (i in 1:nrow(out_tab)){
    gen0_tri <- gen0_triallelic_calls[i]
    gen13_tri <- gen13_triallelic_calls[i]
    if (gen0_tri | gen13_tri){
      out_tab[i,7:19] <- "NA"
    }
  }
  
  jpeg(file=paste(plot1_prefix, tempchr, ".jpg", sep=""))
  
  plot((-1)*log10(out_tab[,6])~out_tab$start_obs, col=as.factor(out_tab$sig_test), pch=19, cex=0.5)
  dev.off()
  
  jpeg(file=paste(plot2_prefix, tempchr, ".jpg", sep=""))
  
  plot(out_tab$end_obs~out_tab$start_obs, col=as.factor(out_tab$sig_test), pch=19, cex=0.5)
  dev.off()
  
  #jpeg(file="Lynch_test_stat_distribution.jpg")
  
  #hist(out_tab[,18], breaks=100)
  #dev.off()
  
  write.table(out_tab, file=output_file_name, quote=FALSE, row.names=FALSE, sep="\t")
  rm(lowestp)
}
  
################################################################

                       