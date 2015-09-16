# setwd("/vlsci/VR0267/pgriffin/hsm/sandra/output_wgs/variant_calls/")
# test_file <- read.csv('MB_C1_noN_reduced_150819.mpileup.sync', sep = "\t", header=FALSE)
# output_file_prefix <- "MB_C1_drift_SNP_test_150819_"
# plot1_prefix <- "Minuslog10_p_vs_starting_freq_MB_C1_"
# plot2_prefix <- "End_freq_vs_starting_freq_MB_C1_"

#library(polysat)
library("bbmle")


setwd("/Users/pgriffin/Documents/Drosophila Selection Experiment/snp_and_gene_lists")
test_file <- read.csv('D5_fet_sig_SNPs_sync.txt', sep = "\t", header=FALSE)
output_file_name <- "MB_D5_lynch_for_fet_sig_snps_150820.txt"
output_file_prefix <- "MB_D5_lynch_for_fet_sig_snps_150820_"

#############################################
#                                           #
# Simulating allele frequencies over 13 gen #
# under drift                               #
#                                           #
#############################################

N=200
#x=1:(2*N)

GenNum=13
Niter=1000
pA_table_nuclear=list()
pA_table_mito=list()
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
  
  pA_table_mito[[paste(pA_start)]] = matrix(NA,nrow=GenNum+1)[,-1]
  
  for (i in 1:Niter) {
    pA=as.numeric(pA_start)
    pA_record=pA
    
    for (g in 1:GenNum) {
      Pop=sample(c("A","B"),N,prob=c(pA,1-pA),replace=T)
      pA=table(Pop)["A"]/(N)
      if (is.na(pA)) pA=0
      pA_record=c(pA_record,pA)
    }
    
    pA_table_mito[[pA_start]]=cbind(pA_table_mito[[pA_start]],pA_record)
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

#For nuclear genome
MLest_nuclear=list()
for (pA_start in Freq_bin) {
	MLest_nuclear[[pA_start]]<- optim(c(mean = .2, var = .1), fn = loglik,y = pA_table_nuclear[[pA_start]][14,], control = list(fnscale = -1,reltol = 1e-16))$par
}

esti=unlist(MLest_nuclear)
esti=cbind(esti[grep("mean",names(esti))],esti[grep("var",names(esti))])
rownames(esti)=gsub(".mean","",rownames(esti))
Out_nuclear=matrix(NA,ncol=5,nrow=length(Freq_bin))
rownames(Out_nuclear)=Freq_bin
colnames(Out_nuclear)=c("val", "lower_est","upper_est","lower_sim","upper_sim")
for (pA_start in Freq_bin) {
  Out_nuclear[pA_start,1]<-round(as.numeric(as.character(pA_start)), digits=2)
  Out_nuclear[pA_start,2]=qlnorm(0.025,esti[pA_start,1],sqrt(esti[pA_start,2]))-1
  Out_nuclear[pA_start,3]=qlnorm(0.975,esti[pA_start,1],sqrt(esti[pA_start,2]))-1
  Out_nuclear[pA_start,4]=quantile(pA_table_nuclear[[pA_start]][14,],.025)
  Out_nuclear[pA_start,5]=quantile(pA_table_nuclear[[pA_start]][14,],.975)
}

#For mitochondrial genome
MLest_mito=list()
for (pA_start in Freq_bin) {
  MLest_mito[[pA_start]]<- optim(c(mean = .2, var = .1), fn = loglik,y = pA_table_mito[[pA_start]][14,], control = list(fnscale = -1,reltol = 1e-16))$par
}

esti=unlist(MLest_mito)
esti=cbind(esti[grep("mean",names(esti))],esti[grep("var",names(esti))])
rownames(esti)=gsub(".mean","",rownames(esti))
Out_mito=matrix(NA,ncol=5,nrow=length(Freq_bin))
rownames(Out_mito)=Freq_bin
colnames(Out_mito)=c("val", "lower_est","upper_est","lower_sim","upper_sim")
for (pA_start in Freq_bin) {
  Out_mito[pA_start,1]<-round(as.numeric(as.character(pA_start)), digits=2)
  Out_mito[pA_start,2]=qlnorm(0.025,esti[pA_start,1],sqrt(esti[pA_start,2]))-1
  Out_mito[pA_start,3]=qlnorm(0.975,esti[pA_start,1],sqrt(esti[pA_start,2]))-1
  Out_mito[pA_start,4]=quantile(pA_table_mito[[pA_start]][14,],.025)
  Out_mito[pA_start,5]=quantile(pA_table_mito[[pA_start]][14,],.975)
}


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

####
####
####

split_by_colon <- function(base_call_object){
  object_split <- strsplit(base_call_object, ":")
  object_as_vector <- as.numeric(unlist(object_split))
  names(object_as_vector) <- c("A", "T", "C", "G", "N", "gap")
  return(object_as_vector)
}

find_maj_and_min <- function(base_count_col){
  base_count_col <- base_count_col[1:4]
  max_location <- base_count_col[which.max(base_count_col)]
  max_name <- names(max_location)
  subtract_max <- base_count_col[-which.max(base_count_col)]
  min_location <- subtract_max[which.max(subtract_max)]
  min_name <- names(min_location)
  return(c(max_name, min_name))
}

extract_counts_matching_letters <- function(base_call_col){
  # This is to be used on a column from a data frame
  # containing counts for A to gap, 
  # then maj_letter and then min_letter
  # e.g. 
  #   
  #   [,1] [,2] [,3] [,4]
  #   A          "30" "0"  "0"  "0" 
  #   T          "5"  "43" "0"  "0" 
  #   C          "0"  "22" "30" "93"
  #   G          "0"  "1"  "2"  "50"
  #   N          "0"  "0"  "0"  "0" 
  #   gap        "0"  "0"  "0"  "0" 
  #   maj_letter "A"  "T"  "C"  "G" 
  #   min_letter "T"  "C"  "G"  "C" 
  #desired_letters <- which(row.names(base_call_col%in%letter_col)
  maj_location <- which(names(base_call_col)==base_call_col['maj_letter'])
  min_location <- which(names(base_call_col)==base_call_col['min_letter'])
  desired_counts <- base_call_col[c(maj_location, min_location)]
  return(as.numeric(desired_counts))
}

find_maj_and_min_for_2_sample_file <- function(base_call_df_1, base_call_df_2){
  split1 <- apply(base_call_df_1, 1, split_by_colon)
  split2 <- apply(base_call_df_2, 1, split_by_colon)
  total_counts <- split1+split2
  overall_maj_min <- apply(total_counts, 2, find_maj_and_min)
  return(list(split1, split2, overall_maj_min))
}

extract_maj_and_min_for_2_sample_file <- function(base_call_df_1, base_call_df_2){
  relevant_info <- find_maj_and_min_for_2_sample_file(base_call_df_1, base_call_df_2)
  rownames(relevant_info[[3]]) <- c("maj_letter", "min_letter")
  base_counts_df_1 <- rbind(relevant_info[[1]], relevant_info[[3]])
  base_counts_df_2 <- rbind(relevant_info[[2]], relevant_info[[3]])  
  allele_calls <- as.data.frame(matrix(NA, ncol=6,
                                       nrow=nrow(base_call_df_1))
  )
  df1_majmin_counts <- apply(base_counts_df_1, 2, extract_counts_matching_letters)
  df2_majmin_counts <- apply(base_counts_df_2, 2, extract_counts_matching_letters)
  allele_calls[,3:4] <- c(df1_majmin_counts[1,], df1_majmin_counts[2,])
  allele_calls[,5:6] <- c(df2_majmin_counts[1,], df2_majmin_counts[2,])
  allele_calls[,1:2] <- c(relevant_info[[3]][1,], relevant_info[[3]][2,])
  colnames(allele_calls) <- c("overall_maj", "overall_min", "df_1_maj_count",
                              "df_1_min_count", "df_2_maj_count", "df_2_min_count")
  rownames(allele_calls) <- 1:nrow(allele_calls)
  return(allele_calls)
}


rearrange_maj_min_and_counts <- function(majmin_call_object){
  #This function rearranges the maj and min calls so that the 'maj' allele
  #is the one that decreased in frequency between sample 1 and sample 2
  #It also outputs the counts correctly
  majmin_call_object <- as.data.frame(majmin_call_object)
  total_df_1_count <- as.numeric(as.character(majmin_call_object$df_1_maj_count))+as.numeric(as.character(majmin_call_object$df_1_min_count))
  df_1_maj_proportion <- as.numeric(as.character(majmin_call_object$df_1_maj_count))/total_df_1_count
  total_df_2_count <- as.numeric(as.character(majmin_call_object$df_2_maj_count))+as.numeric(as.character(majmin_call_object$df_2_min_count))
  df_2_maj_proportion <- as.numeric(as.character(majmin_call_object$df_2_maj_count))/total_df_2_count
  
  maj_has_increased <- df_2_maj_proportion > df_1_maj_proportion
  new_df <- majmin_call_object
  #new_maj <- as.character(majmin_call_object$overall_maj)
  #new_min <- as.character(majmin_call_object$overall_min)
  new_df[which(maj_has_increased), 'overall_maj'] <- as.character(majmin_call_object$overall_min[which(maj_has_increased)])
  new_df[which(maj_has_increased), 'overall_min'] <- as.character(majmin_call_object$overall_maj[which(maj_has_increased)])
  new_df[which(maj_has_increased), 'df_1_maj_count'] <- as.character(majmin_call_object$df_1_min_count[which(maj_has_increased)])
  new_df[which(maj_has_increased), 'df_1_min_count'] <- as.character(majmin_call_object$df_1_maj_count[which(maj_has_increased)])
  new_df[which(maj_has_increased), 'df_2_maj_count'] <- as.character(majmin_call_object$df_2_min_count[which(maj_has_increased)])
  new_df[which(maj_has_increased), 'df_2_min_count'] <- as.character(majmin_call_object$df_2_maj_count[which(maj_has_increased)])
  
  
  return(new_df)  
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
  sum(dbinom((y), size=1, prob=p, log = TRUE)+log(startprobmaj(y, Ny, epsi))+log(startprobmin(y, Ny, epsi)))+
  sum(dbinom((z), size=1, prob=p, log = TRUE)+log(startprobmaj(z, Nz, epsi))+log(startprobmin(z, Nz, epsi)))
  #sum(log(1+dbinom((y), size=1, prob=p, log = FALSE)*(startprobmaj(y, Ny, epsi))*(startprobmin(y, Ny, epsi))))+
  #sum(log(1+dbinom((z), size=1, prob=p, log = FALSE)*(startprobmaj(z, Nz, epsi))*(startprobmin(z, Nz, epsi))))
}


lynch_test <- function(Ny, Nz, epsi, input_data_row){
  # This function implements the method from the Lynch et al. 2014 paper (on p 1215)
  # to compare allele frequencies between two pools
  
  count_y <- sum(input_data_row[1:2])
  count_z <- sum(input_data_row[3:4])
  py <- input_data_row[1]/count_y
  pz <- input_data_row[3]/count_z
  Ny=Ny
  Nz=Nz
  epsi=epsi
  y <- c(rep(0,times=input_data_row[1]), rep(1, times=input_data_row[2]))
  z <- c(rep(0, times=input_data_row[3]), rep(1, times=input_data_row[4]))
  pall <- sum(input_data_row[c(1,3)])/sum(input_data_row)
  
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

# lowestp <- test_file
# rm(test_file)

chrs <- c("2L", "2R", "3L", "3R", "X", "dmel_mitochondrion_genome")
print("starting loop")

for (tempchr in chrs){
  print(paste("About to subset for", tempchr))
  chr_subset <- subset(test_file, test_file[,1]==tempchr)
  print(paste(nrow(chr_subset), "rows in subset"))
  output_file_name <- paste(output_file_prefix, tempchr, '.txt', sep="")
  print(paste("Processing to produce", output_file_name))
    
  temp_base_call_col_1 <- data.frame(base_calls=chr_subset[,4])
  temp_base_call_col_2 <- data.frame(base_calls=chr_subset[,5])
  majmin_calls <- extract_maj_and_min_for_2_sample_file(temp_base_call_col_1,
                                                      temp_base_call_col_2)
  majmin_counts <- rearrange_maj_min_and_counts(majmin_calls)[3:6]
  colnames(majmin_counts) <- c("gen0_maj", "gen0_min", "gen13_maj", "gen13_min")
  majmin_counts <- sapply(majmin_counts, as.numeric)

  
  #############################################
  #                                           #
  #  Perform the Lynch et al. population      #
  #  allele frequency comparison              #
  #                                           #
  #############################################
  
  #Lynch population comparison
  #lynchoutput <- data.frame(matrix(ncol=2, nrow=nrow(majmin_counts)))
  lynch_output <- t(as.data.frame(apply(majmin_counts, 1, lynch_test, Ny=200, Nz=50, epsi=0.01)))
  #for (eachrow in 1:nrow(majmin_counts)){
  #  temprow <- majmin_counts[eachrow,]
  #  templynch <- lynch_test(Ny=200, Nz=50, epsi=0.01, input_data_row=temprow)
  #  lynchoutput[eachrow,] <- templynch
  #  if(eachrow%%10000<1){
  #    print(paste("Lynch MLE estimation performed for" eachrow, "loci"))
  #  }
  #}
  
  colnames(lynch_output) <- c("lynch_chi_val", "lynch_p")

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
  
  limits <- matrix(NA, nrow=nrow(startoutput), ncol=2)
  #end_limit <- c()
  limits[which(direction==TRUE),] <- as.matrix(data.frame(startoutput$start_upper, 
                                                   endoutput$end_lower))[which(direction==TRUE),]
  limits[which(direction==FALSE),] <- as.matrix(data.frame(startoutput$start_lower, 
                                                           endoutput$end_upper))[which(direction==FALSE),]
  limits[which(is.na(as.factor(startoutput$start_lower))),1] <- as.matrix(data.frame(startoutput$start_mle))[which(is.na(as.factor(startoutput$start_lower))),]
  limits[which(is.na(as.factor(endoutput$end_lower))),2] <- as.matrix(data.frame(endoutput$end_mle))[which(is.na(as.factor(endoutput$end_lower))),]

  colnames(limits) <- c('start_limit', 'end_limit')

  final_table <- data.frame(chr_subset, startoutput, endoutput, limits)
  
  # My homemade significance test of whether
  # the allele freq change as calculated by
  # starting at the upper bound of Pop 1 allele freq estimate
  # and ending the lower bound of Pop 2 allele freq CI
  # is in the 2.5% tail of the simulated distribution
  # (calculated way back at the start)
  
  sig_test <- rep("ns", times=nrow(final_table))
  tempstartpre <- final_table$start_limit
  tempstart <- tempstartpre
  tempend <- final_table$end_limit
  tempstart[which(tempstartpre > 0.5)] <- round(1-final_table$start_limit, 2)[which(tempstartpre > 0.5)]
  tempend[which(tempstartpre > 0.5)] <- round(1-final_table$end_limit, 2)[which(tempstartpre > 0.5)]
  #make sure to sample from the right drift simulation
  if(tempchr=="dmel_mitochondrion_genome"){
    temp_test <- Out_mito[match(as.character(tempstart), Out_mito[,1]),]
    #special treatment for starting allele freq of 0
    for(z in which(as.character(tempstart)=='0')){
      temp_test[z,] <- Out_mito[which(Out_mito[,1]=='0.01'),]
    }
  } else {
    temp_test <- Out_nuclear[match(as.character(tempstart), Out_nuclear[,1]),]
    #special treatment for starting allele freq of 0
    for(z in which(as.character(tempstart)=='0')){
      temp_test[z,] <- Out_nuclear[which(Out_nuclear[,1]=='0.01'),]
    }
  }
  sig_test[which(tempend < temp_test[,4] | tempend > temp_test[,5])] <- "sig"

  out_tab <- data.frame(final_table, sig_test, lynch_output)
  #out_tab <- data.frame(tempend, temp_test[,4:5], sig_test) 
  
  write.table(out_tab, file=output_file_name, quote=FALSE, row.names=FALSE, sep="\t")
  #rm(test_file)
}
  
################################################################

                       