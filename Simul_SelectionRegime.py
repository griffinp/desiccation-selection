#! usr/bin/env python

###########################################################################
# 
#              Python script written by Alex Fournier-Level 
#        supporting the results from Griffin et al. (submitted)
#           freely available and modifiable under GNU licence
# requires -pylab
#          -a gene list where each entry is:
# 	chromosome	position_start	position_end	gene_name	number_of_snps
#ex:		3L	21427	27308	mthl8;FBgn0052475	1828
#			3L	52872	55751	CG43149;FBgn0262679	1412
#          -the recombination landscape to run with the Comeron et al. (2012)
#           which can be found here:http://petrov.stanford.edu/RRC_scripts/RRCv2.3.tar.gz
###########################################################################

import sys
from numpy import *
from scipy import stats

Graph=False
Verbose=False

if Graph in sys.argv: Graph=True
if Verbose in sys.argv: Verbose=True

if Graph: from matplotlib import pyplot as plt

#The parameters of the problems
#The heritability of the ONE SNP nearing an effect
h2=.01
nFX=1
#nHaplo have G_rate offsprings each
nHaplo=400 #needs to be an even number
G_rate=10 #growth rate
Sel=0.1 #selection intensity in %, beware that the product G_rate*Sel needs to be kept ~1 otherwise the CPU will die :'(
nGeneration=13
regime=["Selection","Selection","Selection","Selection","Selection","Selection","Selection","Selection","Control","Selection","Control","Selection","Control"]
FET_thresh=3.5e-13

#Pick (1) a gene showing some association and (2) a closely and (3) loosely linked gene
CHR_list=['2L','2R','3L','3R']#,'X']
CHR=random.choice(CHR_list)
GENE_map=genfromtxt(CHR+'_genes_snps_intersection.txt',dtype=str)

x_gene=random.randint(GENE_map.shape[0]-1)
gene=GENE_map[x_gene,:]
x=int(gene[1])
nSNP_X=int(gene[4])
ID_X=gene[3]

gene_St=GENE_map[absolute(x-GENE_map[:,1].astype(int))<10000,:]
gene_St=gene_St[[0,gene_St.shape[0]-1],:]
gene_St=ravel(gene_St[absolute(gene_St[:,1].astype(int)-x)==max(absolute(ravel(gene_St[:,1].astype(int))-x)),:])

nSNP_St=int(gene_St[4])
ID_St=gene_St[3]

gene_Lg=GENE_map[absolute(x-GENE_map[:,1].astype(int))<1e6,:]
gene_Lg=gene_Lg[[0,gene_Lg.shape[0]-1],:]
gene_Lg=ravel(gene_Lg[absolute(gene_Lg[:,1].astype(int)-x)==max(absolute(ravel(gene_Lg[:,1].astype(int))-x)),:])

nSNP_Lg=int(gene_Lg[4])
ID_Lg=gene_Lg[3]

nSNP=array((nSNP_X,nSNP_St,nSNP_Lg),dtype=float)
(nSNP_X,nSNP_St,nSNP_Lg)=tuple((min(nSNP)*nSNP/(3*sum(nSNP))).astype(int))

POS_X=sort(random.randint(x,int(gene[2]),nSNP_X))
POS_St=sort(random.randint(int(gene_St[1]),int(gene_St[2]),nSNP_St))
POS_Lg=sort(random.randint(int(gene_Lg[1]),int(gene_Lg[2]),nSNP_Lg))

#Generating the allele frequency distributon
(a,b)=(0.8294141,0.8291769)
freq=random.beta(a,b,nSNP_X+nSNP_St+nSNP_Lg)

#Generating the genotypes
COV_X1 = random.uniform(0,1,size=(nSNP_X,nSNP_X))
COV_X1 = (COV_X1 + COV_X1.T)/2
COV_X2 = random.uniform(0,1,size=(nSNP_St,nSNP_St))
COV_X2 = (COV_X2 + COV_X2.T)/2
COV_X3 = random.uniform(0,1,size=(nSNP_Lg,nSNP_Lg))
COV_X3 = (COV_X3 + COV_X3.T)/2
COV_St1 = random.uniform(0,.8,size=(nSNP_X,nSNP_St))
COV_St2 = random.uniform(0,.4,size=(nSNP_St,nSNP_Lg))
COV_Lg = random.uniform(0,.2,size=(nSNP_X,nSNP_Lg))

COV=vstack((hstack((COV_X1,COV_St1,COV_Lg)),hstack((COV_St1.T,COV_X2,COV_St2)),hstack((COV_Lg.T,COV_St2.T,COV_X3))))
row_sums = COV.sum(axis=1)
COV = COV / row_sums[:, newaxis]
try:
	X=random.multivariate_normal(mean=freq,cov=COV,size=nHaplo)
except ValueError:
	exit()

X[X<.5]=0
X[X>=.5]=1

X_St=X[:,nSNP_X:(nSNP_X+nSNP_St)]
X_Lg=X[:,(nSNP_X+nSNP_St):]
X=X[:,:nSNP_X]

#And removing monomorphic SNPs
SNPsum=X.sum(axis=0)
mask=all([SNPsum!=0,SNPsum!=nHaplo],axis=0)
X=X[:,mask]
POS_X=POS_X[mask]
nSNP_X=X.shape[1]

SNPsum=X_St.sum(axis=0)
mask=all([SNPsum!=0,SNPsum!=nHaplo],axis=0)
X_St=X_St[:,mask]
POS_St=POS_St[mask]
nSNP_St=X_St.shape[1]

SNPsum=X_Lg.sum(axis=0)
mask=all([SNPsum!=0,SNPsum!=nHaplo],axis=0)
X_Lg=X_Lg[:,mask]
POS_Lg=POS_Lg[mask]
nSNP_Lg=X_Lg.shape[1]

#estimating the allele frequencies in the data
Geno=hstack((X,X_St,X_Lg))
Geno_ini=Geno
Haplo_ini=unique(apply_along_axis(array_str,1,Geno))
nHaplo_ini=Haplo_ini.shape[0]
if Verbose:
	print("The population will be started with "+str(nHaplo_ini)+" distinct haplotypes")
	print(str(nSNP_X)+' SNPs in focal gene/'+str(nSNP_St)+' in close gene/'+str(nSNP_Lg)+' in distant gene')
freq=(Geno[:nHaplo/2,:]+Geno[nHaplo/2:,:]).sum(axis=0)/float(nHaplo)

#Computing the coverage based on a normal distribution truncated in 0
lower=0
upper=inf #inf is part of numpy
if CHR=="2L": (mu,sigma)=(135.82226,33.77574)
if CHR=="2R": (mu,sigma)=(141.64910,34.74283)
if CHR=="3L": (mu,sigma)=(135.20783,34.72228)
if CHR=="3R": (mu,sigma)=(138.68089,34.06114)
Coverage=stats.truncnorm.rvs((lower-mu)/sigma,(upper-mu)/sigma,loc=mu,scale=sigma,size=2*(nSNP_X+nSNP_St+nSNP_Lg)).reshape(2,nSNP_X+nSNP_St+nSNP_Lg)

#drawing genetic effects and assigning associated SNPs
GFX=random.random_sample(nFX)
GFX=GFX/sum(GFX)
AssoSNP=array(random.choice(range(nSNP_X),nFX,replace=False),dtype="int")
if Verbose: print("The associated allele is found at a frequency of "+str(freq[AssoSNP]))

#partionning the variance taking into account the linkage among associated SNPs
Geno=X[:nHaplo/2,:]+X[nHaplo/2:,:]
AssoX=Geno[:,AssoSNP]
COV=corrcoef(AssoX,rowvar=0)
XtX=GFX.reshape(1,nFX)*GFX.reshape(nFX,1)
Vg=sum(COV*XtX)/4
Ve=Vg*(1/h2-1)

#Generate the recombination rates around the Associated SNP
# With Comeron et al. (2012)
REC_map=loadtxt('Comeron_100kb_chr'+CHR+'.txt')
REC_value=REC_map[REC_map[:,0].astype(int)/100000*100000+1==int(x/100000*100000+1),1]/(100*1000000) #since the stat was provided as cM per Mb but we want it in rec per bp

# With Fiston-Lavier & Petrov (2010)
#def REC_value(CHR,pos):
#	pos=pos/1e6
#	if CHR=="X": return((-0.01*pos**3 + 0.30*pos**2 + 1.15*pos -1.87)/(100*1000000))
#	if CHR=="2L": return((-0.01*pos**3 + 0.20*pos**2 + 2.59*pos - 1.56)/(100*1000000))
#	if CHR=="2R": return((-0.007*pos**3 + 0.35*pos**2 - 1.43*pos + 56.91)/(100*1000000))
#	if CHR=="3L": return((-0.006*pos**3 + 0.09*pos**2 + 2.94*pos - 2.90)/(100*1000000))
#	if CHR=="3R": return((-0.004*pos**3 + 0.24*pos**2 - 1.63*pos + 50.26)/(100*1000000))
#
#REC_value=REC_value(CHR,POS_X[AssoSNP])

x_asso=POS_X[AssoSNP]
Rho_X=repeat(REC_value,nSNP_X)*absolute(POS_X-x_asso)
Rho_St=repeat(REC_value,nSNP_St)*absolute(POS_St-x_asso)
Rho_Lg=repeat(REC_value,nSNP_Lg)*absolute(POS_Lg-x_asso)
Rho=concatenate((Rho_X,Rho_St,Rho_Lg))
if Graph:
	plt.plot(Rho,'bo')
	plt.title(ID_X.split(';')[0]+'(focal) / '+ID_St.split(';')[0]+'(close) / '+ID_Lg.split(';')[0]+'(Distant)')
	plt.show()

G=AssoX*GFX
Xb=G.sum(1)
e=random.normal(0,Ve**(0.5),nHaplo/2)
Y=Xb+e
#Use a standard LS fit to the data
#p_SNP=array([])
#for i in range(nSNP_X):
#	slope, intercept, r_value, p_value, std_err = stats.linregress(Geno[:,i],Y)
#	p_SNP=append(p_SNP,p_value)

for g in range(nGeneration):
	#This is where the recombination for the gametes produced by each haplotype will be modeled
	Geno=hstack((X,X_St,X_Lg))
	Geno_grow=repeat(Geno,G_rate,axis=0)
	nHaplo=Geno_grow.shape[0]
	#Random appearance of recombinations
	Geno_rec=Geno_grow
	rec=random.binomial(1,Rho,(nHaplo,nSNP_X+nSNP_St+nSNP_Lg))
	rec=where(rec==1)
	for i in range(len(rec[0])):
		if rec[0][i]<nHaplo/2: 
			Geno_rec[rec[0][i],:]=concatenate((Geno_grow[rec[0][i],:rec[1][i]],Geno_grow[rec[0][i]+nHaplo/2,rec[1][i]:]))
		else :
			Geno_rec[rec[0][i],:]=concatenate((Geno_grow[rec[0][i],:rec[1][i]],Geno_grow[rec[0][i]-nHaplo/2,rec[1][i]:]))
	if Verbose: print("There has been "+str(len(rec[0]))+" recombinations")
	#Random mating of Haplo/gametes
	Mate=random.permutation(nHaplo/2)
	Geno_mate=vstack((Geno_rec[:nHaplo/2,:],Geno_rec[nHaplo/2:,:][Mate,:]))
	#Generating the phenotypes based on the variance components
	Geno=Geno_mate[:nHaplo/2,:]+Geno_mate[nHaplo/2:,:]
	AssoX=Geno[:,AssoSNP]
	G=AssoX*GFX
	Xb=G.sum(1)
	e=random.normal(0,Ve**(0.5),nHaplo/2)
	Y=Xb+e
	#Selecting the best Sel %
	if regime[g]=="Selection":
		Kut=percentile(Y,100-(100.*Sel))
		Geno_mask=[concatenate((Y,Y))>Kut][0]
	else:
		Geno_mask=sort(random.choice(nHaplo/2,Sel*nHaplo/2))
		Geno_mask=concatenate((Geno_mask,Geno_mask+nHaplo/2))
	X=Geno_mate[Geno_mask,:nSNP_X]
	X_St=Geno_mate[Geno_mask,nSNP_X:(nSNP_X+nSNP_St)]
	X_Lg=Geno_mate[Geno_mask,(nSNP_X+nSNP_St):]
	Geno=hstack((X,X_St,X_Lg))
	nHaplo=Geno.shape[0]
	freq_end=Geno.sum(axis=0)/float(nHaplo)
	freq=vstack((freq,freq_end))

Haplo=unique(apply_along_axis(array_str,1,Geno))
if Verbose: 
	print("After "+str(nGeneration)+" generation(s), "+str(Haplo.shape[0])+" distinct haplotypes are present in the population")
	print("The associated allele reached a frequency of "+str(freq[nGeneration-1,AssoSNP]))

if Graph:
	for i in range(nSNP_X):
		plt.plot(freq[:,i],color='r')

	for i in range(nSNP_X,nSNP_X+nSNP_St):
		plt.plot(freq[:,i],color='b')

	for i in range(nSNP_X+nSNP_St,nSNP_X+nSNP_St+nSNP_Lg):
		plt.plot(freq[:,i],color='c')

	plt.plot(freq[:,AssoSNP],color='k')
	plt.title(ID_X.split(';')[0]+'(focal) / '+ID_St.split(';')[0]+'(close) / '+ID_Lg.split(';')[0]+'(Distant)')
	plt.show()

def FET(i):
	table=[[int(sum([Geno_ini[:,i]==0])*Coverage[0,i]/float(nHaplo_ini)),int(sum([Geno_ini[:,i]==1])*Coverage[0,i]/float(nHaplo_ini))],
	       [int(sum([Geno[:,i]==0])*Coverage[1,i]/float(nHaplo)),int(sum([Geno[:,i]==1])*Coverage[1,i]/float(nHaplo))]]
	return(stats.fisher_exact(table,alternative='two-sided')[1])

test=array(())
for i in range(nSNP_X+nSNP_St+nSNP_Lg): test=insert(test,len(test),FET(i))

freq_AssoSNP=freq[0,AssoSNP]
p_AssoSNP=test[AssoSNP]
nAsso_X=sum([test[:nSNP_X]<FET_thresh])

nAsso_St=sum([test[nSNP_X:(nSNP_X+nSNP_St)]<FET_thresh])

nAsso_Lg=sum([test[(nSNP_X+nSNP_St):]<FET_thresh])

OUT=CHR+','+str(x)+','+ID_X+','+str(freq_AssoSNP[0])+','+str(p_AssoSNP[0])+','+str(nSNP_X)+','+str(nAsso_X)+','+gene_St[1]+','+ID_St+','+str(nSNP_St)+','+str(nAsso_St)+','+gene_Lg[1]+','+ID_Lg+','+str(nSNP_Lg)+','+str(nAsso_Lg)

print(OUT)


