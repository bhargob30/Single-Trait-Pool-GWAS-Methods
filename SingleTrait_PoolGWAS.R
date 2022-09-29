# One-dimensional Scheme

#(1) install and load the packages.

# Note that, it is necessary to have 'snpStats' 
# to install 'PhenotypeSimulator'.
# So, if you don't have snpStats, install it using
# the following command.
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("snpStats", force = T)
install.packages('PhenotypeSimulator')

library(PhenotypeSimulator)
library(dplyr)
library(readr)
library(ROCR)
library(glmnet)
library(stabs)

#(2) Generating Phenotypes

# Generating 1-dimensional phenotype (P=1) with 3 components,
# viz, genetic effects, polygenic effects and environmental noise.
# So, no covariates, and no correlated noise.
# Total Genetic heritability (genVar)=60%, Total SNP heritability=50%,
# Note that SNP heritability is controlled by 'genVar * h2s'. 
# N=2000 individuals, 10000 total SNPs(tNrSNP), 10 causal SNPs(cNrSNP).
# for details on other arguments, use ?runSimulation.
phenotype= runSimulation(N=2000, P=1, tNrSNP=10000, cNrSNP=10,
                         SNPfrequencies=c(0.05, 0.1,0.2),mBetaGenetic = 0,
                         sdBetaGenetic =1,
                         pIndependentGenetic = 0, distBetaGenetic = "norm",
                         genVar = 0.6, h2s=0.5/0.6, theta=1, eta=1, phi=1, delta=0,
                         rho=0, alpha=1, gamma=1, seed=10)

# 'Y' below denote phenotype observations for all individuals
Y=as.data.frame(phenotype$phenoComponentsFinal$Y)

# 'genotypes' is the N by S allele frequency matrix.
# S= Total no. of SNPs
genotypes=phenotype[["rawComponents"]][["genotypes"]]
genotypes=as.data.frame(genotypes$genotypes)



#(3) clustering scheme :-division into k clusters based on ranking

# the following function will take no. of reqd clusters (k) as one of the
# arguments and return a list of clusters.
# Y should be a data frame of phenotypes with rownames as IDs and "Traits" in column
# The output will be a list of length k. The i th 
# element of the list contains phenotype observations
# corresponding to the i th cluster.
cls=function(Y, k){
  N=nrow(Y) 
  s=N %% k
  r=N %/% k
  S=c( rep(1,s), rep(0, k-s)      )
  Y_new=arrange(Y, Trait_1)
  l=list()
  for(i in 1:k){
    if(i==1){
      l[[i]] = Y_new[1:( (i*r)+sum(S[1:i]) )  , , drop=F] 
    }
    else{
      l[[i]] = Y_new[(1+((i-1)*r)+sum( S[1:(i-1)]  )):( (i*r)+sum(S[1:i]) )  , , drop=F]
    }
  }
  
  return(l) 
  
}

# making 10 clusters
clusters=cls(Y,10)




#(4) Assigning true genotypes to clusters created

# The following function will create a list of clusters.
# But this time each element of the list will contain
# the genotypes of the individuals belonging to that
# cluster.
# The argument 'clusters' should be clusters object 
# obtained as output of 'cls' function. 'genotype' 
# should be the the allele frequency matrix.

cls_genotype=function(clusters, genotype){
  k=length(clusters)
  l=list()
  for(i in 1:k){
    l[[i]]= genotype %>% filter(row.names(genotype) %in% row.names(clusters[[i]]   )      ) 
    
  } 
  
  return(l)  
}

# clusters_genotype is list of clusters of genotypes
clusters_genotype=cls_genotype(clusters = clusters, genotype=genotypes)



#(5) Getting Pooled allele frequencies (PAF)

# The following function returns the PAF of a paricular SNP
# in all the pools. The output is a vector whose ith element is PAF of 
# the snp in ith pool. Argument 'clusters_genotype' should be the 
# output of the 'cls_genotype' function. The next argument 'snp' should be
# a character string denoting the SNP, for e.g., 'SNP_111', 'SNP_176' etc. 

pool_allele_freq1=function(clusters_genotype, snp){
  k=length(clusters_genotype)
  a=rep(0,k)
  for(i in 1:k){
    a[i]=sum(clusters_genotype[[i]][,snp])/(2*(nrow(clusters_genotype[[i]]))) 
  }
  return(a) 
}

SNPs=phenotype[["setup"]][["id_snps"]] # denote the IDs of the SNPs.

# the following object 'allele_freq1' is pooled allele frequency matrix
# whose rows denote the pools, and columns denote the SNPs. 
allele_freq1=matrix(NA, nrow=length(clusters), ncol=length(SNPs))
for(i in 1:length(SNPs)){
  allele_freq1[,i]=pool_allele_freq1(clusters_genotype = clusters_genotype, snp=SNPs[i])  
}

colnames(allele_freq1)=SNPs


#(6) Regenaration

# Now coming to regenerating the genotype data from 
# these pooled allele frequencies. We compressed our genotype data to
# allele frequencies. Now we will expand again by regeneration.
# The ideal thing will be to regenerate same number of individuals
# per cluster whose genotype we combined to get the allele 
# frequency. 

ind_percluster=sapply(clusters_genotype, nrow) # denotes no. of individuals per cluster

# 'clusters_genotype_re' is a list where i th element of the list
# will be a dataframe with first column all i's (that is, the cluster number)
# and rest columns will denote regenerated allele frequencies for all the SNPs
# in cluster i 
clusters_genotype_re=list()
set.seed(466) # set this seed to have the reproducible regenerated data
for(i in 1:length(clusters)){
  a=ind_percluster[i]
  y=rep(i,a)
  
  d=matrix(NA, nrow=a, ncol=length(SNPs))
  for(j in 1:length(SNPs)){
    d[,j]= rbinom(a,2, allele_freq1[i,j])
    
  }
  
  clusters_genotype_re[[i]]=cbind(y,d)
  
}

data_regenerate=clusters_genotype_re[[1]]
for(i in 2:length(clusters_genotype_re)){
  data_regenerate=rbind(data_regenerate,clusters_genotype_re[[i]] )
}

data_regenerate=as.data.frame(data_regenerate)
colnames(data_regenerate)[2:ncol(data_regenerate)]=SNPs

# the object 'data_regenerate' is the regenerated genotype data from the pooled
# allele frequencies. The first column of this data denotes the cluster numbers.
# So we can use this data for classification purpose, say, random forest. 
# But be sure to shuffle the data before applying to random forest


#(7) Assign the available phenotypes to the regenerated genotypes in an arbitrary order

# 'clusters_genphen_re' is a list with i th element denoting the a dataframe 
# whose first column will be the phenotype observations for all the individuals
# in i th pool, and rest of the columns denoting regenarted allele freq. of all SNPs
# for i th pool. 
clusters_genphen_re=list()
for(i in 1:length(clusters_genotype_re)){
  clusters_genphen_re[[i]]=cbind( clusters[[i]]$Trait_1,  clusters_genotype_re[[i]][, 2:ncol(data_regenerate)]         )
  
}


data_combined_re=clusters_genphen_re[[1]]
for(i in 2:length(clusters_genphen_re)){
  data_combined_re=rbind(data_combined_re,clusters_genphen_re[[i]] )
}

data_combined_re=as.data.frame(data_combined_re)
colnames(data_combined_re)[1:ncol(data_regenerate)]=c('y',SNPs)

# the object 'data_combined_re' is the regenerated genotype data from the pooled
# allele frequencies. It is same as 'data_regenerate', but this time 
# The first column of this data denotes the phenotype observations.
# This data will be used for LASSO.


#(8) LASSO+Stability selection

n0=nrow(data_combined_re) #no. of observations==no. of individuals
data_combined_re=data_combined_re[sample(1:n0, size=n0),] # shuffling the rows
x=data_combined_re[,2:ncol(data_regenerate)] # the feature matrix
y=data_combined_re$y # the response

# perform lasso with stability selection.
# (*) Original stability selection 
# ( For CPSS, just change sampling.type="SS" & assumption="r-concave")

stab.lasso=stabsel(x=x, y=y, intercept=F, fitfun=glmnet.lasso,
                   cutoff=0.8, PFER=0.5, assumption="none",
                   sampling.type="MB")

temp=stab.lasso[["selected"]] # selected variables. These will be just the ID numbers of the SNPs

length(temp) # Total no of selected variables

# find out the true causal snps generated by 'PhenotypeSimulator'
causal_snps=colnames(phenotype[["phenoComponentsIntermediate"]][["genFixed"]][["cov"]])

# find out how many of the causal SNPs are selected by stability selection
sum(parse_number(causal_snps)%in%temp)

# find out the no. of false discoveries:-
length(temp)-sum(parse_number(causal_snps)%in%temp)
