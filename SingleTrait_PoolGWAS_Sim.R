if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("snpStats", force = T)
install.packages('PhenotypeSimulator')

library(PhenotypeSimulator)
library(dplyr)
library(readr)
library(ROCR)
library(glmnet)
library(stabs)

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
###(2)

cls_genotype=function(clusters, genotype){
  k=length(clusters)
  l=list()
  for(i in 1:k){
    l[[i]]= genotype %>% filter(row.names(genotype) %in% row.names(clusters[[i]]   )      ) 
    
  } 
  
  return(l)  
}

###(3)

pool_allele_freq1=function(clusters_genotype, snp){
  k=length(clusters_genotype)
  a=rep(0,k)
  for(i in 1:k){
    a[i]=sum(clusters_genotype[[i]][,snp])/(2*(nrow(clusters_genotype[[i]]))) 
  }
  return(a) 
}

Nsim=100
A=matrix(0, nrow=Nsim, ncol=3)
for(k in 1:Nsim){
  phenotype= runSimulation(N=2000, P=1, tNrSNP=10000, cNrSNP=10,
                           SNPfrequencies=c(0.05, 0.1,0.2),
                           genVar = 0.6, h2s=0.5/0.6, theta=1, eta=1, phi=1,
                           alpha=1, seed=100+k)
  
  
  Y=as.data.frame(phenotype$phenoComponentsFinal$Y)
  X=as.data.frame(phenotype[["rawComponents"]][["genotypes"]][["genotypes"]] ) 
  clusters=cls(Y,10)
  causal_snps=colnames(phenotype[["phenoComponentsIntermediate"]][["genFixed"]][["cov"]])
  clusters_genotype=cls_genotype(clusters = clusters, genotype=X)
  
  SNPs=phenotype[["setup"]][["id_snps"]]
  allele_freq1=matrix(NA, nrow=length(clusters), ncol=length(SNPs))
  for(i in 1:length(SNPs)){
    allele_freq1[,i]=pool_allele_freq1(clusters_genotype = clusters_genotype, snp=SNPs[i])  
  }
  
  colnames(allele_freq1)=SNPs
  
  
  ind_percluster=sapply(clusters_genotype, nrow)
  clusters_genotype_re=list()
  set.seed(seed=300+k)
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
  
  n0=nrow(data_combined_re)
  data_combined_re=data_combined_re[sample(1:n0, size=n0),]
  x=data_combined_re[,2:ncol(data_regenerate)]
  y=data_combined_re$y
  
  stab.lasso=stabsel(x=x, y=y, intercept=F, fitfun=glmnet.lasso,
                     cutoff=0.9, PFER=1, assumption="none",
                     sampling.type="MB")
  
  temp=stab.lasso[["selected"]]
  pred=prediction(stab.lasso[["max"]], as.numeric( colnames(X) %in% causal_snps      ) )
  perf=performance(pred, "tpr", "fpr")
  auc=performance(pred, measure="auc")
  
  A[k,]=c(length(temp), sum(parse_number(causal_snps)%in%temp), auc@y.values[[1]]   )
  ## format =c(Total discovery, causal discovery, AUC)
  
}

A=as.data.frame(A)
colnames(A)=c('TD', 'CD', 'AUC')
write.csv(A, 'MB_0.9_1.csv', row.names = F)
