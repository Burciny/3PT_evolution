

### Site frequency spectra build from population data (i.e, Sample_alignments/count_autosome_wn)

spectra6_build<- function(df,M) {
  polfor_AC<-df[which(df$A!=0 & df$C!=0),]
  polfor_AG<-df[which(df$A!=0 & df$G!=0),]
  polfor_AT<-df[which(df$A!=0 & df$T!=0),]
  polfor_CG<-df[which(df$C!=0 & df$G!=0),]
  polfor_CT<-df[which(df$C!=0 & df$T!=0),]
  polfor_GT<-df[which(df$G!=0 & df$T!=0),]
  monofor_A<-df[which(df$State=="M" & df$A!=0),]
  monofor_T<-df[which(df$State=="M" & df$T!=0),]
  monofor_G<-df[which(df$State=="M" & df$G!=0),]
  monofor_C<-df[which(df$State=="M" & df$C!=0),]
  spectra<- matrix(0,M+1,6) 
  colnames(spectra)<- c("AC","AG","AT","CG","CT","GT")
  rownames(spectra)<- c(0:M)
  # Place monomorphic
  spectra[1,1]<- nrow(monofor_A)
  spectra[1,2]<- nrow(monofor_A)
  spectra[1,3]<- nrow(monofor_A)
  spectra[1,4]<- nrow(monofor_C)
  spectra[1,5]<- nrow(monofor_C)
  spectra[1,6]<- nrow(monofor_G)
  spectra[M+1,1]<- nrow(monofor_C)
  spectra[M+1,2]<- nrow(monofor_G)
  spectra[M+1,3]<- nrow(monofor_T)
  spectra[M+1,4]<- nrow(monofor_G)
  spectra[M+1,5]<- nrow(monofor_T)
  spectra[M+1,6]<- nrow(monofor_T)
  # Place polymorphic
  for (i in 1:(M-1)) {
    spectra[i+1,1]<- sum(polfor_AC$C==i)
    spectra[i+1,2]<- sum(polfor_AG$G==i)
    spectra[i+1,3]<- sum(polfor_AT$T==i)
    spectra[i+1,4]<- sum(polfor_CG$G==i)
    spectra[i+1,5]<- sum(polfor_CT$T==i)
    spectra[i+1,6]<- sum(polfor_GT$T==i)
  }
  
  spectra<- as.data.frame(spectra)
  spectra<- as.data.frame(spectra)
  spectra$Symm=spectra$AT+spectra$CG
  spectra$Asymm=spectra$AC+spectra$AG+
    spectra$CT[(M+1):1]+spectra$GT[(M+1)]
  
  return(spectra)
}

### Log likelihood functions to infer B (strength of gBGC) from site frequency spectrum
logl<- function(B, data, c=3) {
  M=nrow(data)-1
  y=((1+c):(M-1-c))
  C=nrow(data[(2+c):(M-c),2])
  # For GC changing spectra 
  Ly=data[(2+c):(M-c),2]
  
  
  p1<- sum((Ly*B*y)/M)
  p2<- sum(Ly*log(M/(y*(M-y))))
  p3<- sum(Ly)*log(sum(exp(B*y/M)*(M/(y*(M-y)))))
  logl<- p1+p2-p3
  return(-logl)
}

logl_neutral<- function(B, data, c=3) {
  M=nrow(data)-1
  y=((1+c):(M-1-c))
  C=nrow(data[(2+c):(M-c),1])
  # For GC conservative spectra 
  Ly=data[(2+c):(M-c),1]
  
  
  p1<- sum((Ly*B*y)/M)
  p2<- sum(Ly*log(M/(y*(M-y))))
  p3<- sum(Ly)*log(sum(exp(B*y/M)*(M/(y*(M-y)))))
  logl<- p1+p2-p3
  return(-logl)
}

dflogl<- function(B, data, c=3) {
  M=nrow(data)-1
  y=((1+c):(M-1-c))
  C=nrow(data[(2+c):(M-c),2])
  # For GC changing spectra 
  Ly=data[(2+c):(M-c),2]
  
  p1<- sum(Ly*y/M)
  p2<- sum(Ly)
  p3<- sum((y/M)*exp(B*y/M)*(M/(y*(M-y))))/sum(exp(B*y/M)*(M/(y*(M-y))))
  dflogl<- p1-p2*p3
  return(-dflogl)
}

ll_fun<- function(B, data, c=3, s=1) {
  M=nrow(data)-1
  y=((1+c):(M-1-c))
  C=nrow(data[(2+c):(M-c),2])
  # For GC changing spectra 
  Ly=data[(2+c):(M-c),s]  # s=1 For GCcons spectra, s=2 for GCchang spectra
  
  
  p1<- sum((Ly*B*y)/M)
  p2<- sum(Ly*log(M/(y*(M-y))))
  p3<- sum(Ly)*log(sum(exp(B*y/M)*(M/(y*(M-y)))))
  logl<- p1+p2-p3
  return(logl)
}
