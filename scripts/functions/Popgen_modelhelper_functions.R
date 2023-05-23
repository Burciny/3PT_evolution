
source("PG_model_functions.R")


##### Function to estimate Mutation rate matrix from site frequency spectrum 

Logl_EandBPrime <- function(ePrimebPrime, data){
  
  M <- dim(data)[1]-1
  y <- 1:(M - 1)
  # Polymorphics
  L_ij_vec <- data[2:M,1:6]
  L_AC_vec <- L_ij_vec$AC
  L_AG_vec <- L_ij_vec$AG
  L_CT_vec <- L_ij_vec$CT
  L_GT_vec <- L_ij_vec$GT
  
  ePrime <- ePrimebPrime[1]
  bPrime <- ePrimebPrime[2]
  logL_eb <- 
    sum((L_AC_vec + L_GT_vec[M-y])*log(ePrime/(M-y) + (1-bPrime)/y) +
          (L_AG_vec + L_CT_vec[M-y])*log((1-ePrime)/(M-y) + bPrime/y))
  return(-logL_eb)
}


allrate_est<- function(df, region="name", chr="chr") {
  results<- matrix(0,1,6)
  colnames(results)<- c("a","f", "b","d","c","e") # Parameter names see Vogl et al 2020

  M <- dim(df)[1]-1		# number of individuals 
  HM <- sum(1/(1:(M-1)))
  y <- 1:(M-1)
  oneOverYSum <- 1/(M-y) + 1/y
  
  # Monomorphics
  L_A <- df$AC[1]
  L_C <- df$CG[1]
  L_G <- df$GT[1]
  L_T <- df$GT[M+1]
  L_i <- c(L_A, L_C, L_G, L_T)
  # Polymorphics
  L_ij_vec <- df[2:M,1:6]
  L_ij <- colSums(L_ij_vec)
  L_AC_vec <- L_ij_vec$AC
  L_AG_vec <- L_ij_vec$AG
  L_AT_vec <- L_ij_vec$AT
  L_CG_vec <- L_ij_vec$CG
  L_CT_vec <- L_ij_vec$CT
  L_GT_vec <- L_ij_vec$GT
  L_AC <- L_CA <- sum(L_AC_vec)
  L_AG <- L_GA <- sum(L_AG_vec)
  L_AT <- L_TA <- sum(L_AT_vec)
  L_CG <- L_GC <- sum(L_CG_vec)
  L_CT <- L_TC <- sum(L_CT_vec)
  L_GT <- L_TG <- sum(L_GT_vec)
  # Totals
  L_m <- sum(L_i)
  L_p <- sum(L_ij)
  L_total <- L_m + L_p
  
  ###
  L_overATCG <- L_AC + L_AG + L_TC + L_TG
  L_overAT <- L_A + L_T + L_AT
  L_overCG <- L_C + L_G + L_CG
  #
  denominatorAT <- L_overAT + 0.5*L_overATCG
  denominatorCG <- L_overCG + 0.5*L_overATCG
  
  eandbPrime<- optim(c(0.5,0.5), Logl_EandBPrime, data=df, method = "Nelder-Mead",control=list(reltol=1.e-16))$par
  
  ePrime <- eandbPrime[1]								
  bPrime <- eandbPrime[2]	
  
  results[1,1]<- L_AT/denominatorAT/HM # a 
  results[1,2]<- L_CG/denominatorCG/HM # f
  
  results[1,3]<- bPrime*L_overATCG/denominatorCG/2/HM  # b: CT and GA (Transition)
  results[1,4]<- (1-bPrime)*L_overATCG/denominatorCG/2/HM # d: GT and CA
  
  results[1,5]<- (1-ePrime)*L_overATCG/denominatorAT/2/HM # c: AG and TC (Transition)
  results[1,6]<- ePrime*L_overATCG/denominatorAT/2/HM # e: AC and TG

  return(results)
}

######## Function to create joint count matrix from consensus sequence as input (df)

matitoj<- function(df, i=1, j=2) {
  pos<- c(i,j)
  positoj<- filter(df, position %in% pos)
  positojmat<-matrix(0,4,4, dimnames = list(c("A","T","G","C"),c("A","T","G","C")))
  
  pos1A<- which(positoj$Seq=="A" & positoj$position==i)
  positojmat[1,1]<- sum(positoj$Seq[pos1A+1]=="A")
  positojmat[1,2]<- sum(positoj$Seq[pos1A+1]=="T")
  positojmat[1,3]<- sum(positoj$Seq[pos1A+1]=="G")
  positojmat[1,4]<- sum(positoj$Seq[pos1A+1]=="C")
  
  pos1T<- which(positoj$Seq=="T" & positoj$position==i)
  positojmat[2,1]<- sum(positoj$Seq[pos1T+1]=="A")
  positojmat[2,2]<- sum(positoj$Seq[pos1T+1]=="T")
  positojmat[2,3]<- sum(positoj$Seq[pos1T+1]=="G")
  positojmat[2,4]<- sum(positoj$Seq[pos1T+1]=="C")
  
  pos1G<- which(positoj$Seq=="G" & positoj$position==i)
  positojmat[3,1]<- sum(positoj$Seq[pos1G+1]=="A")
  positojmat[3,2]<- sum(positoj$Seq[pos1G+1]=="T")
  positojmat[3,3]<- sum(positoj$Seq[pos1G+1]=="G")
  positojmat[3,4]<- sum(positoj$Seq[pos1G+1]=="C")
  
  pos1C<- which(positoj$Seq=="C" & positoj$position==i)
  positojmat[4,1]<- sum(positoj$Seq[pos1C+1]=="A")
  positojmat[4,2]<- sum(positoj$Seq[pos1C+1]=="T")
  positojmat[4,3]<- sum(positoj$Seq[pos1C+1]=="G")
  positojmat[4,4]<- sum(positoj$Seq[pos1C+1]=="C")
  
  return(positojmat)
}

### Function to calculate base composition per position from consensus sequence (df)

basecount_f<- function(df, l=10) { ## l is the length of the region
  basecount_perposition<-matrix(0,l,4, dimnames = list(c(),c("A","T","G","C")))
  position<- 1:l
  basecount_perposition<- cbind(position, basecount_perposition)
  for (i in 1:l) {
    basecount_perposition[i,2]<-sum(df$position==i & df$Seq=="A") 
    basecount_perposition[i,3]<-sum(df$position==i & df$Seq=="T") 
    basecount_perposition[i,4]<-sum(df$position==i & df$Seq=="G") 
    basecount_perposition[i,5]<-sum(df$position==i & df$Seq=="C") 
  }
  return(basecount_perposition)
}

#### Sample introns 

sample_3PT<- function(df) {  ## input countfile
  intron_count<-length(unique(df$name))
  names<-unique(df$name)
  samp_intron<- sample(names, intron_count, replace = TRUE)
  idx<-c()
  for (i in 1:length(samp_intron)) {
    idx<- append(idx,which(df$name==samp_intron[i]))
  }
  new_samp<- df[idx,]
  return(new_samp)
}


##### Bootstrap gamma under model HIII
## Inputs
# df: consensus sequence, matT0: mutation rate matrix inferred from 5LR
# vgMono, gAG starting gamma values for monomers and AG dimer, respectively 

bootstrap_HIII<- function(df, BS=1000,vgMono, gAG, matT0) {   
  bsres<- matrix(nrow=9*BS, ncol=8)
  colnames(bsres)<- c("rep","Focal_pos","sA","sT","sG", "sC","sAG", "ll")
  bsres<- as.data.frame(bsres)
  techrep<- rep(1:BS, each=9)
  bsres[,1]<- techrep
  for (b in 1:BS) {
    start<-seq(1,9*BS, 9)
    new_samp<- sample_3PT(df)
    pos1to2_3PT<- matitoj(new_samp, i=1, j=2)
    pos2to3_3PT<- matitoj(new_samp, i=2, j=3)
    pos3to4_3PT<- matitoj(new_samp, i=3, j=4)
    pos4to5_3PT<- matitoj(new_samp, i=4, j=5)
    pos5to6_3PT<- matitoj(new_samp, i=5, j=6)
    pos6to7_3PT<- matitoj(new_samp, i=6, j=7)
    pos7to8_3PT<- matitoj(new_samp, i=7, j=8)
    pos8to9_3PT<- matitoj(new_samp, i=8, j=9)
    pos9to10_3PT<- matitoj(new_samp, i=9, j=10)
    
    jointCnt_list<-list(pos1to2_3PT, pos2to3_3PT,pos3to4_3PT,pos4to5_3PT,
                        pos5to6_3PT,pos6to7_3PT,pos7to8_3PT,pos8to9_3PT,
                        pos9to10_3PT)
    
    baseCnt_samp<- basecount_f(new_samp)
    bsres[start[b]:(start[b]+8),2:8]<-Sperposition_III(jointCnt_list, baseCnt_samp, vgMono, g_AG, matT0)[1:9,1:7]
    
    
  }
  return(bsres)
}
