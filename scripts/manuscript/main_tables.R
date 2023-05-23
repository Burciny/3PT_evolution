setwd(dirname(rstudioapi::getSourceEditorContext()$path))



################ ################ ################ ################ ################ 
################ TABLE 1 ################ ################ ################ ################ 

source("../functions/GCbias_functions.R")

SFSauto_mel196 <- read_delim("../../data/SFSauto_mel196", 
                             "\t", escape_double = FALSE, trim_ws = TRUE)
SFSX_mel196 <- read_delim("../../data/SFSX_mel196", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)

spectras_auto=SFSauto_mel196[,7:8]
spectras_X= SFSX_mel196[,7:8]


############ >>> For GC changing sites <<< ############

B_melA_SI<-optim(0.5, logl,data=spectras_auto,c=0, method = "BFGS")$par
B_melA_SI_ll<- -optim(0.5, logl,data=spectras_auto,c=0, method = "BFGS")$value

B_melX_SI<-optim(0.5, logl,data=spectras_X,c=0, method = "BFGS")$par
B_melX_SI_ll<- -optim(0.5, logl,data=spectras_X,c=0, method = "BFGS")$value


############ >>> For GC conservative sites <<< ############

B_melA_SI_GCcons<-optim(0.5, logl_neutral,data=spectras_auto,c=0, method = "BFGS")$par
B_melA_SI_GCcons_ll<- -optim(0.5, logl_neutral,data=spectras_auto,c=0, method = "BFGS")$value

B_melX_SI_GCcons<-optim(0.5, logl_neutral,data=spectras_X,c=0, method = "BFGS")$par
B_melX_SI_GCcons_ll<- -optim(0.5, logl_neutral,data=spectras_X,c=0, method = "BFGS")$value


#### Pool the data

SFSpool_mel196<- SFSauto_mel196+SFSX_mel196
spectras_pool<- SFSpool_mel196[,7:8]

## GC changing 

B_mel_pooled<-optim(0.5, logl,data=spectras_pool,c=0, method = "BFGS")$par
B_mel_pooled_ll<- -optim(0.5, logl,data=spectras_pool,c=0, method = "BFGS")$value

## GC conservative 

B_mel_pooled_GCcons<-optim(0.5, logl_neutral,data=spectras_pool,c=0, method = "BFGS")$par
B_mel_pooled_GCcons_ll<- -optim(0.5, logl_neutral,data=spectras_pool,c=0, method = "BFGS")$value


###### Likelihoods for B=0
### GC changing spectra
## Autosome

ll_mel196_A_GCchang_null0<- ll_fun(B=0, spectras_auto,c=0,s=2)

## X

ll_mel196_X_GCchang_null0<- ll_fun(B=0, spectras_X,c=0,s=2)
## Pool

ll_mel196_pool_GCchang_null0<- ll_fun(B=0, spectras_pool,c=0,s=2)

### GC conservative spectra
## Autosome

ll_mel196_A_GCcons_null0<- ll_fun(B=0, spectras_auto,c=0,s=1)

## X

ll_mel196_X_GCcons_null0<- ll_fun(B=0, spectras_X,c=0,s=1)

## Pool

ll_mel196_pool_GCcons_null0<- ll_fun(B=0, spectras_pool,c=0,s=1)


###########
Spectra<- rep(c("GCchang","GCcons"),3)
Chromosome<- c(rep("Autosome",2),rep("X",2),rep("Pool",2))
B<- c(B_melA_SI,B_melA_SI_GCcons,B_melX_SI, B_melX_SI_GCcons, B_mel_pooled, B_mel_pooled_GCcons)
loglike<- c(B_melA_SI_ll,B_melA_SI_GCcons_ll,B_melX_SI_ll, B_melX_SI_GCcons_ll, B_mel_pooled_ll, B_mel_pooled_GCcons_ll)
nullloglike<- c(ll_mel196_A_GCchang_null0,ll_mel196_A_GCcons_null0,
                ll_mel196_X_GCchang_null0,ll_mel196_X_GCcons_null0,
                ll_mel196_pool_GCchang_null0,ll_mel196_pool_GCcons_null0)

GCbias_summary<- as.data.frame(cbind(Spectra, Chromosome, B, loglike, nullloglike))
sapply(GCbias_summary, class)
GCbias_summary[,3:5]<- lapply(GCbias_summary[,3:5], function(x) {as.numeric(as.character(x))})

################ ################ ################ ################ ################ 
################ TABLE 2 ################ ################ ################ ################ 

nuc_count_5p <- read_delim("../../data/nuc_count_5p", 
                           "\t", escape_double = FALSE, trim_ws = TRUE)

nuc_count_3p <- read_delim("../../data/nuc_count_3p", 
                           "\t", escape_double = FALSE, trim_ws = TRUE)


nuc_count_5p[1,]/sum(nuc_count_5p[1,])
nuc_count_3p[1,]/sum(nuc_count_3p[1,])

### Autosome 
S_TA_5p_A<- (nuc_count_5p$T[2]-nuc_count_5p$A[2])/(nuc_count_5p$T[2]+nuc_count_5p$A[2])
S_TA_3p_A<- (nuc_count_3p$T[2]-nuc_count_3p$A[2])/(nuc_count_3p$T[2]+nuc_count_3p$A[2])

S_CG_5p_A<- (nuc_count_5p$C[2]-nuc_count_5p$G[2])/(nuc_count_5p$C[2]+nuc_count_5p$G[2])
S_CG_3p_A<- (nuc_count_3p$C[2]-nuc_count_3p$G[2])/(nuc_count_3p$C[2]+nuc_count_3p$G[2])



