setwd(dirname(rstudioapi::getSourceEditorContext()$path))

################ ################ ################ ################ ################ 
################ TABLE S1 ################ ################ ################ ################ 


nuc_count_5p <- read_delim("../../data/nuc_count_5p", 
                           "\t", escape_double = FALSE, trim_ws = TRUE)

nuc_count_3p <- read_delim("../../data/nuc_count_3p", 
                           "\t", escape_double = FALSE, trim_ws = TRUE)

nuc_count_5p[2,]/sum(nuc_count_5p[2,])
nuc_count_3p[2,]/sum(nuc_count_3p[2,])

## X 
S_TA_5p_X<- (nuc_count_5p$T[2]-nuc_count_5p$A[2])/(nuc_count_5p$T[2]+nuc_count_5p$A[2])
S_TA_3p_X<- (nuc_count_3p$T[2]-nuc_count_3p$A[2])/(nuc_count_3p$T[2]+nuc_count_3p$A[2])

S_CG_5p_X<- (nuc_count_5p$C[2]-nuc_count_5p$G[2])/(nuc_count_5p$C[2]+nuc_count_5p$G[2])
S_CG_3p_X<- (nuc_count_3p$C[2]-nuc_count_3p$G[2])/(nuc_count_3p$C[2]+nuc_count_3p$G[2])



################ ################ ################ ################ ################ 
################ TABLE S2 ################ ################ ################ ################ 

auto_gammaHI_3PT_selcoeffs <- read_delim("../../data/auto_gammaHI_3PT_selcoeffs", 
                                         delim = "\t", escape_double = FALSE, 
                                         trim_ws = TRUE)
auto_gammaHII_3PT_selcoeffs <- read_delim("../../data/auto_gammaHII_3PT_selcoeffs", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)
auto_gammaHIII_3PT_selcoeffs <- read_delim("../../data/auto_gammaHIII_3PT_selcoeffs", 
                                           delim = "\t", escape_double = FALSE, 
                                           trim_ws = TRUE)

################## Autosome 
####   I (3 params) vs II (1 params)- 2 df

LLR_IvsII<-2*(auto_gammaHI_3PT_selcoeffs$ll-auto_gammaHII_3PT_selcoeffs$ll)

pchisq(LLR_IvsII, df=2, lower.tail = FALSE)<=0.05
pval_IvsII<-pchisq(LLR_IvsII, df=2, lower.tail = FALSE)
pval_IvsII
## I shows better fit

#### III (4 params) vs I (3 params)- 1 df
LLR_IvsIII<-2*(auto_gammaHIII_3PT_selcoeffs$ll- auto_gammaHI_3PT_selcoeffs$ll)

pchisq(LLR_IvsIII, df=1, lower.tail = FALSE)<=0.05
pval_IvsIII<-pchisq(LLR_IvsIII, df=1, lower.tail = FALSE)
pval_IvsIII
## III shows better fit

#### III (4 params) vs II (1 params)- 3 df
LLR_IIvsIII<-2*(auto_gammaHIII_3PT_selcoeffs$ll- auto_gammaHII_3PT_selcoeffs$ll)

pchisq(LLR_IIvsIII, df=3, lower.tail = FALSE)<=0.05
pval_IIvsIII<-pchisq(LLR_IIvsIII, df=3, lower.tail = FALSE)
pval_IIvsIII
## III shows better fit 

################ ################ ################ ################ ################ 
################ TABLE S3 ################ ################ ################ ################ 

X_gammaHI_3PT_selcoeffs <- read_delim("../../data/X_gammaHI_3PT_selcoeffs", 
                                         delim = "\t", escape_double = FALSE, 
                                         trim_ws = TRUE)
X_gammaHII_3PT_selcoeffs <- read_delim("../../data/X_gammaHII_3PT_selcoeffs", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)
X_gammaHIII_3PT_selcoeffs <- read_delim("../../data/X_gammaHIII_3PT_selcoeffs", 
                                           delim = "\t", escape_double = FALSE, 
                                           trim_ws = TRUE)

######### X chromosome

#### I vs II - 2 df

LLR_IvsII<-2*(X_gammaHI_3PT_selcoeffs$ll-X_gammaHII_3PT_selcoeffs$ll)

pchisq(LLR_IvsII, df=2, lower.tail = FALSE)<=0.05
pval_IvsII<-pchisq(LLR_IvsII, df=2, lower.tail = FALSE)
pval_IvsII
## I shows better fit

#### I vs III- 1 df
LLR_IvsIII<-2*(X_gammaHIII_3PT_selcoeffs$ll- X_gammaHI_3PT_selcoeffs$ll)

pchisq(LLR_IvsIII, df=1, lower.tail = FALSE)<=0.05
pval_IvsIII<-pchisq(LLR_IvsIII, df=1, lower.tail = FALSE)
pval_IvsIII

## III shows better fit

#### II vs III- 3 df
LLR_IIvsIII<-2*(X_gammaHIII_3PT_selcoeffs$ll- X_gammaHII_3PT_selcoeffs$ll)

pchisq(LLR_IIvsIII, df=3, lower.tail = FALSE)<=0.05
pval_IIvsIII<-pchisq(LLR_IIvsIII, df=3, lower.tail = FALSE)
pval_IIvsIII
## III shows better fit 



