setwd("~/Documents/Project/3prime_ms/JEB/Repository/scripts/functions/")


eigenSystem <- function(mP) {
  es <- eigen(mP)
  mBackwards <- es$vectors
  N <- length(mP[1, ])
  for(i in 1:length(mP[ ,1])) {
    mBackwards[ ,i] <- N * mBackwards[ ,i] / sum(sqrt(mBackwards[ ,i] * mBackwards[ ,i]))
  }
  mForwards <- solve(mBackwards)
  vLambda <- es$values
  list("vLambda" = vLambda, "mForwards" = mForwards, "mBackwards" = mBackwards)
}

## calculates fixation probability times N
fixp=function(gamma){
  retval=ifelse(gamma==0,1,gamma/(1-exp(-gamma)))
  return(retval)
}

##### One selection coefficient for every position 
Simulate_circular_seq_oneselcoeff<- function(L,tMax, matT0, pi, vgMono, vgAG) {
  vl=sample(x=4, size=L, replace = TRUE, prob = pi) ## sample from stationary distribution of mutation matrix (loci independent)
  vTimes=rep(0,tMax+1) # output times
  matRes=matrix(nrow=tMax+1,ncol=L) # output bases at the L positions
  matProbMutSel=matrix(nrow=L,ncol=4) ## four bases per position
  matRes[1,]=vl
  
  
  for(t in 1:tMax){
    for(l in 1:L){
      prv=ifelse(l>1,l-1,L) # Circular nucleotide sequence
      nxt=ifelse(l<L,l+1,1)
      ## build up selection vector depending on the base at the position
      vSel=-vgMono[vl[l]]+vgMono 
      ## build up selection vector depending on the preceding and following base
      if(vl[prv]==1){ ## If preceding base is A 
        if(vl[l]==3){ ## when focal base is G
          vSel=-vgAG+vSel
        }
        else{
          vSel[3]=vgAG+vSel[3]
        }
      }
      if(vl[nxt]==3){ ## If following base is G
        if(vl[l]==1){ ## when focal base is A
          vSel=-vgAG+vSel
        }
        else{
          vSel[1]=vgAG+vSel[1]
        }
      }
      matProbMutSel[l,]=fixp(vSel)*matT0[vl[l], ]## calculate selection-drift vector
      matProbMutSel[l,vl[l]]=0 ## calculate selection-drift vector 
    }
    totalRate=sum(matProbMutSel)
    vTimes[t]=rexp(1,totalRate) # sample time of change
    l=sample(x=L,size=1,prob=apply(matProbMutSel/totalRate,1,sum)) # sample position
    j=sample(x=4,size=1,prob=matProbMutSel[l,]/sum(matProbMutSel[l,])) # sample base at position
    vl[l]=j
    matRes[t+1,]=vl
  }
  
  res<- matrix(nrow=L, ncol=5)
  colnames(res)<- c("Position","A","T","G", "C")
  res[,1]<- 1:L
  
  for (i in 1:L) {
    tmp=c(sum(vTimes[matRes[,i]==1]),
          sum(vTimes[matRes[,i]==2]),
          sum(vTimes[matRes[,i]==3]),
          sum(vTimes[matRes[,i]==4]))
    res[i,2:5]<- tmp/sum(tmp)
  }
  res<-as.data.frame(res)
  data=list(vTimes=vTimes, matRes=matRes, basefreq=res)
  return(data)
}

##### ##### selection coefficient different for every position 

Simulate_circular_seq<- function(L,tMax, matT0, pi, vgMono, vgAG) {
  vl=sample(x=4, size=L, replace = TRUE, prob = pi) ## sample from stationary distribution of mutation matrix (loci independent)
  vTimes=rep(0,tMax+1) # output times
  matRes=matrix(nrow=tMax+1,ncol=L) # output bases at the L positions
  matProbMutSel=matrix(nrow=L,ncol=4) ## four bases per position
  matRes[1,]=vl
  
  
  for(t in 1:tMax){
    for(l in 1:L){
      prv=ifelse(l>1,l-1,L) # Circular nucleotide sequence
      nxt=ifelse(l<L,l+1,1)
      ## build up selection vector depending on the base at the position
      ## Take respective sel. coeff at the position
      vgMono_tmp<- as.numeric(vgMono[l,])
      vgAG_tmp<- as.numeric(vgAG[l,])
      
      vSel=-vgMono_tmp[vl[l]]+vgMono_tmp 
      ## build up selection vector depending on the preceding and following base
      if(vl[prv]==1){ ## If preceding base is A 
        if(vl[l]==3){ ## when focal base is G
          vSel=-vgAG_tmp+vSel
        }
        else{
          vSel[3]=vgAG_tmp+vSel[3]
        }
      }
      if(vl[nxt]==3){ ## If following base is G
        if(vl[l]==1){ ## when focal base is A
          vSel=-vgAG_tmp+vSel
        }
        else{
          vSel[1]=vgAG_tmp+vSel[1]
        }
      }
      matProbMutSel[l,]=fixp(vSel)*matT0[vl[l], ]## calculate selection-drift vector
      matProbMutSel[l,vl[l]]=0 ## calculate selection-drift vector 
    }
    totalRate=sum(matProbMutSel)
    vTimes[t]=rexp(1,totalRate) # sample time of change
    l=sample(x=L,size=1,prob=apply(matProbMutSel/totalRate,1,sum)) # sample position
    j=sample(x=4,size=1,prob=matProbMutSel[l,]/sum(matProbMutSel[l,])) # sample base at position
    vl[l]=j
    matRes[t+1,]=vl
  }
  
  burnin=1001
  tMaxp1=tMax+1
  
  vTimes_ss<- vTimes[burnin:tMaxp1]
  
  res<- matrix(nrow=L, ncol=5)
  colnames(res)<- c("Position","A","T","G", "C")
  res[,1]<- 1:L
  
  for (i in 1:L) {
    tmp=c(sum(vTimes_ss[matRes[burnin:tMaxp1,i]==1]),
          sum(vTimes_ss[matRes[burnin:tMaxp1,i]==2]),
          sum(vTimes_ss[matRes[burnin:tMaxp1,i]==3]),
          sum(vTimes_ss[matRes[burnin:tMaxp1,i]==4]))
    res[i,2:5]<- tmp/sum(tmp)
  }
  res<-as.data.frame(res)
  data=list(vTimes=vTimes, matRes=matRes, basefreq=res)
  return(data)
}

#### For visualization

mattodf<- function(joint) {
  Pos_i<- rep(c("A","T","G","C"),4)
  Pos_j<- rep(c("A","T","G","C"),each=4)
  Value<-c()
  
  for (i in 1:4) {   # Column  index
    for (j in 1:4) { # Row index
      Value<-append(Value, joint[j,i])
    }
  }
  JointCnt_df<-as.data.frame(cbind(Pos_i, Pos_j,Value))
  JointCnt_df$Value<- as.numeric(as.character(JointCnt_df$Value))
  
  return(JointCnt_df)
}

########

Positoj_alldata<- function(jointCnt_files, k=1, matRes_I, matRes_II, matRes_III, vTimes_I, vTimes_II, vTimes_III) {
  jntCnt<- read.delim(jointCnt_files[k], row.names = 1)
  
  JointCnt_df_itoj<- mattodf(jntCnt)
  JointCnt_df_itoj$Pos_i<- factor(JointCnt_df_itoj$Pos_i, levels = c("A","T","G","C"))
  JointCnt_df_itoj$Pos_j<- factor(JointCnt_df_itoj$Pos_j, levels = c("A","T","G","C"))
  
  burnin=1001
  tMaxp1=  nrow(vTimes_I)
  
  
  ##HI: joint probabilities at position i and j
  vTimes_I_ss<- vTimes_I$X1[burnin:tMaxp1]
  
  mRes_I=matrix(nrow=4,ncol=4)
  for(i in 1:4){
    for(j in 1:4){
      mRes_I[i,j]=sum(vTimes_I_ss[matRes_I[burnin:tMaxp1,k]==i & matRes_I[burnin:tMaxp1,(k+1)]==j])
    }
  }
  mRescount_I<-mRes_I/sum(mRes_I)*sum(JointCnt_df_itoj$Value)  ## to match total empirical counts
  rownames(mRescount_I)=c('A','T','G','C')
  colnames(mRescount_I)=c('A','T','G','C')
  mRescount_I<- as.matrix(mRescount_I)
  
  mRescount_I_df_itoj<- mattodf(mRescount_I)
  mRescount_I_df_itoj$Pos_i<- factor(mRescount_I_df_itoj$Pos_i, levels = c("A","T","G","C"))
  mRescount_I_df_itoj$Pos_j<- factor(mRescount_I_df_itoj$Pos_j, levels = c("A","T","G","C"))
  
  
  ### II
  vTimes_II_ss<- vTimes_II$X1[burnin:tMaxp1]
  
  mRes_II=matrix(nrow=4,ncol=4)
  for(i in 1:4){
    for(j in 1:4){
      mRes_II[i,j]=sum(vTimes_II_ss[matRes_II[burnin:tMaxp1,k]==i & matRes_II[burnin:tMaxp1,(k+1)]==j])
    }
  }
  mRescount_II<-mRes_II/sum(mRes_II)*sum(JointCnt_df_itoj$Value)  ## to match total empirical counts
  rownames(mRescount_II)=c('A','T','G','C')
  colnames(mRescount_II)=c('A','T','G','C')
  mRescount_II<- as.matrix(mRescount_II)
  
  mRescount_II_df_itoj<- mattodf(mRescount_II)
  mRescount_II_df_itoj$Pos_i<- factor(mRescount_II_df_itoj$Pos_i, levels = c("A","T","G","C"))
  mRescount_II_df_itoj$Pos_j<- factor(mRescount_II_df_itoj$Pos_j, levels = c("A","T","G","C"))
  
  ### III
  vTimes_III_ss<- vTimes_III$X1[burnin:tMaxp1]
  
  mRes_III=matrix(nrow=4,ncol=4)
  for(i in 1:4){
    for(j in 1:4){
      mRes_III[i,j]=sum(vTimes_III_ss[matRes_III[burnin:tMaxp1,k]==i & matRes_III[burnin:tMaxp1,(k+1)]==j])
    }
  }
  mRescount_III<-mRes_III/sum(mRes_III)*sum(JointCnt_df_itoj$Value)  ## to match total empirical counts
  rownames(mRescount_III)=c('A','T','G','C')
  colnames(mRescount_III)=c('A','T','G','C')
  mRescount_III<- as.matrix(mRescount_III)
  
  mRescount_III_df_itoj<- mattodf(mRescount_III)
  mRescount_III_df_itoj$Pos_i<- factor(mRescount_III_df_itoj$Pos_i, levels = c("A","T","G","C"))
  mRescount_III_df_itoj$Pos_j<- factor(mRescount_III_df_itoj$Pos_j, levels = c("A","T","G","C"))
  
  
  JointCnt_df_itoj$Valuefreq<- JointCnt_df_itoj$Value/sum(JointCnt_df_itoj$Value)
  mRescount_I_df_itoj$Valuefreq<- mRescount_I_df_itoj$Value/sum(mRescount_I_df_itoj$Value)
  mRescount_II_df_itoj$Valuefreq<- mRescount_II_df_itoj$Value/sum(mRescount_II_df_itoj$Value)
  mRescount_III_df_itoj$Valuefreq<- mRescount_III_df_itoj$Value/sum(mRescount_III_df_itoj$Value)
  
  positoj<-rbind(JointCnt_df_itoj, mRescount_I_df_itoj,mRescount_II_df_itoj,mRescount_III_df_itoj)
  positoj$data<- c(rep("Empirical",16),rep("HI",16),rep("HII",16),rep("HIII",16))
  
  ######
  
  DeviationI<-(positoj$Valuefreq[positoj$data=="Empirical"]-positoj$Valuefreq[positoj$data=="HI"])^2/positoj$Valuefreq[positoj$data=="HI"]
  DeviationII<-(positoj$Valuefreq[positoj$data=="Empirical"]-positoj$Valuefreq[positoj$data=="HII"])^2/positoj$Valuefreq[positoj$data=="HII"]
  DeviationIII<-(positoj$Valuefreq[positoj$data=="Empirical"]-positoj$Valuefreq[positoj$data=="HIII"])^2/positoj$Valuefreq[positoj$data=="HIII"]
  
  Deviation<- c(positoj$Valuefreq[positoj$data=="Empirical"], DeviationI,DeviationII,DeviationIII)
  
  Deviationitoj<-cbind(positoj[,c(1,2,5)], Deviation)
  Deviationitoj_r<-Deviationitoj[17:64,]
  
  
  res<- list(positoj=positoj, deviation=Deviationitoj_r)
  
  return(res)
}