
source("PG_modelhelper_functions.R")

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

###########
### HI: Selection on monomers
###########
# Input 
## par: starting gamma parameters for monomer and AG dimer
## data: data list consisting of mutation rate matrix inferred from 5LR (matT0)
## joint frequency matrix (jointCnt) 
# and base frequency vector in each position (p_ip1)

selFunc_I=function(par,data){ 
  vgMono=par$vgMono
  g_AG=0
  jointCnt=data$jointCnt
  matSAG=matrix(nrow=4,ncol=4) ## build up two selection matrices
  p_im1=rep(0,4)
  for(i in 1:4){
    matSAG[i,i]=0
    p_im1[i]=sum(jointCnt[i,])
  }
  p_im1=p_im1/sum(p_im1)
  ## selection on monomers
  
  for(i in 1:4){
    for(j in 1:4){
      matSAG[i,j]=-vgMono[i]+vgMono[j]    
    }
  }
  ## selection on 'A' before following 'G' (appears like sel on monomer)
  matSAG[,1]= g_AG*data$p_ip1[3]+matSAG[,1]
  matSAG[1,]= -g_AG*data$p_ip1[3]+matSAG[1,]
  matS=matSAG
  ## selection for/against 'G' after 'A'
  matSAG[,3]= g_AG+matSAG[,3] 
  matSAG[3,]= -g_AG+matSAG[3,] 
  matTAG=matrix(nrow=4,ncol=4) ## build up transition matrices
  matT=matTAG
  for(i in 1:4){
    for(j in 1:4){
      matTAG[i,j]=fixp(matSAG[i,j])*data$matT0[i,j] 
      matT[i,j]=fixp(matS[i,j])*data$matT0[i,j]
    }
    matTAG[i,i]=0
    matTAG[i,i]=1-sum(matTAG[i,])
    matT[i,i]=0
    matT[i,i]=1-sum(matT[i,])
  }
  eSmTAG=eigenSystem(matTAG) ## calculate stationary distribution
  piAG=Re(eSmTAG$mForward[1,]/sum(eSmTAG$mForward[1,]))
  eSmT=eigenSystem(matT)
  pi=Re(eSmT$mForward[1,]/sum(eSmT$mForward[1,]))
  matExp=matT
  for(j in 1:4){
    matExp[1,j]=p_im1[1]*piAG[j] ## follows an 'A' => sel against 'G'
  }
  for(i in 2:4){
    for(j in 1:4){
      matExp[i,j]=p_im1[i]*pi[j] ## does no follow an 'A'
    }
  }
  return(matExp)
}

## likelihood function 
ll_I=function(par,matT0,jointCnt,p_ip1){
  parameters=list(vgMono=c(par[1:3],-sum(par[1:3])))
  matExp=selFunc_I(parameters,data=list(matT0=matT0,p_ip1=p_ip1, jointCnt=jointCnt))
  return(-sum(log(matExp)*jointCnt)) 
}


##### 
### jointCount_files: list of files with joint frequency matrices of
## the four bases at position i and i + 1

Sperposition_I<- function(jointCount_files, basefreq, vgMono,g_AG=0, matT0) {
  res<- matrix(nrow=length(jointCount_files), ncol=7)
  colnames(res)<- c("Focal_pos","sA","sT","sG", "sC","sAG", "ll")
  res[,1]<- 1:length(jointCount_files)
  for (i in 1:length(jointCount_files)) {
    jointCnt<- read.delim(jointCount_files[i], row.names = 1)
    
    optim.out=optim(par=c(vgMono[1:3]),ll_I,matT0=matT0,
                    jointCnt=jointCnt,p_ip1=as.numeric(basefreq[(i+1),2:5]/sum(basefreq[(i+1),2:5])))
    
    hat_vgMono=c(optim.out$par[1:3],-sum(optim.out$par[1:3]))
    hat_gAG=g_AG
    
    res[i,2:5]<- hat_vgMono
    res[i,6]<- hat_gAG
    
    res[i,7]<- -optim.out$value
    
    
  }
  res<- as.data.frame(res)
  return(res)
}

###########
### HII: Selection on AG
###########
### vgMono=c(0,0,0,0), gAG=x
selFunc_II=function(par,data){ ## takes parameters and data
  vgMono=rep(0,4)
  g_AG=par$g_AG
  jointCnt=data$jointCnt
  matSAG=matrix(nrow=4,ncol=4) ## build up two selection matrices
  p_im1=rep(0,4)
  for(i in 1:4){
    matSAG[i,i]=0
    p_im1[i]=sum(jointCnt[i,])
  }
  p_im1=p_im1/sum(p_im1)
  ## selection on monomers
  for(i in 1:4){
    for(j in 1:4){
      matSAG[i,j]=-vgMono[i]+vgMono[j]    
    }
  }
  ## selection on 'A' before following 'G' (appears like sel on monomer)
  matSAG[,1]= g_AG*data$p_ip1[3]+matSAG[,1]
  matSAG[1,]= -g_AG*data$p_ip1[3]+matSAG[1,]
  matS=matSAG
  ## selection for/against 'G' after 'A'
  matSAG[,3]= g_AG+matSAG[,3] 
  matSAG[3,]= -g_AG+matSAG[3,] 
  matTAG=matrix(nrow=4,ncol=4) ## build up transition matrices
  matT=matTAG
  for(i in 1:4){
    for(j in 1:4){
      matTAG[i,j]=fixp(matSAG[i,j])*data$matT0[i,j] 
      matT[i,j]=fixp(matS[i,j])*data$matT0[i,j]
    }
    matTAG[i,i]=0
    matTAG[i,i]=1-sum(matTAG[i,])
    matT[i,i]=0
    matT[i,i]=1-sum(matT[i,])
  }
  eSmTAG=eigenSystem(matTAG) ## calculate stationary distribution
  piAG=Re(eSmTAG$mForward[1,]/sum(eSmTAG$mForward[1,]))
  eSmT=eigenSystem(matT)
  pi=Re(eSmT$mForward[1,]/sum(eSmT$mForward[1,]))
  matExp=matT
  for(j in 1:4){
    matExp[1,j]=p_im1[1]*piAG[j] ## follows an 'A' => sel against 'G'
  }
  for(i in 2:4){
    for(j in 1:4){
      matExp[i,j]=p_im1[i]*pi[j] ## does no follow an 'A'
    }
  }
  return(matExp)
}

ll_II=function(par,matT0,jointCnt,p_ip1){
  parameters=list(g_AG=par)
  
  matExp=selFunc_II(parameters,data=list(matT0=matT0,p_ip1=p_ip1, jointCnt=jointCnt))
  return(-sum(log(matExp)*jointCnt)) 
}



##### 

Sperposition_II<- function(jointCount_files, basefreq, vgMono=rep(0,4), g_AG, matT0) {
  res<- matrix(nrow=length(jointCount_files), ncol=7)
  colnames(res)<- c("Focal_pos","sA","sT","sG", "sC","sAG", "ll")
  res[,1]<- 1:length(jointCount_files)
  for (i in 1:length(jointCount_files)) {
    jointCnt<- read.delim(jointCount_files[i], row.names = 1)
    
    optim.out=optim(par=g_AG,ll_II,matT0=matT0,
                    jointCnt=jointCnt,p_ip1=as.numeric(basefreq[(i+1),2:5]/sum(basefreq[(i+1),2:5])),method = "Brent",lower = -10, upper = 10)
    
    hat_gAG=optim.out$par[1]
    res[i,2:5]<- vgMono
    res[i,6]<- hat_gAG
    
    res[i,7]<- -optim.out$value
    
  }
  res<- as.data.frame(res)
  return(res)
}


###########
### HIII: Selection on both AG and monomers
###########
selFunc_III=function(par,data){ ## takes parameters and data
  vgMono=par$vgMono
  g_AG=par$g_AG
  jointCnt=data$jointCnt
  matSAG=matrix(nrow=4,ncol=4) ## build up two selection matrices
  p_im1=rep(0,4)
  for(i in 1:4){
    matSAG[i,i]=0
    p_im1[i]=sum(jointCnt[i,])
  }
  p_im1=p_im1/sum(p_im1)
  ## selection on monomers
  for(i in 1:4){
    for(j in 1:4){
      matSAG[i,j]=-vgMono[i]+vgMono[j]    
    }
  }
  ## selection on 'A' before following 'G' (appears like sel on monomer)
  matSAG[,1]= g_AG*data$p_ip1[3]+matSAG[,1]
  matSAG[1,]= -g_AG*data$p_ip1[3]+matSAG[1,]
  matS=matSAG
  ## selection for/against 'G' after 'A'
  matSAG[,3]= g_AG+matSAG[,3] 
  matSAG[3,]= -g_AG+matSAG[3,] 
  matTAG=matrix(nrow=4,ncol=4) ## build up transition matrices
  matT=matTAG
  for(i in 1:4){
    for(j in 1:4){
      matTAG[i,j]=fixp(matSAG[i,j])*data$matT0[i,j] 
      matT[i,j]=fixp(matS[i,j])*data$matT0[i,j]
    }
    matTAG[i,i]=0
    matTAG[i,i]=1-sum(matTAG[i,])
    matT[i,i]=0
    matT[i,i]=1-sum(matT[i,])
  }
  eSmTAG=eigenSystem(matTAG) ## calculate stationary distribution
  piAG=Re(eSmTAG$mForward[1,]/sum(eSmTAG$mForward[1,]))
  eSmT=eigenSystem(matT)
  pi=Re(eSmT$mForward[1,]/sum(eSmT$mForward[1,]))
  matExp=matT
  for(j in 1:4){
    matExp[1,j]=p_im1[1]*piAG[j] ## follows an 'A' => sel against 'G'
  }
  for(i in 2:4){
    for(j in 1:4){
      matExp[i,j]=p_im1[i]*pi[j] ## does no follow an 'A'
    }
  }
  return(matExp)
}

### 4 parameters 
ll_III=function(par,matT0,jointCnt,p_ip1){
  parameters=list(vgMono=c(par[1:3],-sum(par[1:3])),g_AG=par[4])
  matExp=selFunc(parameters,data=list(matT0=matT0,p_ip1=p_ip1, jointCnt=jointCnt))
  return(-sum(log(matExp)*jointCnt)) 
}

##### 
Sperposition_III<- function(jointCount_files, basefreq, vgMono, g_AG, matT0) {
  res<- matrix(nrow=length(jointCount_files), ncol=7)
  colnames(res)<- c("Focal_pos","sA","sT","sG", "sC","sAG", "ll")
  res[,1]<- 1:length(jointCount_files)
  for (i in 1:length(jointCount_files)) {
    jointCnt<- read.delim(jointCount_files[i], row.names = 1)
    
    
    optim.out=optim(par=c(vgMono[1:3],g_AG),ll,matT0=matT0,
                    jointCnt=jointCnt,p_ip1=as.numeric(basefreq[(i+1),2:5]/sum(basefreq[(i+1),2:5])))
    
    hat_vgMono=c(optim.out$par[1:3],-sum(optim.out$par[1:3]))
    hat_gAG=optim.out$par[4]
    res[i,2:5]<- hat_vgMono
    res[i,6]<- hat_gAG
    
    res[i,7]<- -optim.out$value
    
    
  }
  res<- as.data.frame(res)
  return(res)
}


