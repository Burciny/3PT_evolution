library(gtools)


### #### Get the consensus sequences from alignments 

consensus_find<- function(df) {
  con<- data.frame(matrix(vector(),0,1, dimnames = list(c(),c("Seq"))))
  for (i in 1:nrow(df)) {
    tmp<- unique(df[i,])
    con[i,1]<- tmp[which.max(tabulate(match(df[i,],tmp)))]
  }
  position<- rep(1:10,nrow(df)/10)
  con<- cbind(position,con)
  return(con)
}

### Weightings of oligonucletides
## Functions take input as consensus sequence (df) and
## the length of the sequence (23 for 5LR, 10 for 3PT)

dimer_weight<- function(df, rl=10,) {
  bases<- c("A","T","G","C")
  perm_2<- permutations(n=4, r=2, v=bases, repeats.allowed = TRUE)
  possible<- apply(perm_2, 1, function(x)paste0(x, collapse = ""))

  count<- data.frame(matrix(0,16,1, dimnames = list(c(),c("Count"))))
  weight<- cbind(possible,count)
  df<- as.data.frame(df)
  start_pos<- seq(1,nrow(df),rl)
  for (i in start_pos) {
    for (k in 0:(rl-2)) {
      m<- match(paste0(df[(i+k):(i+k+1),2], collapse = ""), possible)
      weight[m,2]<- weight[m,2]+1 
    }
  }
  return(weight)
}

trimer_weight<- function(df, rl=10) {
  bases<- c("A","T","G","C")
  perm_3<- permutations(n=4, r=3, v=bases, repeats.allowed = TRUE)
  possible_trimer<- apply(perm_3, 1, function(x)paste0(x, collapse = ""))
  
  count<- data.frame(matrix(0,64,1, dimnames = list(c(),c("Count"))))
  weight<- cbind(possible,count)
  df<- as.data.frame(df)
  start_pos<- seq(1,nrow(df),rl)
  for (i in start_pos) {
    for (k in 0:(rl-3)) {
      m<- match(paste0(df[(i+k):(i+k+2),2], collapse = ""), possible)
      weight[m,2]<- weight[m,2]+1 
    }
  }
  return(weight)
}

tetramer_weight<- function(df, rl=10) {
  bases<- c("A","T","G","C")
  perm_4<- permutations(n=4, r=4, v=bases, repeats.allowed = TRUE)
  possible_tetramer<- apply(perm_4, 1, function(x)paste0(x, collapse = ""))

  count<- data.frame(matrix(0,256,1, dimnames = list(c(),c("Count"))))
  weight<- cbind(possible,count)
  df<- as.data.frame(df)
  start_pos<- seq(1,nrow(df),rl)
  for (i in start_pos) {
    for (k in 0:(rl-4)) {
      m<- match(paste0(df[(i+k):(i+k+3),2], collapse = ""), possible)
      weight[m,2]<- weight[m,2]+1 
    }
  }
  return(weight)
}

