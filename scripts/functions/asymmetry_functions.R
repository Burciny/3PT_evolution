
complements <- read_delim("../../data/complements",
                          "\t", escape_double = FALSE, col_names = FALSE,
                          trim_ws = TRUE)
complements_dimer <- read_delim("../../data/complements_dimer",
                          "\t", escape_double = FALSE, col_names = FALSE,
                          trim_ws = TRUE)
complements_trimer <- read_delim("../../data/complements_trimer",
                          "\t", escape_double = FALSE, col_names = FALSE,
                          trim_ws = TRUE)

## Function takes the oligonucletiode counts (df_5LR, df_3PT) and complementary motif files
## to calculate asymmetry scores of both region in one file. 

asymmetry_score_calc <- function(df_5LR,df_3PT, complements) {
  asymm_5LR<- c()
  asymm_3PT<- c()
  
  for (i in 1:nrow(complements)) {
    asymm_5LR[i]<-(df_5LR$Count[df_5LR$possible==complements$X1[i]]-df_5LR$Count[df_5LR$possible==complements$X2[i]])/
      (df_5LR$Count[df_5LR$possible==complements$X1[i]]+df_5LR$Count[df_5LR$possible==complements$X2[i]])
    asymm_3PT[i]<-(df_3PT$Count[df_3PT$possible==complements$X1[i]]-df_3PT$Count[df_3PT$possible==complements$X2[i]])/
      (df_3PT$Count[df_3PT$possible==complements$X1[i]]+df_3PT$Count[df_3PT$possible==complements$X2[i]])
    
  }
  asymm<- as.data.frame(cbind(complements$X1, complements$X2, asymm_5LR,asymm_3PT))
  colnames(asymm)<- c("forward","reverse","5LR","3PT")
  asymm$forward<- factor(asymm$forward, levels = levels(reorder(asymm$forward, asymm$`3PT`)))
  
  return(asymm)
  
}

##### Gamma with asymmetry scores

### Function take oligonucleotide count 
movercomp<- function(counts) {
  moverc<- c()
  for (i in 1:length(counts)) {
    moverc[i]<- counts[i]/(1-counts[i])
  }
  return(moverc)
}

### Function to calculate gamma (Eq. 4 in the article)

gamma_find<- function(counts5,counts3) {    # input: oligonucleotide counts
  sum_trimer5<- sum(counts5$Count)
  sum_trimer3<- sum(counts3$Count)
  obs_count_5p<- counts5$Count
  obs_count_3p<- counts3$Count
  
  obs_ratio_5p<- obs_count_5p/sum_trimer5
  obs_ratio_3p<- obs_count_3p/sum_trimer3
  movercomp_3p<- movercomp(obs_ratio_3p)
  movercomp_5p<- movercomp(obs_ratio_5p)
  gammas<-log(movercomp_3p)-log(movercomp_5p)
  return(gammas)
}


### Bootstrapping functions 

sample_introns<- function(cons, n) {  ## input consensus sequence
  intron_count<-nrow(cons)/n
  intronnames<- 1:intron_count
  introns<-c()
  for (i in intronnames) {
    introns<- append(introns,rep(i, n))
  }
  cons$introns<- introns
  samp_intron<- sample(intronnames, intron_count, replace = TRUE)
  idx<-c()
  for (i in 1:length(samp_intron)) {
    idx<- append(idx,which(cons$introns==samp_intron[i]))
  }
  new_samp<- cons[idx,]
  return(new_samp)
}

bootstrap_oligogamma<- function(cons5, cons3, BS=1000, o_length=2) { # Inputs consensus sequences
  bases<- c("A","T","G","C")
  perm<- permutations(n=4, r=o_length, v=bases, repeats.allowed = TRUE)
  possible<- apply(perm, 1, function(x)paste0(x, collapse = ""))
  
  n5=23
  n3=10
  results<- matrix(0,length(possible),(BS+1))
  results[,1]<- possible
  for (b in 1:BS) {
    samp5<- sample_introns(cons5, n=n5)
    samp3<- sample_introns(cons3, n=n3)
    trimer_counts5<- trimer_weight(samp5,possible,n=n5)
    trimer_counts3<- trimer_weight(samp3,possible,n=n3)
    gammas<- gamma_find(trimer_counts5, trimer_counts3)
    results[,b+1]<- gammas
  }
  return(results)
}



