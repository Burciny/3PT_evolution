setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(readr)
library(gtools)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)

################ ################ ################ ################ ################ 
################ FIGURE S1 A-B ################ ################ ################ ################ 
chitest_trimer_3p_1 <- read_delim("../../data/chitest_trimer_3p_1", 
                                  "\t", escape_double = FALSE, trim_ws = TRUE)

chitest_trimer_3p_1$forward<- factor(chitest_trimer_3p_1$forward, levels = levels(reorder(chitest_trimer_3p_1$forward, chitest_trimer_3p_1$p_val)))
trimers_3p_1<-levels(chitest_trimer_3p_1$forward)

pdf("../../../suppfigures/FigS1A.pdf", width = 8, height = 5)

plot(1:64,chitest_trimer_3p_1$p_val[order(chitest_trimer_3p_1$p_val)], xlab = "Trimers", ylab = "p-values", xaxt="n", col="gray48",pch=16)
axis(1, at=1:64, labels = FALSE)
text(seq(1, 64, by=1), par("usr")[3]-3e-19, labels = trimers_3p_1, cex=0.65, srt = 90, pos = 1, xpd = TRUE)

dev.off()

CI_3p_trimer <- read_delim("../../data//CI_3p_trimer", 
                           "\t", escape_double = FALSE, trim_ws = TRUE)


idx<- c()
for (i in 1:length(trimers_3p_1)) {
  idx[i]<-which(CI_3p_trimer$X1==trimers_3p_1[i])
}
CI_3p_trimer_ro<- CI_3p_trimer[idx,]

pdf("../../../suppfigures/FigS1B.pdf", width = 8, height = 5)

plot(1:64, CI_3p_trimer_ro$Upper, type="l", xlab = "Trimers", ylab = "Confidence Intervals", xaxt="n")
lines(1:64, CI_3p_trimer_ro$Lower)
abline(h=c(1/1.1,1.1), col="red", lty=2)
polygon(x=c(1:64, rev(1:64)), y=c(CI_3p_trimer_ro$Upper, rev(CI_3p_trimer_ro$Lower)), col=adjustcolor("gray48", alpha.f = 0.40))
axis(1, at=1:64, labels = FALSE)
text(seq(1, 64, by=1), par("usr")[3]- 7, labels = trimers_3p_1,cex=0.65, srt = 90, pos = 1, xpd = TRUE)

dev.off()

################ ################ ################ ################ ################ 
################ FIGURE S2 A-B ################ ################ ################ ################ 
SFSauto_mel196 <- read_delim("../../data/SFSauto_mel196", 
                             "\t", escape_double = FALSE, trim_ws = TRUE)
SFSX_mel196 <- read_delim("../../data/SFSX_mel196", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)

M=nrow(SFSauto_mel196)-1
# proportional to theoretical distribution of polymorphic sites (neutral)
y=c(0,1/(1:(M-1))+1/((M-1):1),0) 
HM=c(1/(1:(M-1))+1/((M-1):1))
y=y/sum(y)

postscript("../../../suppfigures/FigS2.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 10, height = 4, colormodel = "cmyk")

par(mfrow=c(1,2))
plot(1:195, SFSauto_mel196$SIasymm[2:M]/sum(SFSauto_mel196$SIasymm[2:M]),pch=4, xlab = "Frequency",ylab = "Probability", main = "A", ylim = c(0,0.22))
lines(1:195, y[2:M]/sum(y[2:M]))

plot(1:195,SFSX_mel196$SIasymm[2:M]/sum(SFSX_mel196$SIasymm[2:M]),pch=4, xlab = "Frequency",ylab = "Probability", main = "B", ylim = c(0,0.22))
lines(1:195, y[2:M]/sum(y[2:M]))
par(mfrow=c(1,1))

dev.off()

################ ################ ################ ################ ################ 
################ FIGURE S3 A-B ################ ################ ################ ################ 

asymmetry_trimers_X <- read_delim("../../data/asymmetry_trimers_X", 
                                           delim = "\t", escape_double = FALSE, 
                                           trim_ws = TRUE)


asymmetry_trimers_X$forward<- factor(asymmetry_trimers_X$forward, levels = levels(reorder(asymmetry_trimers_X$forward, asymmetry_trimers_X$`3'-region`)))

asymmetry_trimers_X_long<- gather(asymmetry_trimers_X[,1:4], key="Region",value="value",3:4)

postscript("../../../suppfigures//FigS3A.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 4, colormodel = "cmyk")
ggplot(asymmetry_trimers_X_long, aes(x=forward, y=value*100, shape=Region))+geom_point()+
  theme(axis.text.x = element_text(size=8, angle=90), legend.position = c(0.07,0.89), legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"))+xlab("Trimers")+ylab("Asymmetry (%)")+
  scale_shape_discrete(labels=c("3PT","5LR"))

dev.off()

CI_Dmel_trimer_gamma_X <- read_delim("../../data/CI_Dmel_trimer_gamma_X", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)


postscript("../../../suppfigures/FigS3B.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 4, colormodel = "cmyk")
ggplot(CI_Dmel_trimer_gamma_X, aes(x=reorder(motif, gamma), y=gamma))+geom_point()+
  theme(axis.text.x = element_text(size=8, angle=90),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"))+xlab("Trimers")+
  geom_errorbar(aes(ymax=Upper, ymin=Lower))

dev.off()
################ ################ ################ ################ ################ 
################ FIGURE S4 A-B-C-D ################ ################ ################ ################ 

asymmetry_dimers <- read_delim("../../data/asymmetry_dimers", 
                                  delim = "\t", escape_double = FALSE, 
                                  trim_ws = TRUE)

asymmetry_dimers$forward<- factor(asymmetry_dimers$forward, levels = levels(reorder(asymmetry_dimers$forward, asymmetry_dimers$`3'-region`)))
asymmetry_dimers_long<- gather(asymmetry_dimers[,1:4], key="Region",value="value",3:4)

postscript("../../../suppfigures//FigS4A.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 4, colormodel = "cmyk")
ggplot(asymmetry_dimers_long, aes(x=forward, y=value*100, shape=Region))+geom_point()+
  theme(axis.text.x = element_text(size=8, angle=90), legend.position = c(0.07,0.89), legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"))+xlab("Dimers")+ylab("Asymmetry (%)")+
  scale_shape_discrete(labels=c("3PT","5LR"))

dev.off()

asymmetry_tetramers <- read_delim("../../data/asymmetry_tetramers", 
                                  delim = "\t", escape_double = FALSE, 
                                  trim_ws = TRUE)

asymmetry_tetramers$forward<- factor(asymmetry_tetramers$forward, levels = levels(reorder(asymmetry_tetramers$forward, asymmetry_tetramers$`3'-region`)))
asymmetry_tetramers_long<- gather(asymmetry_tetramers[,1:4], key="Region",value="value",3:4)

postscript("../../../suppfigures//FigS4B.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 4, colormodel = "cmyk")
ggplot(asymmetry_tetramers_long, aes(x=forward, y=value*100, shape=Region))+geom_point()+
  theme(axis.text.x = element_text(size=8, angle=90), legend.position = c(0.07,0.89), legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"))+xlab("Tetramers")+ylab("Asymmetry (%)")+
  scale_shape_discrete(labels=c("3PT","5LR"))

dev.off()


CI_Dmel_dimers <- read_delim("../../data/CI_Dmel_dimers", 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE)


postscript("../../../suppfigures/FigS4C.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 4, colormodel = "cmyk")
ggplot(CI_Dmel_dimers, aes(x=reorder(motif, gamma), y=gamma))+geom_point()+
  theme(axis.text.x = element_text(size=8, angle=90),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"))+xlab("Dimers")+
  geom_errorbar(aes(ymax=Upper, ymin=Lower))

dev.off()

CI_Dmel_tetramers <- read_delim("../../data/CI_Dmel_tetramers", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)


postscript("../../../suppfigures/FigS4D.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 4, colormodel = "cmyk")
ggplot(CI_Dmel_tetramers, aes(x=reorder(motif, gamma), y=gamma))+geom_point()+
  theme(axis.text.x = element_text(size=8, angle=90),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"))+xlab("Tetramers")+
  geom_errorbar(aes(ymax=Upper, ymin=Lower))

dev.off()

################ ################ ################ ################ ################ 
################ FIGURE S5 ################ ################ ################ ################ 

dimers_ordered <- read_delim("../../data/dimers_ordered", 
                             "\t", escape_double = FALSE, trim_ws = TRUE)

trimers_ordered <- read_delim("../../data/trimers_ordered", 
                              "\t", escape_double = FALSE, trim_ws = TRUE)

tetramers_ordered <- read_delim("../../data/tetramers_ordered", 
                                "\t", escape_double = FALSE, trim_ws = TRUE)

dimers_ordered$gamma[dimers_ordered$motif=="AG"]
postscript("../../../suppfigures/FigS5.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 5, colormodel = "cmyk")
par(mfrow=c(1,3))
hist(dimers_ordered$gamma, xlab = "gamma", main ="Dimers", cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
points(x=dimers_ordered$gamma[dimers_ordered$motif=="AG"], y=1, col="blue", pch=16,cex=2) # AG
points(x=dimers_ordered$gamma[dimers_ordered$motif=="GT"], y=1, col="red", pch=16,cex=2) # GT

hist(trimers_ordered$gamma, xlab = "gamma", main ="Trimers", cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
points(x=trimers_ordered$gamma[trimers_ordered$motif=="AGG"], y=1, col="blue", pch=16,cex=2) # AGG
points(x=trimers_ordered$gamma[trimers_ordered$motif=="TAG"], y=1, col="blue", pch=16,cex=2) # TAG
points(x=trimers_ordered$gamma[trimers_ordered$motif=="CAG"], y=1, col="blue", pch=16,cex=2) # TAG

points(x=trimers_ordered$gamma[trimers_ordered$motif=="GGT"], y=1, col="red", pch=16,cex=2) # GTA
points(x=trimers_ordered$gamma[trimers_ordered$motif=="GTG"], y=1, col="red", pch=16,cex=2) # GTG
points(x=trimers_ordered$gamma[trimers_ordered$motif=="GTA"], y=1, col="red", pch=16,cex=2) # GGT


hist(tetramers_ordered$gamma, xlab = "gamma", main ="Tetramers", cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
points(x=tetramers_ordered$gamma[tetramers_ordered$motif=="TAGG"], y=1, col="blue", pch=16,cex=2) # AGGG
points(x=tetramers_ordered$gamma[tetramers_ordered$motif=="CAGG"], y=1, col="blue", pch=16,cex=2) # TAGG
points(x=tetramers_ordered$gamma[tetramers_ordered$motif=="AGGG"], y=1, col="blue", pch=16,cex=2) # CAGG

points(x=tetramers_ordered$gamma[tetramers_ordered$motif=="GTAA"], y=1, col="red", pch=16,cex=2) # GTAA
points(x=tetramers_ordered$gamma[tetramers_ordered$motif=="GTGA"], y=1, col="red", pch=16,cex=2) # GTGA
points(x=tetramers_ordered$gamma[tetramers_ordered$motif=="GGTA"], y=1, col="red", pch=16,cex=2) # GGTA
points(x=tetramers_ordered$gamma[tetramers_ordered$motif=="GGTG"], y=1, col="red", pch=16,cex=2) # GGTG

dev.off()
par(mfrow=c(1,1))

################ ################ ################ ################ ################ 
################ FIGURE S6 ################ ################ ################ ################ 


CI_Dmel_trimer_gamma <- read_delim("../../data/CI_Dmel_trimer_gamma", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

gammas_trimer_np0_complement <- read_delim("../../data/gammas_trimer_np0_complement", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)
gammas_trimer_allnp0_complement <- read_delim("../../data/gammas_trimer_allnp0_complement", 
                                           delim = "\t", escape_double = FALSE, 
                                           trim_ws = TRUE)
gammas_trimer_phase0_complement <- read_delim("../../data/gammas_trimer_phase0_complement", 
                                              delim = "\t", escape_double = FALSE, 
                                              trim_ws = TRUE)


gammas_trimer_allnp0_complement$motif<- factor(gammas_trimer_allnp0_complement$motif, levels = levels(reorder(gammas_trimer_allnp0_complement$motif, CI_Dmel_trimer_gamma$gamma)))
gammas_trimer_phase0_complement$motif<- factor(gammas_trimer_phase0_complement$motif, levels = levels(reorder(gammas_trimer_phase0_complement$motif, CI_Dmel_trimer_gamma$gamma)))
gammas_trimer_np0_complement$motif<- factor(gammas_trimer_np0_complement$motif, levels = levels(reorder(gammas_trimer_np0_complement$motif, CI_Dmel_trimer_gamma$gamma)))

postscript("../../../suppfigures/FigS6.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 4, colormodel = "cmyk")
ggplot(CI_Dmel_trimer_gamma, aes(x=reorder(motif, gamma), y=gamma))+geom_point()+
  theme(axis.text.x = element_text(size=8, angle=90),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"))+xlab("Trimers")+
  geom_errorbar(aes(ymax=Upper, ymin=Lower))+ylab("gamma")+
  geom_point(data=gammas_trimer_allnp0_complement,aes(x=motif, y=gamma,shape=17))+
  geom_point(data=gammas_trimer_phase0_complement,aes(x=motif, y=gamma, shape=15))+
  geom_point(data=gammas_trimer_np0_complement,size=2,aes(x=motif, y=gamma, shape=18))+
  scale_shape_identity(breaks = c(16, 17, 15,18),
                       labels = c("0","I","II","III"),
                       guide = "legend")+
  theme(legend.position = c(0.045,0.82), legend.title = element_blank())
dev.off()

################ ################ ################ ################ ################ 
################ FIGURE S7 ################ ################ ################ ################ 


CI_Dmel_trimer_gamma <- read_delim("../../data/CI_Dmel_trimer_gamma", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

CI_Dsim_trimers <- read_delim("../../data/CI_Dsim_trimers", 
                                           delim = "\t", escape_double = FALSE, 
                                           trim_ws = TRUE)

postscript("../../../suppfigures/FigS7.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 4, colormodel = "cmyk")
ggplot(CI_Dmel_trimer_gamma, aes(x=reorder(motif, gamma), y=gamma))+geom_point(aes(shape=16))+
  theme(axis.text.x = element_text(size=8, angle=90),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"))+xlab("Trimers")+
  geom_errorbar(data = CI_Dmel_trimer_gamma,aes(ymax=Upper, ymin=Lower))+
  geom_errorbar(data = CI_Dsim_trimers,aes(ymax=Upper, ymin=Lower))+ylab("gamma")+
  geom_point(data=CI_Dsim_trimers,aes(x=motif, y=gamma,shape=17))+
  scale_shape_identity(breaks=c(16,17),
                       labels=c("Mel","Sim"),
                       guide="legend")+
  theme(legend.position = c(0.07,0.89), legend.title = element_blank())

dev.off()


################ ################ ################ ################ ################ 
################ FIGURE S8 ################ ################ ################ ################ 

CI_X_3PT <- read_delim("../../data/X_gammaHIII_3PT_selcoeffs_withBS", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)

X_3PT_selcoeffs <- read_delim("../../data/X_gammaHIII_3PT_selcoeffs", 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE)

postscript("../../../suppfigures/FigS8.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 5, colormodel = "cmyk")

plot(1:9, X_3PT_selcoeffs$sAG[1:9],type = "l", ylim = c(-6,1) ,xlab = "Focal position", ylab = "Gamma")
points(1:9, X_3PT_selcoeffs$sAG[1:9], pch=16, cex=0.5)
arrows(1:9, CI_X_3PT$UB[37:45],1:9, CI_X_3PT$LB[37:45], angle=90, code=3, length=0.06)

lines(1:9, X_3PT_selcoeffs$sA[1:9], col="red")
points(1:9, X_3PT_selcoeffs$sA[1:9], col="red", pch=16,cex=0.5)
arrows(1:9, CI_X_3PT$UB[1:9],1:9, CI_X_3PT$LB[1:9], angle=90, code=3, length=0.06, col="red")

lines(1:9, X_3PT_selcoeffs$sT[1:9], col="purple")
points(1:9, X_3PT_selcoeffs$sT[1:9], col="purple", pch=16,cex=0.5)
arrows(1:9, CI_X_3PT$UB[10:18],1:9, CI_X_3PT$LB[10:18], angle=90, code=3, length=0.06, col="purple")

lines(1:9, X_3PT_selcoeffs$sG[1:9], col="blue")
points(1:9, X_3PT_selcoeffs$sG[1:9], col="blue", pch=16,cex=0.5)
arrows(1:9, CI_X_3PT$UB[19:27],1:9, CI_X_3PT$LB[19:27], angle=90, code=3, length=0.06, col="blue")

lines(1:9, X_3PT_selcoeffs$sC[1:9], col="green")
points(1:9, X_3PT_selcoeffs$sC[1:9], col="green", pch=16,cex=0.5)
arrows(1:9, CI_X_3PT$UB[28:36],1:9, CI_X_3PT$LB[28:36], angle=90, code=3, length=0.06, col="green")

legend("bottomleft", legend = c("AG","A", "T","G","C"), col = c("black","red","purple","blue","green"),lty = c(1,1,1,1,1),cex=0.6)

dev.off()


################ ################ ################ ################ ################ 
################ FIGURE S9 ################ ################ ################ ################ 

AllDeviationsJntCnt_X <- read_delim("../../data/AllDeviationsJntCnt_X", 
                                           delim = "\t", escape_double = FALSE, 
                                           trim_ws = TRUE)

deviation1to2<- AllDeviationsJntCnt_X[AllDeviationsJntCnt_X$pos_i==1 &AllDeviationsJntCnt_X$pos_j==2,]
deviation4to5<- AllDeviationsJntCnt_X[AllDeviationsJntCnt_X$pos_i==4 &AllDeviationsJntCnt_X$pos_j==5,]
deviation8to9<- AllDeviationsJntCnt_X[AllDeviationsJntCnt_X$pos_i==8 &AllDeviationsJntCnt_X$pos_j==9,]

deviation1to2$Base_i<- factor(deviation1to2$Base_i, levels = c("A","T","G","C"))
deviation1to2$Base_j<- factor(deviation1to2$Base_j, levels = c("A","T","G","C"))


postscript("../../../suppfigures/FigS9_1to2.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 9, height = 4, colormodel = "cmyk")

ggplot(deviation1to2,aes(x=Base_j,y=Base_i,fill=Deviation))+
  geom_tile()+
  scale_fill_gradient2(high="red",low="blue",mid="white")+xlab("Position 2")+ylab("Position 1")+
  theme(panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"))+
  facet_grid(.~Model)+geom_text(aes(label = round(Deviation,3)), size=2)

dev.off()

deviation4to5$Base_i<- factor(deviation4to5$Base_i, levels = c("A","T","G","C"))
deviation4to5$Base_j<- factor(deviation4to5$Base_j, levels = c("A","T","G","C"))

postscript("../../../suppfigures/FigS9_4to5.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 9, height = 4, colormodel = "cmyk")

ggplot(deviation4to5,aes(x=Base_j,y=Base_i,fill=Deviation))+
  geom_tile()+
  scale_fill_gradient2(high="red",low="blue",mid="white")+xlab("Position 5")+ylab("Position 4")+
  theme(panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"))+
  facet_grid(.~Model)+geom_text(aes(label = round(Deviation,3)), size=2)

dev.off()

deviation8to9$Base_i<- factor(deviation8to9$Base_i, levels = c("A","T","G","C"))
deviation8to9$Base_j<- factor(deviation8to9$Base_j, levels = c("A","T","G","C"))

postscript("../../../suppfigures/FigS9_8to9.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 9, height = 4, colormodel = "cmyk")

ggplot(deviation8to9,aes(x=Base_j,y=Base_i,fill=Deviation))+
  geom_tile()+
  scale_fill_gradient2(high="red",low="blue",mid="white")+xlab("Position 9")+ylab("Position 1-8")+
  theme(panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"))+
  facet_grid(.~Model)+geom_text(aes(label = round(Deviation,3)), size=2)

dev.off()

################ ################ ################ ################ ################ 
################ FIGURE S10 ################ ################ ################ ################ 
rice <- read_csv("../../data/8Eukaryoutes/lengths_rice_woUTR", col_names = FALSE)
arabidopsis <- read_csv("../../data/8Eukaryoutes/lengths_at_woUTR", col_names = FALSE)
moss <- read_csv("../../data/8Eukaryoutes/lengths_moss_woUTR", col_names = FALSE)
celegans <- read_csv("../../data/8Eukaryoutes/lengths_celegans_woUTR", col_names = FALSE)
human <- read_csv("../../data/8Eukaryoutes/lengths_human_woUTR", col_names = FALSE)
seaurchin <- read_csv("../../data/8Eukaryoutes/lengths_seaurchin_woUTR", col_names = FALSE)
cerevisiae <- read_csv("../../data/8Eukaryoutes/intron_lengths_Scer", col_names = FALSE)
kluyveri <- read_csv("../../data/8Eukaryoutes/intron_lengths_Skluyveri", col_names = FALSE)

postscript("../../../suppfigures/FigS10_1.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 5, colormodel = "cmyk")

par(mfrow=c(3,2))
plot(density(log10(human$X1)), main="Human", xlab="log10(intron length)")
abline(v=c(log10(50),log10(100)), lty=2)
plot(density(log10(seaurchin$X1)), main="Sea Urchin", xlab="log10(intron length)")
abline(v=c(log10(300),log10(500)), lty=2)
plot(density(log10(rice$X1)), main="Rice", xlab="log10(intron length)")
abline(v=c(log10(65),log10(100)), lty=2)
plot(density(log10(arabidopsis$X1)), main="Arabidopsis", xlab="log10(intron length)")
abline(v=c(log10(65),log10(100)), lty=2)
plot(density(log10(moss$X1)), main="Moss", xlab="log10(intron length)")
abline(v=c(log10(90),log10(200)), lty=2)
plot(density(log10(celegans$X1)), main="C. elegans", xlab="log10(intron length)")
abline(v=c(log10(40),log10(65)), lty=2)
par(mfrow=c(1,1))

dev.off()

postscript("../../../suppfigures/FigS10_2.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 3, colormodel = "cmyk")

par(mfrow=c(1,2))
plot(density(log10(cerevisiae$X1)), main="S. cerevisiae", xlab="log10(intron length)")
plot(density(log10(kluyveri$X1)), main="L. kluyveri", xlab="log10(intron length)")
par(mfrow=c(1,1))

dev.off()
################ ################ ################ ################ ################ 
################ FIGURE S11 ################ ################ ################ ################ 
basecomposition_cerevisiae <- read_delim("~/Documents/Project/3prime_ms/JEB/Repository/data/8Eukaryoutes/basecomposition_cerevisiae", 
                                         delim = "\t", escape_double = FALSE, 
                                         trim_ws = TRUE)

postscript("../../../suppfigures/FigS11_Scer.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 5, colormodel = "cmyk")

plot(basecomposition_cerevisiae$positions, basecomposition_cerevisiae$A, type = "l", col="red",ylab="counts",xlab = "positions", main = "S.cerevisiae")
lines(basecomposition_cerevisiae$positions, basecomposition_cerevisiae$T, col="purple")
lines(basecomposition_cerevisiae$positions, basecomposition_cerevisiae$G, col="blue")
lines(basecomposition_cerevisiae$positions, basecomposition_cerevisiae$C, col="green")
legend("top", legend = c("A", "T","G","C"), col = c("red","purple","blue","green"),lty = c(1,1,1,1),cex=0.75)

dev.off()

basecomposition_thermo <- read_delim("~/Documents/Project/3prime_ms/JEB/Repository/data/8Eukaryoutes/basecomposition_thermo", 
                                         delim = "\t", escape_double = FALSE, 
                                         trim_ws = TRUE)

postscript("../../../suppfigures/FigS11_Lther.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 5, colormodel = "cmyk")
plot(basecomposition_thermo$positions, basecomposition_thermo$A, type = "l", col="red",ylab="counts",xlab = "positions", main = "L. thermotolerans")
lines(basecomposition_thermo$positions, basecomposition_thermo$T, col="purple")
lines(basecomposition_thermo$positions, basecomposition_thermo$G, col="blue")
lines(basecomposition_thermo$positions, basecomposition_thermo$C, col="green")
legend("top", legend = c("A", "T","G","C"), col = c("red","purple","blue","green"),lty = c(1,1,1,1),cex=0.75)

dev.off()

basecomposition_l200_moss <- read_delim("~/Documents/Project/3prime_ms/JEB/Repository/data/8Eukaryoutes/basecomposition_l200_moss", 
                                         delim = "\t", escape_double = FALSE, 
                                         trim_ws = TRUE)

postscript("../../../suppfigures/FigS11_moss.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 5, colormodel = "cmyk")
plot(basecomposition_l200_moss$positions, basecomposition_l200_moss$A, type = "l", col="red",ylab="counts",xlab = "positions",main="Moss")
lines(basecomposition_l200_moss$positions, basecomposition_l200_moss$T, col="purple")
lines(basecomposition_l200_moss$positions, basecomposition_l200_moss$G, col="blue")
lines(basecomposition_l200_moss$positions, basecomposition_l200_moss$C, col="green")
legend("top", legend = c("A", "T","G","C"), col = c("red","purple","blue","green"),lty = c(1,1,1,1),cex=0.75)

dev.off()

basecomposition_l500_su <- read_delim("~/Documents/Project/3prime_ms/JEB/Repository/data/8Eukaryoutes/basecomposition_l500_su", 
                                         delim = "\t", escape_double = FALSE, 
                                         trim_ws = TRUE)

postscript("../../../suppfigures/FigS11_su.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 5, colormodel = "cmyk")
plot(basecomposition_l500_su$positions[300:500], basecomposition_l500_su$A[300:500], type = "l", col="red",ylab="counts",xlab = "positions", main = "Sea Urchin")
lines(basecomposition_l500_su$positions[300:500], basecomposition_l500_su$T[300:500], col="purple")
lines(basecomposition_l500_su$positions[300:500], basecomposition_l500_su$G[300:500], col="blue")
lines(basecomposition_l500_su$positions[300:500], basecomposition_l500_su$C[300:500], col="green")
legend("top", legend = c("A", "T","G","C"), col = c("red","purple","blue","green"),lty = c(1,1,1,1),cex=0.75)


dev.off()

basecomposition_l60_celegans <- read_delim("~/Documents/Project/3prime_ms/JEB/Repository/data/8Eukaryoutes/basecomposition_l60_celegans", 
                                         delim = "\t", escape_double = FALSE, 
                                         trim_ws = TRUE)

postscript("../../../suppfigures/FigS11_ce.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 5, colormodel = "cmyk")

plot(basecomposition_l60_celegans$positions, basecomposition_l60_celegans$A, type = "l", col="red",ylab="counts",xlab = "positions", main="C. elegans")
lines(basecomposition_l60_celegans$positions, basecomposition_l60_celegans$T, col="purple")
lines(basecomposition_l60_celegans$positions, basecomposition_l60_celegans$G, col="blue")
lines(basecomposition_l60_celegans$positions, basecomposition_l60_celegans$C, col="green")
legend("top", legend = c("A", "T","G","C"), col = c("red","purple","blue","green"),lty = c(1,1,1,1),cex=0.75)

dev.off()

basecomposition_l90_at <- read_delim("~/Documents/Project/3prime_ms/JEB/Repository/data/8Eukaryoutes/basecomposition_l90_at", 
                                         delim = "\t", escape_double = FALSE, 
                                         trim_ws = TRUE)

postscript("../../../suppfigures/FigS11_at.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 5, colormodel = "cmyk")

plot(basecomposition_l90_at$positions, basecomposition_l90_at$A, type = "l", col="red",ylab="counts",xlab = "positions", main = "Arabidopsis")
lines(basecomposition_l90_at$positions, basecomposition_l90_at$T, col="purple")
lines(basecomposition_l90_at$positions, basecomposition_l90_at$G, col="blue")
lines(basecomposition_l90_at$positions, basecomposition_l90_at$C, col="green")
legend("top", legend = c("A", "T","G","C"), col = c("red","purple","blue","green"),lty = c(1,1,1,1),cex=0.75)

dev.off()

basecomposition_l90_human <- read_delim("~/Documents/Project/3prime_ms/JEB/Repository/data/8Eukaryoutes/basecomposition_l90_human", 
                                         delim = "\t", escape_double = FALSE, 
                                         trim_ws = TRUE)

postscript("../../../suppfigures/FigS11_human.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 5, colormodel = "cmyk")
plot(basecomposition_l90_human$positions, basecomposition_l90_human$A, type = "l", col="red",ylab="counts",xlab = "positions", main="Human")
lines(basecomposition_l90_human$positions, basecomposition_l90_human$T, col="purple")
lines(basecomposition_l90_human$positions, basecomposition_l90_human$G, col="blue")
lines(basecomposition_l90_human$positions, basecomposition_l90_human$C, col="green")
legend("top", legend = c("A", "T","G","C"), col = c("red","purple","blue","green"),lty = c(1,1,1,1),cex=0.75)


dev.off()

basecomposition_l90_rice <- read_delim("~/Documents/Project/3prime_ms/JEB/Repository/data/8Eukaryoutes/basecomposition_l90_rice", 
                                         delim = "\t", escape_double = FALSE, 
                                         trim_ws = TRUE)

postscript("../../../suppfigures/FigS11_rice.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 5, colormodel = "cmyk")

plot(basecomposition_l90_rice$positions, basecomposition_l90_rice$A, type = "l", col="red",ylab="counts",xlab = "positions", main = "Rice")
lines(basecomposition_l90_rice$positions, basecomposition_l90_rice$T, col="purple")
lines(basecomposition_l90_rice$positions, basecomposition_l90_rice$G, col="blue")
lines(basecomposition_l90_rice$positions, basecomposition_l90_rice$C, col="green")
legend("top", legend = c("A", "T","G","C"), col = c("red","purple","blue","green"),lty = c(1,1,1,1),cex=0.75)

dev.off()

################ ################ ################ ################ ################ 
################ FIGURE S12 ################ ################ ################ ################ 


dimer_all <- read_delim("../../data/8Eukaryoutes/dimer_all_8eu", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)

human<-dimer_all[which(dimer_all$species=="Human"),]
su<-dimer_all[which(dimer_all$species=="Sea Urchin"),]
rice<-dimer_all[which(dimer_all$species=="Rice"),]
arabidopsis<-dimer_all[which(dimer_all$species=="Arabidopsis"),]
moss<-dimer_all[which(dimer_all$species=="Moss"),]
celegans<-dimer_all[which(dimer_all$species=="C. elegans"),]
cerevisiae<-dimer_all[which(dimer_all$species=="S. cerevisiae"),]
ther<-dimer_all[which(dimer_all$species=="L. thermotolerans"),]

human$forward<- factor(human$forward, levels = levels(reorder(human$forward, human$Asymmetry)))
su$forward<- factor(su$forward, levels = levels(reorder(su$forward, su$Asymmetry)))
rice$forward<- factor(rice$forward, levels = levels(reorder(rice$forward, rice$Asymmetry)))
arabidopsis$forward<- factor(arabidopsis$forward, levels = levels(reorder(arabidopsis$forward, arabidopsis$Asymmetry)))
moss$forward<- factor(moss$forward, levels = levels(reorder(moss$forward, moss$Asymmetry)))
celegans$forward<- factor(celegans$forward, levels = levels(reorder(celegans$forward, celegans$Asymmetry)))
cerevisiae$forward<- factor(cerevisiae$forward, levels = levels(reorder(cerevisiae$forward, cerevisiae$Asymmetry)))
ther$forward<- factor(ther$forward, levels = levels(reorder(ther$forward, ther$Asymmetry)))



cerevisiaed<- ggplot(cerevisiae, aes(x=forward, y=Asymmetry*100))+geom_point()+ylim(-100,100)+
  theme(axis.text.x = element_text(size=8),axis.text.y = element_text(size=6),plot.title = element_text(size=8),panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"),axis.title.x = element_text(size=5),axis.title.y = element_text(size=5))+xlab("Dimers")+ylab("Asymmetry (%)")+ggtitle("S. cerevisiae")

therd<- ggplot(ther, aes(x=forward, y=Asymmetry*100))+geom_point()+ylim(-100,100)+
  theme(axis.text.x = element_text(size=8),axis.text.y = element_text(size=6),plot.title = element_text(size=8),panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"),axis.title.x = element_text(size=5),axis.title.y = element_text(size=5))+xlab("Dimers")+ylab("Asymmetry (%)")+ggtitle("L. thermotolerans")

humand<-ggplot(human, aes(x=forward, y=Asymmetry*100))+geom_point()+ylim(-100,100)+
  theme(axis.text.x = element_text(size=8),axis.text.y = element_text(size=6),plot.title = element_text(size=8),panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"),axis.title.x = element_text(size=5),axis.title.y = element_text(size=5))+xlab("Dimers")+ylab("Asymmetry (%)")+ggtitle("Human")

arabidopsisd<-ggplot(arabidopsis, aes(x=forward, y=Asymmetry*100))+geom_point()+ylim(-100,100)+
  theme(axis.text.x = element_text(size=8),axis.text.y = element_text(size=6),plot.title = element_text(size=8),panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"),axis.title.x = element_text(size=5),axis.title.y = element_text(size=5))+xlab("Dimers")+ylab("Asymmetry (%)")+ggtitle("Arabidopsis")

celegansd<-ggplot(celegans, aes(x=forward, y=Asymmetry*100))+geom_point()+ylim(-100,100)+
  theme(axis.text.x = element_text(size=8),axis.text.y = element_text(size=6),plot.title = element_text(size=8),panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"),axis.title.x = element_text(size=5),axis.title.y = element_text(size=5))+xlab("Dimers")+ylab("Asymmetry (%)")+ggtitle("C. elegans")

mossd<-ggplot(moss, aes(x=forward, y=Asymmetry*100))+geom_point()+ylim(-100,100)+
  theme(axis.text.x = element_text(size=8),axis.text.y = element_text(size=6),plot.title = element_text(size=8),panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"),axis.title.x = element_text(size=5),axis.title.y = element_text(size=5))+xlab("Dimers")+ylab("Asymmetry (%)")+ggtitle("Moss")

riced<-ggplot(rice, aes(x=forward, y=Asymmetry*100))+geom_point()+ylim(-100,100)+
  theme(axis.text.x = element_text(size=8),axis.text.y = element_text(size=6),plot.title = element_text(size=8),panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"),axis.title.x = element_text(size=5),axis.title.y = element_text(size=5))+xlab("Dimers")+ylab("Asymmetry (%)")+ggtitle("Rice")

sud<-ggplot(su, aes(x=forward, y=Asymmetry*100))+geom_point()+ylim(-100,100)+
  theme(axis.text.x = element_text(size=8),axis.text.y = element_text(size=6),plot.title = element_text(size=8),panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"),axis.title.x = element_text(size=5),axis.title.y = element_text(size=5))+xlab("Dimers")+ylab("Asymmetry (%)")+ggtitle("Sea Urchin")

postscript("../../../suppfigures/FigS12.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 5, colormodel = "cmyk")
grid.arrange(humand, sud, riced, arabidopsisd, mossd, celegansd,cerevisiaed, therd,ncol=2)

dev.off()

################ ################ ################ ################ ################ 
################ FIGURE S13 ################ ################ ################ ################ 

## Trimers

trimer_all <- read_delim("../../data/8Eukaryoutes/trimer_all_8eu", 
                                     "\t", escape_double = FALSE, trim_ws = TRUE)

human<-trimer_all[which(trimer_all$species=="Human"),]
su<-trimer_all[which(trimer_all$species=="Sea Urchin"),]
rice<-trimer_all[which(trimer_all$species=="Rice"),]
arabidopsis<-trimer_all[which(trimer_all$species=="Arabidopsis"),]
moss<-trimer_all[which(trimer_all$species=="Moss"),]
celegans<-trimer_all[which(trimer_all$species=="C. elegans"),]
cerevisiae<-trimer_all[which(trimer_all$species=="S. cerevisiae"),]
ther<-trimer_all[which(trimer_all$species=="L. thermotolerans"),]

human$forward<- factor(human$forward, levels = levels(reorder(human$forward, human$Asymmetry)))
su$forward<- factor(su$forward, levels = levels(reorder(su$forward, su$Asymmetry)))
rice$forward<- factor(rice$forward, levels = levels(reorder(rice$forward, rice$Asymmetry)))
arabidopsis$forward<- factor(arabidopsis$forward, levels = levels(reorder(arabidopsis$forward, arabidopsis$Asymmetry)))
moss$forward<- factor(moss$forward, levels = levels(reorder(moss$forward, moss$Asymmetry)))
celegans$forward<- factor(celegans$forward, levels = levels(reorder(celegans$forward, celegans$Asymmetry)))
cerevisiae$forward<- factor(cerevisiae$forward, levels = levels(reorder(cerevisiae$forward, cerevisiae$Asymmetry)))
ther$forward<- factor(ther$forward, levels = levels(reorder(ther$forward, ther$Asymmetry)))



humant<-ggplot(human, aes(x=forward, y=Asymmetry*100, colour=category))+geom_point(size=1)+ylim(-100,100)+
  theme(axis.text.x = element_text(size=4, angle=90),axis.text.y = element_text(size=4),plot.title = element_text(size=8), panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"), legend.position = "none", axis.title.x = element_text(size=5), axis.title.y = element_text(size=5))+
  xlab("Trimers")+ylab("Asymmetry (%)")+ggtitle("Human")

sut<-ggplot(su, aes(x=forward, y=Asymmetry*100, colour=category))+geom_point(size=1)+ylim(-100,100)+
  theme(axis.text.x = element_text(size=4, angle=90),axis.text.y = element_text(size=4),plot.title = element_text(size=8), panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"), legend.position = "none", axis.title.x = element_text(size=5), axis.title.y = element_text(size=5))+
  xlab("Trimers")+ylab("Asymmetry (%)")+ggtitle("Sea Urchin")

ricet<-ggplot(rice, aes(x=forward, y=Asymmetry*100, colour=category))+geom_point(size=1)+ylim(-100,100)+
  theme(axis.text.x = element_text(size=4, angle=90),axis.text.y = element_text(size=4),plot.title = element_text(size=8), panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"), legend.position = "none", axis.title.x = element_text(size=5), axis.title.y = element_text(size=5))+
  xlab("Trimers")+ylab("Asymmetry (%)")+ggtitle("Rice")

mosst<-ggplot(moss, aes(x=forward, y=Asymmetry*100, colour=category))+geom_point(size=1)+ylim(-100,100)+
  theme(axis.text.x = element_text(size=4, angle=90),axis.text.y = element_text(size=4),plot.title = element_text(size=8), panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"), legend.position = "none", axis.title.x = element_text(size=5), axis.title.y = element_text(size=5))+
  xlab("Trimers")+ylab("Asymmetry (%)")+ggtitle("Moss")

arabidopsist<-ggplot(arabidopsis, aes(x=forward, y=Asymmetry*100, colour=category))+geom_point(size=1)+ylim(-100,100)+
  theme(axis.text.x = element_text(size=4, angle=90),axis.text.y = element_text(size=4),plot.title = element_text(size=8), panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"), legend.position = "none", axis.title.x = element_text(size=5), axis.title.y = element_text(size=5))+
  xlab("Trimers")+ylab("Asymmetry (%)")+ggtitle("Arabidopsis")

celeganst<-ggplot(celegans, aes(x=forward, y=Asymmetry*100, colour=category))+geom_point(size=1)+ylim(-100,100)+
  theme(axis.text.x = element_text(size=4, angle=90),axis.text.y = element_text(size=4),plot.title = element_text(size=8), panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"), legend.position = "none", axis.title.x = element_text(size=5), axis.title.y = element_text(size=5))+
  xlab("Trimers")+ylab("Asymmetry (%)")+ggtitle("C. elegans")

cerevisiaet<-ggplot(cerevisiae, aes(x=forward, y=Asymmetry*100, colour=category))+geom_point(size=1)+ylim(-100,100)+
  theme(axis.text.x = element_text(size=4, angle=90),axis.text.y = element_text(size=4),plot.title = element_text(size=8), panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"), legend.position = "none", axis.title.x = element_text(size=5), axis.title.y = element_text(size=5))+
  xlab("Trimers")+ylab("Asymmetry (%)")+ggtitle("S. cerevisiae")


thermot<-ggplot(ther, aes(x=forward, y=Asymmetry*100, colour=category))+geom_point(size=1)+ylim(-100,100)+
  theme(axis.text.x = element_text(size=4, angle=90),axis.text.y = element_text(size=4),plot.title = element_text(size=8), panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"), legend.position = "none", axis.title.x = element_text(size=5), axis.title.y = element_text(size=5))+
  xlab("Trimers")+ylab("Asymmetry (%)")+ggtitle("L. thermotolerans")

postscript("../../../suppfigures/FigS13.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 6, colormodel = "cmyk")
grid.arrange(humant, sut, ricet, arabidopsist, mosst, celeganst,cerevisiaet,thermot, ncol=2)

dev.off()



