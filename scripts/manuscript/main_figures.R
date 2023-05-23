setwd("~/Documents/Project/3prime_ms/JEB/Repository/scripts/manuscript/")

library(readr)
library(gtools)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)


################ ################ ################ ################ ################ 
################ FIGURE 2A-B ################ ################ ################ ################ 

allcounts_155 <- read_delim("../../data/allcounts_155_A", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)
allcounts_170 <- read_delim("../../data/allcounts_170_A", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)
allcounts_185 <- read_delim("../../data/allcounts_185_A", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)

postscript("../../../mainfigures/Fig2A.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 10, height = 3, colormodel = "cmyk")

par(mfrow=c(1,3))
plot(allcounts_155$position_155, allcounts_155$A, type = "l", col="red", xaxt="n", ylab = "counts", xlab = "position", main = "length 55")
lines(allcounts_155$position_155, allcounts_155$T, col="purple")
lines(allcounts_155$position_155, allcounts_155$G, col="blue")
lines(allcounts_155$position_155, allcounts_155$C, col="green")
legend("topright", legend = c("A", "T","G","C"), col = c("red","purple","blue","green"),lty = c(1,1,1,1),cex=0.75)
axis(1, at=c(-50,-25,0,25,50,55,80,105), labels = c(-50,-25,0,25,50,0,25,50))

plot(allcounts_170$position_170, allcounts_170$A, type = "l", col="red", xaxt="n", ylab = "counts", xlab = "position", main = "length 70")
lines(allcounts_170$position_170, allcounts_170$T, col="purple")
lines(allcounts_170$position_170, allcounts_170$G, col="blue")
lines(allcounts_170$position_170, allcounts_170$C, col="green")
legend("topright", legend = c("A", "T","G","C"), col = c("red","purple","blue","green"),lty = c(1,1,1,1),cex=0.75)
axis(1, at=c(-50,-25,0,25,50,70,95,120), labels = c(-50,-25,0,25,50,0,25,50))

plot(allcounts_185$position_185, allcounts_185$A, type = "l", col="red", xaxt="n", ylab = "counts", xlab = "position", main = "length 85")
lines(allcounts_185$position_185, allcounts_185$T, col="purple")
lines(allcounts_185$position_185, allcounts_185$G, col="blue")
lines(allcounts_185$position_185, allcounts_185$C, col="green")
legend("topright", legend = c("A", "T","G","C"), col = c("red","purple","blue","green"),lty = c(1,1,1,1),cex=0.75)
axis(1, at=c(-50,-25,0,25,50,75,85,110,135), labels = c(-50,-25,0,25,50,75,0,25,50))

par(mfrow=c(1,1))
dev.off()

postscript("../../../mainfigures/Fig2B.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 10, height = 3, colormodel = "cmyk")

par(mfrow=c(1,3))

plot(allcounts_155$position_155[58:105], allcounts_155$CT[58:105], col="red",ylab = "counts",  type = "l", xlab = "position", main = "length 55" ,ylim = c(0,1000))
lines(allcounts_155$position_155[58:105], allcounts_155$AG[58:105],col="blue")
legend("bottomleft", legend = c("CT","AG"), col = c("red","blue"),lty = c(1,1),cex = 0.80)

plot(allcounts_170$position_170[58:120], allcounts_170$CT[58:120],col="red", ylab = "counts",  type = "l", xlab = "position", main = "length 70" ,ylim = c(0,500))
lines(allcounts_170$position_170[58:120], allcounts_170$AG[58:120],col="blue")
legend("bottomleft", legend = c("CT","AG"), col = c("red","blue"),lty = c(1,1),cex = 0.80)

plot(allcounts_185$position_185[58:135], allcounts_185$CT[58:135],col="red", ylab = "counts",  type = "l", xlab = "position", main = "length 85" ,ylim = c(0,110))
lines(allcounts_185$position_185[58:135], allcounts_185$AG[58:135],col="blue")
legend("bottomleft", legend = c("CT","AG"), col = c("red","blue"),lty = c(1,1),cex = 0.80)

par(mfrow=c(1,1))
dev.off()

################ ################ ################ ################ ################ 
################ FIGURE 4A-B ################ ################ ################ ################ 
chitest_trimer_5p_1 <- read_delim("../../data/chitest_trimer_5p_1", 
                                  "\t", escape_double = FALSE, trim_ws = TRUE)

chitest_trimer_5p_1$forward<- factor(chitest_trimer_5p_1$forward, levels = levels(reorder(chitest_trimer_5p_1$forward, chitest_trimer_5p_1$p_val)))
trimers_1<-levels(chitest_trimer_5p_1$forward)

pdf("../../../mainfigures/Fig4A.pdf", width = 8, height = 5)

plot(1:64,chitest_trimer_5p_1$p_val[order(chitest_trimer_5p_1$p_val)], xlab = "Trimers", ylab = "p-values", xaxt="n", col="gray48",pch=16)
axis(1, at=1:64, labels = FALSE)
text(seq(1, 64, by=1), par("usr")[3] - 0.04, labels = trimers_1,cex=0.65, srt = 90, pos = 1, xpd = TRUE)

dev.off()

CI_5p_trimer <- read_delim("../../data//CI_5p_trimer", 
                           "\t", escape_double = FALSE, trim_ws = TRUE)

idx<- c()
for (i in 1:length(trimers_1)) {
  idx[i]<-which(CI_5p_trimer$X1==trimers_1[i])
}
CI_5p_trimer_ro<- CI_5p_trimer[idx,]

pdf("../../../mainfigures/Fig4B.pdf", width = 8, height = 5)

plot(1:64, CI_5p_trimer_ro$Upper, type="l", ylim = c(0.3, 2.1), xlab = "Trimers", ylab = "Confidence Intervals", xaxt="n")
lines(1:64, CI_5p_trimer_ro$Lower)
abline(h=c(1/1.1,1.1), col="red", lty=2)
polygon(x=c(1:64, rev(1:64)), y=c(CI_5p_trimer_ro$Upper, rev(CI_5p_trimer_ro$Lower)), col=adjustcolor("gray48", alpha.f = 0.40))
axis(1, at=1:64, labels = FALSE)
text(seq(1, 64, by=1), par("usr")[3]- 0.08, labels = trimers_1,cex=0.65, srt = 90, pos = 1, xpd = TRUE)

dev.off()

################ ################ ################ ################ ################ 
################ FIGURE 5 ################ ################ ################ ################ 

## 5LR 
S_5LRperposasym <- read_delim("../../data/5LRperposasym", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)

p5<-ggplot(S_5LRperposasym, aes(x=position, y=Scores*100, colour=Asymmetry))+geom_line()+
  xlab("position")+ylab("Asymmetry (%)")+
  theme(legend.position = c(0.07,0.89), legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"))+
  scale_color_manual(labels= c(expression("S"[CG],"S"[TA])),values = c("red", "blue"))+ylim(-100,100)+geom_hline(yintercept = 0, linetype="dashed")
p5  

p5_inset<-ggplot(S_5LRperposasym, aes(x=position, y=Scores*100, colour=Asymmetry))+geom_point()+
  stat_smooth(method = "lm")+xlab("position")+ylab("Asymmetry (%)")+
  theme(legend.position = "none", legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"))+
  scale_color_manual(labels= c(expression("S"[CG],"S"[TA])),values = c("red", "blue"))+ylim(-100,100)+geom_hline(yintercept = 0, linetype="dashed")
p5_inset

## 3PT 

S_3PTperposasym <- read_delim("../../data/3PTperposasym", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)

p3<-ggplot(S_3PTperposasym, aes(x=position, y=Scores*100, colour=Asymmetry))+geom_line()+
  xlab("position")+ylab("Asymmetry (%)")+
  theme(legend.position = c(0.07,0.89), legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"))+
  scale_color_manual(labels= c(expression("S"[CG],"S"[TA])),values = c("red", "blue"))+ylim(-100,100)+geom_hline(yintercept = 0, linetype="dashed")
p3

p3_inset<-ggplot(S_3PTperposasym, aes(x=position, y=Scores*100, colour=Asymmetry))+geom_point()+
  stat_smooth(method = "lm")+xlab("position")+ylab("Asymmetry (%)")+
  theme(legend.position = "none", legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"))+
  scale_color_manual(labels= c(expression("S"[CG],"S"[TA])),values = c("red", "blue"))+ylim(-100,100)+geom_hline(yintercept = 0, linetype="dashed")
p3_inset

## 3' junction

S_3junctionperposasym <- read_delim("../../data/3junctionperposasym", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)


p3_junction<-ggplot(S_3junctionperposasym, aes(x=position, y=Scores*100, colour=Asymmetry))+geom_line()+
  xlab("position")+ylab("Asymmetry (%)")+
  theme(legend.position = c(0.07,0.89), legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"))+
  scale_color_manual(labels= c(expression("S"[CG],"S"[TA])),values = c("red", "blue"))+ylim(-100,100)+geom_hline(yintercept = 0, linetype="dashed")+
  geom_vline(xintercept = 4.5)
p3_junction

p3_all<-p3+annotation_custom(ggplotGrob(p3_inset), xmin = 6, xmax = 10,ymin = -100, ymax = -5)+
  annotate("text",x=8.25,y=-50,label=expression("Slope=0.0420, p-val=1.541e-06"),colour="red")+  annotate("text",x=8.25,y=-60,label=expression("Slope=0.0414, p-val=0.0047"),colour="blue")
p3_all

p5_all<-p5+annotation_custom(ggplotGrob(p5_inset), xmin=12, xmax=23, ymin=-100, ymax = -5)+
  annotate("text",x=18,y=-55,label=expression("Slope=1.843e-05, p-val=0.985"),colour="red")+  annotate("text",x=18,y=-65,label=expression("Slope=-0.000158, p-val=0.855"),colour="blue")
p5_all

postscript("../../../mainfigures/Fig5_5LR.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 4, colormodel = "cmyk")
p5+annotation_custom(ggplotGrob(p5_inset), xmin=11, xmax=23, ymin=-100, ymax = -5)+
  annotate("text",x=18.2,y=-55,label=expression("Slope=1.843e-05, p-val=0.985"),colour="red")+  annotate("text",x=18.2,y=-65,label=expression("Slope=-0.000158, p-val=0.855"),colour="blue")

dev.off()

postscript("../../../mainfigures/Fig5_3PT.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 4, colormodel = "cmyk")
p3+annotation_custom(ggplotGrob(p3_inset), xmin = 5, xmax = 10,ymin = -100, ymax = -5)+
  annotate("text",x=8,y=-50,label=expression("Slope=0.0420, p-val=1.541e-06"),colour="red")+  annotate("text",x=8.2,y=-60,label=expression("Slope=0.0414, p-val=0.0047"),colour="blue")
dev.off()

postscript("../../../mainfigures/Fig5_3junction.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 4, colormodel = "cmyk")
p3_junction
dev.off()

################ ################ ################ ################ ################ 
################ FIGURE 6A-B ################ ################ ################ ################ 
asymmetry_trimers_dro_obsexp <- read_delim("../../data/asymmetry_trimers_dro_obsexp", 
                                           delim = "\t", escape_double = FALSE, 
                                           trim_ws = TRUE)

asymmetry_trimers_dro_obsexp$forward<- factor(asymmetry_trimers_dro_obsexp$forward, levels = levels(reorder(asymmetry_trimers_dro_obsexp$forward, asymmetry_trimers_dro_obsexp$`3'-region`)))

asymmetry_trimers_long<- gather(asymmetry_trimers_dro_obsexp[,1:4], key="Region",value="value",3:4)

postscript("../../../mainfigures/Fig6A.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 4, colormodel = "cmyk")
ggplot(asymmetry_trimers_long, aes(x=forward, y=value*100, colour=Region))+geom_point()+
  theme(axis.text.x = element_text(size=8, angle=90), legend.position = c(0.07,0.89), legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"))+xlab("Trimers")+ylab("Asymmetry (%)")+
  scale_color_manual(values = c("orange", "green3"),labels= c("3PT", "5LR"))
dev.off()

####

observed<- c(asymmetry_trimers_dro_obsexp$`5'-region`, asymmetry_trimers_dro_obsexp$`3'-region`)
expected<- c(asymmetry_trimers_dro_obsexp$`5'-expected`, asymmetry_trimers_dro_obsexp$`3'-expected`)

asymmetry_trimers_df<- as.data.frame(cbind(asymmetry_trimers_long$forward, asymmetry_trimers_long$reverse,asymmetry_trimers_long$Region , observed,expected))
colnames(asymmetry_trimers_df)<- c("forward","reverse","Region","observed","expected")
sapply(asymmetry_trimers_df,class)
asymmetry_trimers_df[,4:5]<- lapply(asymmetry_trimers_df[,4:5], function(x) {as.numeric(as.character(x))})
asymmetry_trimers_df

postscript("../../../mainfigures/Fig6B.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 4, colormodel = "cmyk")
ggplot(asymmetry_trimers_df, aes(x=expected*100, y=observed*100, colour=Region))+geom_point()+
  stat_smooth(method=lm)+theme(legend.position = c(0.07,0.89), legend.title = element_blank(),
                               panel.background = element_rect(fill = "white", colour="black"),
                               panel.grid.minor = element_line(colour = "grey90"),
                               panel.grid.major = element_line(colour = "grey90"))+
  scale_color_manual(values = c("orange","green3"), labels=c("3PT","5LR"))+ annotate("text",x=80,y=-90,label=expression("R"^2*"=0.813"),colour="green3")+
  annotate("text",x=80,y=-80,label=expression("R"^2*"=0.932"),colour="orange")+xlab("Expected asymmetry (%)")+ylab("Observed asymmetry (%)")

dev.off()

################ ################ ################ ################ ################ 
################ FIGURE 7 ################ ################ ################ ################ 

CI_Dmel_trimer_gamma <- read_delim("../../data/CI_Dmel_trimer_gamma", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

postscript("../../../mainfigures/Fig7.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 4, colormodel = "cmyk")
ggplot(CI_Dmel_trimer_gamma, aes(x=reorder(motif, gamma), y=gamma))+geom_point()+
  theme(axis.text.x = element_text(size=8, angle=90),
        panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"))+xlab("Trimers")+
  geom_errorbar(aes(ymax=Upper, ymin=Lower))

dev.off()

################ ################ ################ ################ ################ 
################ FIGURE 8 ################ ################ ################ ################ 

CI_auto_3PT <- read_delim("../../data/auto_gammaHIII_3PT_selcoeffs_withBS", 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE)

auto_3PT_selcoeffs <- read_delim("../../data/auto_gammaHIII_3PT_selcoeffs", 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE)


postscript("../../../mainfigures/Fig8.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 7, height = 5, colormodel = "cmyk")

plot(1:9, auto_3PT_selcoeffs$sAG[1:9],type = "l", ylim = c(-3,1) ,xlab = "Focal position", ylab = "Gamma")
points(1:9, auto_3PT_selcoeffs$sAG[1:9], pch=16, cex=0.5)
arrows(1:9, CI_auto_3PT$UB[37:45],1:9, CI_auto_3PT$LB[37:45], angle=90, code=3, length=0.06)

lines(1:9, auto_3PT_selcoeffs$sA[1:9], col="red")
points(1:9, auto_3PT_selcoeffs$sA[1:9], col="red", pch=16,cex=0.5)
arrows(1:9, CI_auto_3PT$UB[1:9],1:9, CI_auto_3PT$LB[1:9], angle=90, code=3, length=0.06, col="red")

lines(1:9, auto_3PT_selcoeffs$sT[1:9], col="purple")
points(1:9, auto_3PT_selcoeffs$sT[1:9], col="purple", pch=16,cex=0.5)
arrows(1:9, CI_auto_3PT$UB[10:18],1:9, CI_auto_3PT$LB[10:18], angle=90, code=3, length=0.06, col="purple")

lines(1:9, auto_3PT_selcoeffs$sG[1:9], col="blue")
points(1:9, auto_3PT_selcoeffs$sG[1:9], col="blue", pch=16,cex=0.5)
arrows(1:9, CI_auto_3PT$UB[19:27],1:9, CI_auto_3PT$LB[19:27], angle=90, code=3, length=0.06, col="blue")

lines(1:9, auto_3PT_selcoeffs$sC[1:9], col="green")
points(1:9, auto_3PT_selcoeffs$sC[1:9], col="green", pch=16,cex=0.5)
arrows(1:9, CI_auto_3PT$UB[28:36],1:9, CI_auto_3PT$LB[28:36], angle=90, code=3, length=0.06, col="green")

legend("bottomleft", legend = c("AG","A", "T","G","C"), col = c("black","red","purple","blue","green"),lty = c(1,1,1,1,1),cex=0.6)

dev.off()

################ ################ ################ ################ ################ 
################ FIGURE 9 ################ ################ ################ ################ 

AllDeviationsJntCnt_autosome <- read_delim("../../data/AllDeviationsJntCnt_autosome", 
                                           delim = "\t", escape_double = FALSE, 
                                           trim_ws = TRUE)

deviation1to2<- AllDeviationsJntCnt_autosome[AllDeviationsJntCnt_autosome$pos_i==1 &AllDeviationsJntCnt_autosome$pos_j==2,]
deviation4to5<- AllDeviationsJntCnt_autosome[AllDeviationsJntCnt_autosome$pos_i==4 &AllDeviationsJntCnt_autosome$pos_j==5,]
deviation8to9<- AllDeviationsJntCnt_autosome[AllDeviationsJntCnt_autosome$pos_i==8 &AllDeviationsJntCnt_autosome$pos_j==9,]

deviation1to2$Base_i<- factor(deviation1to2$Base_i, levels = c("A","T","G","C"))
deviation1to2$Base_j<- factor(deviation1to2$Base_j, levels = c("A","T","G","C"))

postscript("../../../mainfigures/Fig9_1to2.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 9, height = 4, colormodel = "cmyk")

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

postscript("../../../mainfigures//Fig9_4to5.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 9, height = 4, colormodel = "cmyk")

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

postscript("../../../mainfigures//Fig9_8to9.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 9, height = 4, colormodel = "cmyk")

ggplot(deviation8to9,aes(x=Base_j,y=Base_i,fill=Deviation))+
  geom_tile()+
  scale_fill_gradient2(high="red",low="blue",mid="white")+xlab("Position 9")+ylab("Position 1-8")+
  theme(panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"))+
  facet_grid(.~Model)+geom_text(aes(label = round(Deviation,3)), size=2)

dev.off()


################ ################ ################ ################ ################ 
################ FIGURE 10 ################ ################ ################ ################ 

## Dimers
dimerasymm_PTcelegans <- read_delim("../../data/8Eukaryoutes/dimerasymm_PTcelegans", 
                                    "\t", escape_double = FALSE, trim_ws = TRUE)
dimerasymm_PTarabidopsis <- read_delim("../../data/8Eukaryoutes/dimerasymm_PTarabidopsis", 
                                       "\t", escape_double = FALSE, trim_ws = TRUE)
dimerasymm_PThuman <- read_delim("../../data/8Eukaryoutes/dimerasymm_PThuman", 
                                 "\t", escape_double = FALSE, trim_ws = TRUE)
dimerasymm_PTmoss <- read_delim("../../data/8Eukaryoutes/dimerasymm_PTmoss", 
                                "\t", escape_double = FALSE, trim_ws = TRUE)
dimerasymm_PTrice <- read_delim("../../data/8Eukaryoutes/dimerasymm_PTrice", 
                                "\t", escape_double = FALSE, trim_ws = TRUE)
dimerasymm_PTsu <- read_delim("../../data/8Eukaryoutes/dimerasymm_PTsu", 
                              "\t", escape_double = FALSE, trim_ws = TRUE)
dimerasymm_PTcerevisiae <- read_delim("../../data/8Eukaryoutes/dimerasymm_PTcerevisiae", 
                                      "\t", escape_double = FALSE, trim_ws = TRUE)
dimerasymm_PTthermo <- read_delim("../../data/8Eukaryoutes/dimerasymm_PTthermo", 
                                  "\t", escape_double = FALSE, trim_ws = TRUE)

dimerasymm_PTthermo$forward<- factor(dimerasymm_PTthermo$forward, levels = levels(reorder(dimerasymm_PTthermo$forward, dimerasymm_PTthermo$Asymmetry)))
dimerasymm_PTcerevisiae$forward<- factor(dimerasymm_PTcerevisiae$forward, levels = levels(reorder(dimerasymm_PTcerevisiae$forward, dimerasymm_PTcerevisiae$Asymmetry)))
dimerasymm_PTcelegans$forward<- factor(dimerasymm_PTcelegans$forward, levels = levels(reorder(dimerasymm_PTcelegans$forward, dimerasymm_PTcelegans$Asymmetry)))
dimerasymm_PTarabidopsis$forward<- factor(dimerasymm_PTarabidopsis$forward, levels = levels(reorder(dimerasymm_PTarabidopsis$forward, dimerasymm_PTarabidopsis$Asymmetry)))
dimerasymm_PThuman$forward<- factor(dimerasymm_PThuman$forward, levels = levels(reorder(dimerasymm_PThuman$forward, dimerasymm_PThuman$Asymmetry)))
dimerasymm_PTmoss$forward<- factor(dimerasymm_PTmoss$forward, levels = levels(reorder(dimerasymm_PTmoss$forward, dimerasymm_PTmoss$Asymmetry)))
dimerasymm_PTrice$forward<- factor(dimerasymm_PTrice$forward, levels = levels(reorder(dimerasymm_PTrice$forward, dimerasymm_PTrice$Asymmetry)))
dimerasymm_PTsu$forward<- factor(dimerasymm_PTsu$forward, levels = levels(reorder(dimerasymm_PTsu$forward, dimerasymm_PTsu$Asymmetry)))

species<- c(rep("Human",16),rep("Sea Urchin",16),rep("C. elegans",16),rep("Rice",16), rep("Arabidopsis",16),rep("Moss",16),rep("S. cerevisiae",16), rep("L. thermotolerans",16))
kingdom<- c(rep("Animal",48),rep("Plant",48), rep("Fungi",32))
dimer_all<- rbind(dimerasymm_PThuman, dimerasymm_PTsu,dimerasymm_PTcelegans, dimerasymm_PTrice, dimerasymm_PTarabidopsis, dimerasymm_PTmoss,dimerasymm_PTcerevisiae, dimerasymm_PTthermo)
dimer_all$species<- species
dimer_all$kingdom<- kingdom

dimer_all$species<- factor(dimer_all$species, levels = c("Human","Sea Urchin","Rice","Arabidopsis","Moss","C. elegans","S. cerevisiae","L. thermotolerans"))

dimeroneplot<-ggplot(dimer_all, aes(x=forward, y=Asymmetry*100, shape=species))+geom_point()+ylim(-100,100)+
  theme(panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"))+xlab("Dimers")+ylab("Asymmetry (%)")+
  scale_shape_manual(values=seq(0,8))


## Trimers

trimerasymm_PTcelegans <- read_delim("../../data/8Eukaryoutes/trimerasymm_PTcelegans", 
                                     "\t", escape_double = FALSE, trim_ws = TRUE)
trimerasymm_PTarabidopsis <- read_delim("../../data/8Eukaryoutes/trimerasymm_PTarabidopsis", 
                                        "\t", escape_double = FALSE, trim_ws = TRUE)
trimerasymm_PThuman <- read_delim("../../data/8Eukaryoutes/trimerasymm_PThuman", 
                                  "\t", escape_double = FALSE, trim_ws = TRUE)
trimerasymm_PTmoss <- read_delim("../../data/8Eukaryoutes/trimerasymm_PTmoss", 
                                 "\t", escape_double = FALSE, trim_ws = TRUE)
trimerasymm_PTrice <- read_delim("../../data/8Eukaryoutes/trimerasymm_PTrice", 
                                 "\t", escape_double = FALSE, trim_ws = TRUE)
trimerasymm_PTsu <- read_delim("../../data/8Eukaryoutes/trimerasymm_PTsu", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)
trimerasymm_PTcerevisiae <- read_delim("../../data/8Eukaryoutes/trimerasymm_PTcerevisiae", 
                                       "\t", escape_double = FALSE, trim_ws = TRUE)
trimerasymm_PTthermo <- read_delim("../../data/8Eukaryoutes/trimerasymm_PTthermo", 
                                   "\t", escape_double = FALSE, trim_ws = TRUE)


trimerasymm_PTcelegans$forward<- factor(trimerasymm_PTcelegans$forward, levels = levels(reorder(trimerasymm_PTcelegans$forward, trimerasymm_PTcelegans$Asymmetry)))
trimerasymm_PTarabidopsis$forward<- factor(trimerasymm_PTarabidopsis$forward, levels = levels(reorder(trimerasymm_PTarabidopsis$forward, trimerasymm_PTarabidopsis$Asymmetry)))
trimerasymm_PThuman$forward<- factor(trimerasymm_PThuman$forward, levels = levels(reorder(trimerasymm_PThuman$forward, trimerasymm_PThuman$Asymmetry)))
trimerasymm_PTmoss$forward<- factor(trimerasymm_PTmoss$forward, levels = levels(reorder(trimerasymm_PTmoss$forward, trimerasymm_PTmoss$Asymmetry)))
trimerasymm_PTrice$forward<- factor(trimerasymm_PTrice$forward, levels = levels(reorder(trimerasymm_PTrice$forward, trimerasymm_PTrice$Asymmetry)))
trimerasymm_PTsu$forward<- factor(trimerasymm_PTsu$forward, levels = levels(reorder(trimerasymm_PTsu$forward, trimerasymm_PTsu$Asymmetry)))
trimerasymm_PTcerevisiae$forward<- factor(trimerasymm_PTcerevisiae$forward, levels = levels(reorder(trimerasymm_PTcerevisiae$forward, trimerasymm_PTcerevisiae$Asymmetry)))
trimerasymm_PTthermo$forward<- factor(trimerasymm_PTthermo$forward, levels = levels(reorder(trimerasymm_PTthermo$forward, trimerasymm_PTthermo$Asymmetry)))


species<- c(rep("Human",64),rep("Sea Urchin",64),rep("C. elegans",64),rep("Rice",64), rep("Arabidopsis",64),rep("Moss",64),rep("S. cerevisiae",64), rep("L. thermotolerans",64))
kingdom<- c(rep("Animal",192),rep("Plant",192), rep("Fungi",128))
trimer_all<- rbind(trimerasymm_PThuman, trimerasymm_PTsu,trimerasymm_PTcelegans, trimerasymm_PTrice, trimerasymm_PTarabidopsis, trimerasymm_PTmoss,trimerasymm_PTcerevisiae, trimerasymm_PTthermo)
trimer_all$species<- species
trimer_all$kingdom<- kingdom


AGmotifs<- c("AAG","AGA","AGC","AGG","AGT","CAG","GAG","TAG")
indx<- c()
for (i in 1:8) {
  indx<- append(indx,which(trimer_all$forward==AGmotifs[i]))
}

category<- rep("nonAG",512)
category[indx]="AG"

trimer_all$species<- factor(trimer_all$species, levels = c("Human","Sea Urchin","Rice","Arabidopsis","Moss","C. elegans","S. cerevisiae","L. thermotolerans"))

trimer_all$category<- category

trimeroneplot<-ggplot(trimer_all, aes(x=forward, y=Asymmetry*100, colour=category,shape=species))+geom_point()+ylim(-100,100)+
  theme(axis.text.x = element_text(size=8, angle=90),panel.background = element_rect(fill = "white", colour="black"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"))+xlab("Trimers")+ylab("Asymmetry (%)")+
  scale_shape_manual(values=seq(0,8))


postscript("../../../mainfigures/Fig10.eps",horizontal = FALSE, onefile = FALSE, paper = "special", width = 8, height = 6, colormodel = "cmyk")

ggarrange(dimeroneplot, trimeroneplot, common.legend=TRUE, nrow = 2)

dev.off()

