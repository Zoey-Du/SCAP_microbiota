#A reproducible workflow of the analysis

#Title: Association of the lung microbiota on intensive care unit admission with the clinical improvements in patients with severe community acquired pneumonia 

#Authors: Sisi Du

#Date:2020/5/30

#Bioinformatic analysis
##Merge paired reads
for i in `tail -n+2 sample-metadata1.txt | cut -f 1`;do
vsearch --fastq_mergepairs seq/${i}_1.fastq --reverse seq/${i}_2.fastq \
--fastqout temp/${i}.merged.fastq --relabel ${i}.
done &
cat temp/*.merged.fastq > temp/all1.fastq

##Cut primers and quality filter
time vsearch --fastx_filter temp/all1.fastq \
--fastq_stripleft 17 --fastq_stripright 21 \
--fastq_maxee_rate 0.01 \
--fastaout temp/filtered1.fa

##Find unique read sequences and abundances
time vsearch --derep_fulllength temp/filtered1.fa \
--output temp/uniques.fa --relabel Uni --minuniquesize 8 --sizeout

##Denoise: predict biological sequences and filter chimeras
usearch -unoise3 temp/uniques.fa \
-zotus temp/zotus.fa
sed 's/Zotu/OTU_/g' temp/zotus.fa > temp/otus.fa

##Chimeras were filtered
silva=db/silva_16s_v123.fa
gunzip ${silva}.gz
ls -l ${silva}
grep -c '>' ${silva}
time vsearch --uchime_ref temp/otus.fa \
--db ${silva} \
--nonchimeras temp/otus_silva.fa
gzip ${silva} &
  mv temp/otus_silva.fa result/otus.fa

##Creat ZOTUs table
time vsearch --usearch_global temp/filtered1.fa --db result/otus.fa \
--otutabout result/otutab.txt --id 0.97 --threads 12

##Assign taxonomy
usearch -sintax result/otus.fa -db db/rdp_16s_v16_sp.fa \
-strand both -tabbedout temp/sintax.txt -sintax_cutoff 0.6

#Contaminants identification
library(phyloseq)
library(decontam)
otudec<-read.csv("D:/统计分析2/Dec.csv",header = T,row.names = 1)
meta<-read.csv("D:/统计分析2/samcon.csv",header = T,row.names = 1)
OTUDEC <-otu_table(otudec,taxa_are_rows = F)
META <-sample_data(meta)
phy <-phyloseq(OTUDEC,META)

#with method "frequency"
contamdf.freq <- isContaminant(phy, method="frequency", conc="Qun")
table(contamdf.freq$contaminant)
write.table(contamdf.freq,file='D:/统计分析2/contamdf.freq.txt',quote=FALSE,sep='\t',row.names=T,col.names = T)

#with method "prevalence"
sample_data(phy)$is.neg <- sample_data(phy)$Sam == "neg"
contamdf.prev <- isContaminant(phy, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
write.table(contamdf.prev,file='D:/统计分析2/contamdf.prev.txt',quote=FALSE,sep='\t',row.names=T,col.names = T)

#with method "both"
sample_data(phy)$is.neg <- sample_data(phy)$Sam== "neg"
contamdf.both <- isContaminant(phy, conc="Qun", neg="is.neg", method="both",threshold=c(0.1,0.5))
table(contamdf.both$contaminant)
which(contamdf.both$contaminant)
head(contamdf.both)
write.table(contamdf.both,file='D:/统计分析2/contamdf.both.txt',quote=FALSE,sep='\t',row.names=T,col.names = T)

##Removal of ZOTUs identified as contaminants with decontam package or observed in the controls and whose relative abundance was less than 0.01%

##Normlize by subsample
usearch -otutab_norm result/otutab.txt \
-sample_size 15000 \
-output result/otutab_norm.txt
usearch -otutab_stats result/otutab_norm.txt \
-output result/otutab_norm.stat
cat result/otutab_norm.stat


##Figure 2. Factors associated with the lung microbiota composition
dataall <-read.csv("D:/analysis/otutab.csv",header = T,row.names = 1)
groupall =read.table("D:/analysis/group.txt", header=T, row.names= 1,sep="\t")
library(ggplot2)
library(vegan)
adonis(dataall~AdQtime2,data=groupall,permutations = 999, method="bray") -> au
au

pairwise.adonis <-function(x,factors, sim.method, p.adjust.m)
  
{
  library(vegan)
  
  co = as.matrix(combn(unique(factors),2))
  
  pairs = c()
  
  F.Model =c()
  
  R2 = c()
  
  p.value = c()
  
  
  
  for(elem in 1:ncol(co)){
    
    ad = adonis(x[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem])),] ~
                  
                  factors[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem]))] , method =sim.method);
    
    pairs =c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    
    R2 = c(R2,ad$aov.tab[1,5]);
    
    p.value = c(p.value,ad$aov.tab[1,6])
    
  }
  
  p.adjusted =p.adjust(p.value,method=p.adjust.m)
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  
  return(pairw.res)
  
}

pairwise.adonis(dataall, groupall$AdQtime2, sim.method="bray", p.adjust.m= "fdr")
###Independent association
adonis2(dataall~Bacter+AcuteKid+PCT+MIP_1beta+AdQtime2+MVAP+SmokeStatu+Sputum,data=groupall,permutations = 999, method="bray",by='term') -> au
au

##MIP_1beta PCoA
adonis(dataall~MIP_1beta,data=groupall,permutations = 999, method="bray") -> au
au
Xibaopcoa <- vegdist(dataall,method = 'bray')
Xibaopcoa1 <- as.matrix(Xibaopcoa)
idx =rownames(groupall) %in% colnames(Xibaopcoa1)
sub_design =groupall[idx,]
Xibaopcoa1 =Xibaopcoa1[rownames(sub_design), rownames(sub_design)]
pcoaXibao =cmdscale(Xibaopcoa1, k=2, eig=T)
pointXibao = as.data.frame(pcoaXibao$points) 
eigXibao = pcoaXibao$eig
pointXibao = cbind(pointXibao, sub_design$MIP_1beta)
colnames(pointXibao) = c("PC1", "PC2","MIP_1beta")

p = ggplot(pointXibao, aes(x=PC1, y=PC2, size=MIP_1beta)) + geom_point(alpha=0.8,color="#CC66FF") +
  labs(x=paste("PCoA 1 (", format(100 * eigXibao[1] / sum(eigXibao), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eigXibao[2] / sum(eigXibao), digits=4), "%)", sep=""),
       title="PCoA") + theme_classic()+
  scale_size(range=c(1,25))+
  theme_bw()+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size = 15))+
  theme(legend.text = element_text(size=15),legend.title = element_text(size = 15),legend.key.size = unit(1,'cm'))

##PCT PCoA
adonis(dataall~PCT,data=groupall,permutations = 999, method="bray") -> au
au
PCTpcoa <- vegdist(dataall,method = 'bray')
PCTpcoa1 <- as.matrix(PCTpcoa)
idx =rownames(groupall) %in% colnames(PCTpcoa1)
sub_design =groupall[idx,]
PCTpcoa1 =PCTpcoa1[rownames(sub_design), rownames(sub_design)]
pcoaPCT =cmdscale(PCTpcoa1, k=2, eig=T)
pointPCT = as.data.frame(pcoaPCT$points) 
eigPCT = pcoaPCT$eig
pointPCT = cbind(pointPCT, sub_design$PCT)
colnames(pointPCT) = c("PC1", "PC2","PCT")

p = ggplot(pointPCT, aes(x=PC1, y=PC2, size=PCT)) + geom_point(alpha=0.8,color="green") +
  labs(x=paste("PCoA 1 (", format(100 * eigPCT[1] / sum(eigPCT), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eigPCT[2] / sum(eigPCT), digits=4), "%)", sep=""),
       title="PCoA") + theme_classic()+
  scale_size(range=c(1,10))+
  theme_bw()+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size = 15))+
  theme(legend.text = element_text(size=15),legend.title = element_text(size = 15),legend.key.size = unit(1,'cm'))

##AKI PCoA
groupall$AcuteKid <- factor(groupall$AcuteKid)
adonis(dataall~AcuteKid,data=groupall,permutations = 999, method="bray") -> au
au
Acutpcoa <- vegdist(dataall,method = 'bray')
Acutpcoa1 <- as.matrix(Acutpcoa)
idx =rownames(groupall) %in% colnames(Acutpcoa1)
sub_design =groupall[idx,]
Acutpcoa1 =Acutpcoa1[rownames(sub_design), rownames(sub_design)]
pcoaacutk =cmdscale(Acutpcoa1, k=2, eig=T)
pointacutk = as.data.frame(pcoaacutk$points) 
eigacutk = pcoaacutk$eig
levels(sub_design$AcuteKid)=c("AKI-present","AKI-absent")
pointacutk = cbind(pointacutk, sub_design$AcuteKid)
colnames(pointacutk) = c("PC1", "PC2","AcuteKid")
View(pointacutk)
p = ggplot(pointacutk, aes(x=PC1, y=PC2, color=AcuteKid)) + geom_point(alpha=0.5, size=5) +
  labs(x=paste("PCoA 1 (", format(100 * eigacutk[1] / sum(eigacutk), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eigacutk[2] / sum(eigacutk), digits=4), "%)", sep=""),
       title="PCOA") + theme_classic()+
  scale_color_manual(values = c('red','blue'))+
  theme_bw()+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size = 15))+
  guides(color=F)

p + geom_enterotype(aes(color=AcuteKid,label=AcuteKid),show.legend = T)+
  scale_color_manual(values = c('red','blue'))

##Bacterial detection PCoA
groupall$Bacter<- factor(groupall$Bacter)
adonis(dataall~Bacter,data=groupall,permutations = 999, method="bray") -> au
au
Bacterpcoa <- vegdist(dataall,method = 'bray')
Bacterpcoa1 <- as.matrix(Bacterpcoa)
idx =rownames(groupall) %in% colnames(Stomkpcoa1)
sub_design =groupall[idx,]
Bacterpcoa1 =Bacterpcoa1[rownames(sub_design), rownames(sub_design)]
pcoaBacter =cmdscale(Bacterpcoa1, k=2, eig=T)
pointBacter = as.data.frame(pcoaBacter$points) 
eigBacter = pcoaBacter$eig
levels(sub_design$Bacter)=c("Positive-detection","Negative-detection")
pointBacter = cbind(pointBacter, sub_design$Bacter)
colnames(pointBacter) = c("PC1", "PC2","Bacter")

p = ggplot(pointBacter, aes(x=PC1, y=PC2, color=Bacter)) + geom_point(alpha=0.5, size=5) +
  labs(x=paste("PCoA 1 (", format(100 * eigBacter[1] / sum(eigBacter), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eigBacter[2] / sum(eigBacter), digits=4), "%)", sep=""),
       title="PCOA") + theme_classic()+
  scale_color_manual(values = c('red','blue'))+
  theme_bw()+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size = 15))+
  guides(color=F)
p + geom_enterotype(aes(color=Bacter,label=Bacter), show.legend = T)+
  scale_color_manual(values = c('red','blue'))

##Bronchoscopy time PCoA
groupall$AdQtime2 <- factor(groupall$AdQtime2)
adonis(dataall~AdQtime2,data=groupall,permutations = 999, method="bray") -> au
au
Qipcoa <- vegdist(dataall,method = 'bray')
Qipcoa1 <- as.matrix(Qipcoa)
idx =rownames(groupall) %in% colnames(Qipcoa1)
sub_design =groupall[idx,]
Qipcoa1 =Qipcoa1[rownames(sub_design), rownames(sub_design)]
pcoaQi =cmdscale(Qipcoa1, k=2, eig=T)
pointQi = as.data.frame(pcoaQi$points) 
eigQi = pcoaQi$eig
levels(sub_design$AdQtime2)=c(">24 hours","12-24 hours","<12 hours")
pointQi = cbind(pointQi, sub_design$AdQtime2)
colnames(pointQi) = c("PC1", "PC2","AdQtime2")

p = ggplot(pointQi, aes(x=PC1, y=PC2, color=AdQtime2)) + geom_point(alpha=0.5, size=5) +
  labs(x=paste("PCoA 1 (", format(100 * eigQi[1] / sum(eigQi), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eigQi[2] / sum(eigQi), digits=4), "%)", sep=""),
       title="PCoA") + theme_classic()+
  scale_color_manual(values = c('blue','red','#FFFF00'))+
  theme_bw()+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size = 15))+
  guides(color=F)

p + geom_enterotype(aes(color=AdQtime2,label=AdQtime2), show.legend = T)+
  scale_color_manual(values = c('blue','red','#FFFF00'))


##Co-association between 8 factors
corfactor =read.table("D:/analysis/cor.txt", header=T, row.names= 1,sep="\t")
XX <- cor(corfactor,method = "pearson")
##Kappa test
kappa(XX)
##Corrplot
library(corrplot)
res <- cor.mtest(corfactor,conf.level = 0.95)
corrplot(XX, type = "upper", tl.pos = "d")
corrplot(XX, add = TRUE, type = "lower",p.mat = res$p, low = res$lowCI, upp = res$uppCI, pch.col = "red", sig.level = 0.05, addrect = 3,rect.col = "navy", diag = FALSE, plotCI = "rect",tl.pos = "n", cl.pos = "n")

##Figure 3. The influence of bacteriology results on lung microbiota 

##Random Forest
Family_tab =read.table("D:/analysis/Family_tab.txt", header=T, row.names= 1,sep="\t")
library(randomForest)
set.seed(315)
rf = randomForest(t(Family_tab), groupall$Bacter, importance=TRUE, proximity=TRUE, ntree = 1000)
print(rf)
set.seed(315) 
result = rfcv(t(Family_tab), groupall$Bacter, cv.fold=5)
result$error.cv
with(result, plot(n.var,error.cv,log="x", type="o", lwd=2))
imp= as.data.frame(rf$importance)
imp = imp[order(imp[,1],decreasing = T),]
head(imp,n=10)
write.table(imp,file = "importance_class.txt",quote = F,sep = '\t', row.names = T, col.names = T)

Bactpoint=read.table("D:/analysis/Bactpoint.txt", header=T, row.names= 1,sep="\t")
Bactrect=read.table("D:/analysis/Bactrect.txt", header=T, row.names= 1,sep="\t")

ggplot() + 
  geom_rect(mapping = aes(xmin=Xnmin,xmax=Xnmax,ymin=Meanmin,ymax=Meanmax),fill="blue",alpha=1/2,data = Bactrect)+
  geom_point(Bactpoint,mapping = aes(OR,Poc,color=-log10(RE),size=-log10(RE)),shape=15)+
  scale_color_gradient(low = "red",high = "yellow")+
  scale_size(range = c(10, 3))+
  geom_point(Bactpoint,mapping = aes(OR,Pva),size=2)+
  theme_classic()+
  scale_y_continuous(limits=c(1,2.5),breaks=seq(1,2.5,by=0.1))+
  scale_x_continuous(breaks=Bactrect$Brea,labels = reorder(Bactrect$Family,Bactrect$Xnmin))+
  coord_flip()


##Comprosion of relative abundance
fit3 <- glm(Prevotellaceae~AdQtime,data = Family_tab, family = gaussian())
summary(fit3)
## and
wilcox.test()

###Manhattan plot
otu_stat <- read.delim('D:/analysis/manbac.txt', sep = '\t')
p <- ggplot(otu_stat, aes(OR, -log(Pvalue, 10))) +	#
  geom_point(aes( size = Abun,color = factor(Rich), shape = factor(Rich))) +	
  scale_size(range = c(2, 10))+
  scale_shape_manual(limits = c('rich', 'low','no'), values = c(16, 10,1)) +		
  scale_color_manual(limits = c('rich', 'low','no'),values = c("red","blue","green"))+
  scale_y_continuous(limits=c(0,3))+
  labs(x = NULL, y = '-log10(P)', size = 'relative abundance (%)') +	  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), panel.background = element_rect(fill = 'transparent'), legend.key = element_rect(fill = 'transparent')) +	
  geom_hline(yintercept = -log10(0.05), color = 'gray', linetype = 2, size = 1)	+
  scale_x_continuous(breaks=otu_stat$Break,labels = reorder(otu_stat$Lab,otu_stat$OR))+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

##alpha diversity comparison
Bacteralp =read.table("D:/analysis/Bacteralpha.txt", header=T, row.names= 1,sep="\t")
g1<-ggplot(data =Bacteralp,aes(x=R1,y=Shannon,group=Bacter) )+
  stat_boxplot(geom='errorbar',width=0.6)+
  geom_boxplot(aes(fill=Bacter))+
  geom_jitter(shape=1)+
  scale_x_continuous(limits=c(0,6),breaks=seq(0,6,by=0.5),expand=c(0,0))+
  scale_fill_manual(values=c('blue','red'))+
  theme_classic()

g2<-ggplot(data =Bacteralp,aes(x=R2,y=Richness,group=Bacter) )+
  stat_boxplot(geom='errorbar',width=0.6)+
  geom_boxplot(aes(fill=Bacter))+
  geom_jitter(shape=1)+
  scale_x_continuous(limits=c(0,6),breaks=seq(0,6,by=0.5),expand=c(0,0))+
  scale_fill_manual(values=c('blue','red'))+
  theme_classic()%+replace%
  theme(panel.background = element_rect(fill = "transparent",colour = NA))

ggplot2.two_y_axis(g1, g2)

##Figure 4.Lung microbiota predicting clinical improvements.
library("survival")
library("survminer")
coxall =read.table("D:/analysis/log.txt", header=T, row.names= 1,sep="\t")
###Multiple analysis
res.cox <- coxph(Surv(D14time, D14) ~Season+IL_8R+CURB65+APACHEII+Shock+OxygenIndex+CRR+Bacter, data = coxall)
res.cox
summary(res.cox)
ftest <- cox.zph(res.cox)
ftest
coxall$ActinomycetaceaeN <- factor(coxall$ActinomycetaceaeN)
res.cox <- coxph(Surv(D14time, D14) ~ActinomycetaceaeN, data = coxall)
res.cox
summary(res.cox)
####Plot
fit<- survfit(Surv(D14time, D14) ~RichnessN, data = coxall)
ggsurvplot(fit, pval = TRUE,data = coxall,risk.table = 'abs_pct',palette = c('red','blue','#FFFF00'),xlim=c(0,14), break.time.by=7, ggtheme = theme_bw())

##Supplementary Material
##Figure 1. Factors not associated with the lung microbiota composition

##Age PCoA
adonis(dataall~Age,data=groupall,permutations = 999, method="bray") -> au
au
Agepcoa <- vegdist(dataall,method = 'bray')
Agepcoa1 <- as.matrix(Agepcoa)
idx =rownames(groupall) %in% colnames(Agepcoa1)
sub_design =groupall[idx,]
Agepcoa1 =Agepcoa1[rownames(sub_design), rownames(sub_design)]
pcoaAge =cmdscale(Agepcoa1, k=2, eig=T)
pointAge = as.data.frame(pcoaAge$points) 
eigAge = pcoaAge$eig
pointAge = cbind(pointAge, sub_design$Age)
colnames(pointAge) = c("PC1", "PC2","Age")

p = ggplot(pointAge, aes(x=PC1, y=PC2, size=Age)) + geom_point(alpha=0.8,color="#CC66FF") +
  labs(x=paste("PCoA 1 (", format(100 * eigAge[1] / sum(eigAge), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eigAge[2] / sum(eigAge), digits=4), "%)", sep=""),
       title="PCoA") + theme_classic()+
  scale_size(range=c(1,5))+
  theme_bw()+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size = 15))+
  theme(legend.text = element_text(size=15),legend.title = element_text(size = 15),legend.key.size = unit(1,'cm'))

##Gender PCoA
groupall$Gender <- factor(groupall$Gender)
adonis(dataall~Gender,data=groupall,permutations = 999, method="bray") -> au
au
Genderpcoa <- vegdist(dataall,method = 'bray')
Genderpcoa1 <- as.matrix(Genderpcoa)
idx =rownames(groupall) %in% colnames(Genderpcoa1)
sub_design =groupall[idx,]
Genderpcoa1 =Genderpcoa1[rownames(sub_design), rownames(sub_design)]
pcoaGender =cmdscale(Genderpcoa1, k=2, eig=T)
pointGender = as.data.frame(pcoaGender$points) 
eigGender = pcoaGender$eig
levels(sub_design$Gender)=c("Male","Female")
pointGender = cbind(pointGender, sub_design$Gender)
colnames(pointGender) = c("PC1", "PC2","Gender")

p = ggplot(pointGender, aes(x=PC1, y=PC2, color=Gender)) + geom_point(alpha=0.5, size=5) +
  labs(x=paste("PCoA 1 (", format(100 * eigGender[1] / sum(eigGender), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eigGender[2] / sum(eigGender), digits=4), "%)", sep=""),
       title="PCoA") + theme_classic()+
  scale_color_manual(values = c('blue','red'))+
  theme_bw()+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size = 15))+
  guides(color=F)

p + geom_enterotype(aes(color=Gender,label=Gender), show.legend = T)+
  scale_color_manual(values = c('blue','red'))

##Chronic respiratory disease
groupall$ResDisease <- factor(groupall$ResDisease)
adonis(dataall~ResDisease,data=groupall,permutations = 999, method="bray") -> au
au
ResDiseasepcoa <- vegdist(dataall,method = 'bray')
ResDiseasepcoa1 <- as.matrix(ResDiseasepcoa)
idx =rownames(groupall) %in% colnames(ResDiseasepcoa1)
sub_design =groupall[idx,]
ResDiseasepcoa1 =ResDiseasepcoa1[rownames(sub_design), rownames(sub_design)]
pcoaResDisease=cmdscale(ResDiseasepcoa1, k=2, eig=T)
pointResDisease = as.data.frame(pcoaResDisease$points) 
eigResDisease = pcoaResDisease$eig
levels(sub_design$ResDisease)=c("Diseases-present","Diseases-absent")
pointResDisease = cbind(pointResDisease, sub_design$ResDisease)
colnames(pointResDisease) = c("PC1", "PC2","ResDisease")

p = ggplot(pointResDisease, aes(x=PC1, y=PC2, color=ResDisease)) + geom_point(alpha=0.5, size=5) +
  labs(x=paste("PCoA 1 (", format(100 * eigResDisease[1] / sum(eigResDisease), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eigResDisease[2] / sum(eigResDisease), digits=4), "%)", sep=""),
       title="PCoA") + theme_classic()+
  scale_color_manual(values = c('red','blue'))+
  theme_bw()+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size = 15))+
  guides(color=F)

p + geom_enterotype(aes(color=ResDisease,label=ResDisease), show.legend = T)+
  scale_color_manual(values = c('red','blue'))

##Immunosuppression status PCoA
groupall$Immunosuppression <- factor(groupall$Immunosuppression)
adonis(dataall~Immunosuppression,data=groupall,permutations = 999, method="bray") -> au
au
Immunosuppressionpcoa <- vegdist(dataall,method = 'bray')
Immunosuppressionpcoa1 <- as.matrix(Immunosuppressionpcoa)
idx =rownames(groupall) %in% colnames(Immunosuppressionpcoa1)
sub_design =groupall[idx,]
Immunosuppressionpcoa1 =Immunosuppressionpcoa1[rownames(sub_design), rownames(sub_design)]
pcoaImmunosuppression=cmdscale(Immunosuppressionpcoa1, k=2, eig=T)
pointImmunosuppression = as.data.frame(pcoaImmunosuppression$points) 
eigImmunosuppression = pcoaImmunosuppression$eig
levels(sub_design$Immunosuppression)=c("Immunocompromised-status","Immunocompetent-status")
pointImmunosuppression = cbind(pointImmunosuppression, sub_design$Immunosuppression)
colnames(pointImmunosuppression) = c("PC1", "PC2","Immunosuppression")

p = ggplot(pointImmunosuppression, aes(x=PC1, y=PC2, color=Immunosuppression)) + geom_point(alpha=0.5, size=5) +
  labs(x=paste("PCoA 1 (", format(100 * eigImmunosuppression[1] / sum(eigImmunosuppression), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eigImmunosuppression[2] / sum(eigImmunosuppression), digits=4), "%)", sep=""),
       title="PCoA") + theme_classic()+
  scale_color_manual(values = c('red','blue'))+
  theme_bw()+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size = 15))+
  guides(color=F)

p + geom_enterotype(aes(color=Immunosuppression,label=Immunosuppression), show.legend = T)+
  scale_color_manual(values = c('red','blue'))

####Sample Season PCoA
groupall$Season <- factor(groupall$Season)
adonis(dataall~Season,data=groupall,permutations = 999, method="bray") -> au
au
Seasonpcoa <- vegdist(dataall,method = 'bray')
Seasonpcoa1 <- as.matrix(Seasonpcoa)
idx =rownames(groupall) %in% colnames(Seasonpcoa1)
sub_design =groupall[idx,]
Seasonpcoa1 =Seasonpcoa1[rownames(sub_design), rownames(sub_design)]
pcoaSeason=cmdscale(Seasonpcoa1, k=2, eig=T)
pointSeason = as.data.frame(pcoaSeason$points) 
eigSeason = pcoaSeason$eig
levels(sub_design$Season)=c("Spring","Summer","Autumn","Winter")
pointSeason = cbind(pointSeason, sub_design$Season)
colnames(pointSeason) = c("PC1", "PC2","Season")

p = ggplot(pointSeason, aes(x=PC1, y=PC2, color=Season)) + geom_point(alpha=0.5, size=5) +
  labs(x=paste("PCoA 1 (", format(100 * eigSeason[1] / sum(eigSeason), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eigSeason[2] / sum(eigSeason), digits=4), "%)", sep=""),
       title="PCoA") + theme_classic()+
  scale_color_manual(values = c('#FFFF00','red','#CC66FF','blue'))+
  theme_bw()+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size = 15))+
  guides(color=F)

p + geom_enterotype(aes(color=Season,label=Season), show.legend = T)+
  scale_color_manual(values = c('#FFFF00','red','#CC66FF','blue'))

####PSI Range PCoA
groupall$PSIRange <- factor(groupall$PSIRange)
adonis(dataall~PSIRange,data=groupall,permutations = 999, method="bray") -> au
au
PSIRangepcoa <- vegdist(dataall,method = 'bray')
PSIRangepcoa1 <- as.matrix(PSIRangepcoa)
idx =rownames(groupall) %in% colnames(PSIRangepcoa1)
sub_design =groupall[idx,]
PSIRangepcoa1 =PSIRangepcoa1[rownames(sub_design), rownames(sub_design)]
pcoaPSIRange=cmdscale(PSIRangepcoa1, k=2, eig=T)
pointPSIRange = as.data.frame(pcoaPSIRange$points) 
eigPSIRange = pcoaPSIRange$eig
levels(sub_design$PSIRange)=c("V","IV","I","III","II")
pointPSIRange = cbind(pointPSIRange, sub_design$PSIRange)
colnames(pointPSIRange) = c("PC1", "PC2","PSIRange")

p = ggplot(pointPSIRange, aes(x=PC1, y=PC2, color=PSIRange)) + geom_point(alpha=0.5, size=5) +
  labs(x=paste("PCoA 1 (", format(100 * eigPSIRange[1] / sum(eigPSIRange), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eigPSIRange[2] / sum(eigPSIRange), digits=4), "%)", sep=""),
       title="PCoA") + theme_classic()+
  scale_color_manual(values = c('#CC66FF','#FFFF00','red','blue','green'))+
  theme_bw()+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size = 15))+
  guides(color=F)
p + geom_enterotype(aes(color=PSIRange,label=PSIRange), show.legend = T)+
  scale_color_manual(values = c('#CC66FF','#FFFF00','red','blue','green'))

##CURB-65 PCoA
adonis(dataall~CURB65,data=groupall,permutations = 999, method="bray") -> au
au
CURB65pcoa <- vegdist(dataall,method = 'bray')
CURB65pcoa1 <- as.matrix(CURB65pcoa)
idx =rownames(groupall) %in% colnames(CURB65pcoa1)
sub_design =groupall[idx,]
CURB65pcoa1 =CURB65pcoa1[rownames(sub_design), rownames(sub_design)]
pcoaCURB65 =cmdscale(CURB65pcoa1, k=2, eig=T)
pointCURB65 = as.data.frame(pcoaCURB65$points) 
eigCURB65 = pcoaCURB65$eig
pointCURB65 = cbind(pointCURB65, sub_design$CURB65)
colnames(pointCURB65) = c("PC1", "PC2","CURB65")

p = ggplot(pointCURB65, aes(x=PC1, y=PC2, size=CURB65)) + geom_point(alpha=0.8,color="green") +
  labs(x=paste("PCoA 1 (", format(100 * eigCURB65[1] / sum(eigCURB65), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eigCURB65[2] / sum(eigCURB65), digits=4), "%)", sep=""),
       title="PCoA") + theme_classic()+
  scale_size(range=c(1,6))+
  theme_bw()+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size = 15))+
  theme(legend.text = element_text(size=15),legend.title = element_text(size = 15),legend.key.size = unit(1,'cm'))

##IL_6 PCoA
adonis(dataall~IL_6,data=groupall,permutations = 999, method="bray") -> au
au
IL_6pcoa <- vegdist(dataall,method = 'bray')
IL_6pcoa1 <- as.matrix(IL_6pcoa)
idx =rownames(groupall) %in% colnames(IL_6pcoa1)
sub_design =groupall[idx,]
IL_6pcoa1 =IL_6pcoa1[rownames(sub_design), rownames(sub_design)]
pcoaIL_6 =cmdscale(IL_6pcoa1, k=2, eig=T)
pointIL_6 = as.data.frame(pcoaIL_6$points) 
eigIL_6 = pcoaIL_6$eig
pointIL_6 = cbind(pointIL_6, sub_design$IL_6)
colnames(pointIL_6) = c("PC1", "PC2","IL_6")

p = ggplot(pointIL_6, aes(x=PC1, y=PC2, size=IL_6)) + geom_point(alpha=0.8,color="blue") +
  labs(x=paste("PCoA 1 (", format(100 * eigIL_6[1] / sum(eigIL_6), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eigIL_6[2] / sum(eigIL_6), digits=4), "%)", sep=""),
       title="IL_6 PCoA",color='IL_6') + theme_classic()+
  scale_size(range=c(2,10))+
  theme_bw()+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size = 15))+
  theme(legend.text = element_text(size=15),legend.title = element_text(size = 15),legend.key.size = unit(1,'cm'))

####Antibiotic PCoA
groupall$Antisum <- factor(groupall$Antisum)
adonis(dataall~Antisum,data=groupall,permutations = 999, method="bray") -> au
au
Antisumpcoa <- vegdist(dataall,method = 'bray')
Antisumpcoa1 <- as.matrix(Antisumpcoa)
idx =rownames(groupall) %in% colnames(Antisumpcoa1)
sub_design =groupall[idx,]
Antisumpcoa1 =Antisumpcoa1[rownames(sub_design), rownames(sub_design)]
pcoaAntisum=cmdscale(Antisumpcoa1, k=2, eig=T)
pointAntisum = as.data.frame(pcoaAntisum$points) 
eigAntisum = pcoaAntisum$eig
levels(sub_design$Antisum)=c("Beta-lactams","Beta-lactams+Fluoroquinolones","Carbapenems","Fluoroquinolones")
pointAntisum = cbind(pointAntisum, sub_design$Antisum)
colnames(pointAntisum) = c("PC1", "PC2","Antisum")

p = ggplot(pointAntisum, aes(x=PC1, y=PC2, color=Antisum)) + geom_point(alpha=0.5, size=5) +
  labs(x=paste("PCoA 1 (", format(100 * eigAntisum[1] / sum(eigAntisum), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eigAntisum[2] / sum(eigAntisum), digits=4), "%)", sep=""),
       title="PCoA") + theme_classic()+
  scale_color_manual(values = c('#FFFF00','blue','red','#CC66FF'))+
  theme_bw()+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size = 15))+
  guides(color=F)
p + geom_enterotype(aes(color=Antisum,label=Antisum), show.legend = T)+
  scale_color_manual(values = c('#FFFF00','blue','red','#CC66FF'))

##Virus PCoA
groupall$Virus <- factor(groupall$Virus)
adonis(dataall~Virus,data=groupall,permutations = 999, method="bray") -> au
au
Viruspcoa <- vegdist(dataall,method = 'bray')
Viruspcoa1 <- as.matrix(Viruspcoa)
idx =rownames(groupall) %in% colnames(Viruspcoa1)
sub_design =groupall[idx,]
Viruspcoa1 =Viruspcoa1[rownames(sub_design), rownames(sub_design)]
pcoaVirus=cmdscale(Viruspcoa1, k=2, eig=T)
pointVirus = as.data.frame(pcoaVirus$points) 
eigVirus = pcoaVirus$eig
levels(sub_design$Virus)=c("Positive-detection","Negative-detection")
pointVirus = cbind(pointVirus, sub_design$Virus)
colnames(pointVirus) = c("PC1", "PC2","Virus")

p = ggplot(pointVirus, aes(x=PC1, y=PC2, color=Virus)) + geom_point(alpha=0.5, size=5) +
  labs(x=paste("PCoA 1 (", format(100 * eigVirus[1] / sum(eigVirus), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eigVirus[2] / sum(eigVirus), digits=4), "%)", sep=""),
       title="PCoA") + theme_classic()+
  scale_color_manual(values = c('red','blue'))+
  theme_bw()+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size = 15))+
  guides(color=F)
p + geom_enterotype(aes(color=Virus,label=Virus), show.legend = T)+
  scale_color_manual(values = c('red','blue'))


####ARDS PCoA
groupall$ARDS <- factor(groupall$ARDS)
adonis(dataall~ARDS,data=groupall,permutations = 999, method="bray") -> au
au
ARDSpcoa <- vegdist(dataall,method = 'bray')
ARDSpcoa1 <- as.matrix(ARDSpcoa)
idx =rownames(groupall) %in% colnames(ARDSpcoa1)
sub_design =groupall[idx,]
ARDSpcoa1 =ARDSpcoa1[rownames(sub_design), rownames(sub_design)]
pcoaARDS=cmdscale(ARDSpcoa1, k=2, eig=T)
pointARDS = as.data.frame(pcoaARDS$points) 
eigARDS = pcoaARDS$eig
levels(sub_design$ARDS)=c("ARDS-present","ARDS-absent")
pointARDS = cbind(pointARDS, sub_design$ARDS)
colnames(pointARDS) = c("PC1", "PC2","ARDS")

p = ggplot(pointARDS, aes(x=PC1, y=PC2, color=ARDS)) + geom_point(alpha=0.5, size=5) +
  labs(x=paste("PCoA 1 (", format(100 * eigARDS[1] / sum(eigARDS), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eigARDS[2] / sum(eigARDS), digits=4), "%)", sep=""),
       title="PCoA") + theme_classic()+
  scale_color_manual(values = c('red','blue'))+
  theme_bw()+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size = 15))+
  guides(color=F)
p + geom_enterotype(aes(color=ARDS,label=ARDS), show.legend = T)+
  scale_color_manual(values = c('red','blue'))

####Shock
groupall$Shock <- factor(groupall$Shock)
adonis(dataall~Shock,data=groupall,permutations = 999, method="bray") -> au
au
Shockpcoa <- vegdist(dataall,method = 'bray')
Shockpcoa1 <- as.matrix(Shockpcoa)
idx =rownames(groupall) %in% colnames(Shockpcoa1)
sub_design =groupall[idx,]
Shockpcoa1 =Shockpcoa1[rownames(sub_design), rownames(sub_design)]
pcoaShock=cmdscale(Shockpcoa1, k=2, eig=T)
pointShock = as.data.frame(pcoaShock$points) 
eigShock = pcoaShock$eig
levels(sub_design$Shock)=c("Shock-present","Shock-absent")
pointShock = cbind(pointShock, sub_design$Shock)
colnames(pointShock) = c("PC1", "PC2","Shock")

p = ggplot(pointShock, aes(x=PC1, y=PC2, color=Shock)) + geom_point(alpha=0.5, size=5) +
  labs(x=paste("PCoA 1 (", format(100 * eigShock[1] / sum(eigShock), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eigShock[2] / sum(eigShock), digits=4), "%)", sep=""),
       title="PCoA") + theme_classic()+
  scale_color_manual(values = c('red','blue'))+
  theme_bw()+
  theme(axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  theme(axis.text.x=element_text(size=15),axis.title.x = element_text(size = 15))+
  guides(color=F)
p + geom_enterotype(aes(color=Shock,label=Shock), show.legend = T)+
  scale_color_manual(values = c('red','blue'))



