library(readr)
library(readxl)
library("RColorBrewer")
library("gplots")
library(drc)
library(BSDA)
library(vioplot)
library(lmodel2)
library(data.table)
library('corrplot')
library(plyr)
library(lmodel2)
library("cluster")
library("factoextra")
library("magrittr")
library("NbClust")



Table_S2_eLife <- read_csv("Supplementary File 1b.csv")
Table_S3 <- read_csv("Supplementary File 1c.csv")
Table_S4 <- read_csv("Supplementary File 1d.csv")

###
#fig 2:

dil <- c(0.05,0.1,0.13,0.18,0.3,0.35)
GR_mRNA <- Table_S2_eLife[,c(1:3,46:48,52:66)]

# scale to y-intercept 
GR_mRNA$y_int <- NA
for(i in 1:nrow(GR_mRNA)){
  GR_mRNA[i,"y_int"] <- 10^(summary(lm(unlist(log10(GR_mRNA[i,4:21])) ~ rep(dil, each=3)))$coefficients[1,1])
}
GR_mRNA[,4:21] <- log2(GR_mRNA[,4:21] / GR_mRNA[,"y_int"])
colnames(GR_mRNA)[4:21] <- paste("log2.scaled", colnames(GR_mRNA)[4:21], sep=".")

# NbClust 
Nb_index <- as.data.frame(c("kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew", "friedman", "rubin", 
                            "cindex", "db", "silhouette", "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", 
                            "gap", "frey", "mcclain", "dunn", "sdindex", "sdbw"))
colnames(Nb_index) <- "index"
Nb_index$opt.clusters <- NA
for(i in 1:nrow(Nb_index)){
  j <- as.character(Nb_index[i,1])
  set.seed(101)
  try(Nb_index[i,2] <- NbClust(GR_mRNA[,4:21], distance="euclidean", min.nc=2, max.nc=10, method="kmeans", index=j)$Best.nc[1], silent=T)
}


#Fig 2A  272x330

df_cluster <- as.data.frame(table(Nb_index$opt.clusters))
barplot(df_cluster$Freq, xlab="", ylab="", main="", cex.names=0.8, cex.axis=0.9,
        names.arg=paste("k = ",df_cluster$Var1, sep=""),
        ylim=c(0,15), las=2)
mtext(side=3,"GR experiments", adj=0, line=0)
mtext(side=1, line=4, "Optimal number of \nclusters (k)")
mtext(side=2, line=2,"Frequency among all \nclustering indices")

# visualize k-means clustering with k=3
set.seed(101)
km.res.k3 <- kmeans(GR_mRNA[,4:21], 3, nstart = 50)

GR_mRNA$GR.cluster <- km.res.k3$cluster

## find & plot ksp specific responders
# 3 clusters:
k1 <- GR_mRNA[GR_mRNA$GR.cluster=="1",] #1045
k2 <- GR_mRNA[GR_mRNA$GR.cluster=="2",] #232
k3 <- GR_mRNA[GR_mRNA$GR.cluster=="3",] #1850


#Fig 2B: 272x330
boxplot(matrix(c(rowMeans(k1[,4:6]), rowMeans(k1[,7:9]), rowMeans(k1[,10:12]),
                 rowMeans(k1[,13:15]), rowMeans(k1[,16:18]), rowMeans(k1[,19:21]),rep(1,nrow(k1))),
               ncol=7,byrow=F), #last column of 1's necessary to control boxplot width below
        at=c(dil*17.5,100), width=c(rep(0.5,6),1), outline=F, border="green4", col=rgb(0,0,0,0),
        xlim=c(0,7), ylim=c(0,4.5), xaxt="n", yaxt="n", xlab="", ylab="")
par(new=T)
boxplot(matrix(c(rowMeans(k3[,4:6]), rowMeans(k3[,7:9]), rowMeans(k3[,10:12]),
                 rowMeans(k3[,13:15]), rowMeans(k3[,16:18]), rowMeans(k3[,19:21]),rep(1,nrow(k3))),
               ncol=7,byrow=F), #last column of 1's necessary to control boxplot width below
        at=c(dil*17.5,100)-0.05, width=c(rep(0.5,6),1), outline=F,border="mediumorchid",  col=rgb(0,0,0,0),
        xlim=c(0,7), ylim=c(0,4.5), 
        #xlim=c(0,7), ylim=c(0,3), xaxt="n", yaxt="n", xlab="growth rate", ylab="scaled mRNA expression"
        axes=F)
axis(side=1, at=c(0,.1,.2,.3,.4)*17.5, labels=c(0,.1,.2,.3,.4))
axis(side=2,at=(0:4),labels=2^(0:4), las=1)
mtext(side=2, "scaled mRNA expression", line=2.5)
mtext(side=1, "growth rate", line=2.5)
mtext(side=3,"2,895 out of 3,127 genes", col="black", adj=0, line=0, cex=0.85)
mtext(side=3," GR.cluster.1", col="green4", adj=0, line=-1, cex=0.85)
mtext(side=3," GR.cluster.3", col="mediumorchid", adj=0, line=-2, cex=0.85)

#Fig 2C: 272x330
boxplot(matrix(c(rowMeans(k2[,4:6]), rowMeans(k2[,7:9]), rowMeans(k2[,10:12]),
                 rowMeans(k2[,13:15]), rowMeans(k2[,16:18]), rowMeans(k2[,19:21]),rep(1,nrow(k2))),
               ncol=7,byrow=F), #last column of 1's necessary to control boxplot width below
        at=c(dil*17.5,100), width=c(rep(0.5,6),1), outline=F, #border="snow4",
        xlim=c(0,7), ylim=c(-2.5,2), xaxt="n", yaxt="n", xlab="", ylab="")
axis(side=1, at=c(0,.1,.2,.3,.4)*17.5, labels=c(0,.1,.2,.3,.4))
axis(side=2,at=(-2:2),labels=2^(-2:2), las=1)
mtext(side=2, "scaled mRNA expression", line=3)
mtext(side=1, "growth rate", line=2.5)
mtext(side=3,"232 out of 3,127 genes", col="black", adj=0, line=0, cex=0.85)
mtext(side=3," GR.cluster.2", col="black", adj=0, line=-1, cex=0.85)



GR_mRNA_k2 <- k2
#Fig 2D: write GR_mRNA_k2 into file, plug into SGD's GOslimMapping fuction





#Fig 2E: 
RNApol2 <- read_csv("RNApol2BTC.csv") #67 total
RNApol2 <- RNApol2[RNApol2$annotation=="RNA pol II subunit",] #12 polII subunits
RNApol2 <- merge(RNApol2[,c(1:3)], Table_S2_eLife[,c(1,4:6,10:24)]) #10 detected
RNApol2[RNApol2==0] <- NA
xx <- as.numeric()
for(i in 4:21){
  xx[i-3] <- median(RNApol2[,i], na.rm=T)
}
y_int <- 10^(summary(lm(log10(xx) ~ rep(dil, each=3)))$coefficients[1,1])
xx <- xx/y_int


x <- as.data.frame(matrix(NA, ncol = 4, nrow = 18)) #18 samples
colnames(x) <- c("PolII", "k1", "k2", "k3")
x$PolII <- log2(xx) 
x$k1 <- apply(k1[,4:21],2,median)
x$k2 <- apply(k2[,4:21],2,median)
x$k3 <- apply(k3[,4:21],2,median)

plot(x$PolII, x$k1, ylim=c(-1,3.5), xlim=c(-.2,2.7), xlab="", ylab="", 
     xaxt="n", yaxt="n", col="green4")
par(new=T)
plot(x$PolII, x$k2, ylim=c(-1,3.5), xlim=c(-.2,2.7), xlab="", ylab="", axes=F)
par(new=T)
plot(x$PolII, x$k3, ylim=c(-1,3.5), xlim=c(-.2,2.7), xlab="", ylab="", axes=F, col="mediumorchid")
abline(a=summary(lm(x$k1 ~ x$PolII))$coefficients[1,1],b=summary(lm(x$k1 ~ x$PolII))$coefficients[2,1], col="green4")
abline(a=summary(lm(x$k2 ~ x$PolII))$coefficients[1,1],b=summary(lm(x$k2 ~ x$PolII))$coefficients[2,1])
abline(a=summary(lm(x$k3 ~ x$PolII))$coefficients[1,1],b=summary(lm(x$k3 ~ x$PolII))$coefficients[2,1], col="mediumorchid")
abline(a=0,b=1, lwd=2, col="snow4", lty=5)
axis(side=1, at=c(0:2), labels=2^c(0:2))
axis(side=2,at=(-1:3),labels=2^(-1:3), las=1)
mtext(side=2, "scaled mRNA expression", line=2.5)
mtext(side=1,"scaled protein expression \nof RNA pol II", line=3.5)
mtext(side=3, "GR experiments", line=0, adj=0)
text(2.5,2.9,"y=x", col="snow4", cex=0.9)



remove(RNApol2, x, xx)
remove(k1,k2,k3,km.res.k3,Nb_index)



###


#Fig 3A GR histogram: 272x330
hist(Table_S2_eLife$GR.Spearman.rho, xlim=c(-1,1),ylim=c(0,1000), main="", xlab="", las=1, col="snow2", xaxt="n")
median(Table_S2_eLife$GR.Spearman.rho)
mtext(side=3,line=-1, adj=0,"  GR experiments \n  median = 0.76")
mtext(side=1, line=2.5, expression(paste("Spearman ", rho)))
axis(side=1, at=c(-1,0,1))
axis(side=1, at=c(-0.5,0.5))


#Fig 3B
#first calculate Ksp:

#mine kdP
library(readxl)
kdP_Lahtvee <- read_excel("Protein turnover rate and synthesis rate (Lahtvee Table S5 S6).xlsx", 
                          col_types = c("skip", "skip", "skip", "skip", "text", "numeric", "skip", 
                                        "skip", "skip", "skip"))
colnames(kdP_Lahtvee) <- c("Accession", "kdP (1/h)")

kdP_plus_mu <- merge(Table_S2_eLife[,1:3], kdP_Lahtvee, all.x=TRUE)
kdP_plus_mu[is.na(kdP_plus_mu)] <- median(kdP_Lahtvee$`kdP (1/h)`, na.rm=TRUE)
mu_42 <- c(rep(0.05,3), rep(0.1,6), rep(0.13,3), rep(0.18,3),
           rep(0.3,3), rep(0.35,3), rep(0.1,21))

for(i in 1:42){
  kdP_plus_mu[,(i+4)]<- kdP_plus_mu[,4]+mu_42[i]
}
kdP_plus_mu$`kdP (1/h)`<- NULL
colnames(kdP_plus_mu) <- c("Accession", "Gene", "protein length", paste("kdP+u", 1:42, sep="."))

#calculate ksP
ksP <- Table_S2_eLife[,1:3]
for(i in 4:45){
  for(j in 1:nrow(ksP)){
    ksP[j,i] <- Table_S2_eLife[j,i] / Table_S2_eLife[j,(i+42)] * kdP_plus_mu[j,i]
    #protein (fmol/mgDW) / RNA (fmol/mgDW) * (kdP+u) 
  }
}
colnames(ksP) <- c("Accession", "Gene", "protein length", paste("ksP", 1:42, sep="."))



#Fig 3B:GR ksp 272 x 330 
ksP_lm <- ksP[,1:21]
ksP_lm$y_int <- NA
for(i in 1:nrow(ksP_lm)){
  test <- as.data.frame(unlist(log10(ksP_lm[i,4:21])))
  colnames(test) <- "y"
  test$x <- rep(dil,each=3)
  is.na(test) <- sapply(test, is.infinite)
  ksP_lm[i,"y_int"] <- 10^(summary(lm(test$y ~ test$x))$coefficients[1,1])
}
ksP_lm[,4:21] <- log2(ksP_lm[,4:21] / ksP_lm[,"y_int"])
colnames(ksP_lm)[4:21] <- paste("log2.scaled", colnames(ksP_lm)[4:21], sep=".")

boxplot(matrix(c(rowMeans(ksP_lm[,c(4,6)]), rowMeans(ksP_lm[,7:9]), 
                 rowMeans(ksP_lm[,10:12]), rowMeans(ksP_lm[,13:15]), 
                 rowMeans(ksP_lm[,16:18]), rowMeans(ksP_lm[,19:21]),
                 rep(1,nrow(ksP_lm))),
               ncol=7,byrow=F), #last column of 1's necessary to control boxplot width below
        at=c(dil*17.5,100), width=c(rep(0.5,6),1), outline=F, las=1,
        xlim=c(0,7), ylim=c(-1,3.5), xaxt="n", yaxt="n", xlab="", 
        ylab=expression(paste("scaled ", italic(k[sP])))) 
axis(side=1, at=c(0,.1,.2,.3,.4)*17.5, labels=c(0,.1,.2,.3,.4))
axis(side=2,at=(-1:3),labels=2^(-1:3), las=1)
mtext(side=3,line=0, adj=0,"  GR experiments")
mtext(side=1, line=2.5, "growth rate")

remove(ksP_lm)


#Fig 3C: GR ribosomal fraction 272x330

ribo <- Table_S2_eLife[grep("^RPL|^RPS|^RPP", Table_S2_eLife$Gene), 1:45]

#GR set:
ribo_fraction <- colSums(ribo[,c(4:6,10:24)]) / colSums(Table_S2_eLife[,c(4:6,10:24)])
plot(rep(dil, each=3), ribo_fraction*100, xlim=c(0,0.4), ylim=c(0,50),
     xlab="", ylab="", las=1)
mtext(side=3, line=0, adj=0, "  GR experiments")
mtext("% (mol/mol) ribosomal \nproteins in proteome",side=2, line=2.2)
mtext(side=1, line=2.5, "growth rate")



#######
#Fig 4 AA dataset (NM dataset in eLife resubmission)

AA_mRNA <- Table_S2_eLife[,c(1:3,49:54,67:87)]

# scale by mean of Phe and Ile, N30 
AA_mRNA[,4:30] <- log2(AA_mRNA[,4:30] / rowMeans(AA_mRNA[,c(22:24,28:30)])) 

Nb_index <- as.data.frame(c("kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew", "friedman", "rubin", 
                            "cindex", "db", "silhouette", "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", 
                            "gap", "frey", "mcclain", "dunn", "sdindex", "sdbw"))
colnames(Nb_index) <- "index"
Nb_index$opt.clusters <- NA
for(i in 1:nrow(Nb_index)){
  j <- as.character(Nb_index[i,1])
  set.seed(101)
  try(Nb_index[i,2] <- NbClust(AA_mRNA[,4:30], distance="euclidean", min.nc=2, max.nc=10, method="kmeans", index=j)$Best.nc[1], silent=T)
}

#Fig 4A
df_cluster <- as.data.frame(table(Nb_index$opt.clusters))
barplot(df_cluster$Freq, xlab="", ylab="", main="", cex.names=0.8, cex.axis=0.9,
        names.arg=paste("k = ",df_cluster$Var1, sep=""),
        ylim=c(0,10), las=2)
mtext(side=3," NM experiments", adj=0, line=0)
mtext(side=1, line=4, "Optimal number of \nclusters (k)")
mtext(side=2, line=2,"Frequency among all \nclustering indices")



# visualize k-means clustering with k=2
set.seed(101)
km.res <- kmeans(AA_mRNA[,4:30], 2, nstart = 50)

AA_mRNA$AA.cluster <- km.res$cluster

k1 <- AA_mRNA[AA_mRNA$AA.cluster=="1",] #3011
k2 <- AA_mRNA[AA_mRNA$AA.cluster=="2",] #116


#Fig 4B: 272x330
boxplot(matrix(c(rowMeans(k1[,22:24]), rowMeans(k1[,28:30]), #F(N30), I(N30)
                 rowMeans(k1[,7:9]), rowMeans(k1[,16:18]),  #NH4(N30), Q(N30)
                 rowMeans(k1[,19:21]), rowMeans(k1[,25:27]), #F(std),I(std)
                 rowMeans(k1[,4:6]), rowMeans(k1[,13:15]), #NH4(std),Q-(std)
                 rowMeans(k1[,10:12]),  #Q(+glc), 
                 rep(NA,nrow(k1))),  
               ncol=10,byrow=F), 
        width=c(rep(0.65,9),1),outline=F, xaxt="n", yaxt="n", ylab="",
        at=1:10-0.15, border="maroon1", col=rgb(0,0,0,0),
        xlim=c(0.5,9.5), ylim=c(-1.2,2.2)) 
axis(side=1, at=1:4, labels=c("Phe", "Ile", "NH4", "Gln"), col.axis="tomato", col.ticks="tomato", las=2)
axis(side=1, at=5:9, labels=c("Phe", "Ile", "NH4", "Gln*", "Gln"), col.axis="royalblue", col.ticks="royalblue", las=2)
axis(side=2,at=(-1:2),labels=2^(-1:2), las=1)
mtext(side=2,"scaled mRNA expression", line=2.5)
mtext(side=3,"3,011 out of 3,127 genes", col="black", adj=0, line=0, cex=0.85)
mtext(side=3," NM.cluster.1", col="maroon1", adj=0, line=-1, cex=0.85)

#Fig 4C: 272x330
boxplot(matrix(c(rowMeans(k2[,22:24]), rowMeans(k2[,28:30]), #F(N30), I(N30)
                 rowMeans(k2[,7:9]), rowMeans(k2[,16:18]),  #NH4(N30), Q(N30)
                 rowMeans(k2[,19:21]), rowMeans(k2[,25:27]), #F(std),I(std)
                 rowMeans(k2[,4:6]), rowMeans(k2[,13:15]), #NH4(std),Q-(std)
                 rowMeans(k2[,10:12]),  #Q(+glc), 
                 rep(NA,nrow(k2))),  
               ncol=10,byrow=F), 
        width=c(rep(0.65,9),1),outline=F, xaxt="n", yaxt="n", ylab="scaled mRNA expression",
        at=1:10+0.15, col=rgb(0,0,0,0),
        xlim=c(0.5,9.5), ylim=c(-2.5,1.5)) 
axis(side=1, at=1:4, labels=c("Phe", "Ile", "NH4", "Gln"), col.axis="tomato", col.ticks="tomato", las=2)
axis(side=1, at=5:9, labels=c("Phe", "Ile", "NH4", "Gln*", "Gln"), col.axis="royalblue", col.ticks="royalblue", las=2)
axis(side=2,at=(-2:1),labels=2^(-2:1), las=1)
mtext(side=3,"116 out of 3,127 genes", col="black", adj=0, line=0, cex=0.85)
mtext(side=3," NM.cluster.2", col="black", adj=0, line=-1, cex=0.85)


#Fig 4D: copy the 116 genes in cluster2 into file, plug into SGD's GOslimMapping fuction


remove(clusters, AA_mRNA_k2_GOslim, AA_mRNA_k3)



#Fig 4E

RNApol2 <- read_csv("RNApol2BTC.csv") #67 total
RNApol2 <- RNApol2[RNApol2$annotation=="RNA pol II subunit",] #12 polII subunits
RNApol2 <- merge(RNApol2[,c(1:3)], Table_S2_eLife[,c(1,7:12,25:45)]) #10 detected
xx <- as.numeric()
j <- c(22:24, 28:30, #F(N30), I(N30)
       7:9,16:18, #NH4(N30), Q(N30)
       19:21,25:27, #F(std),I(std)
       4:6,13:15, #NH4(std),Q-(std)
       10:12) #Q(+glc)
for(k in 1:27){
  i <- j[k]
  xx[k] <- median(RNApol2[,i], na.rm=T)
}

xx <- xx/mean(xx[1:6])

x <- as.data.frame(matrix(NA, ncol = 3, nrow = 27)) #27 samples
colnames(x) <- c("PolII", "k1", "k2")
x$PolII <- log2(xx)
x$k1 <- apply(k1[,4:30],2,median)
x$k2 <- apply(k2[,4:30],2,median)

plot(x$PolII, x$k1, ylim=c(-1.5,2), xlab="", ylab="", xaxt="n", yaxt="n",
     xlim=c(-1.2,1),  col="maroon1")
par(new=T)
plot(x$PolII, x$k2, ylim=c(-1.5,2), xlim=c(-1.2,1), xlab="", ylab="", axes=F, col="black")
abline(a=summary(lm(x$k1 ~ x$PolII))$coefficients[1,1],b=summary(lm(x$k1 ~ x$PolII))$coefficients[2,1], col="maroon1", lty=2)
abline(a=summary(lm(x$k2 ~ x$PolII))$coefficients[1,1],b=summary(lm(x$k2 ~ x$PolII))$coefficients[2,1], col="black", lty=2)
axis(side=1, at=(-1:1), labels=2^(-1:1))
axis(side=2,at=(-1:2),labels=2^(-1:2), las=1)
mtext(side=2, "scaled mRNA expression", line=2.5)
mtext(side=1,"scaled protein expression \nof RNA pol II", line=3.5)
mtext(side=3, "NM experiments", line=0, adj=0)
abline(a=0,b=1, lwd=2, col="snow4", lty=5)
text(0.8,1.2,"y=x", col="snow4", cex=0.9)







#Fig 5

AA_DQGS <- as.data.frame(cbind(unlist(colnames(Table_S4[2:28])), unlist(Table_S4[4,2:28]), 
              unlist(Table_S4[5,2:28]), unlist(Table_S4[7,2:28]), unlist(Table_S4[14,2:28])))
colnames(AA_DQGS) <- c("Sample",  "Asp (D)", "Gln (Q)", "Gly (G)", "Ser (S)")
AA_DQGS$`Asp (D)` <- as.numeric(AA_DQGS$`Asp (D)`)
AA_DQGS$`Gln (Q)` <- as.numeric(AA_DQGS$`Gln (Q)`)
AA_DQGS$`Gly (G)` <- as.numeric(AA_DQGS$`Gly (G)`)
AA_DQGS$`Ser (S)` <- as.numeric(AA_DQGS$`Ser (S)`)

#scale to Ile and Phe N30 as before
AA_DQGS[,2] <- log2(AA_DQGS[,2]/mean(unlist(AA_DQGS[c(19:21,25:27),2]), na.rm=T))
AA_DQGS[,3] <- log2(AA_DQGS[,3]/mean(unlist(AA_DQGS[c(19:21,25:27),3]), na.rm=T))
AA_DQGS[,4] <- log2(AA_DQGS[,4]/mean(unlist(AA_DQGS[c(19:21,25:27),4]), na.rm=T))
AA_DQGS[,5] <- log2(AA_DQGS[,5]/mean(unlist(AA_DQGS[c(19:21,25:27),5]), na.rm=T))


#Fig 5C 272x330
plot(AA_DQGS$`Gly (G)`, x$k1, pch="G", ylim=c(-1,2.3), xlim=c(-1.5,3.3), cex=0.7,
     ylab="", xlab="", xaxt="n", yaxt="n", col="maroon1")
par(new=T)
plot(AA_DQGS$`Ser (S)`, x$k1, pch="S", ylim=c(-1,2.3), xlim=c(-1.5,3.3), cex=0.7,axes=F, col="maroon1", xlab="", ylab="")
abline(a=0,b=1, lwd=2, col="snow4", lty=5)
abline(a=lmodel2(c(x$k1, x$k1) ~ c(AA_DQGS$`Gly (G)`, AA_DQGS$`Ser (S)`))$regression.results[3,2],
       b=lmodel2(c(x$k1, x$k1) ~ c(AA_DQGS$`Gly (G)`, AA_DQGS$`Ser (S)`))$regression.results[3,3], col="maroon1")
axis(side=1, at=(-1:3), labels=2^(-1:3))
axis(side=2,at=(-1:3),labels=2^(-1:3), las=1)
text(1.5,2.2,"y=x", col="snow4")
mtext(side=1, "scaled intracellular \nGly (G) and Ser (S)", line=3)
mtext(side=2, "scaled mRNA expression", line=2.5)



#Fig 5D 272x330

plot(AA_DQGS$`Asp (D)`, x$k1, pch="D", ylim=c(-1,2), xlim=c(-1.5,3.4), cex=0.7, yaxt="n", xaxt="n", 
     ylab="", xlab="", col="maroon1")
par(new=T)
plot(unlist(AA_DQGS[c(1:6,16:27), "Gln (Q)"]), unlist(x[c(1:6,16:27),"k1"]), pch="Q",
     ylim=c(-1,2), xlim=c(-1.5,3.4), cex=0.7, axes=F, xlab="", ylab="", col="maroon1")
abline(a=0,b=1, lwd=2, lty=5, col="snow4")
abline(a=lmodel2(x$k1 ~ AA_DQGS$`Asp (D)`)$regression.results[3,2],
       b=lmodel2(x$k1 ~ AA_DQGS$`Asp (D)`)$regression.results[3,3], col="maroon1")
abline(a=lmodel2(unlist(x[c(1:6,16:27),"k1"]) ~ unlist(AA_DQGS[c(1:6,16:27), "Gln (Q)"]) )$regression.results[3,2],
       b=lmodel2(unlist(x[c(1:6,16:27),"k1"]) ~ unlist(AA_DQGS[c(1:6,16:27), "Gln (Q)"]) )$regression.results[3,3], col="maroon1")
axis(side=2, at=(-1:2), labels=2^(-1:2), las=1)
axis(side=1,at=(-1:4),labels=2^(-1:4))
text(2.6,1.9,"y=x", col="snow4")
mtext(side=1, "scaled intracellular \nGln (Q) and Asp (D)", line=3)
mtext(side=2, "scaled mRNA expression", line=2.5)


remove(RNApol2, x, AA_mRNA_k2, GR_mRNA_k2, AA_DQGS)
remove(km.res,k1,k2,k3, Nb_index)
remove(AA_DQGS, DQGS, DQGS_mean)



###

#Fig 6A AA histogram: 272x330
#NB - in eLife Supplemental File replace "NM" with "AA"

hist(Table_S2_eLife$AA.Spearman.rho, xlim=c(-1,1), ylim=c(0,500), main="", xlab="", las=1, col="snow2", xaxt="n")
median(Table_S2_eLife$AA.Spearman.rho)
mtext(side=3,line=-1, adj=0,"  NM experiments \n  median = 0.07")
mtext(side=1, line=2.5, expression(paste("Spearman ", rho)))
axis(side=1, at=c(-1,0,1))
axis(side=1, at=c(-0.5,0.5))
nrow(Table_S2_eLife[Table_S2_eLife$AA.Spearman.rho > 0.5,])
nrow(Table_S2_eLife[Table_S2_eLife$AA.Spearman.rho < -0.5,])

#Fig 6B AA ksp
#scale by Phe and Ile, N30 
ksP_D01 <- ksP[,c(1:3,7:12,25:45)]
ksP_D01[,4:30] <- log2(ksP_D01[,4:30] / rowMeans(ksP_D01[,c(22:24,28:30)]))
is.na(ksP_D01) <- sapply(ksP_D01, is.infinite)

boxplot(matrix(c(rowMeans(ksP_D01[,22:24]), rowMeans(ksP_D01[,28:30]), #F(N30), I(N30)
                 rowMeans(ksP_D01[,16:18]), rowMeans(ksP_D01[,7:9]),   #Q(N30), NH4(N30), 
                 rowMeans(ksP_D01[,19:21]), rowMeans(ksP_D01[,25:27]), #F(std),I(std)
                 rowMeans(ksP_D01[,4:6]),#NH4(std),
                 rowMeans(ksP_D01[,13:15]), #Q-(std)
                 rowMeans(ksP_D01[,11:12]),  #Q(+glc), 
                 rep(NA,nrow(ksP_D01))), 
               ncol=10,byrow=F), #last column of NA's necessary to control boxplot width below
        width=c(rep(0.7,9),1),outline=F, xaxt="n", yaxt="n", ylab=expression(paste("scaled ", italic(k[sP]))),
        col="white",
        xlim=c(0.5,9.5), ylim=c(-3,2)) 
axis(side=1, at=1:4, labels=c("Phe", "Ile", "NH4", "Gln"), col.axis="tomato", col.ticks="tomato", las=2)
axis(side=1, at=5:9, labels=c("Phe", "Ile", "NH4", "Gln*", "Gln"), col.axis="royalblue", col.ticks="royalblue", las=2)
axis(side=2,at=(-3:2),labels=2^(-3:2), las=1)
mtext(side=3,line=0, adj=0,"  NM experiments")



remove(ksP_D01)


#Fig 6C: AA ribosomal fraction 272x330 

ribo_fraction <- colSums(ribo[,c(7:12,25:45)]) / colSums(Table_S2_eLife[,c(7:12,25:45)])
plot(rep(c(8,4,5,9,7,2,1,6,3),each=3), ribo_fraction*100, xlim=c(.5,9.5), ylim=c(0,50),
     xaxt="n", xlab="", ylab="", las=1)
axis(side=1, at=1:4, labels=c("Phe", "Ile", "NH4", "Gln"), col.axis="tomato", col.ticks="tomato", las=2)
axis(side=1, at=5:9, labels=c("Phe", "Ile", "NH4", "Gln*", "Gln"), col.axis="royalblue", col.ticks="royalblue", las=2)
mtext(side=3, line=0, adj=0, "  NM experiments")
mtext("% (mol/mol) ribosomal \nproteins in proteome",side=2, line=2.2)

remove(ribo_fraction)

#Fig 7A 
hist(Table_S2_eLife$all.Spearman.rho, xlim=c(-1,1), ylim=c(0,600), main="", xlab="", las=1, col="snow2", xaxt="n")
median(Table_S2_eLife$all.Spearman.rho)
mtext(side=3,line=-1, adj=0,"  all data (abs-quant) \n  median = 0.22")
mtext(side=1, line=2.5, expression(paste("Spearman ", rho)))
axis(side=1, at=c(-1,0,1))
axis(side=1, at=c(-0.5,0.5))

#Fig 7B:

test <- data.frame(Table_S2_eLife[,c("Gene", "all.Spearman.rho")]) #pull Gene name and all-set correlation
test <- na.omit(test) 
test <- test[order(test[,2]),] #order by correlation, low to high

GO_slim_P_summary <- read_csv("GO_slim_P_summary.csv")
GO_slim_P_map <- read_csv("GO_slim_P_map.csv")

GO_all_Spearman <- GO_slim_P_summary[,1:3] #output file
GO_window_median <- as.numeric() #median Spearman rho of each 200-gene window
for(i in 1:(nrow(test)-199)){ #
  sw <- data.frame(test[i:(199+i), 1]) #grab window = 200 genes
  colnames(sw) <- "Gene"
  sw <- merge(sw, GO_slim_P_map[,c(2,4)], by.x="Gene", by.y="GENE") #merge with GO slim mapper
  tab <- data.frame(table(sw$`GO term`)) #genes in GO in sw
  tab[,3] <- (200 - tab[,2]) #not in GO in sw
  colnames(tab) <- c("Var1", "Var2", "Var3")
  tab <- merge(tab, GO_slim_P_summary[,c(1,4,5)], by.x="Var1", by.y="GO-Slim term") #pull in GO in genome & not in GO in genome
  for(j in 1:nrow(tab)){
    tab[j,6] <- fisher.test(rbind(c(tab[j,4], tab[j,2]),
                                  c(tab[j,5], tab[j,3])),
                            alternative="less")$p.value
  }
  GO_all_Spearman <- merge(GO_all_Spearman, tab[,c(1,6)], by.x="GO-Slim term", by.y="Var1", all.x=TRUE)
  colnames(GO_all_Spearman)[i+3] <- i
  
  GO_window_median <- c(GO_window_median, median(test[i:(199+i), 2]))
}
#code works but will take a while to run 

#draw a heatmap - code works but will take a while to run 
hm <- GO_all_Spearman
hm[is.na(hm)] <- 1
hm <- hm[order(hm$order),]

#for every row, if there is no enrichment across the board, remove it
for (i in 1:nrow(hm)){
  x <- unlist(hm[i,4:ncol(hm)])
  if(min(x)>0.5){
    hm[i,] <- NA
  }
}
hm <- na.omit(hm)

#only keep metab and translation terms
hm <- hm[hm$order=="a" | hm$order=="b",]
rownames(hm) <- 1:nrow(hm)

pdf("GO_all_Spearman_metab_transl.pdf", width=6, height=5)
heatmap.2(data.matrix(hm[,4:ncol(hm)]), 
          dendrogram="none", trace="none", scale="none",
          labRow = hm$`GO-Slim term`, cexRow=0.7, 
          Colv=FALSE, Rowv=FALSE, #RowSideColors=hm$color,
          col=colorRampPalette(c("steelblue4","steelblue1",rep("ghostwhite",3)))(n=6),
          density.info="none", margins=c(5,20), key=F
)
dev.off()


# draw selected metab and translation terms:
# organize the rows manually
hm <- hm[c(2, 5,7, 1, 3, 19, 33:38, 40,41),]

pdf("Fig_7B.pdf", width=6, height=3)
heatmap.2(data.matrix(hm[,4:ncol(hm)]), 
          dendrogram="none", trace="none", scale="none",
          labRow = hm$`GO-Slim term`, cexRow=0.7, 
          labCol = round(GO_window_label,1), cexCol=0.7, 
          Colv=FALSE, Rowv=FALSE, #RowSideColors=hm$color,
          col=colorRampPalette(c("steelblue4","steelblue1",rep("ghostwhite",3)))(n=6),
          density.info="none", margins=c(5,20), key=F
)
dev.off()

remove(test, sw, tab, x, hm)
remove(i,j,k)


# Fig 7C

#read protein fasta:
library("Biostrings")
z <- readAAStringSet("orf_trans.fasta")

x = names(z)
x = sub(" .*", "", x) #remove everything after the first space
y = paste(z)
y = sub("[*]", "", y) #remove the trailing *
prot_fasta <- data.frame(x, y)
colnames(prot_fasta)[2] <- "sequence"

#map Accession & Gene 
z <- as.data.frame(prot_fasta[,1]) 
x <- data.frame(matrix(NA, ncol = 3, nrow = nrow(z)))
for (i in 1:nrow(z)){ 
  if(length(which(rowSums(S288c == as.character(z[i,1]),  na.rm = TRUE) >0)) > 0){
    x[i,1] <- S288c[which(rowSums(S288c == as.character(z[i,1]),  na.rm = TRUE) >0),"Entry"]
    x[i,2] <- S288c[which(rowSums(S288c == as.character(z[i,1]),  na.rm = TRUE) >0),"Gene"]
    x[i,3] <- as.character(z[i,1])
  }
  else{
    x[i,1] <- NA
    x[i,2] <- NA
    x[i,3] <- as.character(z[i,1])
  }
}
colnames(x) <- c("Accession", "Gene", "x")

prot_fasta <- merge(x, prot_fasta)
prot_fasta <- na.omit(prot_fasta[,c("Accession", "Gene", "sequence")])

#tabulate AA content:
AA_count <- data.frame(matrix(NA, ncol=22))
colnames(AA_count) <- c("Accession","Gene", 
                        "A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")

for(i in 1:nrow(prot_fasta)){
  x <- as.character(prot_fasta[i,"sequence"])
  y <- table(unlist(strsplit(x, ""), use.names=FALSE))
  AA_count[i,1] <- prot_fasta[i,"Accession"]
  AA_count[i,2] <- prot_fasta[i,"Gene"]
  for(j in 3:22){
    AA_count[i,j] <- y[colnames(AA_count)[j]]
  }
}

#tabulate AA distribution in the proteome:
AA_in_prot_umol_gDW <- data.frame(matrix(NA, ncol=43, nrow=20))
colnames(AA_in_prot_umol_gDW) <- c("AA", paste("umolAA/gDW.in.prot", 1:42, sep="."))
AA_in_prot_umol_gDW$AA <- colnames(AA_count)[3:22]

for(i in 1:42){
  x <- merge(Table_S2_eLife[,c(1,3+i)],AA_count)
  x[,4:23]<- x[,4:23]*x[,2]
  for(j in 1:nrow(AA_in_prot_umol_gDW)){
    AA_in_prot_umol_gDW[j,i+1] <- sum(x[,which(colnames(x)==AA_in_prot_umol_gDW[j,1])], na.rm=T)/1e6
  }
}

AA_percent <- AA_count
AA_percent[,3:22] <- AA_percent[,3:22] / rowSums(AA_percent[,3:22], na.rm=T)*100

#16 brackets
test <- Table_S2_eLife[,c(1,2,95)]
test <- merge(test, AA_percent[,c(1,3:22)])
test <- test[order(test$all.Spearman.rho),]

i <- 1
test2 <- test[i*200-(199:0),]
test3 <- as.data.frame(colMeans(test2[,4:23], na.rm=T))
colnames(test3)[i] <- i
for(i in 2:15){
  test2 <- test[i*200-(199:0),]
  test3 <- cbind(test3, as.data.frame(colMeans(test2[,4:23], na.rm=T)))
  colnames(test3)[i] <- i
}
test2 <- test[3001:3125,]
test3 <- cbind(test3, as.data.frame(colMeans(test2[,4:23], na.rm=T)))
colnames(test3)[16] <- 16

test3$correlation <- NA
test3$cor.p <- NA
for(i in 1:nrow(test3)){
  test3[i,"correlation"] <- cor(unlist(test3[i,1:16]), 1:16)
  test3[i,"cor.p"] <- cor.test(unlist(test3[i,1:16]), 1:16)$p.value
}

test4 <- test3[,1:16]/rowMeans(test3[,1:16]) 

heatmap.2(data.matrix(log2(test4)), 
          trace="none", scale="none",
          Colv=FALSE, Rowv=T, dendrogram="row", 
          col=colorRampPalette(rev(brewer.pal(5,"RdBu"))),
          density.info="none"
)

# grab just the ones that are changing
test5 <- log2(test4[c("R", "K", "G", "V", "A", "D", "E", "N", "Q", "S"),])
test5 <- -test5 #this is to rig the heatmap so that the increasing ones end up on top (controlled by the dentrogram)

#300 x 350
heatmap.2(data.matrix(test5), 
          trace="none", scale="none",
          Colv=FALSE, Rowv=T, dendrogram="row", 
          col=colorRampPalette((brewer.pal(5,"RdBu"))),
          density.info="none" 
)


remove(test, test2, test3, test4)


# Fig 7D:

library(Biostrings)

z <- readDNAStringSet("orf_coding.fasta")

x = names(z)
x = sub(" .*", "", x) #remove everything after the first space
y = paste(z)
#y = sub("[*]", "", y) #no trailing * in the DNA sequence file 
RNA_fasta <- data.frame(x, y)
colnames(RNA_fasta)[2] <- "sequence"

#map Accession & Gene 
z <- as.data.frame(RNA_fasta[,1]) 
x <- data.frame(matrix(NA, ncol = 3, nrow = nrow(z)))
for (i in 1:nrow(z)){ 
  if(length(which(rowSums(S288c == as.character(z[i,1]),  na.rm = TRUE) >0)) > 0){
    x[i,1] <- S288c[which(rowSums(S288c == as.character(z[i,1]),  na.rm = TRUE) >0),"Entry"]
    x[i,2] <- S288c[which(rowSums(S288c == as.character(z[i,1]),  na.rm = TRUE) >0),"Gene"]
    x[i,3] <- as.character(z[i,1])
  }
  else{
    x[i,1] <- NA
    x[i,2] <- NA
    x[i,3] <- as.character(z[i,1])
  }
}
#x <- na.omit(x)
colnames(x) <- c("Accession", "Gene", "x")

RNA_fasta <- merge(x, RNA_fasta)
RNA_fasta <- na.omit(RNA_fasta[,c("Accession", "Gene", "sequence")])

remove(x,y,z)


codon_key <- read_excel("codon_key.xlsx")
codon_count <- data.frame(matrix(NA, ncol=63))
colnames(codon_count) <- c("Accession", codon_key$`cDNA codon`)
codon_count$Accession <- as.character(codon_count$Accession)

for(i in 1:nrow(RNA_fasta)){
  x <- as.character(RNA_fasta[i,3])
  y <- unlist(str_match_all(x, ".{3}"))
  z <- as.data.frame(table(y[1:(length(y)-1)])) #don't count the stop codon  #bc TGA (opal) codes for W in the yeast mitochondria
  z <- transpose(z)
  colnames(z) <- z[1,]
  z <- as.data.frame(z[2,])
  z$Accession <- as.character(RNA_fasta[i,1])
  
  codon_count <- merge(codon_count, z, all=T)
}

codon_count[,2:65] <- sapply(codon_count[,2:65], as.numeric)



# 16 brackets

test <- Table_S2_eLife[,c(1,2,95)]
test <- merge(test, codon_count)
test <- test[order(test$all.Spearman.rho),]

i <- 1
test2 <- test[i*200-(199:0),]
test3 <- as.data.frame(colMeans(test2[,4:67], na.rm=T))
colnames(test3)[i] <- i
for(i in 2:15){
  test2 <- test[i*200-(199:0),]
  test3 <- cbind(test3, as.data.frame(colMeans(test2[,4:67], na.rm=T)))
  colnames(test3)[i] <- i
}
test2 <- test[3001:3125,]
test3 <- cbind(test3, as.data.frame(colMeans(test2[,4:67], na.rm=T)))
colnames(test3)[16] <- 16


test3$`cDNA codon` <- rownames(test3)
test3 <- merge(codon_key, test3)
test3 <- test3[order(test3$AA),2:20]


test <- c("A","C","D","E","F","G","H","I","K","L","N","P","Q","R","S","T","V","Y") # no M or W bc they're single-codon

pdf("codon usage.pdf", height=9, width=10)
par(mfrow=c(4,5)) #(4 rows x 5 columns)
for(i in 1:length(test)){
  x <- test3[test3$AA==test[i],]
  for(j in 4:19){
    x[,j] <- x[,j] / sum(x[,j])*100
  }
  x[,4:19] <- log2(x[,4:19]/x[,4])
  
  for(k in 1:nrow(x)){
    plot(1:16, x[k,4:19], ylim=c(min(unlist(x[,4:19])), max(unlist(x[,4:19])) ), type="l", xlab="", ylab="", las=1,
         col=c("red", "blue", "cyan", "orange", "green3", "purple")[k])
    text(k*2.5, max(unlist(x[,4:19])), x[k,3], cex=0.6, col=c("red", "blue", "cyan", "orange", "green3", "purple")[k])
    par(new=T)
  }
  plot(0,0,type="n", axes=F, xlab="spearman", ylab="log2 FC", main=x[k,2])
  
}

dev.off()

#hand pick those with big & clear changes
test4 <- test3[test3$AA=="S" | test3$AA=="T" | test3$AA=="I" | test3$AA=="G" | test3$AA=="V" ,]
for(j in 4:19){
  test4[,j] <- test4[,j] / sum(test4[,j])*100
}
test4[,4:19] <- log2(test4[,4:19]/rowMeans(test4[,4:19])) #test4[,4]
rownames(test4) <- 1:21
test4 <- test4[c(17,15,14,16, #Thr
                 22, #break
                 13,11, 10,12,9,8, #Ser
                 22, #break
                 7,6,5, #Ile
                 22,
                 4,2,1,3, #gly
                 22,
                 21,19,18,20 #val
),] 

#300 x 550
heatmap.2(data.matrix(test4[,4:19]), 
          trace="none", scale="none",
          labRow = test4$`RNA codon`, 
          Colv=FALSE, Rowv=F, dendrogram="none",
          col=colorRampPalette(rev(brewer.pal(5,"RdBu"))),
          density.info="none" 
)






#####

# Fig 9

#read Marguerat data:
library(readxl)
Marguerat <- read_excel("Marguerat 2012 Table S4.xlsx", 
                        col_types = c("text", "text", "numeric", 
                                      "numeric", "numeric", "numeric", 
                                      "text", "text"))
Marguerat$x <- rowSums(Marguerat[,3:6], na.rm=T)
Marguerat <- Marguerat[Marguerat$x>0,]
Marguerat <- Marguerat[,1:7]
colnames(Marguerat)[3:6] <- c("prolif.mRNA.counts.per.cell", "quies.mRNA.counts.per.cell", 
                              "prolif.prot.counts.per.cell", "quies.prot.counts.per.cell")


#grab Carbon metab genes from pombebase
library(readr)
Marguerat_Cmetab <- read_delim("Marguerat_related_pombe_GOslim_Cmetab.tsv", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)

# don't adjust for C-metab transcripts 
Marguerat_Cmetab$GR.adj <- 1
Marguerat_Cmetab$N.adj <- 1
Marguerat_predict <- merge(Marguerat, Marguerat_Cmetab[,c(1,4:5)], by.x="Systematic.name", by.y="Systematic ID", all.x=T)

#GR  
x <- GR_mRNA[GR_mRNA$GR.cluster==3,1:3]
x <- merge(x, Table_S2_eLife[,c(1,46:48,61:63)])
Marguerat_predict$GR.adj[is.na(Marguerat_predict$GR.adj)] <- mean(unlist(x[,4:6]))/mean(unlist(x[,7:9]))

#NM
x <- AA_mRNA[AA_mRNA$AA.cluster==1,1:3]
x <- merge(x, absq_fmol_mgDW[,c(1,49:51,52:54)])
Marguerat_predict$N.adj[is.na(Marguerat_predict$N.adj)] <- mean(unlist(x[,7:9]))/mean(unlist(x[,4:6]))
Marguerat_predict$quies.mRNA.predict <- Marguerat_predict$prolif.mRNA.counts.per.cell * Marguerat_predict$GR.adj * Marguerat_predict$N.adj
Marguerat_predict$log2.mRNA.counts.over.predict <- log2(Marguerat_predict$quies.mRNA.counts.per.cell / Marguerat_predict$quies.mRNA.predict)


#Fig 4A: 272x330
plot(log10(Marguerat_predict[Marguerat_predict$log2.mRNA.counts.over.predict > 1, "quies.mRNA.predict"]), 
     log10(Marguerat_predict[Marguerat_predict$log2.mRNA.counts.over.predict > 1,"quies.mRNA.counts.per.cell"]), 
     xlim=c(-3,3), ylim=c(-3,3), col="blue", pch=".", cex=2, xlab="", ylab="", las=1, xaxt="n", yaxt="n") 
par(new=T)
plot(log10(Marguerat_predict[Marguerat_predict$log2.mRNA.counts.over.predict < -1, "quies.mRNA.predict"]), 
     log10(Marguerat_predict[Marguerat_predict$log2.mRNA.counts.over.predict < -1, "quies.mRNA.counts.per.cell"]), 
     xlim=c(-3,3), ylim=c(-3,3), col="blue",  pch=".", cex=2,  xlab="", ylab="",axes=F)
par(new=T)
plot(log10(Marguerat_predict[Marguerat_predict$log2.mRNA.counts.over.predict < 1 & 
                               Marguerat_predict$log2.mRNA.counts.over.predict > -1, "quies.mRNA.predict"]), 
     log10(Marguerat_predict[Marguerat_predict$log2.mRNA.counts.over.predict < 1 & 
                               Marguerat_predict$log2.mRNA.counts.over.predict > -1,"quies.mRNA.counts.per.cell"]), 
     xlim=c(-3,3), ylim=c(-3,3), col="grey", pch=".", cex=2,  xlab="", ylab="",axes=F)
axis(side=1,at=c(-2,0,2))
axis(side=2,at=c(-2,0,2))
mtext(side=1, line=2, "log10 calculated")
mtext(side=2, line=2, "log10 measured")
mtext(side=3, line=0, expression(paste("  quiescent ", italic("S. pombe"))), adj=0)
mtext(side=3, line=-1, "  mRNA", adj=0)

#
## ksP

#grab pombe kdP from Christiano 2012
library(readxl)
Christiano_pombe_kdP <- read_excel("Christiano 2014 Table S2 pombe kdP.xlsx")

Christiano_pombe_kdP <- Christiano_pombe_kdP[,2:3]
Christiano_pombe_kdP[,2] <- Christiano_pombe_kdP[,2]*60
colnames(Christiano_pombe_kdP)[2] <- "kdP (/h)"

Christiano_pombe_kdP <- na.omit(Christiano_pombe_kdP)

Marguerat_predict <- merge(Marguerat_predict, Christiano_pombe_kdP, by.x="Systematic.name", by.y="ENSG", all.x=T)
Marguerat_predict$`kdP (/h)`[is.na(Marguerat_predict$`kdP (/h)`)] <- median(Christiano_pombe_kdP$`kdP (/h)`)

#calculate ksP from proliferating cells
Marguerat_predict$prolif.ksP <- NA
for(i in 1:nrow(Marguerat_predict)){
  Marguerat_predict[i,"prolif.ksP"] <- Marguerat_predict[i,"prolif.prot.counts.per.cell"] / Marguerat_predict[i,"prolif.mRNA.counts.per.cell"] * (0.3+Marguerat_predict[i,"kdP (/h)"])
}

#calculate ksP from quiescent cells
Marguerat_predict$quies.ksP <- NA
for(i in 1:nrow(Marguerat_predict)){
  Marguerat_predict[i,"quies.ksP"] <- Marguerat_predict[i,"quies.prot.counts.per.cell"] / Marguerat_predict[i,"quies.mRNA.counts.per.cell"] * (0.3+Marguerat_predict[i,"kdP (/h)"])
}

Marguerat_predict$quies.ksP.predict <- Marguerat_predict$prolif.ksP * (mean(unlist(ksP[,4:6]))/mean(unlist(ksP[,19:21]))) * (mean(unlist(ksP[,10:12]))/mean(unlist(ksP[,7:9])))
Marguerat_predict$log2.ksP.measured.over.predict <- log2(Marguerat_predict$quies.ksP / Marguerat_predict$quies.ksP.predict)

#Fig 4B: 272x330
plot(log10(Marguerat_predict[Marguerat_predict$log2.ksP.measured.over.predict > 1, "quies.ksP.predict"]), 
     log10(Marguerat_predict[Marguerat_predict$log2.ksP.measured.over.predict > 1,"quies.ksP"]), 
     xlim=c(0,5), ylim=c(1,5), col="blue", pch=".", cex=2, xlab="", ylab="", las=1) 
par(new=T)
plot(log10(Marguerat_predict[Marguerat_predict$log2.ksP.measured.over.predict < -1, "quies.ksP.predict"]), 
     log10(Marguerat_predict[Marguerat_predict$log2.ksP.measured.over.predict < -1, "quies.ksP"]), 
     xlim=c(0,5), ylim=c(1,5), col="blue",  pch=".", cex=2,  xlab="", ylab="",axes=F)
par(new=T)
plot(log10(Marguerat_predict[Marguerat_predict$log2.ksP.measured.over.predict < 1 & 
                               Marguerat_predict$log2.ksP.measured.over.predict > -1, "quies.ksP.predict"]), 
     log10(Marguerat_predict[Marguerat_predict$log2.ksP.measured.over.predict < 1 & 
                               Marguerat_predict$log2.ksP.measured.over.predict > -1,"quies.ksP"]), 
     xlim=c(0,5), ylim=c(1,5), col="grey", pch=".", cex=2,  xlab="", ylab="",axes=F)
mtext(side=1, line=2, "log10 calculated")
mtext(side=2, line=2, "log10 measured")
mtext(side=3, line=0, expression(paste("  quiescent ", italic("S. pombe"))), adj=0)
mtext(side=3, line=-1.2, expression(paste("   ", italic(  k[sP]))), adj=0) #





#Myc data: 
#mine Lin et al nanostring data
Myc_mRNA_Lin_Table_S2 <- read_excel("Lin 2012 Table S2 myc gene expression nanostring.xlsx")
Myc_mRNA_Lin_Table_S2 <- Myc_mRNA_Lin_Table_S2[,c(2,4,5,8,9)]
Myc_mRNA_Lin_Table_S2$sum <- rowSums(Myc_mRNA_Lin_Table_S2[,2:5])
Myc_mRNA_Lin_Table_S2 <- Myc_mRNA_Lin_Table_S2[Myc_mRNA_Lin_Table_S2$sum>0,1:5]

Myc_mRNA_Lin_Table_S2$T0_mean <- rowMeans(Myc_mRNA_Lin_Table_S2[,2:3])
Myc_mRNA_Lin_Table_S2$T24_mean <- rowMeans(Myc_mRNA_Lin_Table_S2[,4:5])


#map Myc_mRNA_Lin_Table_S2 to Accession
# grab human Uniprot ID's & re-map both Myc_prot_Feist and Myc_mRNA_Feist onto the same Accession
library(readr)
Human_Uniprot_ID_map <- read_table2("Human Uniprot ID map.tab")
colnames(Human_Uniprot_ID_map) <- c("Accession", "Gene", "alt.name.1", "alt.name.2", "alt.name.3")

x <- strsplit(unlist(Human_Uniprot_ID_map[,2]), "_")
for(i in 1:nrow(Human_Uniprot_ID_map)){
  y <- unlist(x[i])
  Human_Uniprot_ID_map[i,2] <- y[1]
}


z <- as.data.frame(Myc_mRNA_Lin_Table_S2[,1]) 
x <- data.frame(matrix(NA, ncol = 2, nrow = nrow(z)))
for (i in 1:nrow(z)){ 
  if(length(which(rowSums(Human_Uniprot_ID_map == as.character(z[i,1]),  na.rm = TRUE) >0)) > 0){
    if(length(which(rowSums(Human_Uniprot_ID_map == as.character(z[i,1]),  na.rm = TRUE) >0)) >1){
      try(x[i,1] <-Human_Uniprot_ID_map[which(as.numeric(Human_Uniprot_ID_map$Gene == as.character(z[i,1])) >0),1],
          silent=T)
    }
    else{
      x[i,1] <-Human_Uniprot_ID_map[which(rowSums(Human_Uniprot_ID_map == as.character(z[i,1]),  na.rm = TRUE) >0),1]
    }
    x[i,2] <- as.character(z[i,1])
  }
  else{
    x[i,1] <- NA
    x[i,2] <- as.character(z[i,1])
  }
}
x <- na.omit(x)
colnames(x) <- c("Accession", "Name")
Myc_mRNA_Lin_Table_S2 <- merge(x, Myc_mRNA_Lin_Table_S2)
Myc_mRNA_Lin_Table_S2 <- unique(Myc_mRNA_Lin_Table_S2)


#grab Carbon metab genes from amiGO
library(readr)
Myc_amiGO_Cmetab <- read_delim("Myc_AmiGO_0006092_(288)_clean.txt", 
                               "\t", escape_double = FALSE, col_names = FALSE, 
                               trim_ws = TRUE)
#C-metab genes skip
Myc_amiGO_Cmetab$GR.adj <- 1
Myc_amiGO_Cmetab$N.adj <- 1
Myc_predict <- merge(Myc_mRNA_Lin_Table_S2, Myc_amiGO_Cmetab[,2:4], by.x="Name", by.y="X2", all.x=T)

#GR
x <- GR_mRNA[GR_mRNA$GR.cluster==3,1:3]
x <- merge(x, absq_fmol_mgDW[,c(1,46:48,61:63)])
Myc_predict$GR.adj[is.na(Myc_predict$GR.adj)] <- mean(unlist(x[,7:9]))/mean(unlist(x[,4:6]))

#NM (AA)
x <- AA_mRNA[AA_mRNA$AA.cluster==1,1:3]
x <- merge(x, absq_fmol_mgDW[,c(1,70:72,49:51)])
Myc_predict$N.adj[is.na(Myc_predict$N.adj)] <- mean(unlist(x[,4:6]))/mean(unlist(x[,7:9]))

Myc_predict$T24_predict <- Myc_predict$T0_mean * Myc_predict$GR.adj * Myc_predict$N.adj
Myc_predict$log2.mRNA.counts.over.predict <- log2(Myc_predict$T24_mean / Myc_predict$T24_predict)

#Fig 4C: 272x330
plot(log10(Myc_predict[Myc_predict$log2.mRNA.counts.over.predict > 1, "T24_predict"]), 
     log10(Myc_predict[Myc_predict$log2.mRNA.counts.over.predict > 1,"T24_mean"]), 
     xlim=c(-2.5,4), ylim=c(-2.5,4), col="blue", pch=".", cex=2, xlab="", ylab="", las=1, xaxt="n", yaxt="n") 
par(new=T)
plot(log10(Myc_predict[Myc_predict$log2.mRNA.counts.over.predict < -1, "T24_predict"]), 
     log10(Myc_predict[Myc_predict$log2.mRNA.counts.over.predict < -1, "T24_mean"]), 
     xlim=c(-2.5,4), ylim=c(-2.5,4), col="blue",  pch=".", cex=2,  xlab="", ylab="",axes=F)
par(new=T)
plot(log10(Myc_predict[Myc_predict$log2.mRNA.counts.over.predict < 1 & 
                         Myc_predict$log2.mRNA.counts.over.predict > -1, "T24_predict"]), 
     log10(Myc_predict[Myc_predict$log2.mRNA.counts.over.predict < 1 & 
                         Myc_predict$log2.mRNA.counts.over.predict > -1,"T24_mean"]), 
     xlim=c(-2.5,4), ylim=c(-2.5,4), col="grey", pch=".", cex=2,  xlab="", ylab="",axes=F)
axis(side=1,at=c(-2,0,2,4))
axis(side=2,at=c(-2,0,2,4))
mtext(side=1, line=2, "log10 calculated")
mtext(side=2, line=2, "log10 measured")
mtext(side=3, line=0, expression(paste("  MYC"^"high", " P493-6 cells")), adj=0)
mtext(side=3, line=-1, "  mRNA", adj=0)


#
## ksP
#use Lin 2012 to adjust Feist 2018 TPM into absq:
library(readr)
Myc_mRNA_Feist <- read_csv("Feist_TPM_all_nonzero_13623.csv")

#map Myc_mRNA_Feist to Accession
z <- as.data.frame(Myc_mRNA_Feist[,1]) 
x <- data.frame(matrix(NA, ncol = 2, nrow = nrow(z)))
for (i in 1:nrow(z)){ 
  if(length(which(rowSums(Human_Uniprot_ID_map == as.character(z[i,1]),  na.rm = TRUE) >0)) > 0){
    if(length(which(rowSums(Human_Uniprot_ID_map == as.character(z[i,1]),  na.rm = TRUE) >0)) >1){
      try(x[i,1] <-Human_Uniprot_ID_map[which(as.numeric(Human_Uniprot_ID_map$Gene == as.character(z[i,1])) >0),1],
          silent=T)
    }
    else{
      x[i,1] <-Human_Uniprot_ID_map[which(rowSums(Human_Uniprot_ID_map == as.character(z[i,1]),  na.rm = TRUE) >0),1]
    }
    x[i,2] <- as.character(z[i,1])
  }
  else{
    x[i,1] <- NA
    x[i,2] <- as.character(z[i,1])
  }
}
x <- na.omit(x)
colnames(x) <- c("Accession", "Gene.ID")
Myc_mRNA_Feist <- merge(x, Myc_mRNA_Feist)
Myc_mRNA_Feist <- unique(Myc_mRNA_Feist)

Myc_mRNA_ruler <- merge(Myc_mRNA_Feist, Myc_mRNA_Lin_Table_S2[,c(2,7,8)])
Myc_mRNA_ruler[Myc_mRNA_ruler==0] <- NA
Myc_mRNA_ruler <- na.omit(Myc_mRNA_ruler)

summary(lm(log10(Myc_mRNA_ruler$T0_mean) ~ log10(Myc_mRNA_ruler$mycOFF.mean)))#lm(y~x)
m <- summary(lm(log10(Myc_mRNA_ruler$T0_mean) ~ log10(Myc_mRNA_ruler$mycOFF.mean)))$coefficients[2,1]
b <- summary(lm(log10(Myc_mRNA_ruler$T0_mean) ~ log10(Myc_mRNA_ruler$mycOFF.mean)))$coefficients[1,1]
Myc_mRNA_ruler$mycOFF.abs <- 10^(m * log10(Myc_mRNA_ruler$mycOFF.mean) + b)
Myc_mRNA_Feist$mycOFF.abs <- 10^(m * log10(Myc_mRNA_Feist$mycOFF.mean) + b)

remove(m,b)

summary(lm(log10(Myc_mRNA_ruler$T24_mean) ~ log10(Myc_mRNA_ruler$mycON.mean))) #lm(y~x)
m <- summary(lm(log10(Myc_mRNA_ruler$T24_mean) ~ log10(Myc_mRNA_ruler$mycON.mean)))$coefficients[2,1]
b <- summary(lm(log10(Myc_mRNA_ruler$T24_mean) ~ log10(Myc_mRNA_ruler$mycON.mean)))$coefficients[1,1]
Myc_mRNA_ruler$mycON.abs <- 10^(m * log10(Myc_mRNA_ruler$mycON.mean) + b)
Myc_mRNA_Feist$mycON.abs <- 10^(m * log10(Myc_mRNA_Feist$mycON.mean) + b)

remove(m,b, Myc_mRNA_ruler)


#grab protein half-life from Boisvert 2012 
library(readxl)
Myc_kdP_Boisvert <- read_excel("Boisvert 2012 Table S2 protein turnover HeLa.xlsx")
Myc_kdP_Boisvert <- Myc_kdP_Boisvert[,c(2,5)]
Myc_kdP_Boisvert <- na.omit(Myc_kdP_Boisvert)
# remove those that are 999 (ie infinite) 
Myc_kdP_Boisvert <- Myc_kdP_Boisvert[Myc_kdP_Boisvert$`Whole Half-Life`<999,]
Myc_kdP_Boisvert$kdP <- log(2) / Myc_kdP_Boisvert$`Whole Half-Life`

#remap Myc_kdP_Boisvert to UniProt Accession
z <- as.data.frame(Myc_kdP_Boisvert[,1]) 
x <- data.frame(matrix(NA, ncol = 2, nrow = nrow(z)))
for (i in 1:nrow(z)){ 
  if(length(which(rowSums(Human_Uniprot_ID_map == as.character(z[i,1]),  na.rm = TRUE) >0)) > 0){
    if(length(which(rowSums(Human_Uniprot_ID_map == as.character(z[i,1]),  na.rm = TRUE) >0)) >1 & 
       nrow(Human_Uniprot_ID_map[which(as.numeric(Human_Uniprot_ID_map$Gene == as.character(z[i,1])) >0),])>0){
      x[i,1] <-Human_Uniprot_ID_map[which(as.numeric(Human_Uniprot_ID_map$Gene == as.character(z[i,1])) >0),1]
    }
    else{
      x[i,1] <-Human_Uniprot_ID_map[which(rowSums(Human_Uniprot_ID_map == as.character(z[i,1]),  na.rm = TRUE) >0),1]
    }
    x[i,2] <- as.character(z[i,1])
  }
  else{
    x[i,1] <- NA
    x[i,2] <- as.character(z[i,1])
  }
}
x <- na.omit(x)
colnames(x) <- c("Accession", "Gene")
Myc_kdP_Boisvert <- merge(x, Myc_kdP_Boisvert)
Myc_kdP_Boisvert <- unique(Myc_kdP_Boisvert)


#grab protein from Feist 2018
library(readxl)
Myc_prot_Feist <- read_excel("Feist 2018 myc proteomics.xlsx")

Myc_prot_Feist <- Myc_prot_Feist[,c(1:4,8:10)]
Myc_prot_Feist$low.mean <- rowMeans(Myc_prot_Feist[,2:4])
Myc_prot_Feist$high.mean <- rowMeans(Myc_prot_Feist[,5:7])

#re-map Myc_prot_Feist to UniProt Accession
z <- as.data.frame(Myc_prot_Feist[,1]) 
x <- data.frame(matrix(NA, ncol = 2, nrow = nrow(z)))
for (i in 1:nrow(z)){ 
  if(length(which(rowSums(Human_Uniprot_ID_map == as.character(z[i,1]),  na.rm = TRUE) >0)) > 0){
    if(length(which(rowSums(Human_Uniprot_ID_map == as.character(z[i,1]),  na.rm = TRUE) >0)) >1){
      try(x[i,1] <-Human_Uniprot_ID_map[which(as.numeric(Human_Uniprot_ID_map$Gene == as.character(z[i,1])) >0),1],
          silent=T)
    }
    else{
      x[i,1] <-Human_Uniprot_ID_map[which(rowSums(Human_Uniprot_ID_map == as.character(z[i,1]),  na.rm = TRUE) >0),1]
    }
    x[i,2] <- as.character(z[i,1])
  }
  else{
    x[i,1] <- NA
    x[i,2] <- as.character(z[i,1])
  }
}
x <- na.omit(x)
colnames(x) <- c("Accession", "Protein")

Myc_prot_Feist <- merge(x, Myc_prot_Feist)
Myc_prot_Feist <- unique(Myc_prot_Feist)

#calculate ksP
Myc_predict_protein <- merge(Myc_prot_Feist[,c(1:2,9,10)], Myc_mRNA_Feist[,c(2,5,6)])
colnames(Myc_predict_protein)[3:6] <- c("mycOFF.prot", "mycON.prot", "mycOFF.mRNA", "mycON.mRNA")

Myc_predict_protein <- merge(Myc_predict_protein, Myc_kdP_Boisvert[,c(2,4)], all.x=T)
Myc_predict_protein$kdP[is.na(Myc_predict_protein$kdP)] <- median(Myc_kdP_Boisvert$kdP)
Myc_predict_protein <- unique(Myc_predict_protein)


Myc_predict_protein$mycOFF.ksP <- Myc_predict_protein$mycOFF.prot / Myc_predict_protein$mycOFF.mRNA * (0.01+Myc_predict_protein$kdP)
Myc_predict_protein$mycON.ksP <- Myc_predict_protein$mycON.prot / Myc_predict_protein$mycON.mRNA * (0.13+Myc_predict_protein$kdP)

#use mycOFF.ksp to predict mycON.ksp
Myc_predict_protein$mycON.ksP.pred <- Myc_predict_protein$mycOFF.ksP * (mean(unlist(ksP[,19:21]))/mean(unlist(ksP[,4:6]))) * (mean(unlist(ksP[,28:30]))/mean(unlist(ksP[,7:9])))
Myc_predict_protein$log2.ksP.measured.over.predict <- log2(Myc_predict_protein$mycON.ksP / Myc_predict_protein$mycON.ksP.pred)

#Fig 4D: 272x330
plot(log10(Myc_predict_protein[Myc_predict_protein$log2.ksP.measured.over.predict > 1, "mycON.ksP.pred"]), 
     log10(Myc_predict_protein[Myc_predict_protein$log2.ksP.measured.over.predict > 1,"mycON.ksP"]), 
     xlim=c(-1,4.5), ylim=c(-1,4.5), col="blue", pch=".", cex=2, xlab="", ylab="", las=1) 
par(new=T)
plot(log10(Myc_predict_protein[Myc_predict_protein$log2.ksP.measured.over.predict < -1, "mycON.ksP.pred"]), 
     log10(Myc_predict_protein[Myc_predict_protein$log2.ksP.measured.over.predict < -1, "mycON.ksP"]), 
     xlim=c(-1,4.5), ylim=c(-1,4.5), col="blue",  pch=".", cex=2,  xlab="", ylab="",axes=F)
par(new=T)
plot(log10(Myc_predict_protein[Myc_predict_protein$log2.ksP.measured.over.predict < 1 & 
                                 Myc_predict_protein$log2.ksP.measured.over.predict > -1, "mycON.ksP.pred"]), 
     log10(Myc_predict_protein[Myc_predict_protein$log2.ksP.measured.over.predict < 1 & 
                                 Myc_predict_protein$log2.ksP.measured.over.predict > -1,"mycON.ksP"]), 
     xlim=c(-1,4.5), ylim=c(-1,4.5), col="grey", pch=".", cex=2,  xlab="", ylab="",axes=F)
mtext(side=1, line=2, "log10 calculated")
mtext(side=2, line=2, "log10 measured")
mtext(side=3, line=0, expression(paste("  MYC"^"high", " P493-6 cells")), adj=0)
mtext(side=3, line=-1.2, expression(paste("   ", italic(  k[sP]))), adj=0) 



