### Figure 2a ###

### ggplot2 bargraph with error bars ###
```
no320  <- read.table("/home/xwzhu/transfer/Project_VAC_10942_B01_GRM_WGS.2015-09-18/hg38/prog/strand/09_comb_SR_PE/noMDA.v41.txt", header=F, sep="\t")
no320F1 <- 2*no320[,2]*no320[,5] / (no320[,2]+no320[,5])
no320F2 <- 2*no320[,3]*no320[,6] / (no320[,3]+no320[,6])
no320F3 <- 2*no320[,4]*no320[,7] / (no320[,4]+no320[,7])
no320L1 <- no320[,8]

no316  <- read.table("/home/xwzhu/transfer/Project_VAC_11092_B01_GRM_WGS.2015-09-20/hg38/prog/strand/09_comb_SR_PE/noMDA.v41.txt", header=F, sep="\t")
no316F1 <- 2*no316[,2]*no316[,5] / (no316[,2]+no316[,5])
no316F2 <- 2*no316[,3]*no316[,6] / (no316[,3]+no316[,6])
no316F3 <- 2*no316[,4]*no316[,7] / (no316[,4]+no316[,7])
no316L1 <- no316[,8]

mda316 <- read.table("/home/xwzhu/transfer/Project_VAC_11092_B01_GRM_WGS.2015-09-20/hg38/prog/strand/09_comb_SR_PE/MDA.v41.txt", header=F, sep="\t")
mda316F1 <- 2*mda316[,2]*mda316[,5] / (mda316[,2]+mda316[,5])
mda316F2 <- 2*mda316[,3]*mda316[,6] / (mda316[,3]+mda316[,6])
mda316F3 <- 2*mda316[,4]*mda316[,7] / (mda316[,4]+mda316[,7])
mda316L1 <- mda316[,8]

temp <-c()
boxplot(mda316[,2], mda316[,5], mda316F1, temp, mda316[,3], mda316[,6], mda316F2, temp, mda316[,4], mda316[,7], mda316F3, temp, no316[,2], no316[,5], no316F1, temp, no316[,3], no316[,6], no316F2, temp, no316[,4], no316[,7], no316F3, temp, no320[,2], no320[,5], no320F1, temp, no320[,3], no320[,6], no320F2, temp, no320[,4], no320[,7], no320F3, ylim=c(0,1), xaxt="n", main="Line1 SR/PE support", col=c(2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4))
axis(1,at=c(2,6,10,14,18,22,26,30,34), labels=c("316MDA_PE", "316MDA_SR", "316MDA_both", "316noMDA_PE", "316noMDA_SR", "316noMDA_both","320noMDA_PE", "320noMDA_SR", "320noMDA_both"), las=1, cex.axis=1.2)
legend(32, 0.15, c("recall","precision","Fscore"), lty=c(1,1,1), cex=c(1.4,1.4,1.4), lwd=c(2,2,2),col=c(2,3,4))

### sensitivity ###
library("ggplot2")
sens <- c()
sens316m <- data.frame(rep("316MDA", 5), mda316[,4])
names(sens316m) <- c("subject", "sensitivity")
sens316n <- data.frame(rep("316noMDA", 10), no316[,4])
names(sens316n) <- c("subject", "sensitivity")
sens320n <- data.frame(rep("320noMDA", 53), no320[,4])
names(sens320n) <- c("subject", "sensitivity")
sens <- rbind(sens316m, sens316n, sens320n)
sum_sens <- summarySE(sens, measurevar="sensitivity", groupvars=c("subject"))

ggplot(sum_sens, aes(x=subject, y=sensitivity, fill=subject)) +
        geom_bar(position=position_dodge(), stat="identity",
         colour="black", # Use black outlines,
         size=.3) +      # Thinner lines+
        geom_errorbar(aes(ymin=sensitivity-ci, ymax=sensitivity+ci),
                  size=1,    # Thinner lines
                  width=.2,
                  position=position_dodge(.9)) +
    xlab("Test datasets") +
    ylab("Sensitivity") +
    guides(fill=FALSE) +
    scale_y_continuous(limits = c(0, 1.05), breaks=c(0.2, 0.4, 0.6, 0.8, 1)) +
    theme(axis.text=element_text(size=18),
        axis.title=element_text(size=24,face="bold"))
ggsave("/home/xwzhu/levinson/MEI_plots/L1_sensitvity_v41.pdf")

### precision ###
prec<- c()
prec316m <- data.frame(rep("316MDA", 5), mda316[,7])
names(prec316m) <- c("subject", "precision")
prec316n <- data.frame(rep("316noMDA", 10), no316[,7])
names(prec316n) <- c("subject", "precision")
prec320n <- data.frame(rep("320noMDA", 53), no320[,7])
names(prec320n) <- c("subject", "precision")
prec <- rbind(prec316m, prec316n, prec320n)
sum_prec <- summarySE(prec, measurevar="precision", groupvars=c("subject"))

ggplot(sum_prec, aes(x=subject, y=precision, fill=subject)) +
        geom_bar(position=position_dodge(), stat="identity",
         colour="black", # Use black outlines,
         size=.3) +      # Thinner lines+
        geom_errorbar(aes(ymin=precision-ci, ymax=precision+ci),
                  size=1,    # Thinner lines
                  width=.2,
                  position=position_dodge(.9)) +
    xlab("Test datasets") +
    ylab("Precision") +
    guides(fill=FALSE) +
    scale_y_continuous(limits = c(0, 1.05), breaks=c(0.2, 0.4, 0.6, 0.8, 1)) +
    theme(axis.text=element_text(size=18),
        axis.title=element_text(size=24,face="bold"))
ggsave("/home/xwzhu/levinson/MEI_plots/L1_precision_v41.pdf")

### false positives ###
FP316m <- data.frame(rep("316MDA", 5), mda316L1)
names(FP316m) <- c("subject", "FP")
FP316n <- data.frame(rep("316noMDA", 10), no316L1)
names(FP316n) <- c("subject", "FP")
FP320n <- data.frame(rep("320noMDA", 53), no320L1)
names(FP320n) <- c("subject", "FP")
FP <- rbind(FP316m, FP316n, FP320n)

ggplot(FP, aes(subject, FP, color=subject)) +
   geom_jitter(width = 0.15, ,show.legend=F) +
   theme(axis.text=element_text(size=18),
       axis.title=element_text(size=24,face="bold")) +
   xlab("Test datasets") +
   ylab("False positives") +
   scale_y_continuous(limits = c(-2, 100), breaks=c(0, 5, 25, 50, 75, 100))
ggsave("/home/xwzhu/levinson/MEI_plots/L1_FP_v41.pdf")

### Fscore ###
Fscore <- c()
fs316m <- data.frame(rep("316MDA", 5), 2*mda316[,4]*mda316[,7]/(mda316[,4]+mda316[,7]))
names(fs316m) <- c("subject", "Fscore")
fs316n <- data.frame(rep("316noMDA", 10), 2*no316[,4]*no316[,7]/(no316[,4]+no316[,7]))
names(fs316n) <- c("subject", "Fscore")
fs320n <- data.frame(rep("320noMDA", 53), 2*no320[,4]*no320[,7]/(no320[,4]+no320[,7]))
names(fs320n) <- c("subject", "Fscore")
Fscore <- rbind(fs316m, fs316n, fs320n)
sum_fs <- summarySE(Fscore, measurevar="Fscore", groupvars=c("subject"))

ggplot(sum_fs, aes(x=subject, y=Fscore, fill=subject)) +
        geom_bar(position=position_dodge(), stat="identity",
         colour="black", # Use black outlines,
         size=.3) +      # Thinner lines+
        geom_errorbar(aes(ymin=Fscore-ci, ymax=Fscore+ci),
                  size=1,    # Thinner lines
                  width=.2,
                  position=position_dodge(.9)) +
    xlab("Test datasets") +
    ylab("F score") +
    guides(fill=FALSE) +
    scale_y_continuous(limits = c(0, 1), breaks=c(0.2, 0.4, 0.6, 0.8, 1)) +
    theme(axis.text=element_text(size=18),
        axis.title=element_text(size=24,face="bold"))
ggsave("/home/xwzhu/levinson/MEI_plots/L1_Fscore_v41.pdf")
```
