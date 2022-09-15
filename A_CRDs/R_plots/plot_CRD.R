library(qvalue)

#LOAD MODULES
MOD = as.data.frame(data.table::fread("zcat ~/VitalIT/SGX/V2/data_trees/modules/LCL_ALL.chrALL.module.txt.gz", head=FALSE, stringsAsFactors=FALSE))
colnames(MOD)=c("IDX", "CHILD1", "CHILD2", "UUID", "COUNT", "ACORR", "APVAL", "BCORR", "BPVAL", "CCORR", "CPVAL", "DCORR", "DPVAL", "ECORR", "EPVAL", "PCA1", "DIST", "START", "STOP", "SIZE", "NRE", "ANN1", "ANN2", "ANN3", "CMB1", "CMB2", "CMB3", "CMB4", "CMB5", "CMB6", "CMB7", "CMB8", "BR_RATIO", "DISPER", "COMPL", "LIDX", "RIDX", "N_SIB", "N_REP", "MOD")
MOD = MOD[MOD$NRE > 1 & MOD$MOD == 1, ]

#LOAD PAIRWISE CORRELATION
DATA = as.data.frame(data.table::fread("zcat ~/VitalIT/SGX/V2/data_correlation/cis/modules/LCL_ALL.chrALL.txt.gz", head=FALSE, stringsAsFactors=FALSE))
colnames(DATA) = c("idx1","idx2","pos1","pos2","mid1","mid2","corr","pval")
Q = qvalue(DATA$pval)
DATAs = DATA[Q$qvalue < 0.01, ]

#WRITE SIGNIFICANT HITS
write.table(DATAs, "~/VitalIT/SGX/V2/data_correlation/cis/modules/LCL_ALL.chrALL.significant1FDR.txt", quote=FALSE, row.names=FALSE, col.names=TRUE)

#PLOT1: HISTOGRAM + PI1
pdf("figureC.1.1.pdf",10,5)
par(mfrow=c(1, 2))
hist(DATA$pval, xlab="Nominal P-values", main="", breaks=100)
abline(h=Q$pi0 * nrow(DATA) / 100, col="red")
legend("topright", legend=c(paste("pi1=", signif(100-Q$pi0*100, 3), "%", sep=""), paste("#hits=", nrow(DATAs), " (1% FDR)", sep="")), bty="n")
hist(DATAs$corr, breaks=100, main="", xlab="Correlation coefficient")
dev.off()

#PLOT_OPT: ALL CHROMOSOMES
CORR = as.data.frame(data.table::fread("zcat ~/VitalIT/SGX/V2/data_correlation/cis/chromatin/LSF_ALL.chrALL.subset.txt.gz", head=FALSE, stringsAsFactors=FALSE))
MOD$chr = matrix(unlist(strsplit(MOD$UUID, "_")), ncol=3, byrow=TRUE)[,1]
for (c in 1:22) {
	cat("Processing ", c, "\n")
	CORRs = CORR[CORR$V2 == paste("chr", c, sep=""), ]
	MODs = MOD[MOD$chr == paste("chr", c, sep="") & MOD$COUNT > 10, ]
	png(paste("chromosomes/chr", c, ".png", sep=""), max(CORRs$V5)/2, 250)
	plot((CORRs$V1+CORRs$V5)/2, CORRs$V5 - CORRs$V1, pch=18, col=ifelse(CORRs$V10 > 0, rgb(0,0,1, CORRs$V10^2), rgb(1,0,0, CORRs$V10^2)), xlab = paste("peak index on chr", c, sep=""), cex=0.5)
	for (m in 1:nrow(MODs)) {
		polygon(c(MODs$LIDX[m],(MODs$LIDX[m]+MODs$RIDX[m])/2,MODs$RIDX[m]), c(1,MODs$RIDX[m]-MODs$LIDX[m],1))
	}
	dev.off()
}

#PLOT2: EXAMPLE REGION
CORR = as.data.frame(data.table::fread("zcat ~/VitalIT/SGX/V2/data_correlation/cis/chromatin/LCL_ALL.chr8.subset.txt.gz", head=FALSE, stringsAsFactors=FALSE))
MODS = as.data.frame(data.table::fread("zcat ~/VitalIT/SGX/V2/data_trees/modules/LCL_ALL.chr8.module.txt.gz", head=TRUE, stringsAsFactors=FALSE))
MODSs = MODS[MODS$NRE > 1 & MODS$MOD == 1 & MODS$COUNT > 5, ]
png("figureC.1.2.png", width=3000, height = 1000, res = 200)
col=ifelse(CORR$V10 > 0, rgb(0,0,1,abs(CORR$V10)^2),rgb(1,0,0,abs(CORR$V10)^2))
plot((CORR$V1+CORR$V5)/2, CORR$V5 - CORR$V1, pch=18, col=col, xlim=c(2500, 3500), ylim=c(0, 250), bty="n", ylab="", xlab="Chromatin peak index on chromosome 8", yaxt="n")
abline(h=0)
for (m in 1:nrow(MODSs)) {
	polygon(c(MODSs$LIDX[m], (MODSs$LIDX[m]+MODSs$RIDX[m])/2, MODSs$RIDX[m]), c(0,MODSs$RIDX[m]-MODSs$LIDX[m], 0), border="black", lwd=2)
}
legend("topright", fill=c(rgb(0,0,1), rgb(1,0,0)), legend=c("Positive correlation","Negative correlation"))
dev.off()

#PLOT1: PI1 across tissues as a function of distance
library(RColorBrewer)
COL = brewer.pal(8,"Set1")
LSF = as.data.frame(data.table::fread("zcat ~/VitalIT/SGX/V2/data_correlation/cis/chromatin/LSF_ALL.chrALL.distance.txt.gz", head=FALSE, stringsAsFactors=FALSE))
LSF$TYPE0=substr(LSF$V4, 1,7)
LSF$TYPE1=substr(LSF$V8, 1,7)
HICp=rep(0, nBIN)
HICn=rep(0, nBIN)
CORp=rep(0, nBIN)
CORn=rep(0, nBIN)
LSFs = LSF[LSF$V9 > 0, ]
LSFs$CCORR=cut(abs(LSFs$V10), breaks=quantile(abs(LSFs$V10), probs=seq(0,1,0.02)), labels=1:50, include.lowest=TRUE)
for (b in 1:50) {
	cat("Processing bin ", b, "/50\n")
	LSFss = LSFs[LSFs$CCORR==b,]
	LSFssp = LSFss[LSFss$V10 > 0, ]
	LSFssn = LSFss[LSFss$V10 < 0, ]
	HICp[b] = mean(LSFssp$V9)
	HICn[b] = mean(LSFssn$V9)
	CORp[b] = abs(mean(LSFssp$V10))
	CORn[b] = abs(mean(LSFssn$V10))
}
pdf("figureC.1.3.pdf", 5, 5)
plot(0, 0, main="", xlab="Absolute correlation value", ylab="Mean Hi-C KRnorm intensities", type="n", pch=20, ylim=c(0,200), xlim=c(0,0.4))
points(CORp, HICp, type="b", pch=20, col=COL[2])
points(CORn, HICn, type="b", pch=20, col=COL[1])
legend("topleft", legend=c("Positively correlated", "Negatively correlated"), fill=COL[c(2,1)], bty="n")
dev.off()

