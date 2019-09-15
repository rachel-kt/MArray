#Set the working directory
source("http://bioconductor.org/biocLite.R")

setwd("enter path")
library(simpleaffy)
library(help = "simpleaffy")

dir()

# read .CEL files and create an Affybatch object
AT.raw.data = read.affy("covdesc.txt")

# look up Affybatch help, class and structure of Affybatch
?AffyBatch
class(AT.raw.data)
str(AT.raw.data)

# QC metric
AT.qc = qc(AT.raw.data) # makes a QCstat object
plot(AT.qc) 

# view help files for journalpng function
?journalpng

# to save the diagram in a .png format use the journalpng function
journalpng(file="qc_plot_12_27.png", height=6, width=6)
plot(AT.qc)
dev.off() 

# NUSE plots

AT.plm = fitPLM(AT.raw.data)
pdf("nuse_plot.pdf")

# set the margins to be able to view the sample names. Default margins are (b,l,t,r) = (5.1,4.1,4.1,2.1)

par(mar=c(11.5,4.1,4.1,2)) 
boxplot(AT.plm, main = "NUSE", ylim = c(0.95,1.1), outline = FALSE, col = c("lightgreen","lightgreen","lightgreen","lightblue","lightblue","lightblue"), las=3, at = c(1,2,3,5,6,7), whiskity = 0, staplelty = 0, names = c("12_deg_Rep_1","12_deg_Rep_2","12_deg_Rep_3","17_deg_Rep_1","17_deg_Rep_2","17_deg_Rep_3"))
dev.off()

# RLE plot

pdf("RLE_plot.pdf")
par(mar=c(11.5,4.1,4.1,2))
Mbox(AT.plm, main="RLE", ylim = c(-0.4, 0.4), outline = FALSE, col="mistyrose", las=3, whisklty=0, staplelty=0, names = c("12_deg_Rep_1","12_deg_Rep_2","12_deg_Rep_3","17_deg_Rep_1","17_deg_Rep_2","17_deg_Rep_3"))
dev.off()

# report normalized summary data for each probeset in each array

AT.rma = call.exprs(AT.raw.data, "rma")

# Pairwise comparison of each probeset (grouped)
# compares each gene(probeset) from all the samples and calculates the mean, t-test, fold-change and call(P,M,A)

AT.results = pairwise.comparison(AT.rma, "Temperature", c("tl", "ts"), AT.raw.data)

# write the data in the new slots of the PairComp object into a dataframe to be able to view it

write.table(data.frame(means(AT.results), fc(AT.results), tt(AT.results), call(AT.results)))

# get the annotations for the data
annotation(AT.raw.data)

myAnnot <- data.frame(ACCNUM=sapply(contents(ath1121501ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(ath1121501SYMBOL), paste, collapse=", "), DESC=sapply(contents(ath1121501GENENAME), paste, collapse=", "))
myAnnot <- data.frame(rownames(myAnnot), ACCNUM=sapply(contents(ath1121501ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(ath1121501SYMBOL), paste, collapse=", "), DESC=sapply(contents(ath1121501GENENAME), paste, collapse=", "))

Pairwise_result = data.frame(means(AT.results), fc(AT.results), tt(AT.results), calls(AT.results))
View(Pairwise_result)
idx_n = match(rownames(Pairwise_result), myAnnot$rownames.myAnnot.)
Pairwise_result = cbind(myAnnot[idx_n,"ACCNUM"], Pairwise_result)

> write.table(data.frame(Pairwise_result, myAnnot[idx_n,"DESC"]), file="pairwise_t.xls", sep = "\t", quote=F, col.names = c("PROBEID\tLOCUSID", "12_DEG_MEAN", "27_DEG_MEAN", "FOLD_CHANGE", "T_TEST", "PRESENT_CALL_12_R1", "PRESENT_CALL_12_R2", "PRESENT_CALL_12_R3", "PRESENT_CALL_27_R1", "PRESENT_CALL_27_R2", "PRESENT_CALL_27_R3", "SYMBOL", "GENENAMES"))

write.table(data.frame(N.table$ACCNUM, means(AT.results), fc(AT.results), tt(AT.results), calls(AT.results), N.table$SYMBOL,  N.table$GENENAME), file="pairwise_result_with_map.xls", sep="\t", quote=F,col.names = c("PROBEID\tLOCUSID", "12_DEG_MEAN", "27_DEG_MEAN", "FOLD_CHANGE", "T_TEST", "PRESENT_CALL_12_R1", "PRESENT_CALL_12_R2", "PRESENT_CALL_12_R3", "PRESENT_CALL_27_R1", "PRESENT_CALL_27_R2", "PRESENT_CALL_27_R3", "SYMBOL", "GENENAMES"))
# write all info in a table form into a .xls file
write.table(data.frame(means(AT.results), fc(AT.results), tt(AT.results), calls(AT.results), AT.probe_map$ENTREZID, AT.probe_map$SYMBOL,  AT.probe_map$GENENAME), file="pairwise_comp_with_map.xls", sep="\t", quote=F,col.names = c("PROBEID", "12_DEG_MEAN", "27_DEG_MEAN", "FOLD_CHANGE", "T_TEST", "PRESENT_CALL_12_R1", "PRESENT_CALL_12_R2", "PRESENT_CALL_12_R3", "PRESENT_CALL_27_R1", "PRESENT_CALL_27_R2", "PRESENT_CALL_27_R3\tENTREZID", "SYMBOL", "GENENAMES"))

significant = pairwise.filter(AT.results, fc=log2(4), min.present.no = 6, tt = 0.01, present.by.group = FALSE)

top_probesets = data.frame(means(significant), fc(significant), tt(significant), calls(significant))


lowest_tvalue_probesets = sort(abs(top_probesets$tt), decreasing = F)

# For plotting the volcano plot

P_comp_table = data.frame(means(AT.results), fc(AT.results), tt(AT.results), calls(AT.results))
all_tvalue_probesets = P_comp_table$tt
all_fc_probesets = P_comp_table$fc
lod = -log10(all_tvalue_probesets)
plot(all_fc_probesets, lod, pch = ".", xlab = "fold change", ylab = expression(log[10]~p))
fc = data.frame(fc(significant), -log10(tt(significant)))
points(fc$fc.significant., fc$X.log10.tt.significant.., pch = 1, col = "lightpink")

ls(ath1121501cdf)
contents(ath1121501ACCNUM)
myAnnot <- data.frame(ACCNUM=sapply(contents(ath1121501ACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(ath1121501SYMBOL), paste, collapse=", "), DESC=sapply(contents(ath1121501GENENAME), paste, collapse=", "))

