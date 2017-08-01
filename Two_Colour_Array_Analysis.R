library("limma")
library("CCl4")
dataPath = system.file("extdata", package="CCl4")
dir(dataPath)

#
adf = read.AnnotatedDataFrame("samplesInfo.txt",
                              path=dataPath)
targets = pData(adf)
targets$FileName = row.names(targets)

# Read in the data
RG = read.maimages(targets, path=dataPath, source="genepix")
head(RG$genes)

par_def = par() #save default par() settings

par(mfrow=c(5,1))
imageplot(log2(RG$Rb[,1]), RG$printer, low="white",
          high="red")
imageplot(log2(RG$Gb[,1]), RG$printer, low="white",
          high="green")
imageplot(rank(RG$Rb[,1]), RG$printer, low="white",
          high="red")
imageplot(rank(RG$Gb[,1]), RG$printer, low="white",
          high="green")
imageplot(rank(log(RG$R[,1])+log(RG$G[,1])),
          RG$printer, low="white", high="blue")

# An M A-plot displays the log-ratio of red intensities R and green intensities G on the y-axis versus the overall intensity of each spot on the x-axis.
MA = normalizeWithinArrays(RG, method="none", bc.method="none")

library("geneplotter")
smoothScatter(MA$A[, 1], MA$M[, 1], xlab="A", ylab="M")
abline(h=0, col="red")

# Boxplot
plotformula = log2(RG$G)~col(RG$G)
boxplot(plotformula, ylim=c(5,9), outline=FALSE,
        col="forestgreen", xlab="arrays",
        ylab=expression(log[2]~G), main="boxplot")

# the distributions of the G values on the different arrays with the multidensity function
multidensity(plotformula, xlim=c(5,9),
             main="densities", xlab=expression(log[2]~G))

# Create subset of six arrays which were hybridized with the good RNA 
# sample of the CCl 4 treated hepatocytes
rin = with(MA$targets, ifelse(Cy5=="CCl4", RIN.Cy5,
                              RIN.Cy3))
select = (rin == max(rin))
RGgood = RG[, select]

# created a subset version adfgood of the annotated dataframe adf
adfgood = adf[select, ]

# use the function justvsn from the vsn package to normalize
# the data from these 6 arrays
library("vsn")
ccl4 = justvsn(RGgood, backgroundsubtract=TRUE)

r = assayData(ccl4)$R
g = assayData(ccl4)$G

# scatterplot of standard deviations versus 
#the rank of the mean intensity of each feature.
meanSdPlot(cbind(r, g))

rownames(pData(adfgood)) = sub("\\.gpr$", "",
                               rownames(pData(adfgood)))
pData(adfgood)

# update the channel information in the varMetadata part of adfgood
varMetadata(adfgood)$channel = factor(c("G", "R", "G", "R"),
                                      levels = c("G", "R", "_ALL_"))

phenoData(ccl4) = adfgood
validObject(ccl4)                                      
ccl4AM = ccl4
assayData(ccl4AM) = assayDataNew(A=(r+g)/2, M=r-g)
