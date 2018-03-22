#!/usr/bin/env Rscript


# setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

library("getopt")
options(stringAsFactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec <- matrix(c(
  "help", "h", 0, "logical",
  "countsfile", "c", 1, "character",
  "traitsfile", "t", 1, "character",
  "modulecolors", "m", 1, "character",
  "probes", "p", 1, "character",
  "visant", "v", 1, "character",
  "imageplot", "i", 1, "character",
  "mes", "a", 1, "character"),
  byrow=TRUE, ncol=4)
opt <- getopt(spec)

# if help was asked for print a friendly message
# and exit with a non-zero error code
if (!is.null(opt$help)) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}


library(WGCNA)
library(flashClust)


#setwd(opt$workingdir);



# Uploading data into R and formatting it for WGCNA --------------
# This creates an object called "datExpr" that contains the normalized counts file output from DESeq2
datExpr = read.csv(opt$countsfile)
# "head" the file to preview it
head(datExpr) # You see that genes are listed in a column named "X" and samples are in columns
 

# Manipulate file so it matches the format WGCNA needs
row.names(datExpr) = datExpr$X
datExpr$X = NULL
datExpr = as.data.frame(t(datExpr)) # now samples are rows and genes are columns
dim(datExpr) # 48 samples and 1000 genes (you will have many more genes in reality)
 
 
# Run this to check if there are gene outliers
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK
 
 
 
#Create an object called "datTraits" that contains your trait data
datTraits = read.csv(opt$traitsfile)
head(datTraits)


#form a data frame analogous to expression data that will hold the clinical traits.
rownames(datTraits) = datTraits$Sample
datTraits$Sample = NULL
table(rownames(datTraits)==rownames(datExpr)) #should return TRUE if datasets align correctly, otherwise your names are out of order

pdf(opt$imageplot,width=6,height=4,paper='special') 
 
# You have finished uploading and formatting expression and trait data
# Expression data is in datExpr, corresponding traits are datTraits
 
 
#save(datExpr, datTraits, file="SamplesAndTraits.RData")
#load("SamplesAndTraits.RData")
  
# Cluster samples by expression ----------------------------------------------------------------
A = adjacency(t(datExpr),type="signed") # this calculates the whole network connectivity
k = as.numeric(apply(A,2,sum))-1 # standardized connectivity
Z.k = scale(k)
thresholdZ.k = -2.5 # often -2.5
outlierColor = ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")

 
# Convert traits to a color representation where red indicates high values
traitColors = data.frame(numbers2colors(datTraits,signed=FALSE))
dimnames(traitColors)[[2]] = paste(names(datTraits))
datColors = data.frame(outlier = outlierColor,traitColors)
plotDendroAndColors(sampleTree,groupLabels=names(datColors),
                    colors=datColors,main="Sample Dendrogram and Trait Heatmap")
 
 
# Day "0" outliers have been identified. You could exclude these samples with the code below, but this scientists had a biological reason to NOT exclude these samples. It's up to you. Justify whatever decision you make.
 
 
#Remove outlying samples
#remove.samples = Z.k<thresholdZ.k | is.na(Z.k)
#datExprOut = datExpr[!remove.samples,]
#datTraitsOut = datTraits[!remove.samples,]
#save(datExprOut, datTraitsOut, file="SamplesAndTraits_OutliersRemoved.RData")

enableWGCNAThreads()
softPower = 18
adjacency = adjacency(datExpr, power = softPower, type = "signed") #specify network type
 
 
# Construct Networks- USE A SUPERCOMPUTER IRL -----------------------------
#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM = TOMsimilarity(adjacency, TOMType="signed") # specify network type
dissTOM = 1-TOM
 
# Generate Modules --------------------------------------------------------
 
 
# Generate a clustered gene tree
geneTree = flashClust(as.dist(dissTOM), method="average")
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)
#This sets the minimum number of genes to cluster into a module
minModuleSize = 30
dynamicMods = cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize = minModuleSize)
dynamicColors= labels2colors(dynamicMods)
MEList= moduleEigengenes(datExpr, colors= dynamicColors,softPower = 18)
MEs= MEList$eigengenes
MEDiss= 1-cor(MEs)
METree= flashClust(as.dist(MEDiss), method= "average")
#save(dynamicMods, MEList, MEs, MEDiss, METree, file= "Network_allSamples_signed_RLDfiltered.RData")
 


 
#plots tree showing how the eigengenes cluster together
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")
#set a threhold for merging modules. In this example we are not merging so MEDissThres=0.0
#MEDissThres = 0.0


MEDissThres = 0.2
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
 
 
#plot dendrogram with module colors below it
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

write.table(MEs,opt$mes,sep="\t") 
 
#save(MEs, moduleLabels, moduleColors, geneTree, file= "Network_allSamples_signed_nomerge_RLDfiltered.RData")

# Correlate traits --------------------------------------------------------
 
#Define number of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use= "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
 
 
#Print correlation heatmap between modules and traits
textMatrix= paste(signif(moduleTraitCor, 2), "\n(",
                        signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix)= dim(moduleTraitCor)
par(mar= c(6, 8.5, 3, 3))
 
 
#display the corelation values with a heatmap plot
labeledHeatmap(Matrix= moduleTraitCor,
            xLabels= names(datTraits),
            yLabels= names(MEs),
            ySymbols= names(MEs),
            colorLabels= FALSE,
            colors= blueWhiteRed(50),
            textMatrix= textMatrix,
            setStdMargins= FALSE,
            cex.text= 0.5,
            zlim= c(-1,1),
            main= paste("Module-trait relationships"))

write.table(moduleColors,opt$modulecolors,sep="\t")
TOM = TOMsimilarityFromExpr(datExpr, power = 6);
# Select module
module = "brown";
# Select module probes
probes = names(datExpr)
inModule = (moduleColors==module);
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
file = paste(opt$visant, module, ".txt", sep=""),
weighted = TRUE,
threshold = 0)

write.table(probes,opt$probes,sep="\t")

dev.off()
