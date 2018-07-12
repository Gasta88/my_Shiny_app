# Tutorial from http://bioinformatics.knowledgeblog.org/2011/06/20/analysing-microarray-data-in-bioconductor/
source("https://bioconductor.org/biocLite.R")
biocLite("GEOquery")
library(GEOquery)
getGEOSuppFiles("GSE20986")
untar("GSE20986/GSE20986_RAW.tar", exdir="data")
cels <- list.files("data/", pattern = "[gz]")
sapply(paste("data", cels, sep="/"), gunzip)
cels


biocLite("simpleaffy")
library(simpleaffy)
celfiles <- read.affy(covdesc="phenodata.txt", path="data") 
celfiles.gcrma <- gcrma(celfiles)
celfiles.filtered <- nsFilter(celfiles.gcrma, require.entrez=FALSE, remove.dupEntrez=FALSE)
samples <- celfiles.gcrma$Target
samples <- as.factor(samples)
design <- model.matrix(~0 + samples)
colnames(design) <- c("choroid", "huvec", "iris", "retina")
library(limma)
fit <- lmFit(exprs(celfiles.filtered$eset), design)
contrast.matrix <- makeContrasts(huvec_choroid = huvec - choroid, huvec_retina = huvec - retina, huvec_iris <- huvec - iris, levels=design)
huvec_fits <- contrasts.fit(fit, contrast.matrix)
huvec_ebFit <- eBayes(huvec_fits)

topTable(huvec_ebFit, number=10, coef=1)
probeset.list <- topTable(huvec_ebFit, coef=1, number=10000, lfc=4)
# Annotation
biocLite("hgu133plus2.db")
library(hgu133plus2.db)
library(annotate)
gene.symbols <- getSYMBOL(rownames(probeset.list), "hgu133plus2") 
results <- cbind(probeset.list, gene.symbols)

#write.csv2(results,file='data/input_data.csv',sep = "\t")
pdattab<- exprs(celfiles.filtered$eset)
colnames(pdattab)<-c("iris-1","retina-1","retina-2","iris-2","retina-3","iris-3","choroid-1","choroid-2","choroid-3","huvec-1","huvec-2","huvec-3")
#write.csv2(pdattab,file='data/pdattab_data.csv',sep = "\t")
mean_iris<-apply(pdattab,1,function(x){mean(x[c(1,4,6)])})
