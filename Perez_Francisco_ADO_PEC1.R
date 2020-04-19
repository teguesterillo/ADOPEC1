library(knitr)
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE, 
                      comment = NA, prompt = TRUE, tidy = FALSE, 
                      fig.width = 7, fig.height = 7, fig_caption = TRUE,
                      cache=FALSE)
Sys.setlocale("LC_TIME", "C")

targets <- read.csv2("./data/targets.csv", header = TRUE, sep = ";") 
knitr::kable(
  targets, booktabs = TRUE,
  caption = 'Content of the targets file used for the current analysis')

## ----ReadCELfiles, message=FALSE, results='hide', warning=FALSE------------------------------------
library(oligo)
celFiles <- list.celfiles("./data", full.names = TRUE)
library(Biobase)
my.targets <-read.AnnotatedDataFrame(file.path("./data","targets.csv"), 
                                     header = TRUE, row.names = 1, 
                                     sep=";") 
rawData <- read.celfiles(celFiles, phenoData = my.targets)

## ----ChangeName------------------------------------------------------------------------------------
my.targets@data$ShortName->rownames(pData(rawData))
colnames(rawData) <-rownames(pData(rawData)) 

head(rawData)

## ----QCRaw, message=FALSE, warning=FALSE, eval=FALSE-----------------------------------------------
library(arrayQualityMetrics)
arrayQualityMetrics(rawData)

## Crear PCA3--------------------------------------------------------------------------------------------------
library(ggplot2)
library(ggrepel)
plotPCA3 <- function (datos, labels, factor, title, scale,colores, size = 1.5, glineas = 0.25) {
  data <- prcomp(t(datos),scale=scale)
  # plot adjustments
  dataDf <- data.frame(data$x)
  Group <- factor
  loads <- round(data$sdev^2/sum(data$sdev^2)*100,1)
  # main plot
  p1 <- ggplot(dataDf,aes(x=PC1, y=PC2)) +
    theme_classic() +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_vline(xintercept = 0, color = "gray70") +
    geom_point(aes(color = Group), alpha = 0.55, size = 3) +
    coord_cartesian(xlim = c(min(data$x[,1])-5,max(data$x[,1])+5)) +
    scale_fill_discrete(name = "Group")
  # avoiding labels superposition
  p1 + geom_text_repel(aes(y = PC2 + 0.25, label = labels),segment.size = 0.25, size = size) + 
    labs(x = c(paste("PC1",loads[1],"%")),y=c(paste("PC2",loads[2],"%"))) +  
    ggtitle(paste("PCA de los",title,sep=" "))+ 
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values=colores)
}


## ----PCARaw, message=FALSE, fig.cap="Visualization of the two first Principal Components for raw data"----
plotPCA3(exprs(rawData), labels = targets$ShortName, factor = targets$Group, 
         title="datos crudos", scale = FALSE, size = 3, 
         colores = c("red", "blue", "green", "yellow"))

## ----savePCAraw, echo=TRUE, results='hide'---------------------------------------------------------
tiff("figures/PCA_RawData.tiff", res = 200, width = 4.5, height = 4, units = 'in')
plotPCA3(exprs(rawData), labels = targets$ShortName, factor = targets$Group, 
         title="datos crudos", scale = FALSE, size = 2, 
         colores = c("red", "blue", "green", "yellow"))
dev.off()

#el boxplot quedó raro

tmp1 = exprs(rawData)
tmp1 = as.data.frame(tmp1)
boxplot(tmp1, range=0, las=2, cex.axis=0.7,  col = c(rep("red", 3), rep("blue", 3), rep("green", 3), rep("yellow", 3)))


## ----saveIntensRaw, echo=FALSE, results='hide'-----------------------------------------------------
tiff("figures/Intensity_RawData.tiff", res = 200, width = 4, height = 4, units = 'in')
boxplot(tmp1, range=0, las=2, cex.axis=0.7,  col = c(rep("red", 3), rep("blue", 3), rep("green", 3), rep("yellow", 3)), main="Boxplot de los datos crudos")
dev.off()

## ----Normalization---------------------------------------------------------------------------------
eset_rma <- rma(rawData)

## ----QCNorm, message=FALSE, warning=FALSE, eval=FALSE----------------------------------------------
library(arrayQualityMetrics)
arrayQualityMetrics(eset_rma, outdir = file.path("./results", "QCDir.Norm"), force=TRUE)

## ----PCANorm, message=FALSE, fig.cap="Visualization of first two principal components for normalized data"----
library(ggplot2)
library(ggrepel)
plotPCA3(exprs(eset_rma), labels = targets$ShortName, factor = targets$Group, 
         title="datos normalizados", scale = FALSE, size = 3, 
         colores = c("red", "blue", "green", "yellow"))

## ----savePCAnorm, echo=FALSE, results='hide'-------------------------------------------------------
tiff("figures/PCA_NormData.tiff", res = 150, width = 5, height = 5, units = 'in')
plotPCA3(exprs(eset_rma), labels = targets$ShortName, factor = targets$Group, 
         title="datos normalizados", scale = FALSE, size = 2, 
         colores = c("red", "blue", "green", "yellow"))
dev.off()

tmp = exprs(eset_rma)
tmp = as.data.frame(tmp)
boxplot(tmp, range=0, las=2)

tiff("figures/boxplot_NormData.tiff", res = 150, width = 5, height = 5, units = 'in')
boxplot(tmp, range=0, las=2, cex.axis=0.7,  col = c(rep("red", 3), rep("blue", 3), rep("green", 3), rep("yellow", 3)), main="Boxplot de los datos normalizados")
dev.off()

## ----SDplot, fig.cap="Values of standard deviations allong all samples for all genes ordered from smallest to biggest"----
sds <- apply (exprs(eset_rma), 1, sd)
sdsO<- sort(sds)
plot(1:length(sdsO), sdsO, main="Distribution of variability for all genes",
     sub="Vertical lines represent 90% and 95% percentiles",
     xlab="Gene index (from least to most variable)", ylab="Standard deviation")
abline(v=length(sds)*c(0.9,0.95))

## ----saveSDplot, echo=FALSE, results='hide'--------------------------------------------------------
tiff("figures/SDplot.tiff", res = 150, width = 5, height = 5, units = 'in')
plot(1:length(sdsO), sdsO, main="Distribution of variability for all genes",
     sub="Vertical lines represent 90% and 95% percentiles",
     xlab="Gene index (from least to most variable)", ylab="Standard deviation")
abline(v=length(sds)*c(0.9,0.95))
dev.off()

## ----Filtering1, results='hide', message=FALSE-----------------------------------------------------
library(genefilter)
library(mogene21sttranscriptcluster.db)
annotation(eset_rma) <- "mogene21sttranscriptcluster.db"
filtered <- nsFilter(eset_rma, 
                     require.entrez = TRUE, remove.dupEntrez = TRUE,
                     var.filter=TRUE, var.func=IQR, var.cutoff=0.75, 
                     filterByQuantile=TRUE, feature.exclude = "^AFFX")


## ----FilterResults1, results='hide', echo=FALSE----------------------------------------------------
names(filtered)
class(filtered$eset)


## ----FilterResults2--------------------------------------------------------------------------------
print(filtered$filter.log)
eset_filtered <-filtered$eset

## ----SaveData1, results='hide', message=FALSE------------------------------------------------------
write.csv(exprs(eset_rma), file="./results/normalized.Data.csv")
write.csv(exprs(eset_filtered), file="./results/normalized.Filtered.Data.csv")
save(eset_rma, eset_filtered, file="./results/normalized.Data.Rda")

## ----LoadSavedData---------------------------------------------------------------------------------
if (!exists("eset_filtered")) load (file="./results/normalized.Data.Rda")


## ----DesignMatrix, message=FALSE-------------------------------------------------------------------
library(limma)
designMat<- model.matrix(~0+Group, pData(eset_filtered))
colnames(designMat) <- c("attenuated", "none", "saprophyte", "virulent")
print(designMat)

## ----setContrasts----------------------------------------------------------------------------------
cont.matrix <- makeContrasts (VIRULENTvsNONE = virulent-none,
                              SAPROPHYTEvsNONE = saprophyte-none,
                              ATTENUATEDvsNONE = attenuated-none,
                              levels=designMat)
print(cont.matrix)

## ---- linearmodelfit-------------------------------------------------------------------------------
library(limma)
fit<-lmFit(eset_filtered, designMat)
fit.main<-contrasts.fit(fit, cont.matrix)
fit.main<-eBayes(fit.main)
class(fit.main)

## ---- topTabs1-------------------------------------------------------------------------------------
topTab_VIRULENTvsNONE <- topTable (fit.main, number=nrow(fit.main), coef="VIRULENTvsNONE", adjust="fdr") 
head(topTab_VIRULENTvsNONE)

## ---- topTabs2-------------------------------------------------------------------------------------
topTab_SAPROPHYTEvsNONE <- topTable (fit.main, number=nrow(fit.main), coef="SAPROPHYTEvsNONE", adjust="fdr") 
head(topTab_SAPROPHYTEvsNONE)

## ---- topTabs3-------------------------------------------------------------------------------------
topTab_ATTENUATEDvsNONE <- topTable (fit.main, number=nrow(fit.main), coef="ATTENUATEDvsNONE", adjust="fdr") 
head(topTab_ATTENUATEDvsNONE)

write.csv(topTab_ATTENUATEDvsNONE,"./attenuated.csv", row.names = TRUE)

write.csv(topTab_VIRULENTvsNONE,"./virulent.csv", row.names = TRUE)

write.csv(topTab_SAPROPHYTEvsNONE,"./saprophyte.csv", row.names = TRUE)

## ----GeneAnnotation, message=FALSE, warning=FALSE--------------------------------------------------
annotatedTopTable <- function(topTab, anotPackage)
{
  topTab <- cbind(PROBEID=rownames(topTab), topTab)
  myProbes <- rownames(topTab)
  thePackage <- eval(parse(text = anotPackage))
  geneAnots <- select(thePackage, myProbes, c("SYMBOL", "ENTREZID", "GENENAME"))
  annotatedTopTab<- merge(x=geneAnots, y=topTab, by.x="PROBEID", by.y="PROBEID")
  return(annotatedTopTab)
}

## ----annotateTopTables-----------------------------------------------------------------------------
library(mogene21sttranscriptcluster.db)
topAnnotated_ATTENUATEDvsNONE <- annotatedTopTable(topTab_ATTENUATEDvsNONE,
                                              anotPackage="mogene21sttranscriptcluster.db")
topAnnotated_SAPROPHYTEvsNONE <- annotatedTopTable(topTab_SAPROPHYTEvsNONE,
                                            anotPackage="mogene21sttranscriptcluster.db")
topAnnotated_VIRULENTvsNONE <- annotatedTopTable(topTab_VIRULENTvsNONE,
                                      anotPackage="mogene21sttranscriptcluster.db")
write.csv(topAnnotated_ATTENUATEDvsNONE, file="./results/topAnnotated_ATTENUATEDvsNONE.csv")
write.csv(topAnnotated_SAPROPHYTEvsNONE, file="./results/topAnnotated_SAPROPHYTEvsNONE.csv")
write.csv(topAnnotated_VIRULENTvsNONE, file="./results/topAnnotated_VIRULENTvsNONE.csv")

## ----decideTests.1---------------------------------------------------------------------------------
library(limma)
res<-decideTests(fit.main, method="separate", adjust.method="fdr", p.value=0.1, lfc=1)


## ----resumeDecideTests-----------------------------------------------------------------------------
sum.res.rows<-apply(abs(res),1,sum)
res.selected<-res[sum.res.rows!=0,] 
print(summary(res))


## ---- vennDiagram, fig.cap="Venn diagram showing the genes in common between the three comparisons performed"----
vennDiagram (res.selected[,1:3], cex=0.9)
title("Genes en común entre las tres comparaciones \n Genes seleccionados con FDR < 0.1 y logFC > 1")


## ----vennPlot, echo=FALSE, results='hide'----------------------------------------------------------
tiff("figures/VennPlot.tiff", res = 150, width = 6.5, height = 5.5, units = 'in')
vennDiagram (res.selected[,1:3], cex=0.9)
title("Genes en común entre las tres comparaciones \n Genes seleccionados con FDR < 0.1 y logFC > 1")
dev.off()

## ----selectGenes-----------------------------------------------------------------------------------
listOfTables <- list(VIRULENTvsNONE = topTab_VIRULENTvsNONE, 
                     SAPROPHYTEvsNONE  = topTab_SAPROPHYTEvsNONE, 
                     ATTENUATEDvsNONE = topTab_ATTENUATEDvsNONE)
listOfSelected <- list()
for (i in 1:length(listOfTables)){
  # select the toptable
  topTab <- listOfTables[[i]]
  # select the genes to be included in the analysis
  whichGenes<-topTab["adj.P.Val"]<0.0000005
  selectedIDs <- rownames(topTab)[whichGenes]
  # convert the ID to Entrez
  EntrezIDs<- select(mogene21sttranscriptcluster.db, selectedIDs, c("ENTREZID"))
  EntrezIDs <- EntrezIDs$ENTREZID
  listOfSelected[[i]] <- EntrezIDs
  names(listOfSelected)[i] <- names(listOfTables)[i]
}
sapply(listOfSelected, length)


## --------------------------------------------------------------------------------------------------
mapped_genes2GO <- mappedkeys(org.Mm.egGO)
mapped_genes2KEGG <- mappedkeys(org.Mm.egPATH)
mapped_genes <- union(mapped_genes2GO , mapped_genes2KEGG)

## ----BiologicalSig---------------------------------------------------------------------------------
library(ReactomePA)

listOfData <- listOfSelected[1:3]
comparisonsNames <- names(listOfData)
universe <- mapped_genes

for (i in 1:length(listOfData)){
  genesIn <- listOfData[[i]]
  comparison <- comparisonsNames[i]
  enrich.result <- enrichPathway(gene = genesIn,
                                 pvalueCutoff = 0.05,
                                 readable = T,
                                 pAdjustMethod = "BH",
                                 organism = "mouse",
                                 universe = universe)
  
  cat("##################################")
  cat("\nComparison: ", comparison,"\n")
  print(head(enrich.result))
  
  if (length(rownames(enrich.result@result)) != 0) {
    write.csv(as.data.frame(enrich.result), 
              file =paste0("./results/","ReactomePA.Results.",comparison,".csv"), 
              row.names = FALSE)
    
    pdf(file=paste0("./results/","ReactomePABarplot.",comparison,".pdf"))
    print(barplot(enrich.result, showCategory = 15, font.size = 4, 
                  title = paste0("Reactome Pathway Analysis for ", comparison,". Barplot")))
    dev.off()
    
    pdf(file = paste0("./results/","ReactomePAcnetplot.",comparison,".pdf"))
    print(cnetplot(enrich.result, categorySize = "geneNum", schowCategory = 15, 
                   vertex.label.cex = 0.75))
    dev.off()
  }
}

