out_prefix<-"A_M_Br_Sp"
out_analysis <- "./analysis/"
out_figure <- "./figure/"

## load rna-seq count file and design file
rna.file = "./data/All_Tcell_and_Spleen_GeneCount.txt"
design.file = "./data/design_file.txt"

data.design = read.table(design.file,  sep="\t", head=T, quote="", check.names=F)
data.design$Sample_Name <- gsub(pattern="-", x = data.design$Sample_Name, replacement = ".")

data.raw = read.table(rna.file,  sep="\t", head=T, quote="", check.names=F, row.names=1)
colnames(data.raw)= data.design$Sample_Name

# get rid of YoungMale.Brain and AgedFemale.Brain
library(dplyr)
data.raw <- data.raw %>% select(-contains("YoungMale"), -contains("AgedFemale"))
data.design <- data.design %>% filter(TERM_1=="Aged", TERM_2=="Male") 
all(data.design$Sample_Name == colnames(data.raw)) #check the sample names matches in filtered count files and the design file 


## Idenitfy the number of model terms.
if("TERM_4" %in% colnames(data.design)){
  #Fit a GLM since there are multiple factors.
  #Exact test is only for "1" factor models!
  term.n  = 4
  #Specify the factors and its levels.  
  group = factor(paste(data.design$TERM_1, data.design$TERM_2, data.design$TERM_3, data.design$TERM_4, sep="."))
  #Append group (factor levels) to the design table.
  data.design = cbind(data.design, Group=group)
}else if("TERM_3" %in% colnames(data.design)){
  term.n  = 3    
  group = factor(paste(data.design$TERM_1, data.design$TERM_2, data.design$TERM_3, sep="."))
  data.design = cbind(data.design, Group=group)
}else if("TERM_2" %in% colnames(data.design)){
  term.n  = 2
  group = factor(paste(data.design$TERM_1, data.design$TERM_2, sep="."))
  data.design = cbind(data.design, Group=group)
}else{
  #Run exact test since only one factor.
  term.n  = 1
  group = factor(data.design$TERM_1)
  data.design = cbind(data.design, Group=group)
}

## DGEList, message=FALSE
library(edgeR)
y <- DGEList(data.raw, group=group, genes=row.names(data.raw)) # must specify
options(digits=3)
y$samples

## symbols, message=FALSE
library(org.Mm.eg.db)
y$genes$Symbol <- mapIds(org.Mm.eg.db, rownames(y),
                         keytype="ENSEMBL", column="SYMBOL")  # keytype="ENSEMBL", attach gene Symbol from database to DGElist y 
y$genes$Entrezid <- mapIds(org.Mm.eg.db, rownames(y),
                           keytype="ENSEMBL", column="ENTREZID")
y$genes$Genename <- mapIds(org.Mm.eg.db, rownames(y),
                           keytype="ENSEMBL", column="GENENAME")

### all above is saved for Rdt file for Rmd

### select(org.Hs.eg.db, keys=rownames(resultTable), columns=c("SYMBOL","GENENAME"), keytype="ENTREZID")
head(y$genes)

## dropNAsymbols
y <- y[!is.na(y$genes$Symbol),]
dim(y)

## keep
keep <- rowSums(cpm(y) > 1) >= 2
table(keep)

## ----filter--------------------------------------------------------------
y <- y[keep, , keep.lib.sizes=FALSE]

## ----norm----------------------------------------------------------------
y <- calcNormFactors(y)
y$samples
dim(y)



## ----mdsplot # what's difference to plot log.cpm or the whole object y? 
par(xpd = T, mar = par()$mar + c(0,0,0,7))  # to make legends outside of the plot
plotMDS(y, top = 500, cex = 1, pch=9, dim.plot = c(1,2), ndim = 2, gene.selection = "pairwise", col=as.numeric(group)) # col = as.numeric(group) differentiate colors between groups
legend(2.5, 0.025,levels(group), pch=rep(9, length(group)), col=levels(as.factor(as.numeric(group))))
par(mar=c(5, 4, 4, 7) + 0.1)

## ----mdplot1, fig.cap="MD plot of log2-expression in sample 1 versus the average log2-expression across all other samples. Each point represents a gene, and the red line indicates a log-ratio of zero. The majority of points cluster around the red line."----
plotMD(y, column=1)
abline(h=0, col="red", lty=2, lwd=2)

## ----mdplot11, fig.cap="MD plot of log2-expression in sample 11 versus the average log2-expression across all other samples. The plot shows a number of genes that are both highly expressed and highly up-regulated."----
plotMD(y, column=11)
abline(h=0, col="red", lty=2, lwd=2)

## ----design--------------------------------------------------------------
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
design

## ----estimateDisp--------------------------------------------------------
y <- estimateDisp(y, design, robust=TRUE)

## ----plotBCV, width="3.8in", fig.cap="Scatterplot of the biological coefficient of variation (BCV) against the average abundance of each gene. The plot shows the square-root estimates of the common, trended and tagwise NB dispersions."----
plotBCV(y)

## ----glmQLFit------------------------------------------------------------
fit <- glmQLFit(y, design, robust=TRUE)
head(fit$coefficients)

## ----QLDisp, out.width="3.8in", fig.cap="A plot of the quarter-root QL dispersion against the average abundance of each gene. Estimates are shown for the raw (before EB moderation), trended and squeezed (after EB moderation) dispersions. Note that the QL dispersions and trend shown here are relative to the NB dispersion trend shown in Figure~\ref{fig:plotBCV}."----
plotQLDisp(fit)

## ----df.prior------------------------------------------------------------
summary(fit$df.prior)

## ----MakeContrasts--------------------------------------------------------------
con <- makeContrasts(
  a1_BrvsSp.A.M.NS = Aged.Male.Brain.NoStim - Aged.Male.Spleen.NoStim,
  a2_BrvsSp.A.M.PS = Aged.Male.Brain.PlusStim - Aged.Male.Spleen.PlusStim,
  a3_PSvsNS.A.M.Br = Aged.Male.Brain.PlusStim - Aged.Male.Brain.NoStim,
  a4_PSvsNS.A.M.Sp = Aged.Male.Spleen.PlusStim - Aged.Male.Spleen.NoStim,
  a5_Stim_Tissue_Int.A.M = (Aged.Male.Brain.PlusStim - Aged.Male.Brain.NoStim) - (Aged.Male.Spleen.PlusStim - Aged.Male.Spleen.NoStim),
  levels=design
)

## ----glmQLFTest----------------------------------------------------------
res <- glmQLFTest(fit, contrast=con[,1]) # later comfirmed by PCA that tissue difference (contrast a1 or a2) explain PC1 

## ----topTags-------------------------------------------------------------
topTags(res)

## ----decideTests---------------------------------------------------------
is.de <- decideTestsDGE(res)
summary(is.de)

## ----plotMDfit, fig.cap="MD plot showing the log-fold change and average abundance of each gene. Significantly up and down DE genes are highlighted in red and blue, respectively."----
plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"),
       legend="topright")

## ----treat---------------------------------------------------------------
tr <- glmTreat(fit, contrast=con[,1], lfc=log2(2))
topTags(tr)

## ----treatdecideTests----------------------------------------------------
is.de <- decideTestsDGE(tr)
summary(is.de)

## ----plotMDtreat, fig.cap="MD plot showing the log-fold change and average abundance of each gene. Genes with fold-changes significantly greater than 1.5 are highlighted."----
plotMD(tr, status=is.de, values=c(1,-1), col=c("red","blue"),
       legend="topright")

## ----cpm-----------------------------------------------------------------
logCPM <- cpm(y, prior.count=2, log=TRUE)
logCPM.PCA<-logCPM # save it later for PCA plot
rownames(logCPM) <- y$genes$Symbol
#colnames(logCPM) <- paste(y$samples$group, 1:2, sep="-")
colnames(logCPM) <- data.design$Sample_Name # get it into the y project

## ----order---------------------------------------------------------------
o <- order(tr$table$PValue)
logCPM <- logCPM[o[1:100],]

## ----scale---------------------------------------------------------------
logCPM <- t(scale(t(logCPM)))

## ----heatmap, message=FALSE, fig.width=8, fig.height=12, fig.cap="Heat map across all the samples using the top 100 most DE genes in contrast 1"----
library(gplots)
col.pan <- colorpanel(100, "blue", "white", "red")
heatmap.2(logCPM, col=col.pan, Rowv=TRUE, scale="none", 
          trace="none", dendrogram="both", cexRow=0.5, cexCol=0.7, density.info="none",
          margin=c(10,9), lhei=c(2,10), lwid=c(2,6))

## ---Pincicpal component analysis ---
pca_original = prcomp(t(logCPM.PCA),scale=T, center=T)
pca_x <- pca_original$x
pca_table <- data.frame(data.design, pca_x)
x <- pca_original$sdev^2/sum(pca_original$sdev^2) # Proportion of Variance Explained for all components

## Scree plot
plot(x, xlab="Principal Component", ylab="Proportion of Variance Explained", type="b")
plot(cumsum(x), xlab="Principal Component", ylab="Cumulative Proportion of Variance Explained", type="b")

library(ggplot2)
PCA_plot <- function(pca_table, PC_x, PC_y, color, shape){
  #PC_x,PC_y are type of interger
  #color, shape, are type of string
  g <- ggplot(pca_table, aes_string(x=names(pca_table[PC_x]), y=names(pca_table[PC_y]), color=color, shape=shape)) 
  g <- g + geom_point(alpha=0.7, size=3) 
  # g <- g + labs(color = "Group", shape="Tissue")
  g + labs(x = paste(names(pca_table[PC_x]), scales::percent(x[PC_x]),"variance explained", sep=" "), y=paste(names(pca_table[PC_y]), scales::percent(x[PC_y]),"variance explained", sep=" "))
  #filename <- paste()
  #ggsave(filename, width=7, height=7, units="in")
}

PCA_plot(pca_table, 1, 2, "Group", "TERM_3")
PCA_plot(pca_table, 1, 2, "Group", "TERM_4")
PCA_plot(pca_table, 2, 3, "Group", "TERM_3")
PCA_plot(pca_table, 2, 3, "Group", "TERM_4")
PCA_plot(pca_table, 2, 3, "Group", "TERM_2")


print(con)
as.list(seq(ncol(con)))%>% walk(my_fun,fit=fit, Treat_FC=2)