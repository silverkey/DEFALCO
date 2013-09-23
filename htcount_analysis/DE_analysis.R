#------------------------------------------------------------------------------
# LIMMA ANALYSIS
#------------------------------------------------------------------------------
library(edgeR)
samples = read.table(file='targets.txt',sep='\t',quote='',comment.char='',head=T,stringsAsFactors=F)
rownames(samples) = samples$name
data = read.table('COUNTS.txt',header=T,quote='')
colnames(data) = gsub('.','-',colnames(data),fixed=T)
rownames(data) = data$geneid
data = data[,samples$name]
tot = colSums(data)
lib.size=tot[samples$name]
data = data[grep('ENS',data$geneid),]
notass = data[-grep('ENS',data$geneid),]

comparison = c('FC','FV')
colData = samples[samples$condition %in% comparison,]
countData = data[,colData$name]
dge = DGEList(countData, group=colData$condition,lib.size=lib.size[colData$name])
m = 1e6 * t(t(dge$counts) / dge$samples$lib.size)
ridx = rowSums(m > 1) >= 3
table(ridx)
dge2 = dge[ridx,]
dge2$samples$lib.size <- colSums(dge2$counts)
dge2 = calcNormFactors(dge2)
design = model.matrix(~0+colData$condition)
colnames(design) = gsub('colData$condition','',colnames(design),fixed=T)
v = voom(dge2,design,plot=TRUE)
fit = lmFit(v,design)
contrast = makeContrasts(paste(comparison,sep='',collapse='-'),levels=design)
fit2 = contrasts.fit(fit,contrast)
fit2 = eBayes(fit2)
options(digits=3)
res = topTable(fit2,n=nrow(fit2))
res[res$adj.P.Val<=0.1 & abs(res$logFC)>=log2(1.5),]
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# EDGER ANALYSIS
#------------------------------------------------------------------------------
library(edgeR)
samples = read.table(file='targets.txt',sep='\t',quote='',comment.char='',head=T,stringsAsFactors=F)
rownames(samples) = samples$name
data = read.table('COUNTS.txt',header=T,quote='')
colnames(data) = gsub('.','-',colnames(data),fixed=T)
rownames(data) = data$geneid
data = data[,samples$name]
tot = colSums(data)
lib.size=tot[samples$name]
data = data[grep('ENS',data$geneid),]
notass = data[-grep('ENS',data$geneid),]

comparison = c('FC','FV')
colData = samples[samples$condition %in% comparison,]
countData = data[,colData$name]
dge = DGEList(countData, group=colData$condition,lib.size=lib.size[colData$name])
m = 1e6 * t(t(dge$counts) / dge$samples$lib.size)
ridx = rowSums(m > 1) >= 3
table(ridx)
dge2 = dge[ridx,]
dge2$samples$lib.size <- colSums(dge2$counts)
dge2 = calcNormFactors(dge2)
dge2 = estimateCommonDisp(dge2)
dge2 = estimateTagwiseDisp(dge2)
diff = exactTest(dge2,pair=comparison)
tt = topTags(diff,n=nrow(diff))
res = tt$table
res[res$FDR<=0.1 & abs(res$logFC)>=log2(1.5),]
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# DESEQ2 ANALYSIS
#------------------------------------------------------------------------------
library(DESeq2)
samples = read.table(file='targets.txt',sep='\t',quote='',comment.char='',head=T,stringsAsFactors=F)
rownames(samples) = samples$name
data = read.table('COUNTS.txt',header=T,quote='')
colnames(data) = gsub('.','-',colnames(data),fixed=T)
counts = data[grep('ENS',data$geneid),]
notass = data[-grep('ENS',data$geneid),]
rownames(counts) = counts$geneid
counts = counts[,samples$name]

comparison = c('FC','FV')
colData = samples[samples$condition %in% comparison,]
countData = counts[,colData$name]
dds = DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
colData(dds)$condition = factor(colData(dds)$condition,levels=comparison)
dds = DESeq(dds)
res = results(dds)
res = res[order(res$padj),]
use <- res$baseMean >= 10 & !is.na(res$pvalue)
resFilt <- res[use,]
resFilt$padj <- p.adjust(resFilt$pvalue, method="BH")
sum(res$padj < .1, na.rm=TRUE)
sum(resFilt$padj < .1, na.rm=TRUE)
res2 = as.data.frame(resFilt)
res2[res2$padj<=0.1 & abs(res2$log2FoldChange)>=log2(1.5),]
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# Samples hierarchical clustering using LIMMA
#------------------------------------------------------------------------------
library(edgeR)
samples = read.table(file='targets.txt',sep='\t',quote='',comment.char='',head=T,stringsAsFactors=F)
rownames(samples) = samples$name
data = read.table('COUNTS.txt',header=T,quote='')
colnames(data) = gsub('.','-',colnames(data),fixed=T)
rownames(data) = data$geneid
data = data[,samples$name]
tot = colSums(data)
lib.size=tot[samples$name]
colData = samples
countData = data[,colData$name]
dge = DGEList(countData, group=colData$condition,lib.size=lib.size[colData$name])
dge = calcNormFactors(dge)
design = model.matrix(~0+colData$condition)
colnames(design) = gsub('colData$condition','',colnames(design),fixed=T)
v = voom(dge,design,plot=TRUE)
exprs = v$E
ecTr = dist(t(exprs), method = "euclidean")
hecTr = hclust(ecTr, method = "average")
plot(hecTr, main = "Hierarchical clustering of the Samples", xlab = "", sub = "Average linkage, Euclidean distance for samples")
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
# Samples PCA clustering using DESeq2
#------------------------------------------------------------------------------
library(DESeq2)
samples = read.table(file='targets.txt',sep='\t',quote='',comment.char='',head=T,stringsAsFactors=F)
rownames(samples) = samples$name
data = read.table('COUNTS.txt',header=T,quote='')
colnames(data) = gsub('.','-',colnames(data),fixed=T)
counts = data[grep('ENS',data$geneid),]
notass = data[-grep('ENS',data$geneid),]
rownames(counts) = counts$geneid
counts = counts[,samples$name]
colData = samples
countData = counts[,colData$name]
dds = DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
colData(dds)$condition = factor(colData(dds)$condition,levels=unique(samples$condition))
dds = DESeq(dds)
rld = rlogTransformation(dds)
print(plotPCA(rld[,1:12]))
print(plotPCA(rld[,13:24]))
#------------------------------------------------------------------------------





