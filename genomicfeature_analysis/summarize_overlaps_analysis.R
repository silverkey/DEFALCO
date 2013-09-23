library(parallel)
library(GenomicFeatures)
library(Rsamtools)

# To get all the names of the sequences from from the BAM
#test = scanBam('HC-1.bam')
#test1 = test[[1]]
# We are not interested in NA and chrM
#chr = as.character(unique(test1$rname))[c(1:22,24,25)]
# Once you got it you can save a table and load it
# on demand
chrfile = 'hs19chr'
trdbfile = 'transdb_hsapiens_ensembl_72.sqlite'
mc.cores = 24

options(mc.cores = mc.cores)
chr = as.character(read.table(chrfile)$V1)
transdb = loadDb(file=trdbfile)
genes = exonsBy(transdb,by='gene')

fls = list.files(pattern='*bam$')
bfl = BamFileList(fls)

# Extract a table to build the maximum lengthy of each chromosome
# by taking the maximum end from the transcripts db
t = transcripts(transdb)
a = cbind(as.character(seqnames(t)),end(t))
chrl = tapply(as.numeric(a[,2]),a[,1],max)
names(chrl) = paste('chr',names(chrl),'.fa',sep='')
chrl[chr]
chrl = as.data.frame(chrl[chr])
colnames(chrl) = 'end'
chrl$start = 1
chrl$seqname = rownames(chrl)
chrl = chrl[,c(3,2,1)]

chrr = GRanges(seqnames=Rle(chrl$seqname),ranges=IRanges(chrl$start,(chrl$end+10)),strand=Rle('*'))
param = ScanBamParam(which=chrr)
seqlevels(genes) = paste('chr',seqlevels(genes),'.fa',sep='')

bfls = BamFileList(fls)
genehits = summarizeOverlaps(genes,bfls,ignore.strand=T,param=param,singleEnd=T)
counts = assays(genehits)$counts

#countbam = function(bam) {
#  cat('Analyzing: ')
#  cat(bam)
#  cat("...\n")
#  genehits = summarizeOverlaps(genes,BamFileList(bam),ignore.strand=T,param=param,singleEnd=T)
#  counts = assays(genehits)$counts
#  counts
#}
#counts = sapply(fls,countbam)

write.table(counts,file='counts.txt',sep="\t",row.names=F,quote=F)

#library(GenomicFeatures)
#file = 'Homo_sapiens.GRCh37.72.gtf'
#makeTranscriptDbFromGFF(file,format='gtf',dataSource='ensembl version 72',species='Homo sapiens')

#if(download == 'T') {
#  transdb = makeTranscriptDbFromUCSC(genome=genome,tablename=tablename)
#  saveFeatures(transdb,file=paste(dbdir,'/',genome,'.',tablename,'.','sqlite',sep=''))
#} else {
#  transdb = loadFeatures(file=paste(dbdir,'/',genome,'.',tablename,'.','sqlite',sep=''))
#}

