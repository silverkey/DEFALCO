library(GenomicFeatures)
file = 'Homo_sapiens.GRCh37.72.gtf'
transdb = makeTranscriptDbFromGFF(file,format='gtf',dataSource='ensembl version 72',species='Homo sapiens',exonRankAttributeName='exon_number')
saveDb(transdb,file='transdb_hsapiens_ensembl_72.sqlite')
