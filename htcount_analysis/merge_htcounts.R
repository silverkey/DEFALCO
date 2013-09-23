#------------------------------------------------------------------------------
#  MERGE HTCOUNTS OUTPUT
#------------------------------------------------------------------------------

get.name = function(name) {
  name = gsub('COUNTS_','',name)
  name = gsub('.bam','',name)
  name
}

get.df = function(filename) {
  t = read.table(file=filename)
  name = get.name(filename)
  colnames(t) = c('geneid',name)
  t
}

cfiles = dir()[grep('COUNTS_',dir())]

counts = get.df(cfiles[1])

for(i in 2:length(cfiles)) {
  t = get.df(cfiles[i])
  counts = merge(counts,t,by='geneid',all.x=T,all.y=T,sort=F)
}

write.table(counts,file='COUNTS.txt',sep="\t",quote=F,row.names=F)
#------------------------------------------------------------------------------


