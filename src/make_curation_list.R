options(stringsAsFactors=F)
setwd('~/d/sci/src/drug_target_lof')

htt = read.table('raw_data/HTT-gnomAD-v2.1-ENSG00000197386-2018-11-14-19.51.12.csv',sep=',',header=T,quote='"',comment.char='',colClasses='character')
mapt = read.table('raw_data/MAPT-gnomAD-v2.1-ENSG00000186868-2018-11-15-15.38.03.csv',sep=',',header=T,quote='"',comment.char='',colClasses='character')
snca = read.table('raw_data/SNCA-gnomAD-v2.1-ENSG00000145335-2018-11-15-15.38.32.csv',sep=',',header=T,quote='"',comment.char='',colClasses='character')
sod1 = read.table('raw_data/SOD1-gnomAD-v2.1-ENSG00000142168-2018-11-15-15.39.09.csv',sep=',',header=T,quote='"',comment.char='',colClasses='character')

allvars = rbind(htt,mapt,snca,sod1)

colnames(allvars)[1:5] = c('chrom','pos','rsid','ref','alt')
colnames(allvars)[13] = 'flags'

# set colClasses to 'character' in read.table to avoid alt column only containing T being read as TRUE (logical) rather than 'T' character
# but then that means you need to explicitly cast pos to integer in order for formatC to then work below.
allvars$pos = as.integer(allvars$pos)

allvars$pos_id = paste(allvars$chrom, formatC(allvars$pos, width=9, flag='0'), allvars$ref, allvars$alt, sep='-')

allvars$nonpadded_pos_id = paste(allvars$chrom, allvars$pos, allvars$ref, allvars$alt, sep='-')

write.table(allvars[allvars$flags=='',c('nonpadded_pos_id')],'raw_data/erics-curation-list-2018-12-03.tsv',sep='\t',row.names=F,col.names=F,quote=F)
