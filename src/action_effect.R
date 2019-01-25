# This script is to help me review the drugbank data to decide which drug 
# "actions" count as "positive", "negative", or "unknown/other"
# and then to write out the corresponding gene lists

options(stringsAsFactors = F)
setwd('~/d/sci/src/drug_target_lof/')
library(sqldf)
library(reshape2)
library(venneuler)
drug_gene_action = read.table('data/drug_gene_action_match.tsv',sep='\t',header=T,quote='',comment.char='',as.is=TRUE)
drug_gene_action = drug_gene_action[drug_gene_action$gene != '',]
lookup = read.table('data/hgnc_symbol_lookup_table_2018_09_13.tsv',sep='\t',header=T,quote='',comment.char='')
drug_gene_action$new_symbol = lookup$new_symbol[match(drug_gene_action$gene, lookup$old_symbol)]
drug_gene_action$gene = drug_gene_action$new_symbol
drug_gene_action = drug_gene_action[,c('drug','gene','action')]
gene_action = unique(drug_gene_action[,c('gene','action')])
gene_action = gene_action[gene_action$gene != '',]
# see what-all actions there are:
actions = unique(gene_action$action)
# pulled this into a text editor and added an 'effect' column to group these
# now read back in:
action_effect = read.table('data/action_effect.tsv',sep='\t',header=T)
gene_action$effect = action_effect$effect[match(gene_action$action, action_effect$action)]
gene_action[gene_action$action=='antibody',]
# antibodies could in principle have various effects on target so i labeled as 
# 'other/unknown'. from curating a few it appears most are negative,
# e.g. ERBB2, C5, TNF, ITGAL, CD19, but for example for CD19 the result is
# destruction of cells that express CD19, so a bit different
gene_action$effect[gene_action$action=='antibody' & gene_action$gene %in% c('ERBB5','C5','TNF','ITGAL')] = 'negative'
gene_action = gene_action[with(gene_action, order(gene)),]
table(gene_action$effect)
multi_effect = sqldf("select gene, count(distinct effect) n_effect from gene_action group by 1 having count(*) > 1 order by 1;")
pos_and_neg = sqldf("select gene, count(distinct effect) n_effect from gene_action where effect != 'other/unknown' group by 1 having count(distinct effect) > 1 order by 1;")

write.table(unique(gene_action$gene[gene_action$effect=='negative']),'lists/negative_targets.tsv',sep='\t',row.names=F,quote=F,col.names=T)
write.table(unique(gene_action$gene[gene_action$effect=='positive']),'lists/positive_targets.tsv',sep='\t',row.names=F,quote=F,col.names=T)
write.table(unique(gene_action$gene[gene_action$effect=='other/unknown']),'lists/other_targets.tsv',sep='\t',row.names=F,quote=F,col.names=T)

            
            
            
