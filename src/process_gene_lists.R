setwd('~/d/sci/src/drug_target_lof')
options(stringsAsFactors=F)
library(sqldf)


# HGNC gene symbol mapping table for updating outdated gene symbols
hgnc_lookup = read.table('../drug_target_lof/data/annotation/hgnc_symbol_lookup_table_2018_09_13.tsv',sep='\t',header=T,quote='',comment.char='')

update_symbols = function(old_symbols, remove_na = F) {
  new_symbols = old_symbols
  needs_updating = !(old_symbols %in% hgnc_lookup$new_symbol)
  for (i in which(needs_updating)) {
    updated_symbol = hgnc_lookup$new_symbol[hgnc_lookup$old_symbol==old_symbols[i]]
    if (is.na(updated_symbol)) {
      new_symbols[i] = NA
    } else {
      new_symbols[i] = updated_symbol[i]
    }
  }
  if (remove_na) {
    return (new_symbols[!is.na(new_symbols)])
  } else {
    return(new_symbols)
  }
}


atc = read.table('data/annotation/atc_summary.tsv',sep='\t',header=T)
dg = read.table('data/drugbank/drug_gene_match.tsv',sep='\t',header=T)
dcats = read.table('data/drugbank/drug_categories.tsv',sep='\t',header=T)
dclass = read.table('data/drugbank/drug_classifications.tsv',sep='\t',header=T)

dclass$is_ab = dclass$drug %in% dcats$drug[dcats$category=='Antibodies']
dclass$atc1 = substr(dclass$atc,1,1)

modalities = sqldf("
select   dg.gene, dclass.type, dclass.is_ab
from     dg, dclass
where    dg.drug = dclass.drug
and      dg.gene != ''
group by 1, 2, 3
order by 1, 2, 3
;")

drug_mod_sm = modalities$gene[modalities$type=='small molecule']
drug_mod_ab = modalities$gene[modalities$type=='biotech' & modalities$is_ab]
drug_mod_oth = modalities$gene[modalities$type=='biotech' & !modalities$is_ab]

write.table(drug_mod_sm,'lists/drug_mod_sm.tsv',sep='\t',row.names=F,col.names=F,quote=F)
write.table(drug_mod_ab,'lists/drug_mod_ab.tsv',sep='\t',row.names=F,col.names=F,quote=F)
write.table(drug_mod_oth,'lists/drug_mod_oth.tsv',sep='\t',row.names=F,col.names=F,quote=F)

indications = sqldf("
select   atc.disease_area, dg.gene
from     dg, dclass, atc
where    dg.drug = dclass.drug
and      dclass.atc1 = atc.level1
and      dg.gene != ''
group by 1, 2
order by 1, 2
;")

for (dz in unique(indications$disease_area)) {
  dz_name = gsub(' ','_',dz)
  genes = indications$gene[indications$disease_area==dz]
  write.table(genes,paste('lists/drug_indic_',dz_name,'.tsv',sep=''),sep='\t',row.names=F,col.names=F,quote=F)
}



#### OMIM genes

omim_morbidmap = read.table('../drug_target_lof_full/raw_data/omim/morbidmap.txt',
                            sep='\t',header=F,skip=4,quote='',comment.char='#',
                            col.names=c('phenotype','genes','mim_number','cyto_location'))

omim_unfiltered = sort(unique(update_symbols(unlist(strsplit(omim_morbidmap$genes, split=', ')), remove_na=T)))
length(omim_unfiltered)

# remove questionable and bracketed associations
omim_morbidmap$questionable = grepl("^[\\?\\{\\[]", omim_morbidmap$phenotype)
# remove drug response associations
omim_morbidmap$response = grepl('[Rr]esponse to',omim_morbidmap$phenotype)
# remove somatic associations
omim_morbidmap$somatic = grepl('[Ss]omatic',omim_morbidmap$phenotype)
# the mim_number column sometimes contains the phenotype mim and sometimes the gene mim
# phenotypes without a mim number assigned tend to not be real Mendelian diseases
# need to extract specifically the phenotype mim from the phenotype field and filter on that
omim_morbidmap$phenotype_mim = as.integer(substring(omim_morbidmap$phenotype, regexpr("[0-9]{6}", omim_morbidmap$phenotype), regexpr("[0-9]{6}", omim_morbidmap$phenotype) + 6))

omim_morbidmap_valid = subset(omim_morbidmap, !questionable & !response & !somatic & !is.na(phenotype_mim))

# check that this has removed all the "susceptibility" things
sum(grepl('susceptibility',omim_morbidmap$phenotype))
sum(grepl('susceptibility',omim_morbidmap_valid$phenotype))

omim_genes = sort(unique(update_symbols(unlist(strsplit(omim_morbidmap_valid$genes, split=', ')), remove_na=T)))
length(omim_genes)

write.table(omim_genes,'lists/omim_genes.tsv',sep='\t',row.names=F,col.names=F,quote=F)


##### mouse het lethal genes

hm = read.table('../drug_target_lof_full/raw_data/HMD_HumanPhenotype.txt',sep='\t',header=F)
colnames(hm)[1] = 'hu'
colnames(hm)[5] = 'mo'
mo_het = read.table('../drug_target_lof_full/raw_data/list_mouse_het_lethal_genes.txt',sep='\t',header=F)
mo_het_lethal = hm$hu[hm$mo %in% mo_het$V1]

write.table(mo_het_lethal,'lists/mo_het_lethal.tsv',sep='\t',row.names=F,col.names=F,quote=F)
