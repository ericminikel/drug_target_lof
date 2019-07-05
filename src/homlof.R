# UKBB - DeBoever 2018 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5915386/

# Narasimhan 2016 - see p. 32-33 Table S1 of supplement

options(stringsAsFactors=F)
library(readxl)



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

universe = read.table('lists/universe.tsv',sep='\t',header=F)
homlof_genes = data.frame(gene=universe$V1)

# Sulem 2015 Table S4 Iceland DeCODE data
sulem2015 = read_xlsx('data/homlof/sulem-2015-table-s4.xlsx',skip=3)
colnames(sulem2015) = gsub('[^a-z0-9_]','_',tolower(colnames(sulem2015)))
sulem2015_include = sqldf("
select   gene,
         number_of_compound_heterozygous_carriers_in_a_variant_pair_where_maf_2__for_both_variantsa n_cpdhet,
         observed_number_of_imputed_homozygotes n_hom
from     sulem2015
where    vep_consequence in ('frameshift_variant','splice_acceptor_variant','splice_donor_variant','stop_gained') -- eliminate stop_lost and initiator_codon_variant, which are more missense-like in terms of MAPS
and      sequencing_maf____ < 2 -- rare only - require < 2% MAF
and (    number_of_compound_heterozygous_carriers_in_a_variant_pair_where_maf_2__for_both_variantsa > 0
or       observed_number_of_imputed_homozygotes > 0
)
;")
sulem2015_genes = sort(unique(update_symbols(sulem2015_include$gene, remove_na=T)))
homlof_genes$sulem2015 = homlof_genes$gene %in% sulem2015_genes
sum(homlof_genes$sulem2015)

# Saleheen Table S1 Pakistan PROMIS data
saleheen2017 = read_xlsx('data/homlof/saleheen-2017-table-s1.xlsx')
colnames(saleheen2017) = gsub('[^a-z0-9_]','_',tolower(colnames(saleheen2017)))
# fill in blank gene fields
last_gene = saleheen2017$gene[1]
for (i in 2:nrow(saleheen2017)) {
  if (is.na(saleheen2017$gene[i])) {
    saleheen2017$gene[i] = last_gene
  } else {
    last_gene = saleheen2017$gene[i]
  }
}
# double check there are no zeroes
min(saleheen2017$gene__homozygous_plof_count, na.rm=T)
# check consequences
table(saleheen2017$consequence) # all look ok
# select genes to include
saleheen2017_include = sqldf("
select   gene
from     saleheen2017
where    confident_plof_ = 'Yes'
;")
saleheen2017_genes = sort(unique(update_symbols(saleheen2017_include$gene, remove_na=T)))
homlof_genes$saleheen2017 = homlof_genes$gene %in% saleheen2017_genes
sum(homlof_genes$saleheen2017)
sum(homlof_genes$saleheen2017 | homlof_genes$sulem2015)

# DeBoever 2018 Table S1 UKBB data
deboever2018 = read_xlsx('data/homlof/deboever-2018-table-s1-41467_2018_3910_MOESM4_ESM.xlsx',sheet = 'ko_genes')
colnames(deboever2018) = gsub('[^a-z0-9_]','_',tolower(colnames(deboever2018)))
deboever2018_include = sqldf("
select   x__1 gene
from     deboever2018
where    uk_biobank__maf___1__ -- require rare - MAF < 1%
;")
deboever2018_genes = sort(unique(update_symbols(deboever2018_include$gene, remove_na=T)))
homlof_genes$deboever2018 = homlof_genes$gene %in% deboever2018_genes
sum(homlof_genes$deboever2018)
sum(homlof_genes$deboever2018 | homlof_genes$saleheen2017 | homlof_genes$sulem2015)

# Narasimhan 2016 Table S1 ELGH data - temp until we get updated ELGH data
narasimhan2016 = read.table('data/homlof/narasimhan-2016-table-s1-genes.tsv',sep='\t',header=F)
narasimhan2016_genes = sort(unique(update_symbols(narasimhan2016$V1, remove_na=T)))
homlof_genes$narasimhan2016 = homlof_genes$gene %in% narasimhan2016_genes
sum(homlof_genes$narasimhan2016)
sum(homlof_genes$narasimhan2016 | homlof_genes$deboever2018 | homlof_genes$saleheen2017 | homlof_genes$sulem2015)

# Karcewski 2019 curated gnomAD hom LoF gene set
karcewski2019 = read.table('data/homlof/gnomad_hom_ko_genes_2019-07-01.txt',sep='\t',header=F)
karcewski2019_genes = sort(unique(update_symbols(karcewski2019$V1, remove_na=T)))
homlof_genes$karcewski2019 = homlof_genes$gene %in% karcewski2019_genes
sum(homlof_genes$karcewski2019)
sum(homlof_genes$karcewski2019 | homlof_genes$narasimhan2016 | homlof_genes$deboever2018 | homlof_genes$saleheen2017 | homlof_genes$sulem2015)

homlof_genes$any = homlof_genes$karcewski2019 | homlof_genes$narasimhan2016 | homlof_genes$deboever2018 | homlof_genes$saleheen2017 | homlof_genes$sulem2015

summary(homlof_genes)

write.table(homlof_genes,'data/homlof/homlof_genes_summary.tsv',sep='\t',col.names = T,row.names = F,quote=F)
