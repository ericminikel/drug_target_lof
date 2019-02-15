setwd('~/d/sci/src/drug_target_lof/')
options(stringsAsFactors=FALSE)
library(sqldf)
library(plotrix)

# the original ExAC paper repo still has some useful colors and functions and so on
source('~/d/sci/src/exac_2015/exac_constants.R') # https://github.com/macarthur-lab/exac_2015/blob/master/exac_constants.R

# constants
gnomad = 141456

# generally useful functions
alpha = function(rgb_hexcolor, proportion) {
  hex_proportion = sprintf("%02x",round(proportion*255))
  rgba = paste(rgb_hexcolor,toupper(hex_proportion),sep='')
  return (rgba)
}

percent = function(proportion,digits=2) {
  return ( gsub(' ','',paste(formatC(proportion*100, digits=digits, format='fg'),"%",sep="") ) )
}

expand.range = function(raw_range, by=.5) {
  return ( c(raw_range[1]-by, raw_range[2]+by) )
}

# convert amino acid acronyms
tla_to_ola = function(x) {
  mapping = data.frame(tla=c("ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","TER","THR","TRP","TYR","VAL"),
                       ola=c("A",  "R",  "N",  "D",  "C",  "Q",  "E",  "G",  "H",  "I",  "L",  "K",  "M",  "F",  "P",  "S",  "X",  "T",  "W",  "Y",  "V"))
  for (row in 1:dim(mapping)[1]) {
    x = gsub(mapping$tla[row], mapping$ola[row], toupper(x))
  }
  return (x)
}


### begin reading in datasets

universe = read.table('lists/universe.tsv',sep='\t',header=F)
colnames(universe) = c('symbol')

# gnomad constraint table
if (!('gstraint' %in% ls())) {
  gstraint_all = read.table('data/constraint/constraint.txt.bgz',sep='\t',header=T) # read.table natively handles bgzipped files : ) 
  # save gstraint_all table for transcript-specific lookup later on
  gstraint = gstraint_all
  # for this analysis only use canonical transcripts in the HGNC universe
  gstraint = subset(gstraint, canonical == 'true' & gene %in% universe$symbol)
  # N=43 genes still have >1 transcript marked as canonical. just remove these altogether
  gstraint = subset(gstraint, !(gene %in% gstraint$gene[duplicated(gstraint$gene)]))
  # also remove genes with NA for obs_lof and/or exp_lof -- there are 301 of these
  gstraint = subset(gstraint, !(is.na(obs_lof) | is.na(exp_lof)))
}

# CAF and downsampled CAF
if (!('caf' %in% ls())) {
  cafall = read.table('data/caf/full_lof_metrics_by_transcript_an_adj_by_gene.txt.gz',sep='\t',header=T)
  caf_dsamp = read.table('data/caf/downsampled.subset.txt',sep='\t',header=T) # see src/subset_downsampled_caf.bash for genesis of this file
}

# GTEx expression data
if (!('tissues' %in% ls())) {
  tissues = read.table('data/expression/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz',sep='\t',header=T,comment.char='',skip=2)
  tissues$ensg = gsub('\\..*','',tissues$gene_id)
  tissues$n1tpm = 0
  for (i in 1:nrow(tissues)) {
    tissues$n1tpm[i] = sum(tissues[i,3:55] > 1)
    cat(paste("\r",percent(i/nrow(tissues))," done...",sep=""))
    flush.console()
  }
}

# HGNC-Ensembl mapping for GTEx
hgnc_ensembl = read.table('data/expression/hgnc_ensembl_mapping.tsv',sep='\t',header=T,quote='',comment.char='')
colnames(hgnc_ensembl)[c(2,9)] = c('symbol','ensg')
tissues$symbol = hgnc_ensembl$symbol[match(tissues$ensg,hgnc_ensembl$ensg)]
  
# merge all datasets into one "genes" table to make sure we have the same denominator for everything,
# no confounders due to missing values, etc.
genes = universe
# CAF fields
genes$bp = cafall$cds_length[match(genes$symbol,cafall$gene)]
genes$p = cafall$p[match(genes$symbol,cafall$gene)]
genes$p[is.na(genes$p)] = 0.0 # fill zeroes where p is absent. Konrad confirmed these are where no LoF is observed.
genes$p_fin = cafall$p_fin[match(genes$symbol,cafall$gene)]
genes$p_fin[is.na(genes$p_fin)] = 0.0 # fill zeroes again
# downsampled CAF fields for supplement
caf_df_global = subset(caf_dsamp, pop=='global' & downsampling==10824)
genes$p_dsamp_global = caf_df_global$caf[match(genes$symbol, caf_df_global$gene)]
genes$p_dsamp_global = pmin(genes$p_dsamp_global, 1.0) # CAF has a weird artifact that a small number of genes are >1; truncate the distribution
genes$p_dsamp_global[is.na(genes$p_dsamp_global)] = 0.0
genes$p_dsamp_finn = cafall$classic_caf_fin[match(genes$symbol, cafall$gene)]
genes$p_dsamp_finn = pmin(genes$p_dsamp_finn, 1.0) # CAF has a weird artifact that a small number of genes are >1; truncate the distribution
genes$p_dsamp_finn[is.na(genes$p_dsamp_finn)] = 0.0
# constraint fields
genes$obs_lof = gstraint$obs_lof[match(genes$symbol,gstraint$gene)]
genes$exp_lof = gstraint$exp_lof[match(genes$symbol,gstraint$gene)]
genes$oe_lof = genes$obs_lof / genes$exp_lof
genes$obs_mis = gstraint$obs_mis[match(genes$symbol,gstraint$gene)]
genes$exp_mis = gstraint$exp_mis[match(genes$symbol,gstraint$gene)]
genes$oe_mis = genes$obs_mis / genes$exp_mis
genes$obs_syn = gstraint$obs_syn[match(genes$symbol,gstraint$gene)]
genes$exp_syn = gstraint$exp_syn[match(genes$symbol,gstraint$gene)]
genes$oe_syn = genes$obs_syn / genes$exp_syn
# tissue fields
genes$n1tpm = tissues$n1tpm[match(genes$symbol, tissues$symbol)]

# remove TTN for Figure 1
non_ttn = genes$symbol != 'TTN'

### Begin Figure 1
png('figures/figure1.png',width=2400,height=900,res=300)
par(mfrow=c(1,3))
plot(genes$exp_syn[non_ttn], genes$obs_syn[non_ttn], xaxs='i', yaxs='i', xlim=c(0,1000), ylim=c(0,1000), pch=20, cex=0.5, col=alpha(k_syn,0.2), xlab='expected', ylab='observed', yaxt='n', xaxt='n')
m = lm(obs_syn ~ exp_syn + 0, data=genes[non_ttn,])
abline(m, col=k_syn, lwd=1)
axis(side=1,at=(0:2)*500)
axis(side=2,at=(0:2)*500,las=2)
abline(a=0,b=1)
mtext('a', side=3, cex=2.0, adj = -0.1, line = 0.3)
mtext(side=3,line=1,text='synonymous',font=2,cex=1,col=k_syn)
plot(genes$exp_mis[non_ttn], genes$obs_mis[non_ttn], xaxs='i', yaxs='i', xlim=c(0,2500), ylim=c(0,2500), pch=20, cex=0.5, col=alpha(k_mis,0.2), xlab='expected', ylab='observed', yaxt='n', xaxt='n')
m = lm(obs_mis ~ exp_mis + 0, data=genes[non_ttn,])
abline(m, col=k_mis, lwd=1)
axis(side=1,at=(0:3)*1000)
axis(side=2,at=(0:3)*1000,las=2)
abline(a=0,b=1)
mtext('b', side=3, cex=2.0, adj = -0.1, line = 0.3)
mtext(side=3,line=1,text='missense',font=2,cex=1,col=k_mis)
plot(genes$exp_lof[non_ttn], genes$obs_lof[non_ttn], xaxs='i', yaxs='i', xlim=c(0,150), ylim=c(0,150), pch=20, cex=0.5, col=alpha(k_lof,0.2), xlab='expected', ylab='observed', yaxt='n', xaxt='n')
m = lm(obs_lof ~ exp_lof + 0, data=genes[non_ttn,])
abline(m, col=k_lof, lwd=1)
axis(side=1,at=(0:2)*100)
axis(side=2,at=(0:2)*100,las=2)
abline(a=0,b=1)
mtext(side=3,line=1,text='pLoF',font=2,cex=1,col=k_lof)
mtext('c', side=3, cex=2.0, adj = -0.1, line = 0.3)
dev.off() ### -- End Figure 1


# Read in gene lists for Figure 2
list_metadata = read.table(textConnection("
filename|display
universe|all
olfactory_receptors|olfactory receptors
homozygous_lof_tolerant_twohit|homozygous LoF tolerant
blekhman_ar|autosomal recessive
blekhman_ad|autosomal dominant
CEGv2_subset_universe|essential in culture
clingen_level3_genes_2018_09_13|ClinGen haploinsufficient
drug_targets|drug targets
positive_targets|positive
negative_targets|negative
other_targets|other & unknown
rhodop_gpcr|rhodopsin-like GPCRs
ion_channels|ion channels
nuclear_receptors|nuclear receptors
enzymes|enzymes
gwascatalog|GWAS hits
"),sep='|',header=TRUE)

for (i in 1:nrow(list_metadata)) {
  path = paste('lists/',list_metadata$filename[i],'.tsv',sep='')
  gene_list = read.table(path)
  if (dim(gene_list)[1] < 1) {
    print(paste("No contents in gene list: ",list_metadata$filename[i],sep=''))
    next
  }
  genes[,list_metadata$filename[i]] = genes$symbol %in% gene_list$V1
}


# Make a table with all the forest plot data for Figure 2
lof_oe = list_metadata[1:11,]
lof_oe$mean = 0.0
lof_oe$upper95 = 0.0
lof_oe$lower95 = 0.0
lof_oe$n = 0
i = 1
for (i in 1:dim(lof_oe)[1]) {
  in_list = as.logical(genes[,lof_oe$filename[i]])
  oe_mean = mean(genes$oe_lof[in_list], na.rm=T)
  oe_sd = sd(genes$oe_lof[in_list], na.rm=T)
  oe_n = sum(!is.na(genes$oe_lof[in_list])) # none are NA anyway - sum(is.na(genes$lof_obs_exp)) == 0
  lof_oe$mean[i] = oe_mean
  lof_oe$upper95[i] = oe_mean + 1.96 * oe_sd/sqrt(oe_n)
  lof_oe$lower95[i] = oe_mean - 1.96 * oe_sd/sqrt(oe_n)
  lof_oe$n[i] = oe_n
}

lof_oe$y = -1:(-1*dim(lof_oe)[1])
lof_oe$color = '#777777'
lof_oe$color[lof_oe$filename=='universe'] = '#000000'
lof_oe$color[lof_oe$filename %in% c('drug_targets','positive_targets','negative_targets','other_targets')] = '#2E37FE'



### Begin Figure 2
png('figures/figure2.png',width=2400,height=2000,res=300)
layout(matrix(c(1,1,2,2,2),byrow=T,nrow=5))

# Panel A: histogram
# note it is a histogram even though it is plotted with polygon() so it looks like a density plot
par(mar=c(4,5,4,3))
color_all = '#000000'
color_drug = '#2E37FE'
hist_breaks = c((0:20)*5/100,10) # 0%, 5%, 10%, ... 100% and then 1,000% (~+Inf)
h_all = hist(genes$oe_lof[genes$universe], breaks=hist_breaks,plot=FALSE)
h_drug = hist(genes$oe_lof[genes$drug_targets],plot=FALSE,breaks=hist_breaks)
h_all[['proportion']] = h_all$counts / sum(h_all$counts)
h_drug[['proportion']] = h_drug$counts / sum(h_drug$counts)
h_all_mean = mean(genes$oe_lof[genes$universe],na.rm=TRUE)
h_drug_mean = mean(genes$oe_lof[genes$drug_targets],na.rm=TRUE)
plot(NA, NA, xlim=c(0,1), ylim=c(0,.15), xaxs='i', yaxs='i', axes=FALSE, xlab='', ylab='')
polygon(c(0,h_all$breaks), c(0,h_all$proportion,0), col=alpha(color_all,.5),border=NA)
points(h_all$breaks[1:(length(h_all$breaks)-1)], h_all$proportion, col=color_all, type='l', lwd=5)
polygon(c(0,h_drug$breaks), c(0,h_drug$proportion,0), col=alpha(color_drug,.5),border=NA)
points(h_drug$breaks[1:(length(h_drug$breaks)-1)], h_drug$proportion, col=color_drug, type='l', lwd=5)
abline(h=0,lwd=2)
axis(side=1, at=(0:4)/4, labels=c(percent((0:3)/4),'\u2265100%'), lwd=0, lwd.ticks=1)
axis(side=2, at=(0:4)/20, labels=percent((0:4)/20), lwd=0, lwd.ticks=1, las=2)
abline(v=c(0,1),lwd=2)
mtext(side=1, text='pLoF obs/exp ratio', line=2.5, font=1, cex=1)
mtext(side=2, text='proportion of genes', line=3.0, font=1, cex=1)
segments(x0=c(h_all_mean,h_drug_mean),y0=.10,y1=0,col=c(color_all,color_drug),lwd=2,lty=3)
text(x=c(h_all_mean,h_drug_mean),y=c(.11,.11),pos=c(4,2),font=2,labels=paste(c('all\ngenes\n','drug\ntargets\n'),'mean =',percent(c(h_all_mean,h_drug_mean))),col=c(color_all,color_drug),cex=1.1)
mtext('a', side=3, cex=2, adj = -0.05, line = 0.5)

# Panel B: Forest plot
par(mar=c(4,18,3,3))
plot(NA,NA, xlim=c(0,1), ylim=expand.range(range(lof_oe$y),by=.75), axes=FALSE, xlab='', ylab='')
abline(v=lof_oe$mean[lof_oe$filename=='universe'], lty=3, col='#333333', lwd=2)
par(xpd=T) # allow 95%CI to extend beyond xlims of [0,1]
segments(x0=lof_oe$lower95, x1=lof_oe$upper95, y0=lof_oe$y, col=lof_oe$color, lwd=3)
par(xpd=F)
points(x=lof_oe$mean, y=lof_oe$y, col=lof_oe$color, pch=19, cex=1.5)
axis(side=1, at=(0:4)/4, labels=percent((0:4)/4), lwd=0, lwd.ticks=1)
abline(v=c(0,1))
mtext(side=2, at=lof_oe$y, text=lof_oe$display, las=2, cex = .9)
mtext(side=1, text='mean (95%CI) pLoF obs/exp ratio', cex=1, line = 2.5)
par(xpd=T)
abline(h=lof_oe$y[lof_oe$filename=='universe']+c(.5,-.5), col='#777777', lwd=.5)
abline(h=lof_oe$y[lof_oe$filename=='drug_targets']+c(.5,-.5,-3.5), col='#777777', lwd=.5)
par(xpd=F)
mtext(side=2, at=-4.5, line=13, text='gene lists\nfor comparison\n', cex=0.8, font=2)
mtext(side=2, at=-10.0, line=13, text="by drug\neffect\non target", cex=0.8, font=2)

mtext('b', side=3, cex=2, adj = -0.4, line = 0.3)

dev.off() ### -- End Figure 2



# stats & numbers quoted in text apropos Figure 2
ks.test(genes$oe_lof[genes$drug_targets], genes$oe_lof) # is the obs/exp difference between drug targets & all genes significant?
mean(genes$oe_lof[genes$clingen_level3_genes_2018_09_13], na.rm=T) # mean obs/exp for HI genes
sum(genes$oe_lof < mean(genes$oe_lof[genes$clingen_level3_genes_2018_09_13], na.rm=T), na.rm=T) # how many genes are below this threshold?
mean(genes$oe_lof < mean(genes$oe_lof[genes$clingen_level3_genes_2018_09_13], na.rm=T), na.rm=T) # what proportion of genes are below this threshold?
sum(genes$drug_targets & genes$oe_lof < mean(genes$oe_lof[genes$clingen_level3_genes_2018_09_13], na.rm=T), na.rm=T) # how many drug targets below this threshold?
sum(genes$drug_targets & genes$oe_lof < mean(genes$oe_lof[genes$clingen_level3_genes_2018_09_13], na.rm=T), na.rm=T) / sum(genes$drug_targets & !is.na(genes$oe_lof)) # what proportion of drug targets?
sum(genes$negative_targets & genes$oe_lof < mean(genes$oe_lof[genes$clingen_level3_genes_2018_09_13], na.rm=T), na.rm=T) # how many *negative* drug targets below this threshold?
mean(genes$oe_lof[genes$drug_targets], na.rm=T) # average obs/exp of all drug targets
mean(genes$oe_lof[genes$negative_targets], na.rm=T) # of just negative targets
mean(genes$oe_lof[genes$positive_targets], na.rm=T) # of just positive targets
sum(genes$positive_targets & genes$negative_targets, na.rm=T) # how many genes are both positive and negative targets
ks.test(genes$oe_lof[genes$negative_targets], genes$oe_lof[genes$positive_targets]) # is the difference between positive & negative targets significant?

# For picking exampels for Table 1 - browse some examples of negative drugs and their targets' constraint values
drug_gene_action_match = read.table('data/drugbank/drug_gene_action_match.tsv',sep='\t',header=T)
drug_gene_action_match = subset(drug_gene_action_match, gene != '' & drug != '' & action != '')
drug_gene_action_match$negative = genes$negative_targets[match(drug_gene_action_match$gene, genes$symbol)]
drug_gene_action_match$oe_lof = genes$oe_lof[match(drug_gene_action_match$gene,genes$symbol)]
drug_gene_action_match$obs_lof = genes$obs_lof[match(drug_gene_action_match$gene,genes$symbol)]
drug_gene_action_match$exp_lof = genes$exp_lof[match(drug_gene_action_match$gene,genes$symbol)]
drug_gene_action_match = subset(drug_gene_action_match, !is.na(oe_lof))
subset_to_browse = drug_gene_action_match[drug_gene_action_match$negative & drug_gene_action_match$exp_lof > 5,]
subset_to_browse = subset_to_browse[with(subset_to_browse,order(oe_lof,gene)),]
# View(subset_to_browse)

# Stats and numbers quoted in text apropos Figure 3

# how enriched are the top four "druggable" classes among drug targets?
genes$any_canonically_druggable_class = genes$rhodop_gpcr | genes$enzymes | genes$ion_channels | genes$nuclear_receptors
drug_classes = table(genes[,c("drug_targets","any_canonically_druggable_class")])
drug_classes
f = fisher.test(drug_classes)
f
f$p.value


# check for overlap between the classes
table(genes[,c('rhodop_gpcr','enzymes')]) # none
table(genes[,c('rhodop_gpcr','ion_channels')]) # none
table(genes[,c('rhodop_gpcr','nuclear_receptors')]) # none
table(genes[,c('enzymes','ion_channels')]) # 2 overlap
table(genes[,c('enzymes','nuclear_receptors')]) # none
table(genes[,c('nuclear_receptors','ion_channels')]) # none
# since ion channel is more specific than enzymes, that category should "win" for those 2
# so assign in this order:
genes$family = 'other'
genes$family[genes$rhodop_gpcr] = 'rhodop_gpcr'
genes$family[genes$enzymes] = 'enzymes'
genes$family[genes$ion_channels] = 'ion_channels'
genes$family[genes$nuclear_receptors] = 'nuclear_receptors'

# how does controlling for family affect obs/exp of drug targets vs. all genes?
m = lm(genes$oe_lof ~ genes$drug_targets + genes$family)
summary(m)

# are GWAS hits enriched in drug targets?
drug_gwas = table(genes[,c('gwascatalog','drug_targets')])
drug_gwas
fisher.test(drug_gwas)

# how does controlling for GWAS hit status affect obs/exp of drug targets vs. all genes?
m = lm(genes$oe_lof ~ genes$drug_targets + genes$gwascatalog)
summary(m)


# average number of tissues exprssed for drug targets and non drug targets:
mean(genes$n1tpm[genes$drug_targets], na.rm=T)
mean(genes$n1tpm, na.rm=T)
# is the difference significant?
ks.test(genes$n1tpm[genes$drug_targets], genes$n1tpm)
# look at the full distributions too:
#hist(genes$n1tpm[genes$drug_targets]) 
#hist(genes$n1tpm[!genes$drug_targets])
# is constraint correlated with expression:
cor.test(genes$n1tpm, genes$oe_lof, method='spearman')
cor.test(genes$n1tpm, genes$oe_lof, method='spearman')$p.value
# how does expression across tissues affect obs/exp of drug targets vs. all genes?
m = lm(oe_lof ~ drug_targets + n1tpm, data=genes)
summary(m)


# create forest plot data for Figure 3

forest2 = list_metadata[c(1,8,12:15),]
forest2$mean = 0.0
forest2$upper95 = 0.0
forest2$lower95 = 0.0
forest2$n = 0

forest2$enrichment = 0.0
forest2$enrichment_u95 = 0.0
forest2$enrichment_l95 = 0.0

forest2$targetcount = 0.0

forest2$drug_mean = 0.0
forest2$drug_upper95 = 0.0
forest2$drug_lower95 = 0.0
forest2$drug_n = 0

i = 1
for (i in 1:dim(forest2)[1]) {
  in_list = as.logical(genes[,forest2$filename[i]])
  oe_mean = mean(genes$oe_lof[in_list], na.rm=T)
  oe_sd = sd(genes$oe_lof[in_list], na.rm=T)
  oe_n = sum(!is.na(genes$oe_lof[in_list]))
  forest2$mean[i] = oe_mean
  forest2$upper95[i] = oe_mean + 1.96 * oe_sd/sqrt(oe_n)
  forest2$lower95[i] = oe_mean - 1.96 * oe_sd/sqrt(oe_n)
  forest2$n[i] = oe_n
  in_list = as.logical(genes[,forest2$filename[i]] & genes[,'drug_targets'])
  oe_mean = mean(genes$oe_lof[in_list], na.rm=T)
  oe_sd = sd(genes$oe_lof[in_list], na.rm=T)
  oe_n = sum(!is.na(genes$oe_lof[in_list]))
  forest2$drug_mean[i] = oe_mean
  forest2$drug_upper95[i] = oe_mean + 1.96 * oe_sd/sqrt(oe_n)
  forest2$drug_lower95[i] = oe_mean - 1.96 * oe_sd/sqrt(oe_n)
  forest2$drug_n[i] = oe_n
  if (i > 2) { # don't do this for the 1st two rows
    ctable = table(genes[,c('drug_targets',forest2$filename[i])])
    fisher_result = fisher.test(ctable, alternative = 'two.sided')
    forest2$enrichment[i] = fisher_result$estimate
    forest2$enrichment_l95[i] = fisher_result$conf.int[1]
    forest2$enrichment_u95[i] = fisher_result$conf.int[2]
    
    forest2$targetcount[i] = ctable['TRUE','TRUE']
  }
  
}

forest2$y = -1:(-1*dim(forest2)[1])
forest2$color[1] = '#000000'
forest2$color[2] = '#666666'
forest2$color[3:6] = c('#008837','#EEB422','#2C7BB6','#7B3294')


### Begin Figure 3
png('figures/figure3.png',width=2400,height=1000,res=300)
layout_matrix = matrix(c(1,1,1,2,3),nrow=1,byrow=T)
layout(layout_matrix, heights=c(1,1))

par(mar=c(4,15,3,3))

ylims = expand.range(range(forest2$y),by=.75)

# panel A - forest plot
plot(NA,NA, xlim=c(0,1), ylim=ylims, axes=FALSE, xlab='', ylab='')
abline(v=forest2$mean[forest2$filename=='universe'], lty=3, col='#333333', lwd=2)
segments(x0=forest2$lower95, x1=forest2$upper95, y0=forest2$y, col=forest2$color, lwd=3)
points(x=forest2$mean, y=forest2$y, col=forest2$color, pch=19, cex=1.5)
empty_circle_color = forest2$color[forest2$filename=='drug_targets']
points(x=forest2$drug_mean[3:6], y=forest2$y[3:6], col=empty_circle_color, pch=1, lwd=2, cex=1.5)

axis(side=1, at=(0:4)/4, labels=percent((0:4)/4), lwd=0, lwd.ticks=1)
abline(v=c(0,1))
mtext(side=2, at=forest2$y, text=forest2$display, las=2, cex = .9)
mtext(side=1, text='mean (95%CI) pLoF obs/exp ratio', cex=.7, line = 2.5)
par(xpd=T)
abline(h=forest2$y[forest2$filename=='drug_targets']+c(.5,-.5), col='#777777', lwd=.5)
abline(h=forest2$y[forest2$filename=='enzymes']+c(-.5), col='#777777', lwd=.5)
par(xpd=F)
mtext('a', side=3, cex=2, adj = 0.0, line = 0.3)

par(mar=c(4,2,3,1))

# panel B - barplot of fold enrichment with 95% CI error bars
plot(NA, NA, xlim=c(1,64), ylim=ylims, axes=FALSE, ann=FALSE,log='x',xaxs='i',yaxs='i')
abline(v=1)
for (i in 2:6) {
  segments(x0=1,x1=forest2$enrichment[i],y0=forest2$y[i],y1=forest2$y[i],col=forest2$color[i],lend=1,lwd=25)
  arrows(x0=forest2$enrichment_l95[i],x1=forest2$enrichment_u95[i],y0=forest2$y[i],y1=forest2$y[i],col='black',lwd=1,angle=90,length=0.05,code=3)
}
axis(side=1, at=2^(0:6), labels=NA)
axis(side=1, at=2^(c(0,2,4,6)), labels=2^(c(0,2,4,6)), lwd=0, lwd.ticks=0)
mtext(side=1,line=2.5,text='fold enrichment',cex=.7)
mtext('b', side=3, cex=2, adj = 0.0, line = 0.3)

par(mar=c(4,2,3,1))

# panel C - stacked bar plot of how many drug targets are in each family
target_counts = as.data.frame(table(genes$family[genes$drug_targets]))
colnames(target_counts) = c('family','count')
target_counts$color = forest2$color[match(target_counts$family,forest2$filename)]
target_counts$y = forest2$y[match(target_counts$family,forest2$filename)]
target_counts$color[target_counts$family=='other'] = '#AAAAAA'
target_counts$y[target_counts$family=='other'] = -7
target_counts = target_counts[with(target_counts, order(y)),]
plot(NA,NA,xlim=c(0,1),ylim=c(0,400),axes=FALSE,ann=FALSE,yaxs='i')
running_total = 0
for (i in 1:nrow(target_counts)) {
  rect(ybottom=running_total, ytop=running_total + target_counts$count[i], xleft=0, xright=1, col=target_counts$color[i], border=NA)
  running_total = running_total + target_counts$count[i]
}
text(x=0.5,y=100,labels='other',font=3)
axis(side=1, at=c(-100,100))
axis(side=2,at=(0:4)*100,las=2)
mtext('c', side=3, cex=2, adj = 0.0, line = 0.3)

dev.off() ### -- End Figure 3


### Begin Figure 4 - possibility of ascertaining LoF individuals

png('figures/figure4.png',width=2400,height=3600,res=500)

par(mfrow=c(3,1))

# decide on right histogram breaks for all plots
lowest_p = min(genes$p[genes$p > 0])
lowest_p ^ 2 # 1.58e-11
# so the axis should break at 1e-11 or lower
axis_break = 10^-11.9
# confirm that the leftmost bin (<1e-10) is composed entirely of genes with 0 CAF, even for homs
sum(genes$p^2 < 1e-11) == sum(genes$p == 0)
# but the next bin this is not the case
sum(genes$p^2 < 3.16e-11) == sum(genes$p == 0)

caf_breaks = c(0,10^seq(-12,0,by=.5))
n_bins = length(caf_breaks) - 1

ymax = 6000
label_y = 4500

world_population = 7.66e9
#gnomad = 141456 # defined above

het_k = '#FF9912'
hom_k = '#984EA3'

het_hist = hist(2*genes$p*(1-genes$p), breaks=caf_breaks, plot=FALSE)
hom_hist = hist(genes$p ^ 2, breaks=caf_breaks, plot=FALSE)

par(mar=c(2,5,4,3))
plot(NA, NA, xlim=c(-13,0), ylim=c(0,ymax), ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
axis(side=1, at=c(-13,-9,-6,-3,0), labels=c('','1 in 1 billion','1 in 1 million','1 in 1,000','100%'), tck=-0.03, lwd=1, lwd.ticks=1)
axis(side=1, at=c(-12.3), labels=c('zero'), lwd=0, lwd.ticks=0)
axis(side=1, at=-12:0, labels=NA, lwd=0, lwd.ticks=1, tck=-0.01)
axis.break(axis=1,breakpos=log10(axis_break),style='zigzag')
axis(side=2, at=(0:6)*1000, labels=formatC((0:6)*1000,big.mark=','),las=2)
mtext(side=2, line=3.5, text='number of genes')
rect(xleft=log10(het_hist$mids[1])-.25, xright=log10(het_hist$mids[1])+.25, ybottom=0, ytop=het_hist$counts[1], lwd=3, border=het_k, col=alpha(het_k,.2))
points(log10(het_hist$mids[2:n_bins]), het_hist$counts[2:n_bins], type='l', lwd=3, col=het_k)
polygon(x=c(log10(het_hist$mids[2:n_bins]),1), y=c(het_hist$counts[2:n_bins],0), border=NA, col=alpha(het_k,.2))
rect(xleft=log10(hom_hist$mids[1])-.25, xright=log10(hom_hist$mids[1])+.25, ybottom=0, ytop=hom_hist$counts[1], lwd=3, border=hom_k, col=alpha(hom_k,.2))
points(log10(hom_hist$mids[2:n_bins]), hom_hist$counts[2:n_bins], type='l', lwd=3, col=hom_k)
polygon(x=c(log10(hom_hist$mids[2:n_bins]),1), y=c(hom_hist$counts[2:n_bins],0), border=NA, col=alpha(hom_k,.2))
segments(x0=log10(1/world_population),y0=0,y1=label_y,lwd=0.5)
segments(x0=log10(1/gnomad),y0=0,y1=label_y,lwd=0.5)
text(x=log10(1/world_population),y=label_y,pos=3,labels='1\non\nEarth')
text(x=log10(1/gnomad),y=label_y,pos=3,labels='1\nin\ngnomAD')
par(xpd=T)
legend(x=-4,y=7500,legend=c('heterozygotes','homozygotes &\ncompound hets'),text.col=c(het_k,hom_k),border=c(het_k,hom_k),fill=alpha(c(het_k,hom_k),.2),bty='n',cex=1.2)
par(xpd=F)

mtext('a', side=3, cex=2, adj = 0.0, line = 0.3)



finland_population = 5.50e6
gnomad_finns = 12526

het_hist = hist(2*genes$p_fin*(1-genes$p_fin), breaks=caf_breaks, plot=FALSE)
hom_hist = hist(genes$p_fin ^ 2, breaks=caf_breaks, plot=FALSE)

par(mar=c(3,5,3,3))

plot(NA, NA, xlim=c(-13,0), ylim=c(0,ymax), ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
axis(side=1, at=c(-13,-9,-6,-3,0), labels=c('','1 in 1 billion','1 in 1 million','1 in 1,000','100%'), tck=-0.03, lwd=1, lwd.ticks=1)
axis(side=1, at=c(-12.3), labels=c('zero'), lwd=0, lwd.ticks=0)
axis(side=1, at=-12:0, labels=NA, lwd=0, lwd.ticks=1, tck=-0.01)
axis.break(axis=1,breakpos=log10(axis_break),style='zigzag')
#mtext(side=1, line=2.5, text='expected frequency of pLoF individuals in population')
axis(side=2, at=(0:6)*1000, labels=formatC(c((0:5)*1000,het_hist$counts[1]),big.mark=',',format='fg'),las=2,lwd=1)
axis.break(axis=2,breakpos=5500,style='zigzag')
mtext(side=2, line=3.5, text='number of genes')
rect(xleft=log10(het_hist$mids[1])-.25, xright=log10(het_hist$mids[1])+.25, ybottom=0, ytop=min(het_hist$counts[1],ymax), lwd=3, border=het_k, col=alpha(het_k,.2))
points(log10(het_hist$mids[2:n_bins]), het_hist$counts[2:n_bins], type='l', lwd=3, col=het_k)
polygon(x=c(log10(het_hist$mids[2:n_bins]),1), y=c(het_hist$counts[2:n_bins],0), border=NA, col=alpha(het_k,.2))
rect(xleft=log10(hom_hist$mids[1])-.25, xright=log10(hom_hist$mids[1])+.25, ybottom=0, ytop=min(hom_hist$counts[1],ymax), lwd=3, border=hom_k, col=alpha(hom_k,.2))
points(log10(hom_hist$mids[2:n_bins]), hom_hist$counts[2:n_bins], type='l', lwd=3, col=hom_k)
polygon(x=c(log10(hom_hist$mids[2:n_bins]),1), y=c(hom_hist$counts[2:n_bins],0), border=NA, col=alpha(hom_k,.2))
segments(x0=log10(1/finland_population),y0=0,y1=label_y,lwd=0.5)
segments(x0=log10(1/gnomad_finns),y0=0,y1=label_y,lwd=0.5)
text(x=log10(1/finland_population),y=label_y,pos=3,labels='1\nin\nFinland')
text(x=log10(1/gnomad_finns),y=label_y,pos=3,labels='1\nin\ngnomAD Finns')
par(xpd=T)
par(xpd=F)

mtext('b', side=3, cex=2, adj = 0.0, line = 0.3)



a = 0.05766 # mean % genome autozygous in ELGH according to Hilary Martin's email 2018-12-06

het_hist = hist(2*genes$p*(1-genes$p), breaks=caf_breaks, plot=FALSE)
hom_hist = hist((1-a)*(genes$p ^ 2) + a*genes$p, breaks=caf_breaks, plot=FALSE)

world_consang_population = world_population * .104 # Bittles & Black 2009 PNAS
gnomad_consang = 2912 # SAS w/ F > 0.05 according to Konrad on 2018-12-04 via Slack

par(mar=c(4,5,2,3))

plot(NA, NA, xlim=c(-13,0), ylim=c(0,ymax), ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
axis(side=1, at=c(-13,-9,-6,-3,0), labels=c('','1 in 1 billion','1 in 1 million','1 in 1,000','100%'), tck=-0.03, lwd=1, lwd.ticks=1)
axis(side=1, at=c(-12.3), labels=c('zero'), lwd=0, lwd.ticks=0)
axis(side=1, at=-12:0, labels=NA, lwd=0, lwd.ticks=1, tck=-0.01)
axis.break(axis=1,breakpos=log10(axis_break),style='zigzag')
mtext(side=1, line=2.5, text='expected frequency of pLoF individuals in population')
axis(side=2, at=(0:6)*1000, labels=formatC((0:6)*1000,big.mark=','),las=2)
mtext(side=2, line=3.5, text='number of genes')
rect(xleft=log10(het_hist$mids[1])-.25, xright=log10(het_hist$mids[1])+.25, ybottom=0, ytop=het_hist$counts[1], lwd=3, border=het_k, col=alpha(het_k,.2))
points(log10(het_hist$mids[2:n_bins]), het_hist$counts[2:n_bins], type='l', lwd=3, col=het_k)
polygon(x=c(log10(het_hist$mids[2:n_bins]),1), y=c(het_hist$counts[2:n_bins],0), border=NA, col=alpha(het_k,.2))
rect(xleft=log10(hom_hist$mids[1])-.25, xright=log10(hom_hist$mids[1])+.25, ybottom=0, ytop=hom_hist$counts[1], lwd=3, border=hom_k, col=alpha(hom_k,.2))
points(log10(hom_hist$mids[2:n_bins]), hom_hist$counts[2:n_bins], type='l', lwd=3, col=hom_k)
polygon(x=c(log10(hom_hist$mids[2:n_bins]),1), y=c(hom_hist$counts[2:n_bins],0), border=NA, col=alpha(hom_k,.2))
segments(x0=log10(1/world_consang_population),y0=0,y1=label_y,lwd=0.5)
segments(x0=log10(1/gnomad_consang),y0=0,y1=label_y,lwd=0.5)
text(x=log10(1/world_consang_population),y=label_y,pos=3,labels='1\nconsanguineous\nworldwide')
text(x=log10(1/gnomad_consang),y=label_y,pos=3,labels='1\nconsanguineous\nin gnomAD')
par(xpd=T)
par(xpd=F)

mtext('c', side=3, cex=2, adj = 0.0, line = 0.3)

dev.off() ### -- End Figure 4


# stats for text apropos Figure 4

corrected_threshold = 0.05 / nrow(genes)

# In this sample size, 75% of genes (N=14,340) would still be expected to have <1 double null individual...
sum(!is.na(genes$p) & genes$p ^ 2 * gnomad * 100 < 1)
sum(!is.na(genes$p))
sum(!is.na(genes$p) & genes$p ^ 2 * gnomad * 100 < 1) / sum(!is.na(genes$p))
(sum(!is.na(genes$p)) - sum(!is.na(genes$p) & genes$p ^ 2 * gnomad * 100 < 1))
(sum(!is.na(genes$p)) - sum(!is.na(genes$p) & genes$p ^ 2 * gnomad * 100 < 1)) / sum(!is.na(genes$p))

# and 91% of genes (N=17,546) — including 92% of existing approved drug targets (N=357) —would have sufficiently few expected "knockouts" ...
genes$qbinom_gnomad100 = qbinom(p=corrected_threshold, size=gnomad*100, prob=genes$p^2)
sum(genes$qbinom_gnomad100 == 0) # 17546
sum(genes$qbinom_gnomad100 == 0) / nrow(genes) # 91.4% of genes
(nrow(genes) - sum(genes$qbinom_gnomad100 == 0)) # 1648
(nrow(genes) - sum(genes$qbinom_gnomad100 == 0))/nrow(genes) # 8.6%

# including 92% of existing approved drug targets (N=357) 
sum(!is.na(genes$qbinom_gnomad100) & genes$qbinom_gnomad100 == 0 & genes$drug_targets)
sum(!is.na(genes$qbinom_gnomad100) & genes$drug_targets)
sum(!is.na(genes$qbinom_gnomad100) & genes$drug_targets) - sum(!is.na(genes$qbinom_gnomad100) & genes$qbinom_gnomad100 == 0 & genes$drug_targets)
(sum(!is.na(genes$qbinom_gnomad100) & genes$drug_targets)  - sum(!is.na(genes$qbinom_gnomad100) & genes$qbinom_gnomad100 == 0 & genes$drug_targets))  / sum(!is.na(genes$qbinom_gnomad100) & genes$drug_targets) 

# Indeed, for 38% of genes (N=7,546), even if all humans on Earth were sequenced, observing zero “knockouts” would still not be a statistically significant anomaly
genes$qbinom_world = qbinom(p=corrected_threshold, size=world_population, prob=genes$p^2)
sum(genes$qbinom_world == 0) # 7451
sum(genes$qbinom_world == 0) / nrow(genes) # 38% of genes
sum(genes$p^2 < 1/world_population) / nrow(genes) # how many genes have 0 expected on Earth

# end Figure 4-related stuff



### Begin Figure 5
png('figures/figure5.png',width=2400,height=2700,res=300)

par(mfrow=c(3,1), mar=c(5,4,3,2))

filter_parms = data.frame(disp=c('removed by LOFTEE filter','removed by manual curation','likely true LoF','known or putative GoF'),
                          color=c('#B2DF8A','#A6CEE3','#9D1309','#FF5333'),
                          cat=c('loftee','curated','true','gof'))

# notes from manual curation of variants in these genes:
curation = read.table('data/curation/curation_notes_2018-12-04.csv',sep=',',header=T)

## panel A: HTT
gtf = read.table('data/curation/HTT.distilled.gtf')
colnames(gtf) = c('tx','chrom','start','stop')


intron_plot_length = 50
can = subset(gtf,tx=='ENST00000355072')
can$len = can$stop - can$start
can$xstart = 0
for (i in 2:nrow(can)) {
  can$xstart[i] = can$xstart[i-1] + can$len[i - 1] + intron_plot_length
}
can$xstop = can$xstart + can$len
can$exon_num = as.character(1:67)
can$xmid = (can$xstart + can$xstop)/2


lof = read.table('data/curation/HTT-gnomAD-v2.1-ENSG00000197386-2018-11-14-19.51.12.csv',sep=',',header=T)
colnames(lof)[2] = 'pos'
colnames(lof) = gsub('[^a-z0-9_]','_',tolower(colnames(lof)))
lof$rel_ac = lof$allele_count / max(lof$allele_count[lof$flags==''])
lof$xplot = as.integer(NA)


for (i in 1:nrow(lof)) {
  for (j in 1:nrow(can)) {
    if (lof$pos[i] >= can$start[j] - 10 && lof$pos[i] <= can$stop[j] + 10) { # if this pos falls in this exon
      lof$xplot[i] = lof$pos[i] - can$start[j] + can$xstart[j]
    }
  }
}

lof$pos_id = paste(lof$chrom, formatC(lof$pos, width=9, flag='0'), lof$ref, lof$alt, sep='-')
curation$pos_id = curation$variant_id
lof$curated_verdict = curation$form_verdict[match(lof$pos_id, curation$pos_id)]
lof$form_other_notes = curation$form_other_notes[match(lof$pos_id, curation$pos_id)]
lof$cat = ''
lof$cat[lof$flags!=''] = 'loftee'
lof$cat[lof$curated_verdict %in% c('LoF','likely_LoF')] = 'true'
lof$cat[lof$curated_verdict %in% c('not_LoF','likely_not_LoF','uncertain_LoF','')] = 'curated'
lof$color = filter_parms$color[match(lof$cat,filter_parms$cat)]

# for curated constraint
obs_htt = sum(lof$source=='gnomAD Exomes' & nchar(lof$reference)==1 & nchar(lof$alternate)==1 & lof$cat=='true')
exp_htt = gstraint$exp_lof[gstraint$gene=='HTT']
obs_htt / exp_htt


# for curated CAF
sum(lof$allele_frequency)
sum(lof$allele_count[lof$cat=='true']) / (2*gnomad)

# write out curated data
lof$gene = 'HTT'
htt_supptbl = lof[,c('gene','pos_id','allele_count','cat','flags','curated_verdict','form_other_notes')]
colnames(htt_supptbl) = c('gene','variant','allele_count','filter_status','LOFTEE_flags','curation_result','comments')
#write.table(htt_supptbl,'data/htt_curation.tsv',sep='\t',row.names=F,col.names=T,quote=F)

# AC scale: 0, 1, 10, >100 = 0, 1, 2, 3
ac_scale = data.frame(disp=c('0','1','10','100+'),yval=c(0,1,2,3))
ac_ticks = c(0, log10(c(1:10, (2:10)*10)) + 1)
ac_trunc = 100
lof$rel_ac = log10(pmin(lof$allele_count,ac_trunc)) + 1

exon_ybot = -0.5
exon_ytop = 0


plot(NA, NA, xlim=c(0,max(can$xstop)), ylim=c(exon_ybot,max(ac_scale$yval)), yaxs='i', ann=FALSE, axes=FALSE)
abline(h=mean(c(exon_ybot, exon_ytop)))
rect(xleft=can$xstart, xright=can$xstop, ybottom=rep(exon_ytop,nrow(can)), ytop=rep(exon_ybot,nrow(can)), col='#000000')
axis(side=2, at=ac_ticks, labels=NA, tck=-0.01)
axis(side=2, at=ac_scale$yval, labels=ac_scale$disp, tck=-0.03, lwd=0, lwd.ticks=1, las=2)
mtext(side=2, line=2.5, text='allele count')
mtext(side=1, line=0, at=can$xmid[c(1,nrow(can))], text=can$exon_num[c(1,nrow(can))], font=3, cex=0.8)
mtext(side=1, line=0.5, text='HTT exon structure')

segments(x0=lof$xplot, y0=rep(0,nrow(lof)), y1=lof$rel_ac, lwd=3, lend=1, col=lof$color)

par(xpd=T)
legend(x=max(can$xstop)*.6,y=max(ac_scale$yval)*1.1,text.font=2,legend=filter_parms$disp,col=filter_parms$color,lwd=3,text.col=filter_parms$color,cex=1.5)
par(xpd=F)

mtext(side=3, line=0, text='HTT', font=3, cex=2)

mtext('a', side=3, cex=2, adj = 0.0, line = 0.3)



gtf = read.table('data/curation/MAPT.distilled.gtf')
colnames(gtf) = c('tx','chrom','start','stop')

expressed = read.table('data/expression/ENSG00000186868.11.summary.gz',sep='\t',header=T)
colnames(expressed)
brain_means = grep('^Brain.*mean_TPM',colnames(expressed))
all_means = grep('.*mean_TPM',colnames(expressed))
expressed$brain_mean = rowMeans(expressed[,brain_means])
expressed$brain_rel = expressed$brain_mean / max(expressed$brain_mean)
expressed$all_mean = rowMeans(expressed[,all_means])
expressed$all_rel = expressed$all_mean / max(expressed$all_mean)


intron_plot_length = 20
can = subset(gtf,tx=='ENST00000344290')
can$len = can$stop - can$start
can$xstart = 0
for (i in 2:nrow(can)) {
  can$xstart[i] = can$xstart[i-1] + can$len[i - 1] + intron_plot_length
}
can$xstop = can$xstart + can$len
can$exon_num = c('1', '2', '3', '4', 'PNS-specific', '5', '6', '7', '8', '9', '10', '11', '12', '13')
can$xmid = (can$xstart + can$xstop)/2

expressed$xplot = as.numeric(NA)
expressed$color = '#00000000' # transparent black default

xpr = expressed[expressed$pos >= min(can$start) & expressed$pos <= max(can$stop), c('pos','brain_rel','xplot','color')]
rownames(xpr) = 1:nrow(xpr)

for (i in 1:nrow(xpr)) {
  for (j in 1:nrow(can)) {
    if (xpr$pos[i] >= can$start[j] && xpr$pos[i] <= can$stop[j]) { # if this pos falls in this exon
      xpr$xplot[i] = xpr$pos[i] - can$start[j] + can$xstart[j]
    }
  }
}

xpr$color = alpha('#000000',xpr$brain_rel)


lof = read.table('data/curation/MAPT-gnomAD-v2.1-ENSG00000186868-2018-11-15-15.38.03.csv',sep=',',header=T)
colnames(lof)[2] = 'pos'
colnames(lof) = gsub('[^a-z0-9_]','_',tolower(colnames(lof)))
lof$xplot = as.integer(NA)

for (i in 1:nrow(lof)) {
  for (j in 1:nrow(can)) {
    if (lof$pos[i] >= can$start[j] - 10 && lof$pos[i] <= can$stop[j] + 10) { # if this pos falls in this exon
      lof$xplot[i] = lof$pos[i] - can$start[j] + can$xstart[j]
    }
  }
}


lof$pos_id = paste(lof$chrom, formatC(lof$pos, width=9, flag='0'), lof$ref, lof$alt, sep='-')
curation$pos_id = curation$variant_id
lof$curated_verdict = curation$form_verdict[match(lof$pos_id, curation$pos_id)]
lof$form_other_notes = curation$form_other_notes[match(lof$pos_id, curation$pos_id)]
lof$cat = ''
lof$cat[lof$flags!=''] = 'loftee'
lof$cat[lof$curated_verdict %in% c('LoF','likely_LoF')] = 'true'
lof$cat[lof$curated_verdict %in% c('not_LoF','likely_not_LoF','uncertain_LoF','')] = 'curated'
lof$color = filter_parms$color[match(lof$cat,filter_parms$cat)]

sum(lof$allele_frequency)

# AC scale: 0, 1, 10, >100 = 0, 1, 2, 3
ac_scale = data.frame(disp=c('0','1','10','100+'),yval=c(0,1,2,3))
ac_ticks = c(0, log10(c(1:10, (2:10)*10)) + 1)
ac_trunc = 100
lof$rel_ac = log10(pmin(lof$allele_count,ac_trunc)) + 1

exon_ybot = -0.5
exon_ytop = 0

plot(NA, NA, xlim=c(0,max(can$xstop)), ylim=c(exon_ybot,max(ac_scale$yval)), yaxs='i', ann=FALSE, axes=FALSE)
abline(h=mean(c(exon_ybot, exon_ytop)))
rect(xleft=can$xstart, xright=can$xstop, ybottom=rep(exon_ytop,nrow(can)), ytop=rep(exon_ybot,nrow(can)), col='#FFFFFF' ,border='#000000')
rect(xleft=xpr$xplot,xright=xpr$xplot+1,ybottom=exon_ybot,ytop=exon_ytop,col=xpr$color,border=NA)
axis(side=2, at=ac_ticks, labels=NA, tck=-0.01)
axis(side=2, at=ac_scale$yval, labels=ac_scale$disp, tck=-0.03, lwd=0, lwd.ticks=1, las=2)
mtext(side=2, line=2.5, text='allele count')
mtext(side=1, line=0, at=can$xmid, text=can$exon_num, font=3, cex=0.8)
mtext(side=1, line=1, text='MAPT exon structure')

segments(x0=lof$xplot, y0=rep(0,nrow(lof)), y1=lof$rel_ac, lwd=3, lend=1, col=lof$color)

scale_xleft = 2000
scale_xright = 2500
scale_factor = (scale_xright - scale_xleft) / 100
scale_ybot = 2.5
scale_ytop = 2.8

par(xpd=T)
rect(xleft=scale_xleft,xright=scale_xright,ybottom=scale_ybot,ytop=scale_ytop)
rect(xleft=scale_xleft+(0:100)*scale_factor, xright=scale_xleft+(1:101)*scale_factor, ybottom=rep(scale_ybot,101), ytop=rep(scale_ytop,101), border=NA, col=alpha('#000000',0:100/100))
text(x=(scale_xright + scale_xleft)/2, y=scale_ytop, pos=3, labels='brain expression', cex=1.5)
text(x=c(scale_xleft,scale_xright),y=c(scale_ybot,scale_ybot),pos=1,labels=c('0%','100%'),cex=1.3,font=3)
par(xpd=F)

mtext(side=3, line=0, text='MAPT', font=3, cex=2)

mtext('b', side=3, cex=2, adj = 0.0, line = 0.3)

lof$gene = 'MAPT'
mapt_supptbl = lof[,c('gene','pos_id','allele_count','cat','flags','curated_verdict','form_other_notes')]
colnames(mapt_supptbl) = c('gene','variant','allele_count','filter_status','LOFTEE_flags','curation_result','comments')



# begin PRNP

lof = read.table('data/curation/prnp_lof.tsv',sep='\t',header=T,quote='')


lof$color = filter_parms$color[match(lof$cat,filter_parms$cat)]

# AC scale: 0, 1, 10, >100 = 0, 1, 2, 3
ac_scale = data.frame(disp=c('0','1','10','100+'),yval=c(0,1,2,3))
ac_ticks = c(0, log10(c(1:10, (2:10)*10)) + 1)
ac_trunc = 100
lof$rel_ac = log10(pmin(lof$ac,ac_trunc)) + 1

cleaved_color = '#777777'
prp_color = '#000000'

exon_ybot = -0.5
exon_ytop = 0

plot(NA,NA,xlim=c(1,254),ylim=c(exon_ybot,3),axes=FALSE,xlab='',ylab='')
axis(side=2, at=ac_ticks, labels=NA, tck=-0.01)
axis(side=2, at=ac_scale$yval, labels=ac_scale$disp, tck=-0.03, lwd=0, lwd.ticks=1, las=2)
mtext(side=2, line=2.5, text='allele count')
axis(side=1,at=c(1,(1:4)*50,253),labels=c(1,(1:4)*50,253),lwd=NA,lwd.ticks=1,cex.axis=.8)
rect(xleft=c(1,23,231),xright=c(23,231,254),ybottom=rep(exon_ybot,3),ytop=rep(exon_ytop,3),col=c(cleaved_color,prp_color,cleaved_color),border=NA)

points(lof$codon,lof$rel_ac,col=lof$color,type='h',lwd=3,lend=1)
mtext(side=1, line=2.5, text='PRNP codon number')

label_y = 2.2
segments(x0=15,x1=135,y0=label_y,lwd=3,lend=1)
segments(x0=140,x1=240,y0=label_y,lwd=3,lend=1)
text(x=(15+135)/2,y=label_y,pos=3,labels='gnomAD non-dementia cohorts',cex=1.5)
text(x=(140+240)/2,y=label_y,pos=3,labels='prion or other dementia cohorts',cex=1.5)
mtext(side=1,at=(1+23)/2,line=0.75,text='signal\npeptide',col=cleaved_color,cex=.8,font=2)
mtext(side=1,at=(231+253)/2,line=0.75,text='GPI\nsignal',col=cleaved_color,cex=.8,font=2)


mtext(side=3, line=0, text='PRNP', font=3, cex=2)

mtext('c', side=3, cex=2, adj = 0.0, line = 0.3)


dev.off()



# Begin Table 2 and Table S2 stuff

# extract length and constraint for Table 2
table2_genes = c('HTT','LRRK2','MAPT','PRNP','SNCA','SOD1')
genes[genes$symbol %in% table2_genes, c('symbol','bp','oe_lof')]
# for MAPT, we want only constitutive exons.
# fortuitously, Ensembl contains a transcript with only constitutive brain-expressed exons:
gstraint_all[gstraint_all$transcript=='ENST00000334239',]

# also compute stats for SOD1
lof = read.table('data/curation/SOD1-gnomAD-v2.1-ENSG00000142168-2018-11-15-15.39.09.csv',sep=',',header=T)
colnames(lof)[2] = 'pos'
colnames(lof) = gsub('[^a-z0-9_]','_',tolower(colnames(lof)))
lof$rel_ac = lof$allele_count / max(lof$allele_count[lof$flags==''])
lof$xplot = as.integer(NA)
lof$pos_id = paste(lof$chrom, formatC(lof$pos, width=9, flag='0'), lof$ref, lof$alt, sep='-')
curation$pos_id = curation$variant_id
lof$curated_verdict = curation$form_verdict[match(lof$pos_id, curation$pos_id)]
lof$form_other_notes = curation$form_other_notes[match(lof$pos_id, curation$pos_id)]
lof$cat = ''
lof$cat[lof$flags!=''] = 'loftee'
lof$cat[lof$curated_verdict %in% c('LoF','likely_LoF')] = 'true'
lof$cat[lof$curated_verdict %in% c('not_LoF','likely_not_LoF','uncertain_LoF','')] = 'curated'
lof$color = filter_parms$color[match(lof$cat,filter_parms$cat)]
# for curated constraint:
obs_sod1 = sum(lof$source=='gnomAD Exomes' & nchar(lof$reference)==1 & nchar(lof$alternate)==1 & lof$cat=='true')
exp_sod1 = gstraint$exp_lof[gstraint$gene=='SOD1']
obs_sod1 / exp_sod1
# for curated CAF:
sum(lof$allele_frequency)
sum(lof$allele_count[lof$cat=='true']) / (2*gnomad)
lof$gene = 'SOD1'
sod1_supptbl = lof[,c('gene','pos_id','allele_count','cat','flags','curated_verdict','form_other_notes')]
colnames(sod1_supptbl) = c('gene','variant','allele_count','filter_status','LOFTEE_flags','curation_result','comments')

# and add supp tbl stuff for SNCA and PRNP
lof = read.table('data/curation/SNCA-gnomAD-v2.1-ENSG00000145335-2018-11-15-15.38.32.csv',sep=',',header=T)
colnames(lof)[2] = 'pos'
colnames(lof) = gsub('[^a-z0-9_]','_',tolower(colnames(lof)))
lof$pos_id = paste(lof$chrom, formatC(lof$pos, width=9, flag='0'), lof$ref, lof$alt, sep='-')
lof$gene = 'SNCA'
lof$cat = 'loftee'
lof$curated_verdict = ''
lof$form_other_notes = ''
snca_supptbl = lof[,c('gene','pos_id','allele_count','cat','flags','curated_verdict','form_other_notes')]
colnames(snca_supptbl) = c('gene','variant','allele_count','filter_status','LOFTEE_flags','curation_result','comments')

lof = read.table('data/curation/PRNP-gnomAD_v2.1_ENSG00000171867_2018_12_14_13_55_02.csv',sep=',',header=T)
colnames(lof)[2] = 'pos'
colnames(lof) = gsub('[^a-z0-9_]','_',tolower(colnames(lof)))
lof$pos_id = paste(lof$chrom, formatC(lof$pos, width=9, flag='0'), lof$ref, lof$alt, sep='-')
lof$gene = 'PRNP'
lof$cat = ''
lof$curated_verdict = ''
lof$form_other_notes = ''
prnp_supptbl = lof[,c('gene','pos_id','allele_count','cat','flags','curated_verdict','form_other_notes')]
colnames(prnp_supptbl) = c('gene','variant','allele_count','filter_status','LOFTEE_flags','curation_result','comments')

# calculate constraint obs/exp for PRNP codon 1-144
synth = read.table('data/curation/prnp_synthetic.table',sep='\t',header=T)
colnames(synth) = gsub('[^a-z0-9_]','_',tolower(colnames(synth)))
synth$codon = (synth$cds_position - 1) %/% 3 + 1

hgvsp_split = strsplit(synth$hgvsp,',')
synth$hgvsp1 = mapply('[[',hgvsp_split,1)
hgvsp1_change = strsplit(synth$hgvsp1,'\\.')
synth$hgvsp1_change = mapply('[',hgvsp1_change,3)
synth$hgvsp1_change[grepl('>',synth$hgvsp1_change)] = NA
synth$aa_change = tla_to_ola(synth$hgvsp1_change)
syn = nchar(synth$amino_acids)==1 & !is.na(synth$amino_acids)
synth$aa_change[syn] = paste(synth$amino_acids[syn], synth$codon[syn], synth$amino_acids[syn], sep='')

murates = read.table('data/constraint/fordist_1KG_mutation_rate_table.txt',header=T)
colnames(murates) = c('fromcontext','tocontext','mu_snp')
murates$to_base = substr(murates$tocontext,2,2)

synth_mu = sqldf("
                 select   s.chrom, s.pos, s.ref, s.alt, s.codon, s.context, s.aa_change, m.fromcontext, m.tocontext, m.mu_snp
                 from     synth s, murates m
                 where    s.lof = 'HC'
                 and      s.context = m.fromcontext
                 and      s.alt = m.to_base
                 order by 1, 2, 3
                 ;")

sum_mu_253 = sum(synth_mu$mu_snp)
sum_mu_144 = sum(synth_mu$mu_snp[!is.na(synth_mu$codon) & synth_mu$codon >= 1 & synth_mu$codon < 144])
exp_lof_253 = gstraint$exp_lof[gstraint$gene=='PRNP']
exp_lof_144 = exp_lof_253 * sum_mu_144 / sum_mu_253
obs_lof_144 = 6 # manually curated - R37X, Q41X, Q75X, W81X, W99X, and G131X
obs_lof_144 / exp_lof_144

# rbind everything for Table S2
supptbl = rbind(htt_supptbl, mapt_supptbl, prnp_supptbl, snca_supptbl, sod1_supptbl)

# write out Table S2
write.table(supptbl, 'data/output/table_s2.tsv',sep='\t',na='',row.names=F,col.names=T,quote=F)

# End Table 2 and Table S2 stuff




### Begin Figure S1
png('figures/figures1.png',width=2400,height=3600,res=500)

par(mfrow=c(3,1))

# decide on right histogram breaks
lowest_p = min(genes$p_dsamp_global[genes$p_dsamp_global > 0])
lowest_p ^ 2 # 1.58e-11
# so the axis should break at 1e-11 or lower
axis_break = 10^-11.9
# confirm that the leftmost bin (<1e-10) is composed entirely of genes with 0 CAF, even for homs
sum(genes$p_dsamp_global^2 < 1e-11) == sum(genes$p_dsamp_global == 0)
# but the next bin this is not the case
sum(genes$p_dsamp_global^2 < 3.16e-11) == sum(genes$p_dsamp_global == 0)

caf_breaks = c(0,10^seq(-12,0,by=.5))
n_bins = length(caf_breaks) - 1

ymax = 11000
label_y = 4500

world_population = 7.66e9
#gnomad = 141456 # defined above

het_k = '#FF9912'
hom_k = '#984EA3'

het_hist = hist(2*genes$p_dsamp_global*(1-genes$p_dsamp_global), breaks=caf_breaks, plot=FALSE)
hom_hist = hist(genes$p_dsamp_global ^ 2, breaks=caf_breaks, plot=FALSE)

par(mar=c(2,5,4,3))
plot(NA, NA, xlim=c(-13,0), ylim=c(0,ymax), ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
axis(side=1, at=c(-13,-9,-6,-3,0), labels=c('','1 in 1 billion','1 in 1 million','1 in 1,000','100%'), tck=-0.03, lwd=1, lwd.ticks=1)
axis(side=1, at=c(-12.3), labels=c('zero'), lwd=0, lwd.ticks=0)
axis(side=1, at=-12:0, labels=NA, lwd=0, lwd.ticks=1, tck=-0.01)
axis.break(axis=1,breakpos=log10(axis_break),style='zigzag')
axis(side=2, at=(0:(ymax/1000))*1000, labels=formatC((0:(ymax/1000))*1000,format='d',big.mark=','),las=2)
mtext(side=2, line=3.5, text='number of genes')
rect(xleft=log10(het_hist$mids[1])-.25, xright=log10(het_hist$mids[1])+.25, ybottom=0, ytop=het_hist$counts[1], lwd=3, border=het_k, col=alpha(het_k,.2))
points(log10(het_hist$mids[2:n_bins]), het_hist$counts[2:n_bins], type='l', lwd=3, col=het_k)
polygon(x=c(log10(het_hist$mids[2:n_bins]),1), y=c(het_hist$counts[2:n_bins],0), border=NA, col=alpha(het_k,.2))
rect(xleft=log10(hom_hist$mids[1])-.25, xright=log10(hom_hist$mids[1])+.25, ybottom=0, ytop=hom_hist$counts[1], lwd=3, border=hom_k, col=alpha(hom_k,.2))
points(log10(hom_hist$mids[2:n_bins]), hom_hist$counts[2:n_bins], type='l', lwd=3, col=hom_k)
polygon(x=c(log10(hom_hist$mids[2:n_bins]),1), y=c(hom_hist$counts[2:n_bins],0), border=NA, col=alpha(hom_k,.2))
segments(x0=log10(1/world_population),y0=0,y1=label_y,lwd=0.5)
segments(x0=log10(1/gnomad),y0=0,y1=label_y,lwd=0.5)
text(x=log10(1/world_population),y=label_y,pos=3,labels='1\non\nEarth')
text(x=log10(1/gnomad),y=label_y,pos=3,labels='1\nin\ngnomAD')
par(xpd=T)
legend(x=-4,y=10000,legend=c('heterozygotes','homozygotes &\ncompound hets'),text.col=c(het_k,hom_k),border=c(het_k,hom_k),fill=alpha(c(het_k,hom_k),.2),bty='n',cex=1.2)
par(xpd=F)
mtext('a', side=3, cex=2, adj = 0.0, line = 0.3)

finland_population = 5.50e6
gnomad_finns = 12526

het_hist = hist(2*genes$p_dsamp_finn*(1-genes$p_dsamp_finn), breaks=caf_breaks, plot=FALSE)
hom_hist = hist(genes$p_dsamp_finn ^ 2, breaks=caf_breaks, plot=FALSE)

par(mar=c(3,5,3,3))

plot(NA, NA, xlim=c(-13,0), ylim=c(0,ymax), ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
axis(side=1, at=c(-13,-9,-6,-3,0), labels=c('','1 in 1 billion','1 in 1 million','1 in 1,000','100%'), tck=-0.03, lwd=1, lwd.ticks=1)
axis(side=1, at=c(-12.3), labels=c('zero'), lwd=0, lwd.ticks=0)
axis(side=1, at=-12:0, labels=NA, lwd=0, lwd.ticks=1, tck=-0.01)
axis.break(axis=1,breakpos=log10(axis_break),style='zigzag')
axis(side=2, at=(0:(ymax/1000))*1000, labels=formatC((0:(ymax/1000))*1000,format='d',big.mark=','),las=2)
mtext(side=2, line=3.5, text='number of genes')
rect(xleft=log10(het_hist$mids[1])-.25, xright=log10(het_hist$mids[1])+.25, ybottom=0, ytop=min(het_hist$counts[1],ymax), lwd=3, border=het_k, col=alpha(het_k,.2))
points(log10(het_hist$mids[2:n_bins]), het_hist$counts[2:n_bins], type='l', lwd=3, col=het_k)
polygon(x=c(log10(het_hist$mids[2:n_bins]),1), y=c(het_hist$counts[2:n_bins],0), border=NA, col=alpha(het_k,.2))
rect(xleft=log10(hom_hist$mids[1])-.25, xright=log10(hom_hist$mids[1])+.25, ybottom=0, ytop=min(hom_hist$counts[1],ymax), lwd=3, border=hom_k, col=alpha(hom_k,.2))
points(log10(hom_hist$mids[2:n_bins]), hom_hist$counts[2:n_bins], type='l', lwd=3, col=hom_k)
polygon(x=c(log10(hom_hist$mids[2:n_bins]),1), y=c(hom_hist$counts[2:n_bins],0), border=NA, col=alpha(hom_k,.2))
segments(x0=log10(1/finland_population),y0=0,y1=label_y,lwd=0.5)
segments(x0=log10(1/gnomad_finns),y0=0,y1=label_y,lwd=0.5)
text(x=log10(1/finland_population),y=label_y,pos=3,labels='1\nin\nFinland')
text(x=log10(1/gnomad_finns),y=label_y,pos=3,labels='1\nin\ngnomAD Finns')

mtext('b', side=3, cex=2, adj = 0.0, line = 0.3)

a = 0.05766 # mean % genome autozygous in ELGH according to Hilary Martin's email 2018-12-06

het_hist = hist(2*genes$p_dsamp_global*(1-genes$p_dsamp_global), breaks=caf_breaks, plot=FALSE)
hom_hist = hist((1-a)*(genes$p_dsamp_global ^ 2) + a*genes$p_dsamp_global, breaks=caf_breaks, plot=FALSE)

world_consang_population = world_population * .104 # Bittles & Black 2009 PNAS
gnomad_consang = 2912 # SAS w/ F > 0.05 according to Konrad on 2018-12-04 via Slack

par(mar=c(4,5,2,3))

plot(NA, NA, xlim=c(-13,0), ylim=c(0,ymax), ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
axis(side=1, at=c(-13,-9,-6,-3,0), labels=c('','1 in 1 billion','1 in 1 million','1 in 1,000','100%'), tck=-0.03, lwd=1, lwd.ticks=1)
axis(side=1, at=c(-12.3), labels=c('zero'), lwd=0, lwd.ticks=0)
axis(side=1, at=-12:0, labels=NA, lwd=0, lwd.ticks=1, tck=-0.01)
axis.break(axis=1,breakpos=log10(axis_break),style='zigzag')
mtext(side=1, line=2.5, text='expected frequency of pLoF individuals in population')
axis(side=2, at=(0:(ymax/1000))*1000, labels=formatC((0:(ymax/1000))*1000,format='d',big.mark=','),las=2)
mtext(side=2, line=3.5, text='number of genes')
rect(xleft=log10(het_hist$mids[1])-.25, xright=log10(het_hist$mids[1])+.25, ybottom=0, ytop=het_hist$counts[1], lwd=3, border=het_k, col=alpha(het_k,.2))
points(log10(het_hist$mids[2:n_bins]), het_hist$counts[2:n_bins], type='l', lwd=3, col=het_k)
polygon(x=c(log10(het_hist$mids[2:n_bins]),1), y=c(het_hist$counts[2:n_bins],0), border=NA, col=alpha(het_k,.2))
rect(xleft=log10(hom_hist$mids[1])-.25, xright=log10(hom_hist$mids[1])+.25, ybottom=0, ytop=hom_hist$counts[1], lwd=3, border=hom_k, col=alpha(hom_k,.2))
points(log10(hom_hist$mids[2:n_bins]), hom_hist$counts[2:n_bins], type='l', lwd=3, col=hom_k)
polygon(x=c(log10(hom_hist$mids[2:n_bins]),1), y=c(hom_hist$counts[2:n_bins],0), border=NA, col=alpha(hom_k,.2))
segments(x0=log10(1/world_consang_population),y0=0,y1=label_y,lwd=0.5)
segments(x0=log10(1/gnomad_consang),y0=0,y1=label_y,lwd=0.5)
text(x=log10(1/world_consang_population),y=label_y,pos=3,labels='1\nconsanguineous\nworldwide')
text(x=log10(1/gnomad_consang),y=label_y,pos=3,labels='1\nconsanguineous\nin gnomAD')

mtext('c', side=3, cex=2, adj = 0.0, line = 0.3)

dev.off()


### Table 3

# Table 3 (methods) - get counts of each gene list 
summary(genes)

# also get N for each gene list with non-missing constraint values, for in-line text references
summary(genes[!is.na(genes$oe_lof),])

# For Methods - proportion of gnomAD exomes that is bottlenecked or consanguineous 
asj_exomes = 5040
fin_exomes = 10824
all_exomes = 125748
consang_exomes = 2912
(asj_exomes + fin_exomes) / all_exomes
consang_exomes / all_exomes


