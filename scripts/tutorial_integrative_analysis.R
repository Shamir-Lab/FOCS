##################################################################################################################################
################################### Integrative analysis of E-P links ############################################################
##################################################################################################################################


setwd('../focs/')

# load libraries
bioconductor_libs = c('GenomicRanges','rtracklayer','GenomicInteractions')
for (p in (c(bioconductor_libs))){library(p,character.only=T)}


# Variables
ext.win.size <- 500 #maximum window size upstream/downstream to intersect TSSs/enhancers with eQTLs or ChIA-PET anchors


# load objects and read files
data_directory = '../data/'
script_directory = '../scripts/'
tmp_directory = '../tmp/'
# data types - groseq, fantom, encode, roadmap
data_type <- 'groseq'

# Enhancer genomic positions
enh.bs <- get(load(paste(data_directory,data_type,'.enh.pos.RData',sep='')))
# Promoters positions
prom.bs <- get(load(paste(data_directory,data_type,'.prom.pos.RData',sep='')))
# Sample annotations
sample.annot <- get(load(paste(data_directory,data_type,'.sample.annot.RData',sep='')))
col2type = as.character(sample.annot$cell_type)
# Enhancer-Promoter links
ep_links <- get(load(paste(data_directory,data_type,'.E_P.links.RData',sep='')))
# Gene id map table
gene_id_map <- get(load(paste(data_directory,'gene.id.map.RData',sep='')))
# Get gencodeV10 TSS positions
tss.bs <- get(load(file=paste(data_directory,'tss.bs.genecode10.RData',sep='')))
# Get a map from gencodeV10 TSSs to associated SNPs
snp.inds.gencode <- get(load(file=paste(data_directory,'genecodev10.snp.inds.RData',sep='')))
# Get a map from fantom promoter-TSSs to associated SNPs
snp.inds.fantom <- get(load(file=paste(data_directory,'fantom.snp.inds.RData',sep='')))


#GTEx eQTLs data
sig_pairs.gr <- get(load(file=paste(data_directory,'sig_pairs.gr.RData',sep='')))
var_rs <- get(load(file=paste(data_directory,'var_to_rs.RData',sep='')))


#ChIA-PET data
chia.1.gr <- get(load(file=paste(data_directory,'chia.1.gr.RData',sep='')))
chia.2.gr <- get(load(file=paste(data_directory,'chia.2.gr.RData',sep='')))

#load code
source(paste(script_directory,'external_val_functions.R',sep=''))
source(paste(script_directory,'io_functions.R',sep=''))

##################################################################################################################################
########################################## Part 1: Get E-P links for given gene id ###############################################
##################################################################################################################################

gene_id = "AGRN"
# get the enhancer indices from the E-P links
ep_gene <- getEpLinksForGene(gene_id, gene_type="symbol", ep_links, prom.bs, gene_id_map, data_type=data_type, tss.bs=tss.bs)

# write the genomic positions of the E-P links in bed format - <chr gene>\t<start>\t<end><gene id>\t<strand><chr enh>\t<start>\t<end>\t<strand>
outfile <- paste(tmp_directory,'single_gene_ep_links.bed',sep='')
output_E_P_links_bed(ep_gene,enh.bs,prom.bs,outfile)

##################################################################################################################################
######################################## Part 2: Write genomic positions of enhancers/promoters ##################################
##################################################################################################################################

# write BED12 interaction format for UCSC genome browser
# If you want to write all ep_links then switch ep_gene with ep_links
id_genes <- match(names(ep_gene),names(prom.bs))
id_genes <- match("p_chr1:954780_954930",names(prom.bs))
id_rep <- rep(id_genes,times=sapply(ep_gene,length))
id_rep <- rep(id_genes,times=sapply(ep_links["p_chr1:954780_954930"],length))
prom.ep <- prom.bs[id_rep]
enh.ep <- enh.bs[unlist(ep_gene)]
enh.ep <- enh.bs[unlist(ep_links["p_chr1:954780_954930"])]


outfile <- paste(tmp_directory,'single_gene_ep_links.bed',sep='')

# Export the interactions and get 'GenomicInteractions' object
# Because sometimes enhancer and promoter regions may intersect we limit their region promoter/enhancer widths to 200/100 bp (prom.win.size=100;enh.win.size=100 upstream/downstream)
# Selected promoter regions are by default with respect to the start position (gene.start="start"), otherwise (gene.start="mid") they will be selected w.r.t the center positions
# Enhancer regions (+/-500 bp) are always relative to their center position
gint <- exportBedInteract(prom.ep,enh.ep,names=NA,track_name=paste(data_type," genomic interactions",sep=''),outfile,gene.start="mid",prom.win.size=100,enh.win.size=50,to.export=T)
# export the original region widths of the enhancers and promoters
#gint <- exportBedChIAInteract(prom.ep,enh.ep,names=NA,track_name=paste(data_type," genomic interactions",sep=''),outfile,to.export=T)

##################################################################################################################################
##################################### Part 3: Annotate E-P links with eQTLS/ChIA-PET interactions ################################
##################################################################################################################################

# We show how to annotate E-P links for a given gene ID (from part 1)

####################################################### GRO-seq Only #############################################################

# Annotate with ChIA-PET interactions
# Note that some ChIA-PET anchors may overlap each other. Therefore, the visualization in UCSC will look like as a single black segment

chia_annot <- sapply(names(ep_gene),function(g) find_chia_annot(enh.bs[ep_gene[[g]]],prom.bs[g],chia.1.gr,chia.2.gr,win.size=500))
chia.1.ep <- chia.1.gr[chia_annot[,1]$chia_ind]
chia.2.ep <- chia.2.gr[chia_annot[,1]$chia_ind]
outfile <- paste(tmp_directory,'single_gene_chiapet.bed',sep='')
gint <- exportBedChIAInteract(chia.1.ep,chia.2.ep,names=NA,track_name="ChIA-PET genomic interactions",outfile,to.export=T)

# Annotate with GTEx eQTLs
# Get the eQTL SNP positions that overlap with the enhancer positions
eqtl_annot <- sapply(names(ep_gene),function(g) find_eqtl_annot(enh.bs[ep_gene[[g]]],
									unique(sig_pairs.gr[sig_pairs.gr$ent_id==g]),win.size=500))

# What are the rs_id (dbSNP142) of the overlapping SNPs
# Note that eqtl_annot object is assumed to be a matrix. In the case of a list, the code below should be modified.
enh_ids <- as.integer(rownames(eqtl_annot))
g <- names(ep_gene)
df <- data.frame()
for(i in enh_ids){
	snps <- eqtl_annot[as.character(i),1][[1]]
	var_id <- as.character(snps$variant_id)
	id.match <- match(var_id,var_rs[,1])
	rs_id <- as.character(var_rs[id.match,3])
	enh_pos <- enh.bs[unlist(ep_gene)[i]]
	enh_str <- paste(seqnames(enh_pos),start(enh_pos),end(enh_pos),sep="_")
	df <- rbind(df,data.frame(gene_id=rep(g,length(snps)),enh_pos=enh_str,var_id,rs_id))
	print(paste(rep(g,length(snps)),enh_str,var_id,rs_id,sep=' '))	
}

# Write results
outfile <- paste(tmp_directory,'enh_eQTL_annotation.txt',sep='')
write.table(df,file=outfile,col.names=T,row.names=F,sep="\t",quote=F)


####################################################### FANTOM5 Only #############################################################

# Annotate with ChIA-PET interactions
# Note that some ChIA-PET anchors may overlap each other. Therefore, the visualization in UCSC will look like as a single black segment

chia_annot <- sapply(names(ep_gene),function(g) find_chia_annot(enh.bs[ep_gene[[g]]],prom.bs[g],chia.1.gr,chia.2.gr,win.size=500))
chia.1.ep <- chia.1.gr[unlist(chia_annot[2,])]
chia.2.ep <- chia.2.gr[unlist(chia_annot[2,])]
outfile <- paste(tmp_directory,'single_gene_chiapet.bed',sep='')
gint <- exportBedChIAInteract(chia.1.ep,chia.2.ep,names=NA,track_name="ChIA-PET genomic interactions",outfile,to.export=T)

# Annotate with GTEx eQTLs

# get the reduced promoter genomic positions according to ep_links
prom.bs.ep <- prom.bs[match(names(ep_gene),prom.bs$name)]
gene_ep_list <- prom.bs.ep$gene_id

# Get the eQTL SNP positions that overlap with the enhancer positions
eqtl_annot <- sapply(seq_len(length(ep_gene)),function(i) find_eqtl_annot(enh.bs[ep_gene[[i]]],
									unique(sig_pairs.gr[unlist(snp.inds.fantom[gene_ep_list[[i]]])]),win.size=500))

# What are the rs_id (dbSNP142) of the overlapping SNPs - explained in GRO-seq

################################################## ENCODE/Roadmap data type only ################################################

# Annotate with ChIA-PET interactions
# Note that some ChIA-PET anchors may overlap each other. Therefore, the visualization in UCSC will look like as a single black segment

chia_annot <- sapply(names(ep_gene),function(g) find_chia_annot(enh.bs[ep_gene[[g]]],prom.bs[g],chia.1.gr,chia.2.gr,win.size=500))
chia.1.ep <- chia.1.gr[unlist(chia_annot[2,])]
chia.2.ep <- chia.2.gr[unlist(chia_annot[2,])]
outfile <- paste(tmp_directory,'single_gene_chiapet.bed',sep='')
gint <- exportBedChIAInteract(chia.1.ep,chia.2.ep,names=NA,track_name="ChIA-PET genomic interactions",outfile,to.export=T)

# Annotate with GTEx eQTLs

# get the reduced promoter genomic positions according to ep_links
prom.bs.ep <- prom.bs[match(names(ep_gene),prom.bs$name)]
gene_ep_list <- sapply(seq_len(length(prom.bs.ep)),function(i) get_gene_tss(i,prom.bs.ep[i],tss.bs))

# Get the eQTL SNP positions that overlap with the enhancer positions
eqtl_annot <- sapply(seq_len(length(ep_gene)),function(i) find_eqtl_annot(enh.bs[ep_gene[[i]]],
									unique(sig_pairs.gr[unlist(snp.inds.gencode[gene_ep_list[[i]]$ens_id])]),win.size=500))

# What are the rs_id (dbSNP142) of the overlapping SNPs - explained in GRO-seq







