#### How many E-P pairs are supported by chia-pet interactions ########

get_perm_EP <- function(E_P,rg.cons.ep,dreg.bs,win.size.ep,mc.cores,...){
	perm_E_P <- mclapply(seq_len(length(E_P)),function(i) get_single_perm_EP(i,names(E_P)[i],E_P,rg.cons.ep,dreg.bs,win.size.ep),mc.cores=mc.cores)
	names(perm_E_P) <- names(E_P)
	return(perm_E_P)
}

get_single_perm_EP <- function(i,gene_id,E_P,rg.cons.ep,dreg.bs,win.size.ep,...){
	if(i%%100==0){print(paste("$$$$$$$$$$$$$$$$$$$$$$",i))}
	gene.bs <- rg.cons.ep[gene_id]
	chr <- as.character(seqnames(gene.bs))
	nedges <- length(E_P[[gene_id]])
	orig.enh <- E_P[[gene_id]]
	cand.enh <- which(seqnames(dreg.bs)==chr)
	cand.enh <- setdiff(cand.enh,orig.enh)
	dreg.bs.cand <- c(dreg.bs[orig.enh],dreg.bs[cand.enh])
	shift.size <- (end(dreg.bs.cand)-start(dreg.bs.cand))/2
	dreg.bs.cand <- shift(dreg.bs.cand,shift=shift.size)
	gene.start <- start(gene.bs)
	if(as.character(strand(gene.bs))=='-')
		gene.start <- end(gene.bs)
	X <- as.numeric(gene.start-start(dreg.bs.cand))
	names(X) <- c(orig.enh,cand.enh)
	Tr <- rep(0,times=length(dreg.bs.cand))
	Tr[1:nedges] <- 1
	names(Tr) <- c(orig.enh,cand.enh)
	M <- nedges
	Y <- c(orig.enh,cand.enh)
	rr  <- Match(Y=Y,Tr=Tr,X=X,M=M,version="fast")
	set.id <- unique(rr$index.treated)
	return(as.integer(names(X[sapply(set.id,function(i) sample(rr$index.control[which(rr$index.treated==i)],1))])))
}

count_overlap_with_chia <- function(ind,enh.pos,rg.cons,chia.1.gr,chia.2.gr,win.size=500,...){
	if(ind%%100==0){print(paste("$$$$$$$$$$$$$$$$$$$$$$",ind))}
	shift.size <- (end(enh.pos)-start(enh.pos))/2
	enh.pos <- GenomicRanges::shift(enh.pos,shift=shift.size)
	enh.pos <- promoters(enh.pos,upstream=win.size,downstream=win.size)
	rg.cons <- promoters(rg.cons,upstream=win.size,downstream=win.size)
	hit.within <- findOverlaps(rg.cons,chia.1.gr,ignore.strand=T)
	sub.inds <- c()
	if(length(subjectHits(hit.within))!=0){
		for(sub in subjectHits(hit.within)){
			hit.within.2 <- findOverlaps(chia.2.gr[sub],enh.pos,ignore.strand=T)
			sub.inds <- union(sub.inds,subjectHits(hit.within.2))
		}
	}
	hit.within <- findOverlaps(rg.cons,chia.2.gr,ignore.strand=T)
	if(length(subjectHits(hit.within))!=0){
		for(sub in subjectHits(hit.within)){
			hit.within.2 <- findOverlaps(chia.1.gr[sub],enh.pos,ignore.strand=T)
			sub.inds <- union(sub.inds,subjectHits(hit.within.2))
		}
	}
	return(length(sub.inds))
}

find_overlap_with_eQTL_snp <- function(ind,enh.pos,snp.pos,win.size=500,...){
	if(ind%%100==0){print(paste("$$$$$$$$$$$$$$$$$$$$$$",ind))}
	if(length(snp.pos)==0)
		return(list(total_snp_hit=0,edge_hit=0))
	shift.size <- (end(enh.pos)-start(enh.pos))/2
	enh.pos <- GenomicRanges::shift(enh.pos,shift=shift.size)
	enh.pos <- promoters(enh.pos,upstream=win.size,downstream=win.size)
	hit.within <- findOverlaps(enh.pos,snp.pos,ignore.strand=T)
	s_hit <- subjectHits(hit.within)
	q_hit <- queryHits(hit.within)
	u_q_hit <- unique(q_hit)
	total_snp_hit <- 0
	for(j in u_q_hit){
		match.id <- which(q_hit==j)
		total_snp_hit <- total_snp_hit+length(unique(snp.pos[s_hit[match.id]]))
	}
	return(list(total_snp_hit=total_snp_hit,edge_hit=length(u_q_hit)))
}

################################################################################
## overlap with snps
################################################################################

get_gene_tss <- function(i,rg.cons.prom,tss.bs){
	if(i%%100==0){print(paste("$$$$$$$$$$$$$$$$$$$$$$", i))}
	gene_list <- list()
	flag<-F
	if(!any(is.na(unlist(rg.cons.prom$overlap_TSS_1kb_up)))){
		gene_list <- tss.bs[unlist(rg.cons.prom$overlap_TSS_1kb_up)]
		flag=T
	}
	if(flag)
		return(gene_list)
	else
		return(NA)
}

find_gene_to_snp_ind <- function(i,gene_list,snp.ens){
	if(i%%100==0){print(paste("$$$$$$$$$$$$$$$$$$$$$$",i))}
	ind <- c()
	if(any(is.na(gene_list)))
		return(ind)
	for(gene in gene_list){
		ind <- c(ind,which(snp.ens==gene))
	}
	return(ind)
}

getExtValEmpPval <- function(count_overlap_mat,ep_links,...){
	real_ov <- sum(count_overlap_mat[,1])/length(unlist(ep_links))
	print(real_ov)
	emp.p.v <- 0
	nperms = dim(count_overlap_mat)[2]-1
	for( perm in seq_len(nperms)){
		rand_ov <- sum(count_overlap_mat[,1+perm])/length(unlist(ep_links))
		print(rand_ov)
		if(rand_ov>real_ov)
			emp.p.v <- emp.p.v+1
	}
	print(paste0("Empirical P-value: ",emp.p.v/nperms))
	return(emp.p.v/nperms)
}


############################ Integrative analysis functions #########################

getEpLinksForGene <- function(gene_id, gene_type="entrez", ep_links, prom.bs, gene_id_map, data_type="groseq", tss.bs=NULL,...){
	
	supp_gene_types <- c("entrez","ensembl","symbol")
	if(!is.element(gene_type,supp_gene_types)){
		print("Gene type not supported")
		print("Only entrez, ensembl , and symbol IDs are supported")
		return(NA)
	}
	supp_data_types <- c("groseq","fantom","encode","roadmap")
	if(!is.element(data_type,supp_data_types)){
		print("Data type not supported")
		print("Only groseq, fantom, encode, and roadmap data types are supported")
		return(NA)
	}
	
	if(data_type=="groseq"){
		if(gene_type=="symbol"){
			map <- gene_id_map$ens_sym
			ind <- match(gene_id,map$sym_id, nomatch = NA)
			if(is.na(ind)){ return(NA) }
			gene_id <- as.character(map$ens_id[ind])
			gene_type <- "ensembl"
		}
		if(gene_type=="ensembl"){
			map <- gene_id_map$ens_ent
			ind <- match(gene_id,map$ens_id, nomatch = NA)
			if(is.na(ind)){return(NA)}
			gene_id <- as.character(map$ent_id[ind])
			gene_type=="entrez"
		}
		if(is.element(gene_id,names(ep_links))){
			return(ep_links[gene_id])
		}else{ return(NA) }
	}
	
	if(data_type=="fantom"){
		if(gene_type=="symbol"){
			map <- gene_id_map$ens_sym
			ind <- match(gene_id,map$sym_id, nomatch = NA)
			if(is.na(ind)){ return(NA) }
			gene_id <- as.character(map$ens_id[ind])
			gene_type <- "ensembl"
		}
		if(gene_type=="ensembl"){
			map <- gene_id_map$ens_ent
			ind <- match(gene_id,map$ens_id, nomatch = NA)
			if(is.na(ind)){ return(NA) }
			gene_id <- as.character(map$ent_id[ind])
			gene_type=="entrez"
		}
		
		gene_list <- prom.bs$gene_id
		ind <- sapply(gene_list,function(gl) is.element(gene_id,gl))
		if(all(ind==F)){ return(NA) }
		in.ep <- match(prom.bs$name[ind],names(ep_links),nomatch=NA)
		in.ep <- in.ep[!is.na(in.ep)]
		if(length(in.ep)==0){ return(NA) }
		return(ep_links[in.ep])
	}
	
	if(data_type=="roadmap" | data_type=="encode"){
		if(gene_type=="symbol"){
			map <- gene_id_map$ens_sym
			ind <- match(gene_id,map$sym_id, nomatch = NA)
			if(is.na(ind)){ return(NA) }
			gene_id <- as.character(map$ens_id[ind])
			gene_type <- "ensembl"
		}
		if(gene_type=="entrez"){
			map <- gene_id_map$ens_ent
			ind <- match(gene_id,map$ent_id, nomatch = NA)
			if(is.na(ind)){ return(NA) }
			gene_id <- as.character(map$ens_id[ind])
			gene_type <- "ensembl"
		}
		pos <- which(tss.bs$ens_id==gene_id)
		if(length(pos)==0){ return(NA) }
		gene_list <- prom.bs$overlap_TSS_1kb_up
		ind <- sapply(gene_list,function(gl) if(is.na(gl)){return(F)}else{return(length(intersect(pos,gl))!=0)})
		if(all(ind==F)){ return(NA) }
		in.ep <- match(prom.bs$name[ind],names(ep_links),nomatch=NA)
		in.ep <- in.ep[!is.na(in.ep)]
		if(length(in.ep)==0){ return(NA) }
		return(ep_links[in.ep])		
	}
	return(NA)
}

find_chia_annot <- function(enh.pos,rg.cons,chia.1.gr,chia.2.gr,win.size=500,...){
	shift.size <- (end(enh.pos)-start(enh.pos))/2
	enh.pos <- GenomicRanges::shift(enh.pos,shift=shift.size)
	enh.pos <- promoters(enh.pos,upstream=win.size,downstream=win.size)
	rg.cons <- promoters(rg.cons,upstream=win.size,downstream=win.size)
	hit.within <- findOverlaps(rg.cons,chia.1.gr,ignore.strand=T)
	sub.inds <- c()
	chia.inds <- c()
	if(length(subjectHits(hit.within))!=0){
		for(sub in subjectHits(hit.within)){
			hit.within.2 <- findOverlaps(chia.2.gr[sub],enh.pos,ignore.strand=T)
			sub.inds <- union(sub.inds,subjectHits(hit.within.2))
			if(length(queryHits(hit.within.2))!=0)
				chia.inds <- union(chia.inds,sub)
		}
	}
	hit.within <- findOverlaps(rg.cons,chia.2.gr,ignore.strand=T)
	if(length(subjectHits(hit.within))!=0){
		for(sub in subjectHits(hit.within)){
			hit.within.2 <- findOverlaps(chia.1.gr[sub],enh.pos,ignore.strand=T)
			sub.inds <- union(sub.inds,subjectHits(hit.within.2))
			if(length(queryHits(hit.within.2))!=0)
				chia.inds <- union(chia.inds,sub)
		}
	}
	return(list(cov_ep=length(sub.inds),chia_ind=chia.inds))
}

find_eqtl_annot <- function(enh.pos,snp.pos,win.size=500,...){
	if(length(snp.pos)==0)
		return(list(total_snp_hit=0,edge_hit=0))
	shift.size <- (end(enh.pos)-start(enh.pos))/2
	enh.pos <- GenomicRanges::shift(enh.pos,shift=shift.size)
	enh.pos <- promoters(enh.pos,upstream=win.size,downstream=win.size)
	hit.within <- findOverlaps(enh.pos,snp.pos,ignore.strand=T)
	s_hit <- subjectHits(hit.within)
	q_hit <- queryHits(hit.within)
	u_q_hit <- unique(q_hit)
	total_snp_hit <- 0
	res <- list()
	for(j in u_q_hit){
		match.id <- which(q_hit==j)
		total_snp_hit <- total_snp_hit+length(unique(snp.pos[s_hit[match.id]]))
		res[[as.character(j)]] <- snp.pos[s_hit[match.id]]
	}
	return(res)
}





