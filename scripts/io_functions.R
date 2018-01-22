##################################################################################################################################
########################################## Input/Output functions ################################################################
##################################################################################################################################

output_E_P_links_bed <- function(E_P,enh.bs,gene.bs,outfile,...){
	res <- c()
	gene_names <- names(E_P)
	pos_in_gene_bs <- match(gene_names,names(gene.bs))
	for(i in seq_len(length(pos_in_gene_bs))){
		if(i%%100==0){print(paste("$$$$$$$$$$$$$$$$$$$$$$",i))}
		gname <- gene_names[i]
		gpos <- pos_in_gene_bs[i]
		gstr <- paste(seqnames(gene.bs[gpos]),start(gene.bs[gpos]),end(gene.bs[gpos]),gname,strand(gene.bs[gpos]),sep="\t")
		enh.sub.bs <- E_P[[i]]
		for(j in enh.sub.bs){
			estr <- paste(seqnames(enh.bs[j]),start(enh.bs[j]),end(enh.bs[j]),'*',sep="\t")
			res <- c(res,paste(gstr,estr,sep="\t"))
		}
	}
	con <- file.path(outfile)
	writeLines(res,con=con)

}

output_E_P_links_bed_fantom <- function(E_P,enh.bs,gene.bs,outfile,...){
	res <- c()
	gene_names <- names(E_P)
	pos_in_gene_bs <- match(gene_names,gene.bs$name)
	for(i in seq_len(length(pos_in_gene_bs))){
		if(i%%100==0){print(paste("$$$$$$$$$$$$$$$$$$$$$$",i))}
		gpos <- pos_in_gene_bs[i]
		gname <- paste(unlist(gene.bs[gpos]$gene_id),collapse=";")
		gname_ref <- paste(unlist(gene.bs[gpos]$tx_name),collapse=";")
		gstr <- paste(seqnames(gene.bs[gpos]),start(gene.bs[gpos]),end(gene.bs[gpos]),gname,gname_ref,strand(gene.bs[gpos]),sep="\t")
		enh.sub.bs <- E_P[[i]]
		for(j in enh.sub.bs){
			estr <- paste(seqnames(enh.bs[j]),start(enh.bs[j]),end(enh.bs[j]),'*',sep="\t")
			res <- c(res,paste(gstr,estr,sep="\t"))
		}
	}
	con <- file.path(outfile)
	writeLines(res,con=con)

}

exportBedGraph <- function(bs, filename, ignore.strand=T,...){
	con <- file.path(filename)
	if(!ignore.strand)
		export.bedGraph(bs, con, ignore.strand = ignore.strand, trackLine = new("BasicTrackLine",colorByStrand = col2rgb(c("blue","red"))))
	else
		export.bedGraph(bs, con, ignore.strand = ignore.strand, trackLine = new("BasicTrackLine"))
}

exportWiggle <- function(reads, normCounts=1L, filename, strand="+",...){
	writeWiggle(reads=reads, file=filename, fileType="wig", strand=strand, normCounts=normCounts, reverse=FALSE)
}

exportBedInteract <- function(rg.cons.ep,enh.ep,names=NA,track_name="Genomic interactions",filename,gene.start="start",prom.win.size=100,enh.win.size=50,to.export=T,supp_ep=NA,...){
	
	gint_tmp <- GenomicInteractions(enh.ep, rg.cons.ep, counts=as.integer(seq_len(length(enh.ep))))
	first_line <- paste('track name=\"',track_name,'\" colorByStrand=\'255,0,0 0,0,0\' visibility=3',sep='')
	shift.size <- (end(enh.ep)-start(enh.ep))/2
	enh.ep <- GenomicRanges::shift(enh.ep,shift=shift.size)
	enh.ep <- promoters(enh.ep,upstream=enh.win.size,downstream=enh.win.size)
	
	if(gene.start!="mid")
		rg.cons.ep <- promoters(rg.cons.ep,upstream=prom.win.size,downstream=prom.win.size)
	else{
		shift.size <- (end(rg.cons.ep)-start(rg.cons.ep))/2
		rg.cons.ep <- GenomicRanges::shift(rg.cons.ep,shift=shift.size)
		rg.cons.ep <- promoters(rg.cons.ep,upstream=prom.win.size,downstream=prom.win.size)
	}

	strand(rg.cons.ep) <- '*'
	strand(enh.ep) <- '*'

	gint <- GenomicInteractions(enh.ep, rg.cons.ep, counts=as.integer(seq_len(length(enh.ep))))

	if(!is.na(names))
		gint$name <- names
	bed = asBED(gint, score="counts")
	if(!is.na(supp_ep)){
		id.match = match(supp_ep,bed$score)
		strand(bed[id.match])<-'+'
	}
	bed$score=0
	if(to.export){
#		export.bed12(gint, fn = filename, score = "counts")
		export(bed, filename, format="bed")
		con <- file(filename, "r")
		lines <- readLines(con)
		close(con)
		lines <- c(first_line,lines)
		con <- file(filename, "w")
		writeLines(lines , con = con)
		close(con)
	}
	return(gint)

}

asBED <- function(x, keep.mcols=FALSE, score="score") {
        if (!is.null(names(x)) || !is.null(x$name))
            warning("Names will be dropped during BED12 export")

        x = swapAnchors(x, mode="order")

        is_trans = as.vector(seqnames(regions(x))[x@anchor1] != seqnames(regions(x))[x@anchor2])
        a1_cis = x@anchor1[!is_trans]
        a2_cis = x@anchor2[!is_trans]
        a1_trans = x@anchor1[is_trans]
        a2_trans = x@anchor2[is_trans]

        scores = x$counts

        if (is.null(x$color)) {
            x$color = "#000000"
        }

        names = x$name

        cis_blocks = relist(IRanges(
                start=c(rbind(rep(1L, length(a1_cis)), start(x@regions)[a2_cis] - start(x@regions)[a1_cis] + 1L)),
                width=c(rbind(width(x@regions)[a1_cis], width(x@regions)[a2_cis]))),
            PartitioningByWidth(rep(2, length(a1_cis))))

        output_cis = GRanges(
            seqnames=as.character(seqnames(x@regions)[a1_cis]),
            IRanges(start=start(x@regions)[a1_cis],
                    end=end(x@regions)[a2_cis]),
            name=names[!is_trans],
            score=scores[!is_trans],
            strand=ifelse(
                strand(x@regions)[a1_cis] == strand(x@regions)[a2_cis] &
                       as.vector(strand(x@regions)[a1_cis]) %in% c("+", "-"),
                as.vector(strand(x@regions)[a1_cis]),
                "*"),
            thickStart=start(x@regions)[a1_cis],
            thickEnd=end(x@regions)[a2_cis],
            itemRgb=x$color[!is_trans],
            blocks=cis_blocks
        )

        trans_blocks = relist(IRanges(
                start=rep(1, 2*length(a1_trans)),
                width=c(width(x@regions)[a1_trans], width(x@regions)[a2_trans])),
            PartitioningByWidth(rep(1, 2*length(a1_trans))))

        output_trans = GRanges(
            seqnames=c(as.character(seqnames(x@regions)[a1_trans]),
                    as.character(seqnames(x@regions)[a2_trans])),
            IRanges(start=c(start(x@regions)[a1_trans],
                            start(x@regions)[a2_trans]),
                    end=c(end(x@regions)[a1_trans],
                            end(x@regions)[a2_trans])),
            name=rep(names[is_trans], 2),
            score=rep(scores[is_trans], 2),
            strand=c(as.character(strand(x@regions)[a1_trans]),
                        as.character(strand(x@regions)[a2_trans])),
            thickStart=c(start(x@regions)[a1_trans],
                         start(x@regions)[a2_trans]),
            thickEnd=c(end(x@regions)[a1_trans],
                       end(x@regions)[a2_trans]),
            itemRgb=rep(x$color[is_trans], 2),
            blocks=trans_blocks
        )

        extra_cols = setdiff(colnames(mcols(x)), c("score", "name"))

        if(length(extra_cols) && keep.mcols==TRUE) {
            mcols(output_cis) = cbind(mcols(output_cis), mcols(x)[!is_trans, extra_cols,drop=FALSE])
            mcols(output_trans) =
                cbind(mcols(output_trans),
                      rep(mcols(x)[is_trans, extra_cols,drop=FALSE], 2))
        }

        return(sort(c(output_cis, output_trans)))
}

exportBedChIAInteract <- function(chia.1,chia.2,names=NA,track_name="Genomic interactions",filename,to.export=T,...){
	
	first_line <- paste('track name=\"',track_name,'\" colorByStrand=\'255,0,0 0,0,255\' visibility=3',sep='')

	shift.size <- (end(rg.cons.ep)-start(rg.cons.ep))/2
	rg.cons.ep <- GenomicRanges::shift(rg.cons.ep,shift=shift.size)
	rg.cons.ep <- promoters(rg.cons.ep,upstream=win.size,downstream=win.size)
	
	gint <- GenomicInteractions(chia.1,chia.2, counts=as.integer(rep(1,times=length(chia.1))))
	if(!is.na(names))
		gint$name <- names
	if(to.export){
		export.bed12(gint, fn = filename)
		con <- file(filename, "r")
		lines <- readLines(con)
		close(con)
		lines <- c(first_line,lines)
		con <- file(filename, "w")
		writeLines(lines , con = con)
		close(con)
	}
	return(gint)

} 
