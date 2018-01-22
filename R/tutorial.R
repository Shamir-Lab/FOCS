##################################################################################################################################
########################################## Tutorial for running FOCS analyses ####################################################
##################################################################################################################################

########################################## Note: Only works via R unix platform ##################################################

setwd('../focs/')

# load libraries
cran_libs = c('pscl','MASS','parallel','AUC','glmnet','Matching','RColorBrewer')
bioconductor_libs = c('GenomicRanges')
for (p in (c(bioconductor_libs,cran_libs))){library(p,character.only=T)}


# Variables
k <- 10 #k closest enhancers to each gene
win.size <- 5*10**5 #maximum window size upstream/downstream to select genes with at least 10 enhancers
# Number of parallel cores to run on - use detectCores(all.tests = FALSE, logical = TRUE) to identify the number of cores
mc.cores <- 40
# FDR correction method
method='BY'
# FDR threshold
fdr_th = 0.1
# Regression type (1:ZINB, 2:OLS, 3:GLM.NB)
reg_type <- 2
#External validation variables
ext.win.size <- 500 #maximum window size upstream/downstream to intersect TSSs/enhancers with eQTLs or ChIA-PET anchors 
nperms <- 100 #number of permutation tests to validate the true number of supported E-P links by eQTLs/ChIA-PET versus random E-P links

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
# Enhancer count matrix:
Me = as.matrix(get(load(paste(data_directory,data_type,'.enh.count.RData',sep=''))))
#rpkm values of Me
Me_rpkm <- as.matrix(get(load(paste(data_directory,data_type,'.enh.rpkm.RData',sep=''))))
# Promoter count matrix:
Mg = as.matrix(get(load(paste(data_directory,data_type,'.prom.count.RData',sep=''))))
#rpkm values of Mg
Mg_rpkm <- as.matrix(get(load(paste(data_directory,data_type,'.prom.rpkm.RData',sep=''))))
# Sample annotations
sample.annot <- get(load(paste(data_directory,data_type,'.sample.annot.RData',sep='')))
lib_sizes = sample.annot$lib_size
col2type = as.character(sample.annot$cell_type)
# ln of the library size for zinb and GLM-NB2 regression methods
offs=log(lib_sizes)
# Promoter to closest enhnacers in fixed window:
g_to_e = get(load(paste(data_directory,data_type,'.cand.enh.RData',sep='')))
# Enhancer-Promoter links
ep_links <- get(load(paste(data_directory,data_type,'.E_P.links.RData',sep='')))

#GTEx eQTLs data
sig_pairs.gr <- get(load(file=paste(data_directory,'sig_pairs.gr.RData',sep='')))

#ChIA-PET data
chia.1.gr <- get(load(file=paste(data_directory,'chia.1.gr.RData',sep='')))
chia.2.gr <- get(load(file=paste(data_directory,'chia.2.gr.RData',sep='')))

#load code
source(paste(script_directory,'reg_functions.R',sep=''))
source(paste(script_directory,'external_val_functions.R',sep=''))

#regression functions
reg_func <- c(run_zeroinf_test, run_lm_test, run_glmnb_test)

##################################################################################################################################
########################################## Part 1: Regression analysis ###########################################################
##################################################################################################################################

# count per million mapped reads (CPM) matrices
Mg_normalized = t(t(Mg)/lib_sizes)
Me_normalized = t(t(Me)/lib_sizes)

# Run all tests - Works only on UNIX platform 
# list of arguments
arg_list <- list(Mg,Mg_normalized,Me,Me_normalized,g_to_e,offs,Mg_rpkm,col2type)

# Run leave cell-type out cross validation for each regression type (OLS,ZINB,GLM.NB2)
res <- mclapply(seq_len(nrow(Mg)),function(j) apply_cv_test_by_type(j,arg_list,k),mc.cores = mc.cores)
auroc_m <- do.call(rbind,sapply(res,function(x) x[1]))
corr_m <- do.call(rbind,sapply(res,function(x) x[2]))
corr_pvals_m <- do.call(rbind,sapply(res,function(x) x[3]))
preds_m[[1]] <- do.call(rbind,sapply(res,function(x) x[[4]][1]))
preds_m[[2]] <- do.call(rbind,sapply(res,function(x) x[[4]][2]))
preds_m[[3]] <- do.call(rbind,sapply(res,function(x) x[[4]][3]))

# Add the rownames to the matrices
rownames(auroc_m) = rownames(Mg)
rownames(corr_m) = rownames(Mg)
rownames(corr_pvals_m) = rownames(Mg)
for(j in 1:length(preds_m)){
	rownames(preds_m[[j]]) = rownames(Mg)
}

save(auroc_m,file=paste(tmp_directory,data_type,'.auroc_m.RData',sep=''))
save(corr_m,file=paste(tmp_directory,data_type,'.corr_m.RData',sep=''))
save(corr_pvals_m,file=paste(tmp_directory,data_type,'.corr_pvals_m.RData',sep=''))
save(preds_m,file=paste(tmp_directory,data_type,'.preds_m.RData',sep=''))

# Get AUROC scores p-values
auroc_pvals_m = c()
for(i in 1:nrow(Mg)){
	if(i%%100==0){print(paste("$$$$$$$$$$$$$$$$$$$$$$",i))}
	y_rpkm = Mg_fpkm[i,]
	# check if we have at least k enhancers and at least 3 expressed samples > 1 RPKM
	if(length(g_to_e[[i]]) < k || sum(y_rpkm > 1) < 3){
		res = c(NA,NA,NA)
		auroc_pvals_m = rbind(auroc_pvals_m,res)
		next
	}
	if(sum(y_rpkm <= 1) ==0){
		res = c(NA,NA,NA)
		auroc_pvals_m = rbind(auroc_pvals_m,res)
		next
	}
	pos_inds = y_rpkm > 1
	neg_inds = !pos_inds
	curr_pvals = c()
	for (j in 1:length(preds_m)){
		curr_preds = preds_m[[j]][i,]
		x1 = curr_preds[pos_inds]
		x2 = curr_preds[neg_inds]
		curr_pvals[j] = wilcox.test(x1,x2)$p.value
	}
	auroc_pvals_m = rbind(auroc_pvals_m,curr_pvals)			
}
rownames(auroc_pvals_m) = rownames(Mg)
save(auroc_pvals_m,file=paste(tmp_directory,data_type,'.auroc_pvals_m.RData',sep=''))

# get the OLS R.squared values based on all samples without cross-validation
r_square_lm <- c()
for(j in 1:nrow(Mg_normalized)){
	if(j%%100==0){print(paste("$$$$$$$$$$$$$$$$$$$$$$",j))}
	if(length(g_to_e[[j]]) < k){
		r_square_lm <- c(r_square_lm,NA)
		next
	}
	regr_data_normalized = fetch_regr_data(j,Mg_normalized,Me_normalized,g_to_e)
	x_n = regr_data_normalized$x;y_n=regr_data_normalized$y
	res <- run_lm_test(x_n,y_n)
	r_square_lm <- c(r_square_lm,res$r.squared)

}

save(r_square_lm,file=paste(tmp_directory,data_type,'.r_square_lm.RData',sep=''))

##################################################################################################################################
########################################## Part 2: Analyze regression results ####################################################
##################################################################################################################################

# load previous part results
auroc_m <- get(load(file=paste(tmp_directory,data_type,'.auroc_m.RData',sep='')))
auroc_pvals_m <- get(load(file=paste(tmp_directory,data_type,'.auroc_pvals_m.RData',sep='')))
corr_m <- get(load(file=paste(tmp_directory,data_type,'.corr_m.RData',sep='')))
corr_pvals_m <- get(load(file=paste(tmp_directory,data_type,'.corr_pvals_m.RData',sep='')))
preds_m <- get(load(file=paste(tmp_directory,data_type,'.preds_m.RData',sep='')))
r_square_lm <- get(load(file=paste(tmp_directory,data_type,'.r_square_lm.RData',sep='')))

gene2num_positive_samples = apply(Mg_rpkm>1,1,sum,na.rm=T)

# Look at the OLS results: obtain estimation for the FDR and power of using R.squared values
binary_validation_qvals = p.adjust(auroc_pvals_m[,reg_type],method=method)
binary_validation_qvals[is.na(binary_validation_qvals)] = 1
binary_validation_qvals = set_min_greater_than_zero(binary_validation_qvals)
level_validation_qvals = p.adjust(corr_pvals_m[,reg_type],method=method)
level_validation_qvals[is.na(level_validation_qvals)] = 1
level_validation_qvals = set_min_greater_than_zero(level_validation_qvals)

# Divide promoter models into groups based on the FDR threshold
gene_val_groups  = list()
# Promoter models that passed both 'binary' and 'expression level' validations
gene_val_groups[["Both"]] = names(which(binary_validation_qvals<=fdr_th & level_validation_qvals<=fdr_th))
# Promoter models that passed only the binary validation
gene_val_groups[["Binary only"]] = names(which(binary_validation_qvals<=fdr_th & !level_validation_qvals<=fdr_th))
# Promoter models that passed only the expression level validation
gene_val_groups[["Level only"]] = names(which(!binary_validation_qvals<=fdr_th & level_validation_qvals<=fdr_th))
# Promoter models that did not pass any validation
gene_val_groups[["None"]] = names(which(!binary_validation_qvals<=fdr_th & !level_validation_qvals<=fdr_th))
save(gene_val_groups,file=paste(tmp_directory,data_type,'.gene_val_groups.RData',sep=''))


################################################ Performance plots ##########################################################

# Web figure - a/b plots
# Plot percentiles of gene/promoter models vs. q-values in each regression method: ZINB, OLS, GLM.NB
# You may want to change some parameters, e.g. ylim/width/height, to adjust the plots

png(paste(tmp_directory,data_type,'.reg.png',sep=''),width=1400,height=750)
par(mfrow=c(1,2),mex=1.15,cex=1.8)
na_rows = apply(is.na(auroc_pvals_m),1,all)
roc_scores_p = auroc_pvals_m[!na_rows,]
roc_scores_q = apply(roc_scores_p,2,p.adjust,method=method)
cor(roc_scores_p)

probs <- seq(0.2,0.9,0.05)
# Binarized expression validation plot
dat <- apply(-log10(roc_scores_q),2,quantile,probs=probs )
plot(probs ,dat[,1],type="l",main="Binarized expression validation",xlab="Percentile", ylab=expression('-log'[10]*'(q-value)'),cex.lab=1.6,cex.axis=1.6,col="green",lwd = 2,lty=2,ylim=c(0,11))
lines(probs ,dat[,2], col="red",lwd = 2,lty=1)
lines(probs ,dat[,3], col ="blue",lwd = 2,lty=4)
legend("topleft",legend=c("ZINB","OLS","GLM.NB"),col=c("green","red","blue"),lwd = 2,bty="n",lty=c(2,1,4))

na_rows = apply(is.na(corr_pvals_m),1,any)
spearman_pvals = corr_pvals_m[!na_rows,]
spearman_qvals = apply(spearman_pvals,2,p.adjust,method=method)

# Expression level validation plot
dat <- apply(-log10(spearman_qvals),2,quantile,probs=probs )
plot(probs ,dat[,1],type="l",main="Expression level validation",xlab="Percentile", ylab=expression('-log'[10]*'(q-value)'),cex.lab=1.6,cex.axis=1.6,col="green",lwd = 2,lty=2,ylim=c(0,4))
lines(probs ,dat[,2], col="red",lwd = 2,lty=1)
lines(probs ,dat[,3], col ="blue",lwd = 2,lty=4)
legend("topleft",legend=c("ZINB","OLS","GLM.NB"),col=c("green","red","blue"),lwd = 2,bty="n",lty=c(2,1,4))

dev.off()

# Web figure - plot c
# Pie chart of the number of models passed each group (Both, Binary only, Level only)
# Boxplot of the number of positive samples each model contains across the groups
# You may want to change some parameters, e.g. ylim/width/height, to adjust the plots

png(paste(tmp_directory,data_type,'.pie_box.png',sep=''),width=1000,height=1200)
par(mfrow=c(2,1),mar=c(6,6,1,1),cex=2.0,mex=1.35)
marker = list(color = brewer.pal(4, "Pastel1"))
cols = marker$color
group_sizes <- sapply(gene_val_groups,length)
names(cols) = names(group_sizes)
labls = paste("\n\n",names(group_sizes),group_sizes,"\n",sep="\n")
labls[1] = paste(names(group_sizes)[1],group_sizes[1],"\n",sep="\n")
labls[2] = paste(names(group_sizes)[2],group_sizes[2],"\n",sep="\n")
labls[3] = paste(names(group_sizes)[3],group_sizes[3],"\n",sep="\n")
pie(group_sizes,labels=labls,col=cols,lwd=2,cex=1,radius=0.9)
expr_by_group = lapply(gene_val_groups,function(x,y)y[x],y=gene2num_positive_samples)
boxplot(expr_by_group,ylab = paste0("Number of positive samples","\n","(RPKM>1)"),las=2,col=cols,cex.lab=1.2,cex.names=1.2)
dev.off()


# Web figure - plot d
# R.squared without CV vs. R.squared with CV
# You may want to change some parameters, e.g. ylim/width/height, to adjust the plo

r_cv <- sapply(seq_len(dim(Mg)[1]),function(i) getRsquared(i,Mg_normalized[i,],preds_m[[reg_type]][i,]))
png(paste(tmp_directory,data_type,'.r_cv_vs_r2.png',sep=''),width=1000,height=1000)
par(mar=c(6,6,1,1),cex=1.8)
cols = rep("black",length(r2))
cols[r_square_lm>=0.5] = "red"
cols[r_square_lm>=0.5 & r_cv >= 0.25] = "blue"
plot(r_square_lm,r_cv,xlab = expression("R"^2), ylab = expression("R"[CV]^2), cex.lab=1.6,cex.axis=1.6,col=cols,mex=1.15)
table(cols)
sum(cols=="red")/sum(cols!="black")
dev.off()


# Plot the spearman coefficient vs. R.squared of the OLS regression
# You may want to change some parameters, e.g. ylim/width/height, to adjust the plo

png(paste(tmp_directory,data_type,'.sp_vs_r2.png',sep=''),width=1000,height=1000)
par(mar=c(6,6,1,1),cex=1.8)
spearman_rs = corr_m[,reg_type]
cols = rep("black",length(r_square_lm))
cols[r_square_lm>=0.5] = "red"
cols[r_square_lm>=0.5 & spearman_rs >= 0.5] = "blue"
plot(r_square_lm,spearman_rs,xlab = expression("R"^2), ylab = expression("Spearman "*symbol("\162")), cex.lab=1.6,cex.axis=1.6,col=cols,mex=1.15)
# Percentage of models with R.squared>0.5 and Spearman coefficient <0.5 - This may suggest problematic models (e.g., models with high percentages of zero samples)
sum(cols=="red")/sum(cols!="black")
dev.off()



##################################################################################################################################
########################################## Part 3: Promoter model shrinkage using elastic-net ####################################
##################################################################################################################################

# Note: this part only works with OLS regression type

# load previous part results
gene_val_groups <- get(load(file=paste(tmp_directory,data_type,'.gene_val_groups.RData',sep='')))

# Which promoter/gene group models to perform elastic-net - in GRO-seq we selected "Both", "Binary only" and "Level only"
gene_list <- union(union(gene_val_groups[["Both"]],gene_val_groups[["Binary only"]]),gene_val_groups[["Level only"]])
# In ENCODE/Roadmap/FANTOM5 we selected "Both" and "Level only" validation groups
#gene_list <- union(gene_val_groups[["Both"]],gene_val_groups[["Level only"]])
inds <- match(gene_list,rownames(Mg));

g_to_e_lmfit <- matrix(0,length(inds),k)
func <- reg_func[reg_type] #default OLS
rownames(g_to_e_lmfit) <- gene_list

for(i in seq_len(length(inds))){
	regr_data_normalized = fetch_regr_data(inds[i],Mg_normalized,Me_normalized,g_to_e)
	x_n = regr_data_normalized$x;y_n=regr_data_normalized$y
	m <- func(x_n,y_n)
	coef_pv <- m$coeffs[2:11,4]
	g_to_e_lmfit[i,] <- coef_pv
	if(i%%100==0){print(paste("$$$$$$$$$$$$$$$$$$$$$$", i))}
}

# Step 1 - eBY filtering
# In this step we select enhancers that should appear in the shrunken model 
# Collect all enhancer p.values from all models to a single vector 
# Correct enhancer p.values by correction method - default "BY" and FDR 0.01
g_to_e_lmfit_correction_obj = correct_p_matrix(g_to_e_lmfit,0.01,"BY")
# maximum enhancer p.value that passed the FDR threshold
g_to_e_lmfit_correction_obj[[1]]

# Step 2 - elastic-net
# Perform elastic-net with alpha mixing parameter 0.5
g_to_e_cvfit <- matrix(0,length(inds),k)
rownames(g_to_e_cvfit) <- gene_list

for(i in seq_len(length(inds))){
	regr_data_normalized = fetch_regr_data(inds[i],Mg_normalized,Me_normalized,g_to_e)
	x_n = regr_data_normalized$x;y_n=regr_data_normalized$y
	gfit <- glmnet(x_n,y_n,alpha = 0.5)
	# Take the current model survived enhancers and make sure that they appear in the shrunken model
	coef_lm <- g_to_e_lmfit_correction_obj[[2]][i,]
	pos.id <- which(coef_lm==T)
	for(j in 1:ncol(gfit$beta)){
		coef_gfit <- gfit$beta[,j]
		# if no enhancer survived step 1 then find the maximum lambda with at least one non-zero enhancer coefficient
		if(all(coef_lm==F)){
			if(any(coef_gfit!=0))
				break
			else
				next
		}
		else{
			# Check if the enhancers that survived step 1 have non-zero coefficients in the shrunken model
			if(all(coef_gfit[pos.id]!=0))
				break
		}

	}
	g_to_e_cvfit[i,] <- coef_gfit!=0
	if(i%%100==0){print(paste("$$$$$$$$$$$$$$$$$$$$$$", i))}
}

#which enhancer ids were selected per gene model after glmnet
g_to_e_indx <- apply(g_to_e_cvfit,1,function(row) which(row!=0))
#how mnay enhancers were mapped to each gene - stats
summary(sapply(g_to_e_indx,function(x) length(x)))
table(sapply(g_to_e_indx,function(x) length(x)))
#how many gene models had at least 1 enhancer selected after glmnet - should be the same as length(inds)
length(which(sapply(g_to_e_indx,function(x) length(x))!=0))

# Plot the enhancer rank by distance to their closest gene's TSS - web figure: plot a
png(paste(tmp_directory,data_type,'.ranks.png',sep=''),width=1100,height=1000)
par(mar=c(6,6,2,1),cex=2.3)
barplot(table(unlist(g_to_e_indx))/length(inds),names.arg = 1:10,xlab="Enhancer (ranked by proximity)",ylim=c(0,0.6), ylab="Model inclusion frequency", cex.lab=1.6,cex.axis=1.6,cex.names=1.6,mex=1.15)
dev.off()

# Plot the number of optimal models with 1-10 enhancers - web figure: plot b
png(paste(tmp_directory,data_type,'.totEnhPerModel.png',sep=''),width=1100,height=1000)
par(mar=c(6,6,4,2)+0.1,cex=2.0)
barplot(table(sapply(g_to_e_indx,function(x) length(x))),names.arg = 1:10, ylab="Number of models",ylim=c(0,1400), cex.lab=1.6,cex.axis=1.6,cex.names=1.4)
mtext(paste0("Number of enhancers","\n","in optimal model"),side=1,line=4,cex=2.0*1.6)
dev.off()

# create the E-P network as list
ep_links <- sapply(seq_len(length(inds)),function(i) g_to_e[[inds[i]]][g_to_e_indx[[i]]])
#how many E-P links we have?
length(unlist(ep_links))
#how many unique enhancers were mapped?
length(unique(unlist(ep_links)))
names(ep_links) <- gene_list
save(ep_links,file=paste(tmp_directory,data_type,'.E_P.links.RData',sep=''))

##################################################################################################################################
################################### Part 4: External validation with ChIA-PET and GTEx eQTLs #####################################
##################################################################################################################################

################################################## GRO-seq data type only ########################################################
# get E-P links
ep_links <- get(load(file=paste(tmp_directory,data_type,'.E_P.links.RData',sep='')))

# get the reduced promoter genomic positions according to ep_links
prom.bs.ep <- prom.bs[match(names(ep_links),names(prom.bs))]

# The matrices mark for each promoter/gene the number of E-P links supported by at least one SNPs/ChIA-PET interaction
count_overlap_mat_chia <- matrix(0,length(ep_links),1+nperms)
count_overlap_mat_snp <- matrix(0,length(ep_links),1+nperms)

# check overlap with ChIA-PET anchors
res <- mclapply(seq_len(length(ep_links)),function(i) count_overlap_with_chia(i,enh.bs[ep_links[[i]]],prom.bs.ep[i],chia.1.gr,chia.2.gr,ext.win.size),mc.cores=mc.cores)
sum(unlist(res))/length(unlist(ep_links))

count_overlap_mat_chia[,1] <- unlist(res)

# check overlap with eQTLs - gene's associated SNPs that fall within linked enhancers
res <- mclapply(seq_len(length(ep_links)), function(i) find_overlap_with_eQTL_snp(i,enh.bs[ep_links[[i]]],
									unique(sig_pairs.gr[which(sig_pairs.gr$ent_id==names(ep_links[i]))]),ext.win.size),mc.cores=mc.cores)
edge_count <- unlist(sapply(res,function(x) x$edge_hit),use.names=F)
sum(edge_count)/length(unlist(ep_links))

count_overlap_mat_snp[,1] <- edge_count

# perfrom the random permutation test
# In each permutation we preserve the promoter degree (#enhancers) and select different enhancers with
# similar enhancer distance from TSS distribution as the original enhancers using Matching R package
for( perm in seq_len(nperms)){
	perm_E_P <- get_perm_EP(ep_links,prom.bs.ep,enh.bs,win.size,mc.cores)
	res <- mclapply(seq_len(length(perm_E_P)),function(i) count_overlap_with_chia(i,enh.bs[perm_E_P[[i]]],prom.bs.ep[i],chia.1.gr,chia.2.gr,ext.win.size),mc.cores=mc.cores)
	count_overlap_mat_chia[,1+perm] <- unlist(res)
	print(sum(unlist(res))/length(unlist(ep_links)))
	res <- mclapply(seq_len(length(perm_E_P)), function(i) find_overlap_with_eQTL_snp(i,enh.bs[perm_E_P[[i]]],
										unique(sig_pairs.gr[which(sig_pairs.gr$ent_id==names(perm_E_P[i]))]),ext.win.size),mc.cores=mc.cores)
	edge_count <- unlist(sapply(res,function(x) x$edge_hit),use.names=F)
	print(sum(edge_count)/length(unlist(ep_links)))
	count_overlap_mat_snp[,1+perm] <- edge_count
	print(paste("Finished permutation number:   ",perm))
}

#Get empirical P-value from random permutation test
chia_emp_pval <- getExtValEmpPval(count_overlap_mat_chia,ep_links)
chia_emp_pval
eqtls_emp_pval <- getExtValEmpPval(count_overlap_mat_snp,ep_links)
eqtls_emp_pval

################################################## FANTOM5 data type only ########################################################

# Get a map from fantom promoter-TSSs to associated SNPs
snp.inds <- get(load(file=paste(data_directory,data_type,'.snp.inds.RData',sep='')))

# get E-P links
ep_links <- get(load(file=paste(tmp_directory,data_type,'.E_P.links.RData',sep='')))

# get the reduced promoter genomic positions according to ep_links
prom.bs.ep <- prom.bs[match(names(ep_links),prom.bs$name)]


# The matrices mark for each promoter/gene the number of E-P links supported by at least one SNPs/ChIA-PET interaction
count_overlap_mat_chia <- matrix(0,length(ep_links),1+nperms)
count_overlap_mat_snp <- matrix(0,length(ep_links),1+nperms)

# check overlap with ChIA-PET anchors
res <- mclapply(seq_len(length(ep_links)),function(i) count_overlap_with_chia(i,enh.bs[ep_links[[i]]],prom.bs.ep[i],chia.1.gr,chia.2.gr,ext.win.size),mc.cores=mc.cores)
sum(unlist(res))/length(unlist(ep_links))

count_overlap_mat_chia[,1] <- unlist(res)

# check overlap with eQTLs - gene's associated SNPs that fall within linked enhancers
gene_ep_list <- prom.bs.ep$gene_id

res <- mclapply(seq_len(length(ep_links)), function(i) find_overlap_with_eQTL_snp(i,enh.bs[ep_links[[i]]],
									unique(sig_pairs.gr[unlist(snp.inds[gene_ep_list[[i]]])]),ext.win.size),mc.cores=mc.cores)

#find only promoter-TSSs annotated with known gene
id <- which(sapply(prom.bs.ep$gene_id,function(x) any(is.na(x)==T))==F)
edge_count <- unlist(sapply(res,function(x) x$edge_hit),use.names=F)
sum(edge_count)/length(unlist(ep_links[prom.bs.ep$name[id]]))
count_overlap_mat_snp[,1] <- edge_count

# perfrom the random permutation test
# In each permutation we preserve the promoter degree (#enhancers) and select different enhancers with
# similar enhancer distance from TSS distribution as the original enhancers using Matching R package
for( perm in seq_len(nperms)){
	perm_E_P <- get_perm_EP(ep_links,prom.bs.ep,enh.bs,win.size,mc.cores)
	res <- mclapply(seq_len(length(perm_E_P)),function(i) count_overlap_with_chia(i,enh.bs[perm_E_P[[i]]],prom.bs.ep[i],chia.1.gr,chia.2.gr,ext.win.size),mc.cores=mc.cores)
	count_overlap_mat_chia[,1+perm] <- unlist(res)
	print(sum(unlist(res))/length(unlist(ep_links)))
	res <- mclapply(seq_len(length(perm_E_P)), function(i) find_overlap_with_eQTL_snp(i,enh.bs[perm_E_P[[i]]],
										unique(sig_pairs.gr[unlist(snp.inds[gene_ep_list[[i]]])]),ext.win.size),mc.cores=mc.cores)
	edge_count <- unlist(sapply(res,function(x) x$edge_hit),use.names=F)
	count_overlap_mat_snp[,1+perm] <- edge_count
	print(paste("Finished permutation number:   ",perm))
}

#Get empirical P-value from random permutation test
chia_emp_pval <- getExtValEmpPval(count_overlap_mat_chia,ep_links[prom.bs.ep$name[id]])
chia_emp_pval
eqtls_emp_pval <- getExtValEmpPval(count_overlap_mat_snp,ep_links[prom.bs.ep$name[id]])
eqtls_emp_pval

################################################## ENCODE/Roadmap data type only ########################################################

# Get gencodeV10 TSS positions
tss.bs <- get(load(file=paste(data_directory,'tss.bs.genecode10.RData',sep='')))

# Get a map from gencodeV10 TSSs to associated SNPs
snp.inds <- get(load(file=paste(data_directory,'genecodev10.snp.inds.RData',sep='')))

# get E-P links
ep_links <- get(load(file=paste(tmp_directory,data_type,'.E_P.links.RData',sep='')))

# get the reduced promoter/enhancer genomic positions according to ep_links
prom.bs.ep <- prom.bs[match(names(ep_links),prom.bs$name)]


# The matrices mark for each promoter/gene the number of E-P links supported by at least one SNPs/ChIA-PET interaction
count_overlap_mat_chia <- matrix(0,length(ep_links),1+nperms)
count_overlap_mat_snp <- matrix(0,length(ep_links),1+nperms)

# check overlap with ChIA-PET anchors
res <- mclapply(seq_len(length(ep_links)),function(i) count_overlap_with_chia(i,enh.bs[ep_links[[i]]],prom.bs.ep[i],chia.1.gr,chia.2.gr,ext.win.size),mc.cores=mc.cores)
sum(unlist(res))/length(unlist(ep_links))

count_overlap_mat_chia[,1] <- unlist(res)

gene_ep_list <- mclapply(seq_len(length(prom.bs.ep)),function(i) get_gene_tss(i,prom.bs.ep[i],tss.bs),mc.cores=mc.cores)

res <- mclapply(seq_len(length(ep_links)), function(i) find_overlap_with_eQTL_snp(i,enh.bs[ep_links[[i]]],unique(sig_pairs.gr[unlist(snp.inds[gene_ep_list[[i]]$ens_id])]),ext.win.size),mc.cores=mc.cores)

#find only promoter-DHSs annotated with known gene
id <- which(sapply(prom.bs.ep$overlap_TSS_1kb_up,function(x) any(is.na(x)==T))==F)
edge_count <- unlist(sapply(res,function(x) x$edge_hit),use.names=F)
sum(edge_count)/length(unlist(ep_links[prom.bs.ep$name[id]]))
count_overlap_mat_snp[,1] <- edge_count

# perfrom the random permutation test
# In each permutation we preserve the promoter degree (#enhancers) and select different enhancers with
# similar enhancer distance from TSS distribution as the original enhancers using Matching R package
for( perm in seq_len(nperms)){
	perm_E_P <- get_perm_EP(ep_links,prom.bs.ep,enh.bs,win.size,mc.cores)
	res <- mclapply(seq_len(length(perm_E_P)),function(i) count_overlap_with_chia(i,enh.bs[perm_E_P[[i]]],prom.bs.ep[i],chia.1.gr,chia.2.gr,ext.win.size),mc.cores=mc.cores)
	count_overlap_mat_chia[,1+perm] <- unlist(res)
	print(sum(unlist(res))/length(unlist(ep_links)))
	res <- mclapply(seq_len(length(perm_E_P)), function(i) find_overlap_with_eQTL_snp(i,enh.bs[perm_E_P[[i]]],
										unique(sig_pairs.gr[unlist(snp.inds[gene_ep_list[[i]]$ens_id])]),ext.win.size),mc.cores=mc.cores)
	edge_count <- unlist(sapply(res,function(x) x$edge_hit),use.names=F)
	count_overlap_mat_snp[,1+perm] <- edge_count
	print(paste("Finished permutation number:   ",perm))
}

#Get empirical P-value from random permutation test
chia_emp_pval <- getExtValEmpPval(count_overlap_mat_chia,ep_links[prom.bs.ep$name[id]])
chia_emp_pval
eqtls_emp_pval <- getExtValEmpPval(count_overlap_mat_snp,ep_links[prom.bs.ep$name[id]])
eqtls_emp_pval



##################################################################################################################################
################################### Part 5: Integrative analysis with ChIA-PET and GTEx eQTLs ####################################
##################################################################################################################################

# This part is in 'tutorial_integrative_analysis.R' script



