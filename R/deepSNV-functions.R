# functions of the deepSNV package
# 
# Author: Moritz Gerstung
###############################################################################

#' The actual function for the deepSNV test
#' @noRd
.deepSNV <- function(test, control, nucleotides, dirichlet.prior, alternative, model, over.dispersion, combine.method, pseudo.count) {
	if(length(nucleotides)==10)
		CV <- consensusSequence(control[,1:5]+control[,6:10], vector=TRUE)
	else
		CV <- consensusSequence(control[,1:5], vector=TRUE)	
	
	res0 <- .deepSNVsingle(test = test[,1:5], control = control[,1:5], dirichlet.prior = dirichlet.prior, alternative = alternative, CV = CV, model=model, alpha=over.dispersion, pseudo.count=pseudo.count)
	log.lik = res0$log.lik 
	
	#Combine P-values from both strands.
	if(length(nucleotides) == 10){
		res1 <- .deepSNVsingle(test = test[,6:10], control = control[,6:10], dirichlet.prior = dirichlet.prior[,COMPLEMENT[colnames(dirichlet.prior)]], alternative = alternative, CV = COMPLEMENT[CV], strand=1, model=model, alpha=over.dispersion, pseudo.count = pseudo.count)
		p.val <- p.combine(res0$p.val, res1$p.val, method=combine.method)
		log.lik = log.lik + res1$log.lik
	}
	else{
		p.val <- res0$p.val
	}
	for (i in 1:nrow(p.val)) p.val[i, CV[i]] <- NA
	return(list(p.val = p.val, log.lik = log.lik))
}

#' The actual workhorse for the deepSNV test
#' @noRd
.deepSNVsingle <- function(test, control, dirichlet.prior, alternative, CV, strand=0, model, alpha=NULL, pseudo.count) {	
	sum.test <- rowSums(test)
	sum.control <- rowSums(control)
	colnames(test) <- colnames(control) <- c("A","T","C","G","-")
	
	f.test <- RF(test[,1:5]+pseudo.count)
	f.control <- RF(control[,1:5]+dirichlet.prior[CV,]) 
	f.both <- RF(test[,1:5]+control[,1:5]+dirichlet.prior[CV,])
	
	if(model == "bin"){
		## H0: Both the same:
		L0 <- (test+control) * log(f.both) + (sum.test - test + sum.control - control) * log(1-f.both)
		
		## H1: Not the same:
		L1 <- test * log(f.test) + (sum.test - test) * log(1-f.test) + control * log(f.control) + (sum.control - control) * log(1-f.control)
	}
	if(model == "betabin"){
		beta.test <- alpha * (1-f.test)/f.test
		beta.control <- alpha * (1-f.control)/f.control
		beta.both <- alpha * (1-f.both)/f.both
		L0 <- lbeta(test + alpha, sum.test - test + beta.both) + lbeta(control + alpha, sum.control - control + beta.both) - 2 * lbeta(alpha, beta.both)
		L1 <- lbeta(test + alpha, sum.test - test + beta.test) + lbeta(control + alpha, sum.control - control + beta.control) - lbeta(alpha, beta.test) - lbeta(alpha, beta.control)
	}
	
	L1[L1 >=0] <- .Machine$double.xmin
	L0[L0 >=0] <- .Machine$double.xmin
	D <- -2 * (L0 - L1 )
	
	p.val <- pchisq(D, 1, lower.tail=FALSE)
	colnames(p.val) <- colnames(test)
	
	log.choose <- lgamma(sum.test +1) - lgamma(test[,1:5] +1) - lgamma(sum.test-test[,1:5]+1) +lgamma(sum.control +1) - lgamma(control[,1:5] +1) - lgamma(sum.control-control[,1:5]+1) 
	CV.ind <- rep(1:5, each=nrow(control)) == apply(control, 1, which.max)
	log.lik <- sum(L0[!CV.ind] + log.choose[!CV.ind])
	
	## Without prior
	f.test <- RF(test[,1:5])
	f.control <- RF(control[,1:5]) 
	
	## If smaller, take H0
	if (alternative == 'greater'){
		f.se.g <- f.test < f.control
		#f[f.se.g] <- g[f.se.g] <- fg[f.se.g]
		p.val <- p.val / 2
		p.val[f.se.g] <- 1 # Can be set to one, eq. to L0/L1 = 1 ..
		log.lik <- sum(L0[f.se.g & !CV.ind] + log.choose[f.se.g & !CV.ind])
	}  
	
	else if (alternative == 'less'){
		g.se.f <- f.test > f.control
		#f[f.le.g] <- g[f.le.g] <- fg[f.le.g]
		p.val <- p.val / 2
		p.val[g.se.f] <- 1 # Can be set to one, eq. to L0/L1 = 1 ..
		log.lik <- sum(L0[g.se.f & !CV.ind] + log.choose[g.se.f & !CV.ind])
	}
	
	return(list(p.val=p.val, log.lik = log.lik))
}

#' Parameter estimation for the beta-binomial model
#' This function estimates the parameter alpha of the beta-binomial model that is shared across sites. Since the mean is estimated for each locus separately,
#' alpha effectively controls the variance.
#' @param test The counts of the test experiment.
#' @param control The counts of the control experiment.
#' @param dirichlet.prior The prior.
#' @param CV A vector of the consensus.
#' @return alpha The parameter alpha of the beta-binomial model.
#' @author Moritz Gerstung
#' @noRd
.estimateDispersion <- function(test, control, dirichlet.prior = NULL, CV=NULL, alternative = c("two.sided", "greater", "less"), interval=c(0,1000)) {
	if(is.null(CV)) CV = consensusSequence(control[,1:5], vector=T)
	if(is.null(dirichlet.prior)) dirichlet.prior = matrix(1, nrow=5, ncol=5, dimnames = list(colnames(control)[1:5], colnames(control)[1:5]))
	alternative = match.arg(alternative)
	
	if(alternative == "less") alt.ind = RF(test, total=T) > RF(control, total=T) 
	else if(alternative == "greater") alt.ind = RF(test, total=T) < RF(control, total=T) 
	else alt.ind = T
	
	totalLoglik <- function(alpha, test, control, dirichlet.prior, CV, CV.ind){
		
		l.tot <- 0
		for(s in 1:(ncol(test)/5) - 1){
			sum.test <- rowSums(test[,1:5 + 5*s])
			sum.control <- rowSums(control[,1:5+ 5*s])
			
			#f.test <- RF(test[,1:5+ 5*s]+dirichlet.prior[CV,ifelse(s==0, 1:5, COMPLEMENT[colnames(dirichlet.prior)])])
			#f.control <- RF(control[,1:5+ 5*s]+dirichlet.prior[CV,ifelse(s==0, 1:5, COMPLEMENT[colnames(dirichlet.prior)])])
			f.both <- RF(test[,1:5+ 5*s]+control[,1:5+ 5*s]+dirichlet.prior[CV,ifelse(s==0, 1:5, COMPLEMENT[colnames(dirichlet.prior)])])
			
			#beta.test <- alpha * (1-f.test)/f.test
			#beta.control <- alpha * (1-f.control)/f.control
			beta.both <- alpha * (1-f.both)/f.both
			
			#log.choose <- lgamma(sum.test +1) - lgamma(test[,1:5+ 5*s] +1) - lgamma(sum.test-test[,1:5+ 5*s]+1) +lgamma(sum.control +1) - lgamma(control[,1:5+ 5*s] +1) - lgamma(sum.control-control[,1:5+ 5*s]+1) 
			
			l0 <- lbeta(test[,1:5+ 5*s] + alpha, sum.test - test[,1:5+ 5*s] + beta.both) + lbeta(control[,1:5+ 5*s] + alpha, sum.control - control[,1:5+ 5*s] + beta.both) - 2 * lbeta(alpha, beta.both) #+ log.choose
			#l1 <- lbeta(test[,1:5+ 5*s] + alpha, sum.test - test[,1:5+ 5*s] + beta.test) + lbeta(control[,1:5+ 5*s] + alpha, sum.control - control[,1:5+ 5*s] + beta.control) - lbeta(alpha, beta.test) - lbeta(alpha, beta.control) + log.choose
			
			l.tot <- l.tot + sum(l0[!CV.ind & alt.ind])
		}
		return(l.tot)
	}
	
	CV.ind <- sapply(colnames(control[,1:5]), function(nt) CV == nt)
	
	opt <- optimize(f=totalLoglik, interval=interval, test, control, dirichlet.prior, CV, CV.ind, maximum=TRUE)
	alpha <- opt$maximum
	return(opt)
}

#' Relative frequencies.
#' 
#' Convenience function to compute the relative frequencies from a matrix with absolute counts.
#' @param freq A matrix with nucleotide counts.
#' @param total If the nucleotide counts have columns for forward and reverse direction, return each strand sepratatelu (FALSE), or add the two (TRUE). 
#' @return A matrix with the relative frequencies.
#' @author Moritz Gerstung
#' @export
#' @examples data(HIVmix)
#' RF(test(HIVmix))[1:10,]
#' RF(test(HIVmix), total=TRUE)[1:10,]
RF <- function(freq, total = FALSE){
	if(total){
		freq <- matrix(freq[,1:5]+freq[,6:10], ncol=5, dimnames=list(NULL, colnames(freq)[1:5]))
	}
	if(ncol(freq)==10){
		return(cbind(freq[,1:5]/rowSums(freq[,1:5]),freq[,6:10]/rowSums(freq[,6:10])))
	}
	else if(ncol(freq)==5){
		return(freq/rowSums(freq))
	}
}



#' Minor allele sequence
#' @noRd
.MASequence <- function(freq, vector=FALSE){
	tfreq <- t(freq)
	cons <- factor(apply(tfreq,2,function(x) order(x, decreasing=TRUE)[2]), levels=1:5)
	levels(cons) <- rownames(tfreq)
	seq <- DNAString(paste(as.character(cons), collapse=""))
	return(seq)
}

#' Base substitution matrix
#' @noRd
.BSM <- function(freq, CV, complement=F){
	bases <- c("A","T","C","G")
	if(complement) {
		cfreq <- freq[,c(COMPLEMENT[bases], setdiff(colnames(freq), bases))]
		colnames(cfreq) <- colnames(freq)
		freq <- rbind(freq,cfreq)
	}
	t(sapply(bases, function(i) colMeans(RF(freq[CV == i,]))))
}

#' Utility for making unstructured VCF headers (e.g. ##date)
#' @noRd
.makeVCFheader <- function(name, value){
	return(as(splitAsList(DataFrame(Value=value, row.names=name), name), "SimpleDataFrameList"))
}

#' Significant SNVs.
#' @noRd
.significantSNV <- function(deepSNV, sig.level = 0.05, adjust.method = "bonferroni", fold.change=1, value="data.frame"){
	if(is.null(adjust.method)) q = deepSNV@p.val
	else q = p.adjust(deepSNV@p.val, method=adjust.method)
	CV.control <- consensusSequence(deepSNV, vector=TRUE)
	f <- RF(deepSNV@test, total=T)
	g <- RF(deepSNV@control, total=T)
	d <- f - g
	fc <- pmax(f,g)/pmin(f,g)
	var.d <- f/rowSums(deepSNV@test) + g/rowSums(deepSNV@control)
	
	if(is.null(fold.change))
		fccond <- TRUE
	else
		fccond <- fc > fold.change 
		
	cond = q <= sig.level & !is.na(q) & fccond 
	cond.idx <- which(cond, arr.ind=TRUE)
	cond.rows <- cond.idx[,1]
	cond.cols <- cond.idx[,2]
	
#	strand.bias.test <- (RF(deepSNV@test[,1:5]+1)/RF(deepSNV@test[,6:10]+1))[cond]
#	strand.bias.control <- (RF(deepSNV@control[,1:5]+1)/RF(deepSNV@control[,6:10]+1))[cond]
	
	lengths = deepSNV@regions$stop-deepSNV@regions$start +1
	beg = cumsum(c(1,lengths[-length(lengths)]))
	end = cumsum(lengths)
	reg_ix = as.numeric(sapply(cond.rows, function(i) which(end >= i)[1]))
	pos = deepSNV@regions$start[reg_ix] + cond.rows - beg[reg_ix]
	
	table <- data.frame(
					chr = deepSNV@regions$chr[reg_ix],
					pos = pos,
					ref = CV.control[cond.rows],
					var = factor(cond.cols, levels=1:5, labels = colnames(f)),
					p.val = q[cond],
					freq.var = d[cond],
					sigma2.freq.var = var.d[cond],
#					count.test = testCounts(deepSNV)[cond],
#					count.control = controlCounts(deepSNV)[cond],
#					fold.change = fc[cond],
#					strand.bias.test = log2(strand.bias.test),
#					strand.bias.control = log2(strand.bias.control),
					n.tst.fw = deepSNV@test[,1:5][cond],
					cov.tst.fw = rowSums(deepSNV@test[cond.rows,1:5, drop=FALSE]),
					n.tst.bw = deepSNV@test[,6:10][cond],
					cov.tst.bw = rowSums(deepSNV@test[cond.rows,6:10, drop=FALSE]),
					n.ctrl.fw = deepSNV@control[,1:5][cond],
					cov.ctrl.fw = rowSums(deepSNV@control[cond.rows,1:5, drop=FALSE]),
					n.ctrl.bw = deepSNV@control[,6:10][cond],
					cov.ctrl.bw = rowSums(deepSNV@control[cond.rows,6:10, drop=FALSE])
			)
	if(ncol(deepSNV@regions) > 3){
		table <- cbind(table, deepSNV@regions[reg_ix,-c(1,2,3)])
		colnames(table)[15 + 1:(ncol(deepSNV@regions)-3)] = colnames(deepSNV@regions)[-c(1,2,3)]
	}
	if(!is.null(adjust.method)){
		table$raw.p.val <- deepSNV@p.val[cond]
	}
	rownames(table) <- NULL
	
	if(value == "data.frame"){
		o <- order(table$chr, table$pos)
		table <- table[o,]
		return(table)
	}
	else{
		if(nrow(table)==0)
			return(NULL)
		isDel <- table$var == "-"
		isIns <- table$ref == "-"
		if(length(deepSNV@files)==2)
			samples <- sub("/.+/","",deepSNV@files)
		else
			samples <- c("test","control")
		headerTemplate = scanVcfHeader(system.file("extdata", "deepSNV.vcf", package="deepSNV"))
		v = VCF(
				rowRanges=GRanges(table$chr, 
						IRanges(table$pos - (isDel | isIns), width=1 + (isDel | isIns)), 
				),
				fixed = DataFrame(
						REF = DNAStringSet(paste(ifelse(isDel | isIns, as.character(CV.control[cond.rows - 1]),""), sub("-", "", table$ref), sep="")),
						#ALT = DNAStringSet(paste(ifelse(isDel | isIns, as.character(CV.control[cond.rows - 1]),""), sub("-", "", table$var), sep="")),
						ALT = do.call(DNAStringSetList,as.list(paste(ifelse(isDel | isIns, as.character(CV.control[cond.rows - 1]),""), sub("-", "", table$var), sep=""))),
						QUAL = round(-10*log10(table$raw.p.val)),
						FILTER = "PASS"
				),
				info = DataFrame(
						VF = table$freq.var,
						VFV = table$sigma2.freq.var),
				geno = SimpleList(
						FW = cbind(table$n.tst.fw, table$n.ctrl.fw),
						BW = cbind(table$n.tst.bw, table$n.ctrl.bw),
						DFW = cbind(table$cov.tst.fw, table$cov.ctrl.fw),
						DBW = cbind(table$cov.tst.bw, table$cov.ctrl.bw)),
				exptData = list(header = VCFHeader(
						reference = reference(headerTemplate),
						samples = as.character(samples),
						header = append(header(headerTemplate), .makeVCFheader("date", paste(Sys.time()))))),
				colData = DataFrame(samples=1:length(samples), row.names=samples),
				collapsed=TRUE
		)
		
		return(sort(v))
	}
}

#' Manhattan plot.
#' 
#' This functions performs a Manhattan plot of the p-values of a deepSNV test against the position
#' @param x An \code{\link{deepSNV}} object.
#' @param col An optional vector of colors for the nucleotides.
#' @return NULL.
#' @author Moritz Gerstung
#' @export
#' @examples data(HIVmix)
#' manhattanPlot(HIVmix)
manhattanPlot <- function(x, col=nt.col){
	if(is.null(col)) col <- 1:ncol(x@p.val)
	plot((1:length(x@p.val) -1) %% nrow(x@p.val) +1 ,
				1/x@p.val,
				pch=1, 
				col=col[(1:length(x@p.val) -1) %/% nrow(x@p.val) +1  ], log="y", xlab="index", ylab="1/p-value")
	legend("topleft",1, 
			colnames(slot(x,"test")),
			pch=1, 
			pt.cex=sqrt(2), 
			col=col, 
			bty='n',
			title = "nt")
}

#' Read nucleotide counts from a .bam file
#' 
#' This function uses a C interface to read the nucleotide counts on each position of a .bam alignment. The counts of both strands are reported separately 
#' and nucleotides below a quality cutoff are masked. It is called by \code{\link{deepSNV}} to parse the alignments of the test and control experiments,
#' respectively. 
#' 
#' @param file The name of the .bam file as a string.
#' @param chr The chromosome as a string.
#' @param start The start position (1-indexed).
#' @param stop The end position (1-indexed).
#' @param q An optional cutoff for the nucleotide Phred quality. Default q = 25. Nucleotides with Q < q will be masked by 'N'.
#' @param mq An optional cutoff for the read mapping quality. Default mq = 0 (no filter). reads with MQ < mq will be discarded.
#' @param s Optional choice of the strand. Defaults to s = 2 (both).
#' @param head.clip Should n nucleotides from the head of reads be clipped? Default 0.
#' @param max.depth The maximal depth for the pileup command. Default 1,000,000.
#' @param verbose Boolean. Set to TRUE if you want to get additional output.
#' @param mask Integer indicating which flags to filter. Default 0 (no mask). Try 3844  (UNMAP|SECONDARY|QCFAIL|DUP|SUPPLEMENTARY).
#' @param keepflag Integer indicating which flags to keep. Default 0 (no mask). Try 3  (PAIRED|PROPERLY_PAIRED).
#' @param max.mismatches Integer indicating maximum NM value to allow in a read. Default NULL (no filter).
#' @return A named \code{\link{matrix}} with rows corresponding to genomic positions and columns for the nucleotide counts (A, T, C, G, -), masked nucleotides (N), (INS)ertions, (DEL)etions, (HEAD)s and (TAIL)s that count how often a read begins and ends at the given position, respectively, 
#' and the sum of alignment (QUAL)ities, which can be indicative of alignment problems. 
#' Counts from matches on the reference strand (s=0) are uppercase, counts on the complement (s=1) are lowercase. The returned matrix has 11 * 2 (strands) = 22 columns and (stop - start + 1) rows.
#' @author Moritz Gerstung
#' @export bam2R
#' @examples ## Simple example:
#' counts <- bam2R(file = system.file("extdata", "test.bam", package="deepSNV"), chr="B.FR.83.HXB2_LAI_IIIB_BRU_K034", start = 3120, stop=3140, q = 10, mask = 3844)
#' show(counts)
#' ## Not run: Requires an internet connection, but try yourself.
#' # bam <- bam2R(file = "http://www.bsse.ethz.ch/cbg/software/deepSNV/data/test.bam", chr="B.FR.83.HXB2_LAI_IIIB_BRU_K034", start = 2074, stop=3585, q=10, mask = 3844)
#' # head(bam)
bam2R = function(file, chr, start, stop, q=25, mq=0, s=2, head.clip = 0, max.depth=1000000, verbose=FALSE, mask=0, keepflag=0, max.mismatches=NULL){
	if(is.null(max.mismatches)) max.mismatches <- -1
	region = paste(chr,":",start,"-",stop, sep="")
	result = .C("bam2R",
			as.character(file),
			as.character(chr),
			as.integer(start),
			as.integer(stop),
			vector("integer",(stop-start+1)*11*2),
			as.integer(q),
			as.integer(mq),
			as.integer(s),
			as.integer(head.clip),
			as.integer(max.depth),
			as.integer(verbose),
			as.integer(mask),
      as.integer(keepflag),
      as.integer(max.mismatches),
			PACKAGE="deepSNV"
	)[[5]]
	return(matrix(result, nrow=stop-start+1, 
					dimnames = list(NULL,c('A','T','C','G','-','N','INS','DEL','HEAD','TAIL','QUAL','a','t','c','g','_','n','ins','del','head','tail','qual'))))
}

#' Combine two p-values
#' 
#' This function combines two P-values into a single one using a statistic defined by method. 
#' "fisher" uses the product of the two, in this case the logarithm of the product is
#' \eqn{\chi^2_4} distributed. If the method = "max", the resulting P-value is \eqn{\max\{P_1,P_2\}^2}. 
#' For method = "average" the mean is used, yielding a P-value of \eqn{2 x^2}{2 x^2} if \eqn{x=(P_1+P_2)/2 < .5}{x=(P_1+P_2)/2 < .5}
#' and  \eqn{1-2 x^2}{1-2 x^2} otherwise. "negfisher" is the negative of Fisher's method using $1-F(1-P_1, 1-P_2)$, where $F$ is the combination 
#' function of Fisher's method; for small $P_1,P_2$, the result is very similar to method="average". Fisher's method behaves a bit like a logical AND
#' of the joint null-hypothesis, whereas negative Fisher is like an OR.
#' @param p1 P-value 1
#' @param p2 P-value 2
#' @param method One of "fisher" (default), "max" or "average"
#' @return p-values
#' @examples p1 <- runif(1000)
#' p2 <- runif(1000)
#' hist(p1)
#' p.avg = p.combine(p1,p2, method="average")
#' hist(p.avg)
#' p.fish = p.combine(p1,p2, method="fisher")
#' hist(p.fish)
#' p.max = p.combine(p1,p2, method="max")
#' hist(p.max)
#' pairs(data.frame(p1,p2,p.fish,p.max,p.avg))
#' @author Moritz Gerstung
#' @export
p.combine <- function(p1,p2, method=c("fisher", "max", "average", "prod", "negfisher")){
	method <- match.arg(method)
	if(method == "prod")
		p <- ifelse(p1*p2 !=0, pgamma(-log(p1 * p2), shape = 2,scale=1, lower.tail=F), 0)
	if(method == "fisher")
		p <- ifelse(p1*p2 !=0, p1 * p2 * (1 - log(p1 * p2)), 0)
	if(method == "negfisher"){
		p1 <- 1-p1
		p2 <- 1-p2
		p <- 1 - ifelse(p1*p2 !=0, p1 * p2 * (1 - log(p1 * p2)), 0)
	}
	if(method == "max")
		p <- pmax(p1, p2)^2
	if(method == "average"){
		ptriangle <- function(x) ifelse(x<0.5, 2*x^2, 1-2*(1-x)^2)
		p <- ptriangle(0.5*(p1 + p2))
	}
	return(p)
}
