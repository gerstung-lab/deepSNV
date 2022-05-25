# Functions of the Bayesian beta-binomial test (aka shearwater) for subclonal mutation in large cohorts
# 
# Author: mg14
###############################################################################


#' Function to load all data from a list of bam files
#' 
#' This function uses the parallel package and the bam2R interface to load all nucleotide counts from a list of bam files and a set of regions into a large array.
#' @param files A character vector with the paths to all bam files
#' @param regions Either a GRanges or data.frame with the coordinates of interest
#' @param ... Arguments passed to bam2R
#' @param mc.cores Number of cores used for loading, default = 1
#' @return counts
#' 
#' @author mg14
#' @note Experimental code, subject to changes
#' @export
loadAllData = function(files, regions, ..., mc.cores=1){
	if(class(regions) == "GRanges"){
		regions = as.data.frame(regions)[,1:3]
		colnames(regions) = c("chr", "start", "stop")
	}
	nucleotides = c("A","T","C","G","-","a","t","c","g","_")
	lengths = regions$stop-regions$start +1
	rows = sum(lengths)
	beg = cumsum(c(1,lengths[-length(lengths)]))
	end = cumsum(lengths)
	coordinates = data.frame(chr=unlist(sapply(1:nrow(regions), function(i) rep(regions$chr[i], lengths[i]))),
			pos = unlist(sapply(1:nrow(regions), function(i) regions$start[i]:regions$stop[i]))
	)
	c = mclapply(files, function(f){
				test.matrix <- matrix(0, ncol=length(nucleotides), nrow=rows)
				for (j in 1:nrow(regions)){
					test.matrix[beg[j]:end[j],] =  bam2R(f, regions$chr[j], regions$start[j], regions$stop[j], ...)[,nucleotides]
				}
				mode(test.matrix) <- "integer"
				test.matrix
			}, mc.cores=mc.cores)
	counts = array(0, c(length(files), rows, length(nucleotides)))
	mode(counts) <- "integer"
	for(i in 1:length(files))
		counts[i,,] = c[[i]]
	return(counts)
}


#' Little helper function to split the count objects into a smaller digestible chunks and run function FUN on each subset
#' @param FUN The function to call on each chunk
#' @param X The object to be subsetted using [,i,]
#' @param split The size of each chunk
#' @param mc.cores The number of cores to use
#' @param ... Additional arguments passed to FUN
#' @return The value of FUN
#' @export
#' 
#' @author mg14
#' @note Experimental code, subject to changes

mcChunk = function(FUN, X, split=250, mc.cores=1, ...){
	f = match.fun(FUN)
	ix = rep(1:ceiling(ncol(X)/split), each=split)[1:ncol(X)]
	tmp = f(X[,1:10,])
	p.val = array(dim=c(nrow(X), ncol(X), dim(tmp)[3]))
	res = mclapply(1:ceiling(ncol(X)/split), function(i) f(X[,ix==i,], ...), mc.cores=mc.cores)
	for(i in 1:length(res)){
		p.val[,ix==i,] = res[[i]]
		gc()
	}
	return(p.val)
}


#' Function to create a \code{\link{VCF}} object with variant calls from an array of Bayes factors.
#' 
#' This function thresholds the Bayes factors computed by the shearwater algorithm and creates a \code{\link{VCF}} object as output.
#' @param BF array of Bayes factors from \code{\link{bbb}}.
#' @param counts array of counts from \code{\link{loadAllData}}.
#' @param regions \code{\link{GRanges}} with the regions corresponding to counts and BF.
#' @param samples vector of samples names.
#' @param cutoff Cutoff for the posterior artifact probability below which a variant is considered to be true (default = 0.05)
#' @param prior matrix of prior probabilities for finding a true call, typically from \code{\link{makePrior}}. Alternatively a single fixed number.
#' @param mvcf boolean flag, if TRUE compute a large VCF with as many genotype columns as samples. Default TRUE. Otherwise use duplicate rows and only one genotype column. The sample is then provided by the info:PD field. Can be inefficient for large sample sizes.
#' @param err Optional matrix of error rates, otherwise recomputed from counts.
#' @param mu Optional matrix of relative frequencies, otherwise recomputed from counts.
#' @return A \code{\link{VCF}} object
#' 
#' @author mg14
#' @note Experimental code, subject to changes
#' @export
## TODO: check for cases with zero result
bf2Vcf <- function(BF, counts, regions, samples = 1:nrow(counts), err = NULL, mu = NULL, cutoff = 0.05, prior = 0.5, mvcf=TRUE){
	coordinates <- regions2Coordinates(regions)
	prior = array(rep(prior, each = length(BF)/length(prior)), dim=dim(BF))
	odds = prior/(1-prior)
	posterior = BF / (BF + odds)
	w = which(posterior < cutoff, arr.ind=TRUE)
	w = w[order(w[,2], w[,3], w[,1]),,drop=FALSE]
	
	## If no variants found, select first possible and set to NULL later (avoids a lot of errors in the following)
	if(dim(w)[1]==0) {
		w = matrix(1,ncol = 3)
		isNull = TRUE
	}else{
		isNull = FALSE
	}
	
	totCounts = counts[,,1:5] + counts[,,6:10]
	
	if(is.null(err)){
		err = colSums(totCounts)
		err = err / (rowSums(err) + .Machine$double.eps)
	}
	if(is.null(mu))
		mu = (totCounts)/rep(rowSums(counts, dims=2)+.Machine$double.eps, 5)
	
	
	#coordinates = regions2Coordinates(regions.GR)
	ref = factor(apply(err,1,which.max), levels = 1:5, labels=c("A","T","C","G","-"))
#	gr = GRanges(paste(w[,1], w[,3]==5, sep="."), IRanges(w[,2], width=1), alt=w[,3], AF = select(w,subcl))
#	o = order(gr)
#	rd = reduce(gr)
#	m = match(gr[o], rd)
#	alt = lapply(split(c("A","T","C","G","")[w[o,3]], m), paste, collapse="")
	if(!mvcf){
		headerTemplate = scanVcfHeader(system.file("extdata", "shearwater.vcf", package="deepSNV"))
		v = VCF(
				rowRanges=GRanges(coordinates$chr[w[,2]], 
						IRanges(coordinates$pos[w[,2]] - (w[,3]==5), width=1 + (w[,3]==5)), ## If del make one longer..
				),
				fixed = DataFrame(
						REF = DNAStringSet(paste(ifelse(w[,3]==5,as.character(ref[w[,2]-1]),""), ref[w[,2]], sep="")),
						ALT = do.call(DNAStringSetList,as.list(paste(ifelse(w[,3]==5,as.character(ref[w[,2]-1]),""), c("A","T","C","G","")[w[,3]], sep=""))),
						#ALT = DNAStringSet(paste(ifelse(w[,3]==5,as.character(ref[w[,2]-1]),""), c("A","T","C","G","")[w[u,3]], sep="")),
						QUAL = round(-10*log10(select(w,posterior))),
						FILTER = "PASS"
				),
				info = DataFrame(
						#paramRangeID=NA,
						PD = samples[w[,1]], 
						AF = select(w, mu),
						ER = select(w[,-1], err),
						FW = select(w, counts[,,1:5]),
						BW = select(w, counts[,,6:10]),	
						DP = select(w[,-3], rowSums(totCounts, dims=2)),
						BF = select(w, BF),
						PI = select(w, prior),
						LEN = 1),
				exptData = list(header = VCFHeader(
						reference = reference(headerTemplate),
						samples = as.character(samples),
						header = append(header(headerTemplate), .makeVCFheader("date", paste(Sys.time()))))),
				collapsed = FALSE
		)}else{
		u = !duplicated(w[,-1, drop=FALSE])
		wu = w[u,,drop=FALSE]
		pp = mapply(function(i,j){ posterior[,i,j]}, wu[,2],wu[,3])
		
		geno = SimpleList(
				GT = t(pp < cutoff)+0,
				GQ = t(round(-10*log10(pp))),
				BF = signif(t(mapply(function(i,j){ BF[,i,j]}, wu[,2],wu[,3])),3),
				VF = signif(t(mapply(function(i,j){ mu[,i,j]}, wu[,2],wu[,3])),3),
				FW = t(mapply(function(i,j){ counts[,i,j]}, wu[,2],wu[,3])),
				BW = t(mapply(function(i,j){ counts[,i,j+5]}, wu[,2],wu[,3])),
#				DP = t(mapply(function(i,j){ rowSums(totCounts[,i,])}, wu[,2],wu[,3]))
				FD = t(mapply(function(i,j){ rowSums(counts[,i,1:5])}, wu[,2],wu[,3])),
				BD = t(mapply(function(i,j){ rowSums(counts[,i,6:10])}, wu[,2],wu[,3]))
		)
		
		#rownames(w) = samples[w[,1]]
		headerTemplate = scanVcfHeader(system.file("extdata", "shearwater2.vcf", package="deepSNV"))
		v = VCF(
				rowRanges=GRanges(coordinates$chr[wu[,2]], 
						IRanges(coordinates$pos[wu[,2]] - (wu[,3]==5), width=1 + (wu[,3]==5)), 
				),
				fixed = DataFrame(
						REF = DNAStringSet(paste(ifelse(wu[,3]==5,as.character(ref[wu[,2]-1]),""), ref[wu[,2]], sep="")), ##TODO: Warning if w[u,3][1]==5 & w[u,2]==1...
						ALT = do.call(DNAStringSetList,as.list(paste(ifelse(wu[,3]==5,as.character(ref[wu[,2]-1]),""), c("A","T","C","G","")[wu[,3]], sep=""))),
						#ALT = DNAStringSet(paste(ifelse(w[u,3]==5,as.character(ref[w[u,2]-1]),""), c("A","T","C","G","")[w[u,3]], sep="")),
						QUAL = 1,
						FILTER = "PASS"
				),
				info = DataFrame(
						#paramRangeID=NA,
						ER = select(wu[,-1], err),
						PI = select(wu, prior),
						AF = rowMeans(geno$GT),
						LEN = 1),
				geno = geno,
				exptData = list(header = VCFHeader(
						reference = reference(headerTemplate),
						samples = as.character(samples),
						header = append(header(headerTemplate), .makeVCFheader("date", paste(Sys.time()))))),
				colData = DataFrame(samples=1:length(samples), row.names=samples),
				collapsed = TRUE
		)
		colnames(v) = samples
	}
	
	## If no variants found set to zero..
	if(isNull)
		v <- v[NULL]
	
	return(v)
}

#' Select elements using an array index
#' @param arr.ind 
#' @param array 
#' @return values
#' @noRd 
#' 
#' @author mg14
#' @note Experimental code, subject to changes
select = function(arr.ind, array){
	ix = (arr.ind-1) %*% cumprod(c(1, dim(array)[-length(dim(array))])) + 1
	array[ix]
}

#' Little helper to compute a data.frame with coordinates
#' @param regions.GR 
#' @return a data.frame with chrom and position for each position within regions
#' @noRd 
#' 
#' @author mg14
#' @note Experimental code, subject to changes
regions2Coordinates <- function(regions.GR) {
	regions = as.data.frame(regions.GR)
	data.frame(chr=unlist(sapply(1:nrow(regions), function(i) rep(regions$seqnames[i], regions$width[i]))),
			pos = unlist(sapply(1:nrow(regions), function(i) regions$start[i]:regions$end[i]))
	)
}

#' Compute a prior from a COSMIC VCF object
#' 
#' This function computes the prior probability of detecting a true variant from a variation data base. It assumes a
#' VCF file with a CNT slot for the count of a given base substitution. Such a VCF file can be downloaded at 
#' ftp://ngs.sanger.ac.uk/production/cosmic/. The prior probability
#' is simply defined as pi.mut * CNT[i]/sum(CNT). On sites with no count, a background probability of pi0 is used.
#' @param COSMIC A VCF object from COSMIC VCF export. 
#' @param regions A GRanges object with the regions (gene) of interest.
#' @param pi.gene Probability that a gene is mutated
#' @param pi.backgr Background probability of a locus being mutated. Default 1e-4, corresponding to an expected value of 1 SNV per 1e4 bases.
#' @return A vector of prior values with length given by the length of the regions GRanges object.
#' @examples ## Make prior (not run)
#' #COSMIC <- readVcf("PATHTO/CosmicCodingMuts_v64_02042013_noLimit.vcf.gz", genome="GChr37")
#' #prior <- makePrior(COSMIC[info(COSMIC)$GENE=="TP53"], regions=GRanges(17, IRanges(7571720,7578811)))
#' #plot(prior[,1], type="h")
#' 
#' @author mg14
#' @note Experimental code, subject to changes
#' @export
makePrior = function(COSMIC, regions, pi.gene = 0.1, pi.backgr = 1e-4){
	c = COSMIC[width(COSMIC) == 1 & !grepl("ins", info(COSMIC)$CDS)]
	cnt = info(c)$CNT
	coordinates = regions2Coordinates(regions)
	var = match(sub("(.+)(\\>|del|ins)(.+)","\\3",info(c)$CDS, perl=TRUE), c("A","T","C","G"))
	var = ifelse(info(c)$STRAND=="+",var,c(2,1,4,3)[var])
	pos = match(start(c),coordinates$pos)
	prior = matrix(0,nrow = nrow(coordinates), ncol=5)
	for(i in 1:nrow(c)){
		prior[pos[i], var[i]] = cnt[i] + prior[pos[i], var[i]]
	}
	prior = prior/sum(prior) * pi.gene + pi.backgr
	prior[prior>=1] = 1 - .Machine$double.eps
	return(prior)
}

#' Helper function for estimating the dispersion factor rho
#' 
#' It uses a method of moments approximation to estimate rho from the variances of the relative frequencies nu across samples
#' @param x counts
#' @param mu relative frequency across all samples
#' @param ix index indicating the set of samples to use (typically indicating those with relative frequency smaller than 0.1).
#' @param pseudo.rho a pseudo count added to each sample to avoid problems with zeros. Default = .Machine$double.eps
#' @return rho
#' 
#' @author mg14
#' @note Experimental code, subject to changes
estimateRho = function(x, mu, ix, pseudo.rho=.Machine$double.eps){
	n = array(rep(rowSums(x, dims=2), dim(x)[3]), dim=dim(x)[1:2])
	Xix = colSums(x * ix, dims=1)
	nix = array(rep(n, dim(x)[3]), dim=dim(x))
	nix[!ix | nix == 0] = NA
	N = colSums(ix)
	nu = (Xix + pseudo.rho) /  (rep(rowSums(colSums(x, dims=1)) + dim(x)[3]* pseudo.rho, dim(x)[3]))
	
	s2 =   N * colSums(nix  * (mu - rep(nu, each = nrow(x)))^2, na.rm=TRUE) / ( (N-1) * colSums(nix, na.rm=TRUE))
	rho = (N * (s2/nu/(1-nu)) - colSums(1/nix, na.rm=TRUE)) / (N -  colSums(1/nix, na.rm=TRUE))
	rho = pmin(pmax(rho, 0), 1)
	return(rho)
}


#' Logarithmic beta binomial without binomial coefficients.
#' @param x Counts
#' @param n Coverage
#' @param mu Relative frequency (scaled by dispersion)
#' @param disp Dispersion factor rho/(1-rho)
#' @return The logarithmic beta binomial
#' @noRd 
#' 
#' @author mg14
#' @note Experimental code, subject to changes
logbb <- function(x, n, mu, disp) {
	lbeta(x + mu, n - x - mu + disp) - lbeta(mu, disp - mu)
}

#' Bayesian beta-binomal test, codename shearwater
#' 
#' This is the workhorse of the shearwater test. It computes the Bayes factor for each sample, nucleotide and position of the null-model vs. the alternative of a real variant.
#' @param counts An \code{\link{array}} of nucleotide counts (samples x positions x 10 nucleotides in forward and reverse orientation), typically from \code{\link{loadAllData}}
#' @param truncate The model uses a compound control sample which is the sum of all samples with a relative nucleotide frequency below truncate at this locus. Default = 0.1.
#' @param alternative The alternative. Currently only "greater" is implemented.
#' @param rho Disperision factor. If NULL, estimated from the data.
#' @param rho.min Lower bound for the method of moment estimate of the dispersion factor rho.
#' @param rho.max Upper bound for the method of moment estimate of the dispersion factor rho.
#' @param mu.min Minimum of the error rate mu.
#' @param mu.max Maximal error rate mu.
#' @param pseudo A pseudo count to be added to the counts to avoid problems with zeros.
#' @param return.value Return value. Either "BF" for Bayes Factor of "P0" for the posterior probability (assuming a prior of 0.5).
#' @param model The null model to use. For "OR" it requires the alternative model to be violated on either of the strands, for "AND" the null is specified such that the error rates of the sample 
#' of interest and the compound control sample are identical on both strands. "AND" typically yield many more calls. The most recent addition is "adaptive", which switches from "OR" to "AND", if the coverage 
#' is less than min.cov, or if the odds of forward and reverse coverage is greater than max.odds. Default = "OR".
#' @param min.cov Minimal coverage to swith from OR to AND, if model is "adaptive"
#' @param max.odds Maximal odds before switching from OR to AND if model is "adaptive" and min.cov=NULL.
#' @return An \code{\link{array}} of Bayes factors
#' @example inst/example/shearwater-example.R
#' 
#' @author mg14
#' @rdname shearwater
#' @aliases shearwater
#' @note Experimental code, subject to changes
#' @export
bbb <- function(counts, rho = NULL, alternative="greater", truncate=0.1, rho.min = 1e-4, rho.max = 0.1, pseudo = .Machine$double.eps, return.value=c("BF","P0", "err"), model=c("OR","AND", "adaptive"), min.cov=NULL, max.odds=10, mu.min = 1e-6, mu.max = 1 - mu.min) {
	pseudo.rho = .Machine$double.eps
	## minum value for rho
	
	model = match.arg(model)
	return.value = match.arg(return.value)
	
	ncol = dim(counts)[3]/2
	
	x.fw = counts[,,1:ncol, drop=FALSE]
	x.bw = counts[,,1:ncol + ncol, drop=FALSE]
	
	n.fw = rep(rowSums(x.fw, dims=2), dim(x.fw)[3])
	n.bw = rep(rowSums(x.bw, dims=2), dim(x.bw)[3])
	
	
	x <- x.fw+x.bw
	n = array(n.fw + n.bw, dim=dim(x)[1:2])
	mu = (x + pseudo.rho) / (rep(n + ncol*pseudo.rho, dim(x)[3]) )
	ix = (mu < truncate)
	if(is.null(rho)){
		rho = estimateRho(x, mu, ix)
		rho = pmin(pmax(rho, rho.min), rho.max)
		rho[is.na(rho)] = rho.min
		
	}
	X =  colSums(x, dims=1)
	#rm(x)
	
	bound = function(x, xmin, xmax){
		x = pmax(x, xmin)
		x = pmin(x, xmax)
		return(x)
	}
	
	disp = (1-rho)/rho 
	rdisp <- rep(disp, each=nrow(counts))
	mu = (x + pseudo) / (rep(n + ncol*pseudo, dim(x)[3]) ) ## sample rate forward+bwackward (true allele frequency)
	mu = bound(mu, mu.min, mu.max) * rdisp
	tr.fw = x.fw * ix
	
	X.fw = rep(colSums(tr.fw, dims=1), each = nrow(counts)) - tr.fw ## control samples
	N.fw = rep(colSums(n.fw * ix), each = nrow(counts)) - n.fw * ix
	#nu0.fw <- array(rep((colSums(tr.fw) +pseudo) / (colSums(n.fw * ix) + ncol*pseudo) * disp, each = nrow(tr.fw)), dim = dim(tr.fw)) * mumax
	nu0.fw <- (X.fw + x.fw + pseudo)/(N.fw + n.fw + ncol*pseudo)
	nu0.fw <- bound(nu0.fw, mu.min, mu.max)* rdisp
	mu0.bw <- (x.bw+pseudo) / (n.bw + ncol*pseudo) 
	mu0.bw <- bound(mu0.bw, mu.min, mu.max) * rdisp
	nu.fw <- (X.fw+pseudo) / (N.fw + ncol*pseudo)
	nu.fw <- bound(nu.fw, mu.min, mu.max) * rdisp
	
	tr.bw = x.bw * ix
	X.bw = rep(colSums(tr.bw, dims=1), each = nrow(counts)) - tr.bw 
	N.bw = rep(colSums(n.bw * ix), each = nrow(counts)) - n.bw * ix
	#nu0.bw <- array(rep((colSums(tr.bw) + pseudo) / (colSums(n.bw * ix) + ncol*pseudo) * disp, each = nrow(tr.bw)), dim = dim(tr.bw)) * mumax
	nu0.bw <- (X.bw + x.bw + pseudo)/(N.bw + n.bw + ncol*pseudo)
	nu0.bw <- bound(nu0.bw, mu.min, mu.max) * rdisp
	mu0.fw <- (x.fw+pseudo) / (n.fw +  ncol*pseudo)
	mu0.fw <- bound(mu0.fw, mu.min, mu.max) * rdisp
	nu.bw <- (X.bw+pseudo) / (N.bw + ncol*pseudo) 
	nu.bw <- bound(nu.bw, mu.min, mu.max) * rdisp
	
	## Return rates only
	if(return.value == "err"){ 
		nu0 <- (X.bw + tr.fw + X.fw + tr.bw+ pseudo)/(N.bw + n.bw +N.fw +n.fw+ ncol*pseudo)
		nu0 <- bound(nu0, mu.min, mu.max)
		return(list(nu = nu0[1,,], nu.fw=(nu0.fw/rdisp)[1,,], nu.bw=(nu0.bw/rdisp)[1,,], rho=rho))
	}
	rm(tr.fw)
	rm(tr.bw)
	
	## Enforce mu > nu
	mu = pmax(mu, nu0.fw)
	mu = pmax(mu, nu0.bw)
	
	mu0.fw = pmax(mu0.fw, nu0.fw)
	mu0.bw = pmax(mu0.bw, nu0.bw)
	
	if(model %in% c("OR","adaptive")){
		## Bayes factor forward
		Bf.fw <- logbb(x.fw, n.fw, nu0.fw, rdisp) + logbb(x.bw, n.bw, mu0.bw, rdisp) + logbb(X.fw, N.fw, nu0.fw, rdisp) - logbb(x.fw, n.fw, mu, rdisp) - logbb(x.bw, n.bw, mu, rdisp) - logbb(X.fw, N.fw, nu.fw, rdisp)
		Bf.fw = exp(Bf.fw)
		
		Bf.both = logbb(x.fw, n.fw, nu0.fw, rdisp) + logbb(X.fw, N.fw, nu0.fw, rdisp) - logbb(x.fw, n.fw, mu, rdisp) - logbb(X.fw, N.fw, nu.fw, rdisp)
		
		rm(X.fw, N.fw, mu0.bw, nu.fw)
		
		## Bayes factor reverse
		Bf.bw <- logbb(x.fw, n.fw, mu0.fw, rdisp) + logbb(x.bw, n.bw, nu0.bw, rdisp) + logbb(X.bw, N.bw, nu0.bw, rdisp) - logbb(x.fw, n.fw, mu, rdisp) - logbb(x.bw, n.bw, mu, rdisp) - logbb(X.bw, N.bw, nu.bw, rdisp)
		Bf.bw = exp(Bf.bw)
		
		Bf.both = Bf.both + logbb(x.bw, n.bw, nu0.bw, rdisp) + logbb(X.bw, N.bw, nu0.bw, rdisp) - logbb(x.bw, n.bw, mu, rdisp) - logbb(X.bw, N.bw, nu.bw, rdisp)
		Bf.both = exp(Bf.both)
		
		rm(X.bw, N.bw, mu0.fw, nu.bw)
		
		#f.se.g <- mu < 0.5*(nu0.fw + nu0.bw)	
		
		rm(mu, nu0.fw, nu0.bw)	
		
		Bf = Bf.fw + Bf.bw - Bf.both + .Machine$double.xmin
	}else{
		Bf.both = logbb(x.fw, n.fw, nu0.fw, rdisp) + logbb(X.fw, N.fw, nu0.fw, rdisp) - logbb(x.fw, n.fw, mu, rdisp) - logbb(X.fw, N.fw, nu.fw, rdisp) + logbb(x.bw, n.bw, nu0.bw, rdisp) + logbb(X.bw, N.bw, nu0.bw, rdisp) - logbb(x.bw, n.bw, mu, rdisp) - logbb(X.bw, N.bw, nu.bw, rdisp)
		Bf = exp(Bf.both)
	}
	
	if(model=="adaptive"){
		if(!is.null(min.cov))
			ix <- n.fw < min.cov | n.bw < min.cov
		else
			ix <- na.omit(abs(log10(n.fw/n.bw)) > log10(max.odds))
		Bf[ix] <- Bf.both[ix]
	}
	
	#Bf[f.se.g] = Inf
	## If smaller, take H0
	cons= apply(X,1, which.max)
	for(i in 1:ncol(Bf))
		Bf[,i,cons[i]] = NA
	
	if(return.value=="P0"){
		return(Bf/(1+Bf))
	}else{
		return(Bf)
	}
}


#' Get mutation IDs
#' @param vcf 
#' @return character
#' @noRd
#' 
#' @author mg14
mutID = function(vcf){
	alt <- as.character(unlist(alt(vcf)))
	ref <- as.character(ref(vcf))
	isDel <- width(ref) > 1
	ifelse(!isDel, 
			paste(seqnames(vcf), start(vcf), alt, sep=":"),
			paste(seqnames(vcf), start(vcf)+1, paste("del",substring(ref,2),sep=""), sep=":"))
}

