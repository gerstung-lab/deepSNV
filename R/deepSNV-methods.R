# Special methods of the deepSNV package
# 
# Author: Moritz Gerstung
###############################################################################


#' Test two matched deep sequencing experiments for low-frequency SNVs.
#' 
#' This generic function can handle different types of inputs for the test and control experiments. It either reads from two .bam files,
#' uses two matrices of nucleotide counts, or re-evaluates the test results from a \code{\link{deepSNV-class}} object. The actual test is a
#' likelihood ratio test of a (beta-)binomial model for the individual nucleotide counts on each position under the hypothesis that both experiments share the same parameter,
#' and the alternative that the parameters differ. Because the difference in degrees of freedom is 1, the test statistic \eqn{D = -2 \log \max{L_0}/\max{L_1}}
#' is asymptotically distributed as \eqn{\chi_1^2}. The statistic may be tuned by a nucleotide specific Dirichlet prior that is learned across all genomic sites, 
#' see \code{\link{estimateDirichlet}}. If the model is beta-binomial, a global dispersion parameter is used for all sites. It can be learned with \code{\link{estimateDispersion}}.
#' @param test The test experiment. Either a .bam file, or a matrix with nucleotide counts, or  a \code{\link{deepSNV-class}} object.
#' @param control The control experiment. Must be of the same type as test, or missing if test is  a \code{\link{deepSNV-class}} object.
#' @param alternative The alternative to be tested. One of greater, less, or two.sided.
#' @param model Which model to use. Either "bin", or "betabin". Default "bin".
#' @param dirichlet.prior A base-sepecific Dirichlet prior specified as a matrix. Default NULL.
#' @param pseudo.count If dirichlet.prior=NULL, a pseudocount can be used to define a flat prior.
#' @param over.dispersion A numeric factor for the over.dispersion, if the model is beta-binomial. Default 100.
#' @param combine.method The method to combine p-values. One of "fisher" (default), "max", or "average". See \code{\link{p.combine}} for details.
#' @param regions The regions to be parsed if test and control are .bam files. Either a \code{\link{data.frame}} with columns "chr" (chromosome), 
#' "start", "stop", or a \code{\link{GRanges}} object. If multiple regions are specified, the appropriate slots of the returned object are concatenated by row.
#' @param q The quality arguement  passed to \code{\link{bam2R}} if the experiments are .bam files.
#' @param s The strand argument  passed to \code{\link{bam2R}} if the experiments are .bam files.
#' @param head.clip The head.clip argument  passed to \code{\link{bam2R}} if the experiments are .bam files.
#' @param ... Additional arguments.
#' @return A \code{\linkS4class{deepSNV}} object
#' @author Moritz Gerstung
#' @example inst/example/deepSNV-example.R
#' @exportMethod deepSNV
#' @rdname deepSNV-methods
#' @docType methods
setGeneric("deepSNV",
		function(test, control, ...) {
			standardGeneric("deepSNV")
		})

#' @rdname deepSNV-methods
#' @aliases deepSNV,matrix,matrix-method
#' @usage \S4method{deepSNV}{matrix,matrix}(test,control, alternative = c('greater', 'less', 'two.sided'), dirichlet.prior = NULL, pseudo.count=1, combine.method = c("fisher", "max", "average"), over.dispersion = 100, model = c("bin", "betabin"), ...)
setMethod("deepSNV",
		signature = signature(test="matrix", control="matrix"),
		function(test,control, alternative = c('greater', 'less', 'two.sided'), dirichlet.prior = NULL, pseudo.count=1, combine.method = c("fisher", "max", "average"), over.dispersion = 100, model = c("bin", "betabin"), ...){
			alternative = match.arg(alternative)
			combine.method = match.arg(combine.method)
			model = match.arg(model)
			
			stopifnot(all(colnames(test) == colnames(control)), all(dim(test) == dim(control)), !(is.null(colnames(test)) & is.null(colnames(control))))
			nucleotides <- colnames(test)
			test <- as.matrix(test)
			control <- as.matrix(control)
			if (is.null(dirichlet.prior)){
				dirichlet.prior <- matrix(pseudo.count, nrow = 5, ncol = 5)
				rownames(dirichlet.prior) <- colnames(dirichlet.prior) <- c("A","T","C","G","-")
			}
			
			#Test
			res = .deepSNV(test = test, control = control, nucleotides = nucleotides, dirichlet.prior = dirichlet.prior, alternative = alternative, model = model, over.dispersion = over.dispersion, combine.method = combine.method, pseudo.count = pseudo.count)
			
			#Collate results
			result = new("deepSNV", 
					control=control, 
					test=test, 
					p.val = res$p.val, 
					dirichlet.prior=dirichlet.prior,
					alternative=alternative,
					nucleotides=nucleotides,
					combine.method=combine.method,
					model = model,
					over.dispersion=over.dispersion,
					log.lik = res$log.lik,
					pseudo.count = pseudo.count,
					...	)
			result@call = match.call()
			return(result)
		}
)

#' @rdname deepSNV-methods
#' @aliases deepSNV,deepSNV,missing-method
#' @usage \S4method{deepSNV}{deepSNV,missing}(test, control, ...)
setMethod("deepSNV",
		signature = signature(test="deepSNV", control="missing"),
		function(test, control, ...) {
			newobject = initialize(test, ...)
			res = .deepSNV(test = newobject@test, 
					control = newobject@control, 
					nucleotides = newobject@nucleotides,
					alternative = newobject@alternative, 
					dirichlet.prior = newobject@dirichlet.prior, 
					combine.method=newobject@combine.method, 
					model=newobject@model, 
					over.dispersion=newobject@over.dispersion,
					pseudo.count = newobject@pseudo.count)
			newobject@p.val = res$p.val
			newobject@log.lik = res$log.lik
			newobject@call = match.call()
			return(newobject)
		}
)

#' @rdname deepSNV-methods
#' @aliases deepSNV,character,character-method
#' @usage \S4method{deepSNV}{character,character}(test, control, regions, q=25, s=2, head.clip=0, ...)
setMethod("deepSNV",
		signature = signature(test="character", control="character"),
		function(test, control, regions, q=25, s=2, head.clip=0, ...) {
			stopifnot(class(regions) %in% c("data.frame","GRanges"))
			if(class(regions) == "GRanges"){
				regions = as.data.frame(regions)[,1:3]
				colnames(regions) = c("chr", "start", "stop")
			}
			nucleotides = c("A","T","C","G","-","a","t","c","g","_")
			lengths = regions$stop-regions$start +1
			rows = sum(lengths)
			test.matrix <- control.matrix <- matrix(0, ncol=length(nucleotides), nrow=rows)
			beg = cumsum(c(1,lengths[-length(lengths)]))
			end = cumsum(lengths)
			coordinates = data.frame(chr=unlist(sapply(1:nrow(regions), function(i) rep(regions$chr[i], lengths[i]))),
					pos = unlist(sapply(1:nrow(regions), function(i) regions$start[i]:regions$stop[i]))
			)
			for (j in 1:nrow(regions)){
				test.matrix[beg[j]:end[j],] =  bam2R(test, regions$chr[j], regions$start[j], regions$stop[j], q=q, s=s, head.clip=head.clip)[,nucleotides]
			}				
			for (j in 1:nrow(regions)){
				control.matrix[beg[j]:end[j],] = bam2R(control, regions$chr[j], regions$start[j], regions$stop[j], q=q, s=s, head.clip=head.clip)[,nucleotides]
			}
			colnames(test.matrix) <- colnames(control.matrix) <- nucleotides
			result = deepSNV(test.matrix, control.matrix, coordinates=coordinates, regions=regions, ...)
			result@call = match.call()
			result@files = list(test=test, control=control)
			return(result)
		})

#' @rdname deepSNV-methods
#' @aliases deepSNV,matrix,character-method
#' @usage \S4method{deepSNV}{matrix,character}(test, control, regions, q=25, s=2, ...)
setMethod("deepSNV",
		signature = signature(test="matrix", control="character"),
		function(test, control, regions, q=25, s=2, ...) {
			stopifnot(class(regions) %in% c("data.frame","GRanges"))
			if(class(regions) == "GRanges"){
				regions = as.data.frame(regions)[,1:3]
				colnames(regions) = c("chr", "start", "stop")
			}			
			nucleotides = c("A","T","C","G","-","a","t","c","g","_")
			lengths = regions$stop-regions$start +1
			rows = sum(lengths)
			control.matrix <- matrix(0, ncol=length(nucleotides), nrow=rows)
			beg = cumsum(c(1,lengths[-length(lengths)]))
			end = cumsum(lengths)
			coordinates = data.frame(chr=unlist(sapply(1:nrow(regions), function(i) rep(regions$chr[i], lengths[i]))),
					pos = unlist(sapply(1:nrow(regions), function(i) regions$start[i]:regions$stop[i]))
			)
			test.matrix <- test
			for (j in 1:nrow(regions)){
				control.matrix[beg[j]:end[j],] = bam2R(control, regions$chr[j], regions$start[j], regions$stop[j], q=q, s=s)[,nucleotides]
			}
			colnames(control.matrix) <- nucleotides
			result = deepSNV(test, control.matrix, coordinates=coordinates, regions=regions, ...) 
			result@files$control = control
			result@call = match.call()			
			return(result)
		})

#' @rdname deepSNV-methods
#' @aliases deepSNV,character,matrix-method
#' @usage \S4method{deepSNV}{character,matrix}(test, control, regions, q=25, s=2, ...)
setMethod("deepSNV",
		signature = signature(test="character", control="matrix"),
		function(test, control, regions, q=25, s=2, ...) {
			stopifnot(class(regions) %in% c("data.frame","GRanges"))
			if(class(regions) == "GRanges"){
				regions = as.data.frame(regions)[,1:3]
				colnames(regions) = c("chr", "start", "stop")
			}			
			nucleotides = c("A","T","C","G","-","a","t","c","g","_")
			lengths = regions$stop-regions$start +1
			rows = sum(lengths)
			test.matrix <- matrix(0, ncol=length(nucleotides), nrow=rows)
			beg = cumsum(c(1,lengths[-length(lengths)]))
			end = cumsum(lengths)
			coordinates = data.frame(chr=unlist(sapply(1:nrow(regions), function(i) rep(regions$chr[i], lengths[i]))),
					pos = unlist(sapply(1:nrow(regions), function(i) regions$start[i]:regions$stop[i]))
			)
			for (j in 1:nrow(regions)){
				test.matrix[beg[j]:end[j],] =  bam2R(test, regions$chr[j], regions$start[j], regions$stop[j], q=q, s=s)[,nucleotides]
			}
			control.matrix <- control
			colnames(test.matrix) <- nucleotides
			result = deepSNV(test.matrix, control, coordinates=coordinates, regions=regions, ...) 
			result@files$test = test
			result@call = match.call()
			return(result)
		})

#' Learn a base-specific Dirichlet prior.
#' 
#' The prior learns the parameters of a Dirichlet distribution seperately for each consensus base. The expected value of the Dirichlet distributions 
#' is the base-substitution matrix, where rows correspond to the initial nucleotide and columns to the substituted nucleotide. The absolute values determine 
#' the higher moments of the Dirichlet distributions. After having learned the prior the \code{\link{deepSNV-class}} test is recomputed.
#' @param control Either a matrix with nucleotide counts or a \code{\link{deepSNV-class}} object.  
#' @return An \code{\link{deepSNV-class}} object.
#' @author Moritz Gerstung
#' @exportMethod estimateDirichlet
#' @rdname estimateDirichlet-methods
#' @docType methods
#' @examples data(phiX)
#' estimateDirichlet(phiX)
setGeneric("estimateDirichlet", function(control) standardGeneric("estimateDirichlet"))

#' @rdname estimateDirichlet-methods
#' @aliases estimateDirichlet,matrix-method
setMethod("estimateDirichlet", 
		signature = signature(control="matrix"),
		function(control){
			y <- c("A","T","C","G","-")
			CV.control <- consensusSequence(control[,1:5]+control[,6:10], vector=TRUE)
			prior = matrix(1, nrow = length(y), ncol=length(y))
			colnames(prior) <- rownames(prior) <- y
			for(nt in unique(CV.control)) {
				if(sum(CV.control==nt) > 1){ #Need at least two data points!
					data <- data.matrix(rbind(
									RF(control[CV.control==nt,1:5]+1),
									RF(control[CV.control==COMPLEMENT[nt],c(7,6,9,8,10)]+1)
							))
					fit <- vglm(data ~ 1, dirichlet, crit="c")
					prior[nt,] = exp(coefficients(fit))
				}
			}
			return(prior)}
)

#' @rdname estimateDirichlet-methods
#' @aliases estimateDirichlet,deepSNV-method
setMethod("estimateDirichlet", 
		signature = signature(control="deepSNV"),
		function(control) {
			prior = estimateDirichlet(control@control)
			deepSNV(initialize(control, dirichlet.prior = prior))
		})




#' Mask homopolymeric repeats.
#' 
#' This function masks homopolymeric repeats longer than a given width. These are hot-spots of sequencing error and can confound the analysis.
#' @param x An object. Either a \code{\link{deepSNV-class}} object or a \code{\link{DNAString}} with the nucleotide sequence.
#' @param flank Boolean. Indicates whether the sites adjacent to the repeat should also be masked.
#' @param w Integer. The minimal length at which repeats should be masked. Default \code{w=0}.
#' @param ... Additional param passed to specific methods
#' @return A boolean vector where TRUE indicates a non-homopolymeric region.
#' @author Moritz Gerstung
#' @exportMethod repeatMask
#' @rdname repeatMask-methods
#' @docType methods
#' @examples data(HIVmix)
#' which(repeatMask(HIVmix))
setGeneric("repeatMask", function(x, ...) standardGeneric("repeatMask"))

#' @rdname repeatMask-methods
#' @aliases repeatMask,DNAString-method
#' @usage \S4method{repeatMask}{DNAString}(x, w=5, flank=TRUE)
setMethod("repeatMask",
		signature = signature(x="DNAString"),
		function(x, w=5, flank=TRUE){
			l <- length(x)
			unq <- sapply(w:l, function(i) length(uniqueLetters(x[(i-w + 1):i])) == 1)
			ll <- length(unq)
			mask <- rep(FALSE,l)
			for (i in 1:w)
				mask[i:(l-w+i)] <- mask[i:(l-w+i)] | unq
			if (flank){
				mask[1:(l-w)] <- mask[1:(l-w)] | unq[-1]
				mask[(w+1):l] <- mask[(w+1):l] | unq[1:(l-w)]
			}
			return(!mask)
		}
)

#' @rdname repeatMask-methods
#' @aliases repeatMask,deepSNV-method
#' @usage \S4method{repeatMask}{deepSNV}(x, w=5, flank=TRUE)
setMethod("repeatMask",
		signature = signature(x="deepSNV"),
		function(x, w=5, flank=TRUE) repeatMask(consensusSequence(control(x, total=TRUE)), w, flank)
)

#' Calculate the consensus sequence.
#' 
#' This function computes the consensus sequence from a matrix of nucleotide counts, or the control slot of a deepSNV object.
#' @param x An object. Either an \code{\link{deepSNV-class}} object, or a named matrix with nucleotide counts.
#' @param vector Boolean where TRUE indicates that a character vector should be returned.
#' @param haploid Should the consensus be called for a haploid control? Otherwise, also all bases larger than het.cut are rerported. Default haploid = TRUE.
#' @param het.cut Heterozygous cutoff. If haploid = FALSE, report all nucleotides with relative frequency larger than het.cut. Default = 0.333.
#' @param ... Additional arguments passed to methods.
#' @return A \code{\link{DNAString}} with the consensus sequence, or if vector = TRUE, a character vector.
#' @author Moritz Gerstung
#' @exportMethod consensusSequence
#' @rdname consensusSequence-methods
#' @docType methods
#' @examples data(HIVmix)
#' seq = consensusSequence(HIVmix)
#' consensusSequence(HIVmix, vector=TRUE)[1:10]
setGeneric("consensusSequence", 
		function(x, ...) standardGeneric("consensusSequence"))

#' @rdname consensusSequence-methods
#' @aliases consensusSequence,matrix-method
#' @usage \S4method{consensusSequence}{matrix}(x, vector=FALSE, haploid=TRUE, het.cut = .333)
setMethod("consensusSequence", 
		signature = signature(x="matrix"), 
		function(x, vector=FALSE, haploid=TRUE, het.cut = .333){
			if(haploid){
				cons <- factor(apply(x,1,which.max), levels = 1:ncol(x))
				levels(cons) <- colnames(x)
				if(vector)  return(cons)
				else  return(DNAString(paste(as.character(cons), collapse="")))
		}else{
			cons <- apply(x,1, function(c) {i = c > het.cut; sum(i * 2^(0:4))})
			u <- sort(unique(cons))
			cons <- factor(cons, levels=u, labels = sapply(u, function(i){paste(colnames(x)[(i %/% 2^(0:4)) %% 2 == 1], collapse="/") }))
			return(cons)
		}
		})

#' @rdname consensusSequence-methods
#' @aliases consensusSequence,deepSNV-method
#' @usage \S4method{consensusSequence}{deepSNV}(x, vector=FALSE, haploid=TRUE, het.cut = .333)
setMethod("consensusSequence", signature = signature(x="deepSNV"), function(x, vector=FALSE, haploid=TRUE, het.cut = .333) consensusSequence(control(x, total=TRUE), vector, haploid, het.cut))


#' Estimate the Dispersion factor in a beta-binomial model.
#' 
#' This function estimates the dispersion factor in a beta-binomial model of the nucleotide counts. This model assumes that the count for nucleotide j at position i is
#' distributed after a beta-binomial \eqn{X_{i,j}\sim \mathrm{BB}(n_i; \alpha, \beta_{ij})}{X_ib ~ BB(n_i; alpha, beta_ij)}, where \eqn{n_i}{n_i} is the coverage.
#' The base and nucleotide specific parameter \eqn{\beta_{ij}}{beta_ij} is estimated from the local mean by the method-of-moments estimate, \eqn{\alpha}{alpha} is a shared
#' overdispersion parameter. It is estimated via a numerical optimization of the likelihood under the null-hypothesis. 
#' @param test Either a deepSNV object, or a matrix with the test counts.
#' @param control Missing if test is a deepSNV object, otherwise missing.
#' @param alternative The alternative to be tested. One of "greater", "less", "two-sided" (default). If test is a deepSNV object, automatically taken from the corresponding slot
#' if unspecified. 
#' @param interval The interval to be screened for the overdispersion factor. Default (0,1000).
#' @param ... Additional param passed to specific methods
#' @return A \code{\link{deepSNV-class}} object if the input was a deepSNV object. Otherwise the loglikelihood and the estimated parameter.
#' @author Moritz Gerstung
#' @examples data("RCC", package="deepSNV")
#' plot(RCC)
#' summary(RCC)[,1:6]
#' RCC.bb = estimateDispersion(RCC, alternative = "two.sided")
#' summary(RCC.bb)
#' @exportMethod estimateDispersion
#' @rdname estimateDispersion-methods
#' @docType methods
setGeneric("estimateDispersion", function(test, control, ...) standardGeneric("estimateDispersion"))

#' @rdname estimateDispersion-methods
#' @aliases estimateDispersion,deepSNV,missing-method
#' @usage \S4method{estimateDispersion}{deepSNV,missing}(test, control, alternative = NULL, interval = c(0,1000))
setMethod("estimateDispersion", 
		signature = signature(test = "deepSNV", control = "missing"),
		function(test, control, alternative = NULL, interval = c(0,1000)){
			if(test@model =="bin") cat("Note: The initial object used a binomial model. Will be changed to beta-binomial.\n")
			res = .estimateDispersion(test = test@test, control = test@control, dirichlet.prior = test@dirichlet.prior, alternative = ifelse(is.null(alternative),test@alternative,alternative), interval = interval)
			cat(paste("Estimated dispersion factor", res$maximum, "\n"))
			deepSNV(initialize(test), over.dispersion = res$maximum, model="betabin")
		}
)

#' @rdname estimateDispersion-methods
#' @aliases estimateDispersion,matrix,matrix-method
#' @usage \S4method{estimateDispersion}{matrix,matrix}(test, control, alternative = NULL, interval = c(0,1000))
setMethod("estimateDispersion", 
		signature = signature(test = "matrix", control = "matrix"),
		function(test,control, alternative = NULL, interval = c(0,1000)) .estimateDispersion(test=test, control=control, alternative=alternative, interval=interval)
)
		
		
#' Normalize nucleotide counts.
#' 
#' This functions performs a \code{\link{loess}} normalization of the nucleotide. This experimental feature can 
#' be used to compare experiments from different libraries or sequencing runs that may have differing noise characteristics.
#' @note This feature is somewhat experimental and the results should be treated with care. Sometimes it can be better to leave the data unnormalized and use a model with greater dispersion instead.
#' @param test Either an \code{\link{deepSNV-class}} object or a named matrix with nucleotide counts.
#' @param control Missing if test is an \code{link{deepSNV-class}} object, otherwise a matrix with nucleotide counts.
#' @param round Logical. Should normalized counts be rounded to integers? Default=TRUE 
#' @param ... Parameters passed to \code{\link{loess}}.
#' @return A \code{\link{deepSNV-class}} object.
#' @author Moritz Gerstung
#' @examples data(phiX, package = "deepSNV")
#' plot(phiX)
#' phiN <- normalize(phiX, round = TRUE)
#' plot(phiN)
#' @rdname normalize-methods
#' @exportMethod normalize
		setGeneric("normalize",
				function(test, control, ...) {
					standardGeneric("normalize")
				})
		
#' Normalize nucleotide counts.
#' @rdname normalize-methods
#' @aliases normalize,matrix,matrix-method
#' @usage \S4method{normalize}{matrix,matrix}(test, control, round=TRUE, ...)
setMethod("normalize",
		signature = signature(test="matrix", control="matrix"),
		function(test, control, round=TRUE, ...) {
			CV <- consensusSequence(control[,1:5]+control[,6:10], vector=TRUE)
			control.rf <- RF(control+.5)
			control.sum <- rowSums(control +.5)
			test.rf <- RF(test+.5)
			test.sum <- rowSums(test+.5)
			M <- 0.5*(log10(test.rf) + log10(control.rf))
			A <- 0.5*(log10(test.rf) - log10(control.rf))
			#plot(M,A)
			fits <- sapply(colnames(control), function(i){
						fit <- numeric(nrow(control))
						for (nt in unique(CV)){
							if(sum(CV==nt) > 1){ #Need at least two data points!
									fit[CV==nt] <- loess(
											A[CV==nt,i] ~ M[CV==nt,i], 
											weights=sqrt(control.rf[CV==nt,i]/control.sum[CV==nt] + test.rf[CV==nt,i]/test.sum[CV==nt]) , 
											...
									)$fitted
							}
						}
						return(fit)
					})
			
			test.norm <- 10^(M+A-fits)
			test.norm[,1:5] <- test.norm[,1:5]*rowSums(test[,1:5])
			test.norm[,5+1:5] <- test.norm[,5+1:5]*rowSums(test[,5+1:5])			
			control.norm <- 10^(M-A+fits)
			control.norm[,1:5] <- control.norm[,1:5]*rowSums(control[,1:5])
			control.norm[,5+1:5] <- control.norm[,5+1:5]*rowSums(control[,5+1:5])	
			if(round)
			{
				test.norm <- round(test.norm)
				control.norm <- round(control.norm)
			}
			list(test=test.norm, control=control.norm)
		})
		
#' Normalize nucleotide counts.
#' @rdname normalize-methods
#' @aliases normalize,deepSNV,missing-method
#' @usage \S4method{normalize}{deepSNV,missing}(test, control,  ...)
setMethod("normalize",
		signature = signature(test="deepSNV", control="missing"),
		function(test,control,...){
			fit = normalize(test@test,test@control,...)
			obj = initialize(test)
			obj@test = fit$test
			obj@control = fit$control
			deepSNV(obj)
		}
)


#' Get control counts
#' 
#' Convenience function to obtain the control counts from a deepSNV object.
#' @param deepSNV a \code{\link{deepSNV-class}} object
#' @param total Logical. If true the sum of both strands is returned
#' @param ... Additional param passed to specific methods
#' @return A matrix with the absolute frequencies summed over both strands.
#' @examples data(HIVmix)
#' control(HIVmix)[1:10,]
#' control(HIVmix, total=TRUE)[1:10,]
#' @exportMethod control
#' @rdname control-methods
#' @docType methods
setGeneric("control",
		function(deepSNV, ...) {
			standardGeneric("control")
		})

#' @rdname control-methods
#' @aliases control,deepSNV-method
#' @usage \S4method{control}{deepSNV}(deepSNV, total = FALSE)
setMethod("control",
		signature = signature(deepSNV="deepSNV"),
		function(deepSNV, total = FALSE) {
			count <- slot(deepSNV, "control")
			if(total)
				matrix(count[,1:5] + count[,6:10], ncol=5, dimnames=list(NULL, deepSNV@nucleotides[1:5]))
			else
				count
		}
)

#' Get test counts
#' 
#' Convenience function to obtain the test counts from a deepSNV object.
#' @param deepSNV a \code{\link{deepSNV-class}} object
#' @param total Logical. If true the sum of both strands is returned
#' @param ... Additional param passed to specific methods
#' @return A matrix with the absolute frequencies summed over both strands.
#' @examples data(HIVmix)
#' test(HIVmix)[1:10,]
#' test(HIVmix, total=TRUE)[1:10,]
#' @exportMethod test
#' @rdname test-methods
#' @docType methods
setGeneric("test",
		function(deepSNV, ...) {
			standardGeneric("test")
		})

#' @rdname test-methods
#' @aliases test,deepSNV-method
#' @usage \S4method{test}{deepSNV}(deepSNV, total = FALSE)
setMethod("test",
		signature = signature(deepSNV="deepSNV"),
		function(deepSNV, total = FALSE) {
			count <- slot(deepSNV, "test")
			if(total)
				matrix(count[,1:5] + count[,6:10], ncol=5, dimnames=list(NULL, deepSNV@nucleotides[1:5]))
			else
				count
		}
)

#' Get p-values
#' 
#' Convenience function to get the p-values from a deepSNV object.
#' @param deepSNV a \code{\link{deepSNV-class}} object
#' @param ... Additional param passed to specific methods
#' @return A matrix with the p-values.
#' @examples data(HIVmix)
#' p.val(HIVmix)[1:10,]
#' @exportMethod p.val
#' @rdname p.val-methods
#' @docType methods
setGeneric("p.val",
		function(deepSNV, ...) {
			standardGeneric("p.val")
		})

#' @rdname p.val-methods
#' @aliases p.val,deepSNV-method
#' @usage \S4method{p.val}{deepSNV}(deepSNV)
setMethod("p.val",
		signature = signature(deepSNV="deepSNV"),
		function(deepSNV) {
			slot(deepSNV, "p.val")
		}
)

#' Get coordinates
#' 
#' Convenience function to get the coordinates from a deepSNV object.
#' @param deepSNV a \code{\link{deepSNV-class}} object
#' @param ... Additional param passed to specific methods
#' @return A \code{\link{data.frame}} with columns "chrom(osome)" and "pos(ition)".
#' @examples data(HIVmix)
#' coordinates(HIVmix)[1:10,]
#' @exportMethod coordinates
#' @rdname coordinates-methods
#' @docType methods
setGeneric("coordinates",
		function(deepSNV, ...) {
			standardGeneric("coordinates")
		})

#' @rdname coordinates-methods
#' @aliases coordinates,deepSNV-method
#' @usage \S4method{coordinates}{deepSNV}(deepSNV)
setMethod("coordinates",
		signature = signature(deepSNV="deepSNV"),
		function(deepSNV) {
			slot(deepSNV, "coordinates")
		}
)