usethis::use_pipe(export =TRUE)

#' Convert pearson R to LOD score
#' 
#' @param n sample size 
#' @param r pearson r
#' @return LOD score
#' @export
RtoLOD=function(n, r) {(-n*log(1-r^2))/(2*log(10)) }

#' Convert LOD to p-value
#' 
#' Convert LOD to LRS, and then assume LRS is χ2 distributed with 1 d.f. 
#' @param x LOD score
#' @return p-value
#' @export
LODToPval =  function(x){  pchisq(x*(2*log(10)),df=1,lower.tail=FALSE)/2 }

#' Convert p-value to LOD score
#' 
#' Convert p-value to LOD score, inverse of LODToPval
#' @param x p-value
#' @return LOD score
#' @export
PvalToLOD =  function(x){ l= qchisq(x*2,df=1,lower.tail=FALSE)/(2*log(10)) 
                          l[is.na(l)]=0
                          return(l) }

#' Do loess regression of p1/(p1+p2) per block, in user defined blocks
#' 
#' @param p1 counts for parent 1 at each position
#' @param p2 counts for parent 2 at each position
#' @param pos positions for each marker 
#' @param bin.width collapse counts to bins every 500bp (default) or user specificied physical marker positions (optional)
#' @return imputed fitted allele frequency delta at each pos
#' @export
doBlockLoess=function(p1, p2, pos, bin.width=500) {
    loess.results=calcBlockLoess(p1,p2,pos, bin.width=bin.width, DEG=2)
    imputed.afd.delta=approxfun(loess.results$x0,loess.results$value)(pos)
    return(imputed.afd.delta) 
}

#SEE BRM https://github.com/huanglikun/BRM and
#https://academic.oup.com/bioinformatics/article/36/7/2150/5631910
# meta value of block
#https://academic.oup.com/bioinformatics/article/36/7/2150/5631910
def_block_pos <- function(chr, size){
	# block number
	n <- as.integer(2*chr/size)+2;
	if( n%%2 != 0 ) n <- n+1;
	n <- as.integer(n/2);
	# block index and the middle position of each block
	i <- c(1:n);
	pos <- (i-1)*size+floor(size/2); # middle position
	if(pos[n]>chr) pos[n]<-chr;
	return(pos);
}

cal_block_meta <- function(loc, val, chr, size, depth, MIN){
	# input: location vector, value vector
	# input: chr length, block size, location depth vector
	pos <- def_block_pos(chr, size);
	idx <- as.integer(0.5+loc/size)+1;
	#
	avg <- c();
	blockDepth <- c();
	for(i in 1:length(pos)){
		k  <- which(idx==i);
		no <- length(k);
		a <- NA;
		n <- 0;
		if (no > 0) {n <- sum(depth[k])};
		if( n >= MIN ) a <- sum(val[k])/n;
		avg <- c(avg, a);
	}
	return( list(pos=pos, avg=avg) );
}

calcBlockLoess=function(p1,p2,pos,bin.width=1000,min.depth=10,DEG=2,doSE=T) {
    opt.span <- function(model, criterion = c("aicc", "gcv"), span.range = c(0.05, 0.95)) {
	    as.crit <- function(x) {
	        span   <- x$pars$span;
	        traceL <- x$trace.hat;
	        sigma2 <- sum(x$residuals^2)/(x$n - 1);
	        aicc   <- log(sigma2) + 1 + 2 * (2 * (traceL + 1))/(x$n - traceL - 2);
	        gcv    <- x$n * sigma2/(x$n - traceL)^2;
	        result <- list(span = span, aicc = aicc, gcv = gcv);
	        return(result);
	    }
	    criterion <- match.arg(criterion);
	    fn <- function(span) {
	        mod <- update(model, span = span);
	        as.crit(mod)[[criterion]];
	    }
	    result <- optimize(fn, span.range);
	    return(list(span = result$minimum, criterion = result$objective));
	}
	# return

    #pmap binning 
    val=p1
    chr=max(pos,na.rm=T)+bin.width
    depth=p1+p2
    MIN=min.depth

    
    block <- cal_block_meta(pos, val,chr, bin.width, depth,MIN) 
    xx0 <- as.numeric(block$pos);
 #   print(xx0)
    yy <- as.numeric(block$avg);
    #print(x0)
    jdx <- which(!is.na(yy));
#    print(jdx)
 #   print(xx0[jdx])
    #print(jdx)
    fit0  <- loess(yy[jdx]~xx0[jdx], degree=DEG);

    span1 <- opt.span(fit0, criterion="aicc")$span;
    print(span1)
 #   print(span1)
   fit1  <- loess(yy[jdx]~xx0[jdx], degree=DEG, span=span1);
    plo <- predict(fit1, xx0, se=doSE);
    if(doSE){
        value <- plo$fit;
        return(list(x0=xx0, value=value, se=plo$se.fit))
    } else {
        value=plo
        return(list(x0=xx0, value=value))
        
    }
}

calcEffectiveNumberOfTests=function(Lg=4900,Lp=1.2e4, el=8.4, nchr=16){
    #Lg=4900
    #Lp=sum(sapply(gmaps[['A']], function(x) max(x$ppos)))/1e3
    #Lp=1.2e4
    rM=Lp/Lg # ~3.5kb per cM
    #el=8.4
    eff.tests=16+Lp/(rM*el)
    return(eff.tests)
}


#Calculate effect size given ratio of alleles
#K is selection threshold
#OR is observed P2/P1
# calculate x given OR and K
calcX.R = function (K, OR ) {
    T=qnorm(K, lower.tail = FALSE)
    fx=function(x,T,OR) { exp((pnorm(T, -x/2,sd=1, lower.tail=FALSE,log.p=T)
                              -pnorm(T,  x/2,sd=1, lower.tail=FALSE,log.p=T)))-OR }
    abs(uniroot ( 
        fx,
        c(-50,50), 
        T=T,
        OR=OR)$root)
}

# Leonid's clever approximation
calcX.LK = function(K,OR) {
    T=qnorm(K, lower.tail = FALSE)
    i=dnorm(T)/K
    (log(OR)/i)
}

#VE=x^2/4
#kvals = c(.1, .05, .01, .001, .0001, .00001, .000001)

getMeff_Li_and_Ji=function(cor.mat) {
    evals = eigen(cor.mat,symmetric=T)$values
    M = length(evals)
    L = M-1
    # Equation 5 from Li 
    intevals=ifelse(evals>=1, 1, 0)
    # modification for negative eigenvalues JB
    nonintevals=c(evals-floor(evals)) #[evals>0]
    Meff.li=sum(intevals+nonintevals)
    print(Meff.li)
    return(Meff.li)
}



#Lg=sum(sapply(gmaps[['A']], function(x) max(x$map)))
#eff.tests=calcEffectiveNumberOfTests()


#' Esimate sample sizes either given total pop (unselected) sample size and selection strengths
#' or effective depths if sample size is constrained by sequencing depth
#'
#' @param results.data data frame of results
#' @param sample.size size of unselected population
#' @param sel.high fraction of population in high tail
#' @param sel.low fraction of population in low tail
#' @param eff.length (default=600) effective number of tests across the genome
#' @return list of effective sample size estimates
#' @export 
getSampleSizes=function(results.data, sample.size, sel.high, sel.low, eff.length=600) {
    #get sample size based on 
    n1=sample.size*sel.high
    n2=sample.size*sel.low
    
    effective.sample.n.low=n2
    effective.sample.n.high=n1

    effective.sample.n.two.tailed=1/((n1+n2)/(n1*n2))
  
    print(paste("sample size based on n for low tail = ", round(effective.sample.n.low)))
    print(paste("sample size based on n for high tail = ", round(effective.sample.n.high)))
    print(paste("sample size based on n both tails = ", round(effective.sample.n.two.tailed)))

    #sample.n.tail=sample.size*sel.frac
    #print(paste('actual sample size= ', sample.n.tail, ' |  estimated sample size', sum(r+a)/eff.length ))
    depth.sample.size.low=sum(results.data$p1.low+results.data$p2.low)/eff.length 
    depth.sample.size.high=sum(results.data$p1.high+results.data$p2.high)/eff.length 
    depth.both=1/((depth.sample.size.low+depth.sample.size.high)/(depth.sample.size.low*depth.sample.size.high))
    
    print(paste("sample size based on depth for low tail = ", round(depth.sample.size.low)))
    print(paste("sample size based on depth for high tail = ", round(depth.sample.size.high)))
    print(paste("sample size based on depth for both tails = ", round(depth.both)))

    n.low=min(effective.sample.n.low,    depth.sample.size.low)
    n.high=min(effective.sample.n.high,  depth.sample.size.high)
    n.both=1/((n.low+n.high)/(n.low*n.high))
    return(list(n.low=n.low, n.high=n.high,n.both=n.both))
}


