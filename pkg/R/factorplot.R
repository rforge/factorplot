factorplot <- function(obj,...){
	UseMethod("factorplot")
}

factorplot.glm <-function(obj, factor.variable=NULL, order="natural", adjust.method="none", mcp.args = list(), ...){
	require(multcomp)
	tmp.classes <- attr(terms(obj), "dataClasses")
	tmp.classes <- tmp.classes[tmp.classes == "factor"]
	tmp.levs <- NULL
	if(level %in% names(...)){pval <- 1-level}
	else{pval <- .05}
	for(i in 1:length(tmp.classes)){
	    tmp.levs <- c(tmp.levs, length(levels(obj$model[[names(tmp.classes)[i]]])))
	}
	tmp.class <- tmp.classes[tmp.levs > 2]
	if(is.null(factor.variable)){
	{if(length(tmp.classes) > 1){
	    myvar <- names(tmp.classes)[menu(names(tmp.classes))]
	}
	else{ 
	    myvar <- names(tmp.classes[1])
	}}
	}
	else{
	    myvar <- factor.variable
	}
	mcp.args[[myvar]] <- 'Tukey'
	glht.obj <- do.call(glht, list(model=obj, linfct = do.call(mcp, mcp.arg.list)))
	s.glht.obj <- summary(glht.obj, type=adjust.method, level=1-pval, ...) 
	levs <- mod$xlevels[[myvar]]
	if(!(order %in% c("alph", "natural", "size")))stop("Order must be one of 'size', 'natural', 'alph'")
	varind <- which(attr(terms(obj), "term.labels") == myvar)
	bcols <- which(attr(model.matrix(obj), "assign") == varind)
	bcoef <- c(0, obj$coef[bcols])
	names(bcoef) <- levs
	tmp.ord <- switch(order, 
		alph = order(names(bcoef)), 
		size = order(bcoef), 
		natural = 1:length(levs))

	b.diff <- b.sd <- b.p <- array(NA, dim=c(length(levs), length(levs)))
	colnames(b.diff) <- rownames(b.diff) <- colnames(b.sd) <- rownames(b.sd) <- colnames(b.p) <- rownames(b.p) <- levs[tmp.ord]
	proc.names <- do.call(rbind, strsplit(names(s.glht.obj$test$coef), split=" - "))
	rc <- apply(proc.names, c(1,2), function(x)which(colnames(b.diff) == x))[,c(2,1)]
	b.diff[rc] <- s.glht.obj$test$coef
	b.sd[rc] <- s.glht.obj$test$sigma
	b.p[rc] <- s.glht.obj$test$pvalues
	b.diff <- b.diff[-nrow(b.diff),-1]
	b.sd <- b.sd[-nrow(b.sd), -1]
	b.p <- b.p[-nrow(b.p), -1]
	rns <- rownames(b.diff)
	cns <- colnames(b.diff)
	rns.split <- strsplit(rns, split=" ")
	rns.out <- sapply(rns.split, paste, collapse="_")
	cns.split <- strsplit(cns, split=" ")
	cns.out <- sapply(cns.split, paste, collapse="_")
	rownames(b.p) <- rownames(b.diff) <- rownames(b.sd) <- rns.out
	colnames(b.p) <- colnames(b.diff) <- colnames(b.sd) <- cns.out
	ret <- list(b.diff=b.diff, b.sd=b.sd, pval = b.p, p=pval)
	class(ret) <- c("factorplot", "list")
	return(ret)
}

factorplot.glht <-function(obj, adjust.method="none", ...){
	require(multcomp)
	if(level %in% names(...)){pval <- 1-level}
	else{pval <- .05}
	s.glht.obj <- summary(obj, type=adjust.method, ...) 
	proc.names <- do.call(rbind, strsplit(names(s.glht.obj$test$coef), split=" - "))
	levs <- unique(c(proc.names[,c(2,1)]))
	b.diff <- b.sd <- b.p <- array(NA, dim=c(length(levs), length(levs)))
	colnames(b.diff) <- rownames(b.diff) <- colnames(b.sd) <- rownames(b.sd) <- colnames(b.p) <- rownames(b.p) <- levs
	proc.names <- do.call(rbind, strsplit(names(s.glht.obj$test$coef), split=" - "))
	rc <- apply(proc.names, c(1,2), function(x)which(colnames(b.diff) == x))[,c(2,1)]
	b.diff[rc] <- s.glht.obj$test$coef
	b.sd[rc] <- s.glht.obj$test$sigma
	b.p[rc] <- s.glht.obj$test$pvalues
	b.diff <- b.diff[-nrow(b.diff),-1]
	b.sd <- b.sd[-nrow(b.sd), -1]
	b.p <- b.p[-nrow(b.p), -1]
	rns <- rownames(b.diff)
	cns <- colnames(b.diff)
	rns.split <- strsplit(rns, split=" ")
	rns.out <- sapply(rns.split, paste, collapse="_")
	cns.split <- strsplit(cns, split=" ")
	cns.out <- sapply(cns.split, paste, collapse="_")
	rownames(b.p) <- rownames(b.diff) <- rownames(b.sd) <- rns.out
	colnames(b.p) <- colnames(b.diff) <- colnames(b.sd) <- cns.out
	ret <- list(b.diff=b.diff, b.sd=b.sd, pval = b.p, p=pval)
	class(ret) <- c("factorplot", "list")
	return(ret)
}

factorplot.default <-function(obj, var, resdf, pval=0.05, two.sided=TRUE, order="natural", adjust.method="none"){
	require(multcomp)
	b <- obj
	if(!is.matrix(var)){
		v <- diag(var)
	} else{
		v <- var
	}
	if(is.null(names(b))){
		names(b) <- as.character(1:length(b))
	}
	levs <- colnames(v) <- rownames(v) <- names(b)
	if(!(order %in% c("alph", "natural", "size")))stop("Order must be one of 'size', 'natural', 'alph'")
	tmp.ord <- switch(order, 
		alph = order(names(b)), 
		size = order(b), 
		natural = 1:length(levs))
	b <- b[tmp.ord]
	v <- v[tmp.ord,tmp.ord]
	cmbn <- t(combn(length(b), 2))
	diffs <- matrix(0, nrow=length(b), ncol=nrow(cmbn))
	diffs[cbind(cmbn[,1], 1:ncol(diffs))] <- -1
	diffs[cbind(cmbn[,2], 1:ncol(diffs))] <- 1

	b.diff <- b.sd <- matrix(NA, ncol=length(b), nrow=length(b))
	b.diff[cmbn] <- b %*% diffs
	b.sd[cmbn] <- sqrt(diag(t(diffs) %*% v %*% diffs))
	colnames(b.diff) <- rownames(b.diff) <- colnames(b.sd) <- rownames(b.sd) <- names(b)
	b.diff <- b.diff[-nrow(b.diff),-1]
	b.sd <- b.sd[-nrow(b.sd),-1]

	b.t <- b.diff/b.sd
	rns <- rownames(b.t)
	cns <- colnames(b.t)
	rns.split <- strsplit(rns, split=" ")
	rns.out <- sapply(rns.split, paste, collapse="_")
	cns.split <- strsplit(cns, split=" ")
	cns.out <- sapply(cns.split, paste, collapse="_")
	rownames(b.t) <- rownames(b.diff) <- rownames(b.sd) <- rns.out
	colnames(b.t) <- colnames(b.diff) <- colnames(b.sd) <- cns.out
	b.p <- (2^(as.numeric(two.sided)))*pt(abs(b.t), resdf,lower.tail=FALSE)
	b.bp <- array(p.adjust(b.p, method=adjust.method), dim=dim(b.p))
	ret <- list(b.diff=b.diff, b.sd=b.sd, pval = b.bp, p = pval)
	class(ret) <- c("factorplot", "list")
	return(ret)
}



factorplot.multinom <- function(obj, variable = NULL, pval = .05, two.sided=T, order="natural", adjust.method="none"){
	v <- vcov(obj)
	b <- c(t(coef(obj)))
	names(b) <- rownames(v)
	inds <- which(gsub("^[^\\:]+\\:\\s*", "", names(b))	 == variable)
	v <- v[inds,inds]
	b <- b[inds]
	b <- c(0,b)
	v <- rbind(0, cbind(0, v))
	levs <- names(b) <- obj$lev
	resdf <- nrow(obj$residuals) - length(c(coef(obj)))
	if(!(order %in% c("alph", "natural", "size")))stop("Order must be one of 'size', 'natural', 'alph'")
	tmp.ord <- switch(order, 
		alph = order(names(b)), 
		size = order(b), 
		natural = 1:length(levs))
	b <- b[tmp.ord]
	v <- v[tmp.ord,tmp.ord]
	cmbn <- t(combn(length(b), 2))
	diffs <- matrix(0, nrow=length(b), ncol=nrow(cmbn))
	diffs[cbind(cmbn[,1], 1:ncol(diffs))] <- -1
	diffs[cbind(cmbn[,2], 1:ncol(diffs))] <- 1

	b.diff <- b.sd <- matrix(NA, ncol=length(b), nrow=length(b))
	b.diff[cmbn] <- b %*% diffs
	b.sd[cmbn] <- sqrt(diag(t(diffs) %*% v %*% diffs))
	colnames(b.diff) <- rownames(b.diff) <- colnames(b.sd) <- rownames(b.sd) <- names(b)
	b.diff <- b.diff[-nrow(b.diff),-1]
	b.sd <- b.sd[-nrow(b.sd),-1]

	b.t <- b.diff/b.sd
	rns <- rownames(b.t)
	cns <- colnames(b.t)
	rns.split <- strsplit(rns, split=" ")
	rns.out <- sapply(rns.split, paste, collapse="_")
	cns.split <- strsplit(cns, split=" ")
	cns.out <- sapply(cns.split, paste, collapse="_")
	rownames(b.t) <- rownames(b.diff) <- rownames(b.sd) <- rns.out
	colnames(b.t) <- colnames(b.diff) <- colnames(b.sd) <- cns.out
	b.p <- (2^(as.numeric(two.sided)))*pt(abs(b.t), resdf,lower.tail=FALSE)
	b.bp <- array(p.adjust(b.p, method=adjust.method), dim=dim(b.p))
	ret <- list(b.diff=b.diff, b.sd=b.sd, pval = b.bp,  p = pval)
	class(ret) <- c("factorplot", "list")
	return(ret)
}



squares <- function(ll, width=1,col){ 
    poly.x <- c(ll[1], ll[1]+width, ll[1]+width, ll[1], ll[1])
    poly.y <- c(ll[2], ll[2], ll[2]+width, ll[2]+width, ll[2])
    polygon(poly.x, poly.y, col=col)
}

plot.factorplot <- function(x, ..., abbrev.char=10, polycol=NULL, textcol = NULL, trans=NULL, 
	print.sig.leg=TRUE, print.square.leg=TRUE, scale.text=1, space.text=1){
r.bdiff <- x$b.diff[rev(1:nrow(x$b.diff)), ]
r.bsd <- x$b.sd[rev(1:nrow(x$b.sd)), ]
use.pval <- x$adj.pval
cns.out <- abbreviate(colnames(x$b.diff), abbrev.char)
rns.out <- abbreviate(rownames(x$b.diff), abbrev.char)

ymarg <- max(strwidth(rns.out, units="inches"))
tmarg <- max(strwidth(cns.out, units="inches"))
par(mai=c(0,ymarg,tmarg,0), oma=c(0,0,1,0))
plot(c(1,nrow(x$b.diff)+1), 
    c(1, nrow(x$b.diff)+1), type="n", main="", xlab="", ylab="", axes=FALSE)
axis(3, at=seq(from=1.5, to=nrow(x$b.diff)+.5, by=1), labels=gsub("_", " ", cns.out, fixed=T), 
    tick=FALSE, lwd=0, line=-1, las=2)
axis(2, at=seq(from=1.5, to=nrow(x$b.diff)+.5, by=1), labels=rev(gsub("_", " ", rns.out, fixed=T)), 
    tick=FALSE, lwd=0, line=-1, las=1)
rseq <- rev(1:nrow(x$b.diff))

if(is.null(polycol)){
	colvec <- c("gray80", "white", "gray40")
} else{
	colvec <- polycol
}
if(is.null(textcol)){
	text.col <- c("black", "black", "white")
} else{
	text.col <- textcol
}
if(!is.null(trans)){
	r.bdiff <- do.call(trans, list(r.bdiff))
}
m <- 1
for(i in rseq){ 
    for(j in m:nrow(x$b.diff)){
        if(use.pval[m,j] < x$p & x$b.diff[m,j] < 0){
            col.ind <- 1
            }
            else if(use.pval[m,j] < x$p & x$b.diff[m,j] > 0){
                col.ind <- 3
                }
                else{
                col.ind <- 2
                }
 	squares(c(j, i), col = colvec[col.ind])
    text(j+.5, i+.5+((.05*log(nrow(x$b.diff)))*space.text), sprintf("%.2f", r.bdiff[i,j]), font=2, 
        cex=(1-(.0275*(nrow(x$b.diff)-2)))*scale.text, col=text.col[col.ind])
    text(j+.5, i+.5-((.05*log(nrow(x$b.diff)))*space.text), sprintf("%.2f", r.bsd[i,j]), font=3, 
       cex=1-(.0275*(nrow(x$b.diff)-2)), col=text.col[col.ind])
    }
m <- m+1
}
leg <- legend(1,1, c("Significantly < 0", "Not Significant", "Significantly > 0"), fill=colvec, 
    bty="n", xjust=0, yjust=0, cex=ifelse(nrow(x$b.diff) == 2, .75, 1), plot=print.sig.leg)
legend(1+leg$rect$w*as.numeric(print.sig.leg), 1, c(expression(bold("bold = ")~b[row]-b[col]), 
	expression(italic("ital = ")~SE(b[row]-b[col]))), xjust=0, yjust=0, bty="n",
	cex=ifelse(nrow(x$b.diff) == 2, .75, 1), plot=print.square.leg)
}




print.factorplot <- function(x, ..., digits=3, sig=FALSE, trans=NULL){
	eg <- expand.grid(rownames(x$b.diff), colnames(x$b.diff))
	mc <- apply(eg, 2, function(x)max(nchar(x)))
	strnames <- paste(sprintf("%*s", mc[1], eg[,1]), sprintf("%*s", mc[1], eg[,2]), 
		sep = " - ")
	if(!is.null(trans)){
		x$b.diff <- do.call(trans, list(x$b.diff))
	}
	tmp <- cbind(c(x$b.diff), c(x$b.sd), c(x$pval))
	tmp <- round(tmp, 3)
	rownames(tmp) <- strnames
	colnames(tmp) <- c("Difference", "SE", "p.val")
	tmp <- as.data.frame(tmp)
	tmp <- na.omit(tmp)
	if(sig == TRUE){
		if(any(tmp$p.val > x$p)){
		tmp <- tmp[-which(tmp$p.val > x$p), ]
		}
	}
	t1c <- do.call(rbind, strsplit(as.character(tmp[,1]), split=".", fixed=T))
	t1mc <- apply(t1c, 2, function(x)max(nchar(x)))
	t2c <- do.call(rbind, strsplit(as.character(tmp[,2]), split=".", fixed=T))
	t2mc <- apply(t2c, 2, function(x)max(nchar(x)))

	tmp[,1] <- sprintf(paste("%*.", digits, "f", sep=""),t1mc[1], tmp[,1])
	tmp[,2] <- sprintf(paste("%*.", digits, "f", sep=""),t2mc[1], tmp[,2])
	tmp[,3] <- sprintf(paste("%1.", digits, "f", sep=""),tmp[,3])
	tmp
}

summary.factorplot <- function(object, ...){
	tmp <- object$b.diff
	tmp <- cbind(NA, tmp)
	tmp <- rbind(tmp, NA)
	rownames(tmp)[nrow(tmp)] <- colnames(tmp)[ncol(tmp)]
	colnames(tmp)[1] <- rownames(tmp)[1]
	tmp[lower.tri(tmp)] <- -t(tmp)[lower.tri(t(tmp))]
	tmp.p <- object$pval
	tmp.p <- cbind(NA, tmp.p)
	tmp.p <- rbind(tmp.p, NA)
	tmp.p[lower.tri(tmp.p)] <- t(tmp.p)[lower.tri(t(tmp.p))]
	rownames(tmp.p)[nrow(tmp.p)] <- colnames(tmp.p)[ncol(tmp.p)]
	colnames(tmp.p)[1] <- rownames(tmp.p)[1]
	diag(tmp.p) <- 1
	tmp1 <- tmp
	tmp1[which(tmp.p > object$p, arr.ind=T)] <- 0
	sig.plus <- apply(tmp1, 1, function(object)sum(object > 0))
	sig.minus <- apply(tmp1, 1, function(object)sum(object < 0))
	insig <- (nrow(tmp) -1) - (sig.plus + sig.minus)
    out <- cbind(sig.plus, sig.minus, insig)
    colnames(out) <- c("sig+", "sig-", "insig")
    out
}

