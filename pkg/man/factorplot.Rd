\encoding{UTF-8}
\name{factorplot}
\alias{factorplot}
\title{Exact Factor Contrast Differences}
\description{
This function calculates all pairwise difference in the coefficients for factor covariates in linear models and GLMs.  
}
\usage{
factorplot(obj, factor.variable=NULL, pval = 0.05, 
	two.sided=TRUE, order="natural")
}
\arguments{
  \item{obj}{An object of class lm or glm}
  \item{factor.variable}{String containing the name of the factor to be plotted}
  \item{pval}{The (uncorrected) Type I error probability required, defaults to 0.05}
  \item{two.sided}{Logical argument indicating whether the hypothesis test should be against a two-sided alternative if TRUE (default) or a one-sided alternative if FALSE}
  \item{order}{One of \sQuote{natural}, \sQuote{alph}, or \sQuote{size} indicating how the levels of the factor should be ordered for presentation.  The \sQuote{natural} option (the default) leaves the levels as they are in the factor contrasts.  \sQuote{alph} sorts the levels alphabetically and \sQuote{size} sorts the levels by size of coefficient.}
}
\details{
This function calculates all pairwise differences in the coefficients for factor covariates in linear models and GLMs.  Specifically, it calculates:\cr\cr
 \deqn{\frac{b_{row}-b_{col}}{\sqrt{V(b_{row}-b_{col})}}}{(b[row]-b[col])/V(b[row]-b[col])^(1/2)}\cr where \deqn{V(b_{row}-b_{col}) = V(b_{row}) + V(b_{col}) + 2(1)(-1)V(b_{row},b_{col})}{V(b[row]-b[col]) = V(b[row]) + V(b[col]) + 2(1)(-1)V(b[row],b[col])}\cr\cr  This function is an alternative to using floating (quasi-) variances to present pairwise differences in coefficients (or floating absolute risk). \cr\cr

While the factorplot function and its print and summary methods work equally well regardless of the number of levels in the \code{factor.variable}, the plot function automatically scales the resulting graph to the appropriate size, but will be less useful as the number of contrasts gets large (e.g., > 30).  If more than one factor covariate is present and the \code{factor.variable} option is NULL, the function generates a text-based menu in the R GUI that will allow the users to pick the term for which they want to calculate the results.    
}
\value{
\item{b.diff}{An upper-triangular matrix of pairwise differences between row and column levels of the factor}
\item{b.sd}{An upper-triangular matrix of standard errors of the pairwise differences represented in b.diff}
\item{pval}{An upper-triangular matrix of uncorrected (one-sided) p-values corresponding to the entries of b.diff}
\item{bon.pval}{An upper-triangular matrix of bonferroni-corrected p-values corresponding to the entries of b.diff}
\item{p}{The p-value specified in the command}
}
\references{
Easton, D.F., J. Peto and G.A.G. Babiker. 1991. Floating absolute risk: kAn alternative to relative risk in survival and case control analysis avoiding an arbitrary reference group.  \emph{Statistics in Medicine} \bold{10}: 1025--1035.\cr
Firth, David and Renee X. de Menzes.  2004.  Quasi-variances.  \emph{Biometrika} \bold{91.1}: 65--80.\cr
Plummer, M. 2004. Improved estimates of floating absolute risk.  \emph{Statistics in Medicine} \bold{23}: 93--104.\cr}
\author{Dave Armstrong (Department of Political Science, UW-Milwaukee)}
\examples{
x <- as.factor(round(runif(1000, .5,20.5)))
levels(x) <- paste("lab", 1:20, sep="")
X <- model.matrix(~x)
Y <- X \%*\% rnorm(ncol(X),0,4) + rnorm(1000)
mod <- lm(Y ~ x)
fp <- factorplot(mod, "x", pval = 0.05, order="alph")
}
