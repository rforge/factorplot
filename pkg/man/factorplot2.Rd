\encoding{UTF-8}
\name{factorplot2}
\alias{factorplot2}
\title{Factor Contrast Differences Based on Published Estimates}
\description{
This function calculates all pairwise difference in the coefficients for factor covariates using published estimates (e.g., floating absolute risk) and floating variances.
}
\usage{
factorplot2(est, var, resdf, pval = 0.05, 
	two.sided=TRUE, order="natural")
}
\arguments{
  \item{est}{A named vector of estimates point estimates among differences are to be calculated and tested}
  \item{var}{A vector or matrix containing (perhaps floating) variances to be used in the hypothesis testing.  If a matrix is given, the matrix is used, otherwise the vector is made into a matrix with the \code{diag} command. }
  \item{resdf}{The residual degrees of freedom for the model.  For a normal distribution rather than t, use \code{Inf}}
  \item{pval}{The (uncorrected) Type I error probability required, defaults to 0.05}
  \item{two.sided}{Logical argument indicating whether the hypothesis test should be against a two-sided alternative if TRUE (default) or a one-sided alternative if FALSE}
  \item{order}{One of \sQuote{natural}, \sQuote{alph}, or \sQuote{size} indicating how the levels of the factor should be ordered for presentation.  The \sQuote{natural} option (the default) leaves the levels as they are in the factor contrasts.  \sQuote{alph} sorts the levels alphabetically and \sQuote{size} sorts the levels by size of coefficient.}
}
\details{
This function calculates all pairwise differences in the coefficients for factor covariates from provided estimates and (perhaps floating) variances.  Specifically, it calculates:\cr\cr
 \deqn{\frac{b_{row}-b_{col}}{\sqrt{V(b_{row}-b_{col})}}}{(b[row]-b[col])/V(b[row]-b[col])^(1/2)}\cr where \deqn{V(b_{row}-b_{col}) = V(b_{row}) + V(b_{col}) + 2(1)(-1)V(b_{row},b_{col})}{V(b[row]-b[col]) = V(b[row]) + V(b[col]) + 2(1)(-1)V(b[row],b[col])}\cr
If the variances passed to the function are organized in a matrix, then they are used as is.  Otherwise, they are made into a matrix using the \code{diag} function.    
}
\value{
\item{b.diff}{An upper-triangular matrix of pairwise differences between row and column levels of the factor}
\item{b.sd}{An upper-triangular matrix of standard errors of the pairwise differences represented in b.diff}
\item{pval}{An upper-triangular matrix of uncorrected (one-sided) p-values corresponding to the entries of b.diff}
\item{bon.pval}{An upper-triangular matrix of bonferroni-corrected p-values corresponding to the entries of b.diff}
\item{p}{The p-value specified in the command}
}
\references{
Easton, D.F., J. Peto and G.A.G. Babiker. 1991. Floating absolute risk: kAn alternative to relative risk in survival and case control analysis avoiding an arbitrary reference 
group.  \emph{Statistics in Medicine} \bold{10}: 1025--1035.\cr
Plummer, M. 2004. Improved estimates of floating absolute risk.  \emph{Statistics in Medicine} \bold{23}: 93--104.\cr
Plummer, M. et al. 2007  \emph{Helicobacter pylori} cytotoxin-associated genotype and gastric precancerous lesions.  \emph{Journal of the National Cancer Institute} \bold{99}, 1328--1334. \cr}
\author{Dave Armstrong (Department of Political Science, UW-Milwaukee)}
\examples{
### Example from Plummer et al. (2007)
est1 <- log(c(1.00,2.12,1.44,1.31,1.44,
     1.46,0.90))
var1 <- c(0.242,0.096,0.156,0.140,
     0.380,0.484,0.375)^2
names(est1) <- c(
     "Normal & superficial gastritis", 
     "Chronic gastritis", 
     "Chronic atrophic gastritits", 
     "Intestinal metaplasia I", 
     "Intestinal metaplasia II", 
     "Intestinal metaplasia III", 
     "Dysplasia")

plummer_fp1 <- factorplot2(est1, var1, Inf)
plot(plummer_fp1, bonferroni=FALSE, 
	abbrev.char = 100)}
