\name{factorplot-package}
\alias{factorplot-package}
\alias{factorplot-package}
\docType{package}
\title{factorplot}
\description{Factorplot is a way to summarize and plot information from categorical predictors from linear models and GLMs.  It creates all simple contrasts and analytical standard errors for those contrasts.}
\details{
\tabular{ll}{
Package: \tab factorplot\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-07-20\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}
After a linear model or GLM has been estimated, the factorplot command creates all pairwise differences among the levels (including the reference category) of the indicated factor as well as their associated standard errors to facilitate hypothesis testing directly.  The print method prints the pairwise difference, standard error, p-value and Bonferroni-corrected p-value.  The summary method prints the number of significant positive/negative pairwise differences.  The plot method makes something akin to an upper-triangular levelplot that indicates whether differences are positive/negative and statistically significant. 
}
\author{
Dave Armstrong
Maintainer: Dave Armstrong <davearmstrong.ps@gmail.com>
}
\references{
~~ Literature or other references for background information ~~
}
\keyword{package}
\examples{
x <- as.factor(round(runif(1000, .5,20.5)))
levels(x) <- paste("lab", 1:20, sep="")
X <- model.matrix(~x)
b <- rnorm(ncol(X),0,4)
Y.hat <- X \%*\% b 
Y <- Y.hat  + rnorm(1000)
mod <- lm(Y ~ x)
fp <- factorplot(mod, "x", pval=0.05, order="alph")
print(fp)
summary(fp)
plot(fp)
}
