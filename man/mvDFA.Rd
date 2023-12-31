% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mvDFA.R
\name{mvDFA}
\alias{mvDFA}
\title{Analyze multivariate correlated time series and estimate long memory by the extension of the using univariate Detrended Fluctuations Analysis (DFA; Peng et al., 1995) to multivariate time series: mvDFA}
\usage{
mvDFA(
  X,
  steps = 50,
  degree = 1,
  verbose = FALSE,
  cores = 1,
  covlist = FALSE,
  brownian = FALSE
)
}
\arguments{
\item{X}{Matrix or data.frame containing the time series in long format.}

\item{steps}{Maximum number of window sizes. These are spread logarithmically. If time series is short and steps is large, fewer window sizes are drawn. Default to \code{50}. The dimensions (\code{ncol(X)}) and the \code{degree} influence the smallest possible window size.}

\item{degree}{The maximum order of the detrending polynomial in the segments. This influences the smallest window size \code{minS} such that \code{minS} = \code{d + degree + 2}, where \code{d} is the dimension of the time series.}

\item{verbose}{Indicator whether additional info should be printed. Default to \code{TRUE}.}

\item{cores}{Number of cores used in computation. Default to \code{1}.}

\item{covlist}{Indicator whether covariance of the time series per window size should be saved in a list.}

\item{brownian}{Indicator whether time series are assumed to be brownian (i.e. variance increases proportional to time)}
}
\value{
An object of class \code{mvDFA} containing long memory coefficients (Hurst exponents) and corresponding further informations:

\item{Ltot}{ the estimated long memory coefficient for the multivariate time series using the total variance approach}
\item{Lgen}{the generalized approach}
\item{Lfull}{the average covariance approach}
\item{LmeanUni}{average Hurst exponent across all time series}
\item{univariate_DFA}{univariate Hurst exponents}
\item{R2tot}{R-squared of total variance approach in regression of log10(RMS) vs log10(S)}
\item{R2gen}{R-squared of generalized variance approach in regression of log10(RMS) vs log10(S)}
\item{R2full}{R-squared of  covariance approach in regression of log10(RMS) vs log10(S)}
\item{R2meanUni}{average R-squared across all time series in regression of log10(RMS) vs log10(S)}
\item{R2univariate_DFA}{R-squares of single time series approach in regression of log10(RMS) vs log10(S)}
\item{RMS_tot}{a list of Root Mean Squares per window size corresponding to the total variance approach}
\item{RMS_gen}{a list of Root Mean Squares per window size corresponding to the total generalized approach}
\item{Cov_RMS_s}{a list of Root Mean Squares per window size corresponding to the covariance approach}
\item{S}{window sizes used}
\item{CovRMS_list}{a list of covariance matrices per \code{S} may be returned}
}
\description{
Analyze multivariate correlated time series and estimate long memory by the extension of the using univariate Detrended Fluctuations Analysis (DFA; Peng et al., 1995) to multivariate time series: mvDFA
}
\examples{
Sigma <- matrix(.5, 3, 3); diag(Sigma) <- 1
# generate correlated Gaussian white noise (i.i.d. multivariate normal variables)
X <- mvtnorm::rmvnorm(n = 500, sigma = Sigma)
mvDFA(X = X, steps = 5) # steps = 5 is only for demonstration,
                        # use many steps instead, e.g. steps = 50!
}
\references{
Peng, C. K., Havlin, S., Stanley, H. E., & Goldberger, A. L. (1995). Quantification of scaling exponents and crossover phenomena in nonstationary heartbeat time-series. Chaos, 5, 82–87. <doi:10.1063/1.166141>
}
