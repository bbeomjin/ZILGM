\name{find_lammax}
\alias{find_lammax}
\title{
Compute the max value of \code{lambda}
}
\description{
Compute the max value of \code{lambda} from Karush-Kuhn-Tucker condition
}
\usage{find_lammax(X)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  A \emph{n} x \emph{p} data matrix.
}
}

\value{
  a numeric value
}
\references{Choi, H., J. Gim, S. Won, Y. J. Kim, S. Kwon, and C. Park, 2017: Network analysis for count data with excess zeros. \emph{BMC genetics}, \bold{18}, no. 1, 1-10. \cr
			Park, B., H. Choi, C. Park, 2021: Negative binomial graphical model with excess zeros.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
require(ZILGM)
set.seed(1)
n = 100; p = 10; prob = 2 / p;
A = generate_network(p, prob, type = "random")
simul_dat = zilgm_sim(A = A, n = n, p = p, zlvs = 0.1, 
					  family = "negbin", signal = 1.5, theta = 0.5, noise = 0.0)
find_lammax(simul_dat$X)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
