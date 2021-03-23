\name{classic3}
\alias{classic3}
\title{
Classic3 dataset
}
\description{
The classic3 dataset consists of CISI, MED, and CRAN document collections. CISI, MED, and CRAN contain 1460, 1032, and 1398 abstracts from information science, medical, and aerospace engineering journals, respectively.
}

\usage{
data(classic3)
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
