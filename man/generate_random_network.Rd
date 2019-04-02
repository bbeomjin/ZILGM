\name{zilgm_sim}
\alias{zilgm_sim}

\title{
Generate random graph structure
}
\description{
Generate random graph structure
}
\usage{
generate_random_network(p, prob)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
  The number of nodes to generate graph structure.
  }
  \item{prob}{
  The probability that the edges of the random graph appear.
  }
}
\value{
A \emph{p} x \emph{p} adjacency matrix
}


\references{Choi, H., Gim, J., Won, S., Kim, Y. J., Kwon, S., and Park, C. (2017). Network analysis for count data with excess zeros. \emph{BMC genetics}, \bold{18}(1), 93.}

\author{
Park Beomjin,
}

\seealso{
\link{zilgm}
}
\examples{
require(ZILGM)
n = 100; p = 30; prob = 2 / p;
A = generate_random_network(p, prob)
simul_dat = zilgm_sim(A = A, n = n, p = p, zlvs = 0.1, 
					  family = "negbin", signal = 1.5, theta = 1, noise = 0.1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
