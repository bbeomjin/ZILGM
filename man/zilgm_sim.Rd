\name{zilgm_sim}
\alias{zilgm_sim}

\title{
Generate simulated data for ZILGM models
}
\description{
Generate data from zero-inflated multivariate Poisson or negative binomial distribution with random graph structure.
}
\usage{
zilgm_sim(A, n, p, zlvs, family = c("poisson", "negbin"),
          signal, theta = NULL, noise)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prob}{
  The probability that the edges of the random graph appear.
  }
  \item{n}{
  The number of observation to generate simulated data.
  }
  \item{p}{
  The number of variables or nodes to generate simulated data.
  }
  \item{zlvs}{
  The probability of structural zero for each variable or nodes.
  }
  \item{family}{
  The distribution to draw samples.
  }
  \item{signal}{
  The location parameter for distribution. 
  }
  \item{theta}{
  The dispersion parameter for negative binomial distribution.
  }
  \item{noise}{
  The location parameter for noise distribution.
  }
}
\value{
A list of two elements:
  \item{A}{\emph{p} x \emph{p} adjacency matrix of the random graph.}
  \item{X}{\emph{n} x \emph{p} simulated data.}
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
simul_dat = zilgm_sim(prob = prob, n = n, p = p, zlvs = 0.1, 
					  family = "negbin", signal = 1.5, theta = 1, noise = 0.1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
