\name{zilgm_sim}
\alias{zilgm_sim}

\title{
Generate simulation data for zilgm
}
\description{
Generate simulation data for zilgm
}
\usage{
zilgm_sum(A, n, p, zlvs, family = c("poisson", "negbin"),
          signal, theta = NULL, noise, is.symm = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
  adjacency matrix
  }
  \item{n}{
  The number of observation to generate data
  }
  \item{p}{
  The number of nodes to generate
  }
  \item{zlvs}{
  zero probability
  }
  \item{family}{
  family of distribution 
  }
  \item{signal}{
  signal for distributio (the mean)
  }
  \item{theta}{
  dispersion parameter for negative binomial distribution
  }
  \item{noise}{
  noise parameter
  }
  \item{is.symm}{
  for symmetric matrix
  }
}
\value{
  \item{X}{The simulated data}
}

\details{
For simulation data
}

\references{
  Network analysis for excess zeros
}
\author{
Park Beomjin, Choi hosik
}

\seealso{
\link{zilgm}
}
\examples{
require(ZILGM)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
