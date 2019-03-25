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
  The number of observation to generate simulated data.
  }
  \item{p}{
  The number of variables or nodes to generate simulated data.
  }
  \item{zlvs}{
  The probability of structural zero for each variable.
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
  \item{is.symm}{
  A logical. wheter to generate symmetric matrix.
  }
}
\value{
  \item{X}{$n \times p$ simulated data.}
}

\details{
For simulated data
}

\references{
  Network analysis for excess zeros
}
\author{
Park Beomjin,
}

\seealso{
\link{zilgm}
}
\examples{
require(ZILGM)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
