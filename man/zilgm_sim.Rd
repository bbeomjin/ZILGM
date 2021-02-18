\name{zilgm_sim}
\alias{zilgm_sim}

\title{
Generate simulated data under the graph structure for ZILGM models
}
\description{
Generate data from zero-inflated multivariate Poisson or negative binomial distribution under given graph structure.
}
\usage{
zilgm_sim(A, n, p, zlvs, family = c("poisson", "negbin"),
          signal, theta = NULL, noise)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
  A \emph{p} x \emph{p} adjacency matrix following graph structure.
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


\references{Choi, H., Gim, J., Won, S., Kim, Y. J., Kwon, S., and Park, C. (2017). Network analysis for count data with excess zeros. \emph{BMC genetics}, \bold{18}(1), 93. \\
			Park, B., Choi, H., Park, C. (2021). Negative binomial graphical model with excess zeros.}

\author{
Beomjin, Park.
}

\seealso{
\link{zilgm}
}
\examples{
require(ZILGM)
set.seed(1)
n = 100; p = 10; prob = 2 / p;
A = generate_network(p, prob, type = "random")
simul_dat = zilgm_sim(A = A, n = n, p = p, zlvs = 0.1, 
					  family = "negbin", signal = 1, theta = 0.5, noise = 0.0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
