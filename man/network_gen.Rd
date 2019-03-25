\name{network_gen}
\alias{network_gen}

\title{
Generate random graph structure
}
\description{
Generate random graph structure
}
\usage{
generate_network(node, prob)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{node}{
  The number of nodes to generate graph.
  }
  \item{prob}{
  The probability of appearing edges. 
  }
}
\value{
  \item{networkmat}{Adjacency matrix for random graph.}
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
