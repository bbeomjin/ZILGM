\name{zilgm_sim}
\alias{zilgm_sim}

\title{
Generate graph structure
}
\description{
Generate graph structure.
}
\usage{
generate_network(p, prob, NofHub = 3, type = c("scale-free", "hub", "random"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
  The number of nodes to generate graph structure.
  }
  \item{prob}{
  The probability that the edges of the random graph appear.
  }
  \item{NofHub}{  
  The number of hub nodes in hub graph.
  }
  \item{type}{
  Types of graph structure. Either a character string representing one of the graph structures.
  }
}
\value{
A \emph{p} x \emph{p} adjacency matrix
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
n = 100; p = 30; prob = 2 / p;

# Random graph
random_graph = generate_network(p, prob, type = "random")

# Scale-free graph
scale_free_graph = generate_network(p, prob, type = "random")

# Hub graph
hub_graph = generate_network(p, prob, NofHub = 3, type = "hub")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
