\name{zilgm}
\alias{zilgm}
\title{
Zero-inflated Local Graphical Model
}
\description{
Zero-inflated Local Graphical Model
}
\usage{
zilgm(X, lambda = NULL, nlambda = 50, family = c("Poisson", "NBI", "NBII"), update_type = c("IRLS", "MM"), 
                sym = c("AND", "OR"), thresh = 1e-6, weights_mat = NULL, penalty_mat = NULL,
                do_boot = FALSE, boot_num = 10, beta = 0.05, lambda_min_ratio = 1e-4,
                init_select = FALSE, nCores, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  A $n \times p$ data matrix, where $n$ is the number of observations and $p$ is the number of variables or nodes.
  }
  \item{lambda}{
  A sequence of regularization parameter to control a level of $l_1$-penalty.
  }
  \item{nlambda}{
  The number of regularization parameter.
  }
  \item{family}{
  Types of node-conditional distribution to be assumed among zero-inflated distributions.
  }
  \item{update_type}{
  Algorithm for estimating edge coefficients.
  }
  \item{sym}{
  Symmetrize the output graphs. If sym = "AND", the edge between node $i$ and node $j$ is selected only when both node $i$ and node $j$ are selected as neighbors for each other. If sym = "OR", the edge is selected when either node $i$ or node $j$ is selected as the neighbor for each other.
  }
  \item{thresh}{
  Threshold value for the estimated edge coefficients.
  }
  \item{weights_mat}{
  A $n \times p$ matrix containing weights for observations in each node
  }
  \item{penalty_mat}{
  A $p \times p$ matrix containing weights for each edge coefficients
  }
  \item{do_boot}{
  A logical. Whether to use Stability Apprach to Regularization Selection (StARS).
  }
  \item{boot_num}{
  The number of iteration for StARS.
  }
  \item{beta}{
  Threshold value on sparsity of the network.
  }
  \item{lambda_min_ratio}{
  The smallest value for lambda, as a fraction of the $\lambda^{\text{max}}$ of the regularization parameter.
  }
  \item{init_select}{
  A logical. Whether to use initial selection.
  }
  \item{nCores}{
  The number of cores to use for parallel computing.
  }
  \item{...}{
  Other arguments that can be passed to neighborhood selection function(zilgm_poisson, zilgm_negbin, and zilgm_negbin2)
  }
}
\details{
Zero-inflated local graphical model
}
\value{
  An S3 object with the following slots
  \item{network}{a list of $p \times p$ matrices of estimated networks along the regularization path.}
  \item{coef_network}{a array of $p \times p \times \text{nlambda}$ of estimated edge coefficients matrix along the regularization path.}
  \item{lambda}{vector used for regularization path.}
  \item{v}{vector of network variability measured for each regularization level.}
  \item{opt_lambda}{The lambda that gives the optimal network.}
}
\references{Network analysis for excess zero}
\author{
Park Beomjin, 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
Not yet
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
