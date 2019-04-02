\name{zilgm}
\alias{zilgm}
\title{
Undirected Local Graphical Model for Count data with excess zeros and overdispersion
}
\description{
Network estimation following local markov property for count data with zero inflated and overdispersion by specified node-conditional distribution
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
  A \emph{n} x \emph{p} data matrix, where \emph{n} is the number of observations and \emph{p} is the number of variables or nodes.
  }
  \item{lambda}{
  A sequence of regularization parameter to control a level of \emph{l_1}-penalty.
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
  Symmetrize the output graphs. If \code{sym = "AND"}, the edge between node \emph{i} and node \emph{j} is selected only when both node \emph{i} and node \emph{j} are selected as neighbors for each other. If \code{sym = "OR"}, 
  the edge is selected when either node \emph{i} or node \emph{j} is selected as the neighbor for each other.
  }
  \item{thresh}{
  Threshold value for the estimated edge coefficients.
  }
  \item{weights_mat}{
  A \emph{n} x \emph{p} matrix containing weights for observations in each node
  }
  \item{penalty_mat}{
  A \emph{p} x \emph{p} matrix containing weights for each edge coefficients
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
  The smallest value for lambda, as a fraction of the \emph{lambda^{max}} of the regularization parameter.
  }
  \item{init_select}{
  A logical. Whether to use initial selection.
  }
  \item{nCores}{
  The number of cores to use for parallel computing.
  }
  \item{...}{
  Other arguments that can be passed to neighborhood selection function(\code{zilgm_poisson, zilgm_negbin, and zilgm_negbin2})
  }
}
\details{
Zero-inflated local graphical model
}
\value{
  An S3 object with the following slots
  \item{network}{a list of \emph{p} x \emph{p} matrices of estimated networks along the regularization path.}
  \item{coef_network}{a array of \emph{p} x \emph{p} x nlambda of estimated edge coefficients matrix along the regularization path.}
  \item{lambda}{vector used for regularization path.}
  \item{v}{vector of network variability measured for each regularization level.}
  \item{opt_lambda}{The lambda that gives the optimal network.}
}
\references{Choi, H., Gim, J., Won, S., Kim, Y. J., Kwon, S., and Park, C. (2017). Network analysis for count data with excess zeros. \emph{BMC genetics}, \bold{18}(1), 93.}
\author{
Park Beomjin, 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
require(ZILGM)
n = 100; p = 30; prob = 2 / p;
A = generate_random_network(p, prob)
simul_dat = zilgm_sim(A = A, n = n, p = p, zlvs = 0.1, 
					  family = "negbin", signal = 1.5, theta = 1, noise = 0.1)
lambda_max = find_lammax(simul_dat$X)
lambda_min = 1e-4 * lambda_max
lambs = exp(seq(log(lambda_max), log(lambda_min), length.out = 50))
simul_nb2_result = zilgm(X = simul_dat$X, lambda = lambs, family = "NBII", thresh = 1e-6, EM_tol = 1e-6,
                        EM_iter = 1000, nCores = 30, update_type = "IRLS", weights_mat = NULL, do_boot = TRUE,
                        boot_num = 30, init_select = FALSE, beta = 0.05, sym = "OR", tol = 1e-6, maxit = 1e+3)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
