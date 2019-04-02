\name{zilgm_negbin2}
\alias{zilgm_negbin2}
\title{
Neiborhood selection under the zero-inflated negative binomial II distribution for Zero-inflated Local Graphical Model
}
\description{
Zero-inflated negative binomial regression with \emph{l_1}-regularization.
}
\usage{
zilgm_negbin2(y, x, lambda, weights = NULL, update_type = c("IRLS", "MM"), penalty.factor = NULL,
              thresh = 1e-6, EM_tol = 1e-6, EM_iter = 500, tol = 1e-6, maxit = 1e+3, fixed_sigma = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
  A response \emph{y}.
}
  \item{x}{
  A design matrix \emph{x}.
  }
  \item{lambda}{
  A regularization parameter to control a level of \emph{l_1}-penalty.
  }
  \item{weights}{
  Weights vector for observations
  }
  \item{update_type}{
  Algorithm for estimating coefficients
  }
  \item{penalty.factor}{
  Weights vector for coefficients of each variable.
  }
  \item{thresh}{
  Threshold value for the estimated coefficients.
  }
  \item{EM_tol}{
  Convergence tolerance for EM algorithm.
  }
  \item{EM_iter}{
  Maximum number of EM algorithm iterations.
  }
  \item{tol}{
  Convergence tolerance for coordinate descent.
  }
  \item{maxit}{
  Maximum number of coordinate descent iterations.
  }
  \item{fixed_sigma}{
  A logical. Whether to estimate dispersion parameter sigma.
  }
}
\details{
Zero-inflated local graphical model
}
\value{
  An S3 object with the following slots
  \item{beta}{Estimated coefficients vector.}
  \item{sigma}{Estimated dispersion parameter sigma.}
  \item{prob}{Estimated probability of structural zero.}
  \item{pos_zero}{Indices of zero values.}
  \item{iteration}{Iteration numbers until convergence.}
  \item{loglik}{\emph{l_1}-penalized negative log-likelihood value.}
  \item{call}{The mathced call.}
}
\references{Choi, H., Gim, J., Won, S., Kim, Y. J., Kwon, S., and Park, C. (2017). Network analysis for count data with excess zeros. \emph{BMC genetics}, \bold{18}(1), 93.}
\author{
Park Beomjin, 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# Zero-inflated negative binomial2 regression with l_1 penalty
require(ZILGM)
n = 100; p = 30; prob = 2 / p;
A = generate_random_network(p, prob)
simul_dat = zilgm_sim(A = A, n = n, p = p, zlvs = 0.1, 
					  family = "negbin", signal = 1.5, theta = 1, noise = 0.1)
y = simul_dat$X[, 1]
X = simul_dat$X[, -1]
negbin2_fit = zilgm_negbin2(y = y, x = X, lambda = 1, update_type = "IRLS", fixed_sigma = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
