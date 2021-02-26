\name{zilgm_negbin}
\alias{zilgm_negbin}
\title{
Neiborhood selection under the zero-inflated negative binomial I distribution for Zero-inflated Local Graphical Model
}
\description{
Neiborhood selection under the zero-inflated negative binomial distribution I by zero-inflated negative binomial regression with \emph{l_1}-regularization.
}
\usage{
zilgm_negbin(y, x, lambda, weights = NULL, update_type = c("IRLS", "MM"), 
             penalty.factor = NULL, thresh = 1e-6, EM_tol = 1e-5, EM_iter = 3e+2,
 			 tol = 1e-6, maxit = 3e+2, fixed_theta = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
  A numeric vector. A response variable \emph{y}.
}
  \item{x}{
  A design matrix \emph{x}.
  }
  \item{lambda}{
  A regularization parameter to control a level of \emph{l_1}-penalty.
  }
  \item{weights}{
  Weights vector for observations. A default value for each observation is 1.
  }
  \item{update_type}{
  Types of algorithm for estimating coefficients.
  }
  \item{penalty.factor}{
  Weights vector for coefficients of each variable. A default value for each variable is 1. 
  }
  \item{thresh}{
  Threshold value for the estimated coefficients.
  }
  \item{EM_tol}{
  Convergence tolerance for EM algorithm.
  }
  \item{EM_iter}{
  A integer value. Maximum number of EM algorithm iterations.
  }
  \item{tol}{
  Convergence tolerance for coordinate descent algorithm.
  }
  \item{maxit}{
  A integer value. Maximum number of coordinate descent iterations.
  }
  \item{fixed_theta}{
  A logical. Whether to estimate dispersion parameter theta.
  }
}
\value{
  An S3 object with the following slots
  \item{beta}{Estimated coefficients vector.}
  \item{theta}{Estimated dispersion parameter theta.}
  \item{prob}{Estimated probability of structural zero.}
  \item{pos_zero}{Indices of zero values.}
  \item{iteration}{Iteration numbers until convergence.}
  \item{loglik}{\emph{l_1}-penalized negative log-likelihood value.}
  \item{call}{The matched call.}
}

\references{Wang, Z., S. Ma, M. Zappitelli, C. Parikh, C.-Y. Wang, and P. Devarajan, 2016: Penalized count data regression with application to hospital stay after pediatric cardiac surgery. \emph{Stat. Methods Med. Res.}, \bold{25}, no.6, 2685-2703. \cr
			Choi, H., J. Gim, S. Won, Y. J. Kim, S. Kwon, and C. Park, 2017: Network analysis for count data with excess zeros. \emph{BMC genetics}, \bold{18}, no. 1, 1-10. \cr
			Park, B., H. Choi, C. Park, 2021: Negative binomial graphical model with excess zeros.}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# Zero-inflated negative binomial I regression with l_1 penalty
require(ZILGM)
set.seed(1)
n = 100; p = 10; prob = 2 / p;
A = generate_network(p, prob, type = "random")
simul_dat = zilgm_sim(A = A, n = n, p = p, zlvs = 0.1, 
					  family = "negbin", signal = 1.5, theta = 0.5, noise = 0.0)
y = simul_dat$X[, 1]
X = simul_dat$X[, -1]
negbin_fit = zilgm_negbin(y = y, x = X, lambda = 1, update_type = "IRLS", fixed_theta = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
