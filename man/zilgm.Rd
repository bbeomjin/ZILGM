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
  Matrix
}
  \item{lambda}{
  tuning parameter
  }
  \item{nlambda}{
  The number of lambda
  }
  \item{family}{
  family
  }
  \item{update_type}{
  Algorithm for estimating coefficients
  }
  \item{sym}{
  estimated matrix 
  }
  \item{thresh}{
  thresholding
  }
  \item{weights_mat}{
  Matrix containing weights for observations in each node
  }
  \item{penalty_mat}{
  Matrix containing penalty for each node
  }
  \item{do_boot}{
  logical. whether use bootstrap
  }
  \item{boot_num}{
  The number of bootstrap.
  }
  \item{beta}{
  stability selection
  }
  \item{lambda_min_ratio}{
  lambda min ratio
  }
  \item{init_select}{
  whether use initial selection
  }
  \item{nCores}{
  The number of using cores
  }
  \item{...}{
  pass argument for coord zilgm
  }
}
\details{
Zero-inflated local graphical model
}
\value{
  An S3 object with the following slots
  \item{network}{Estimated network.}
  \item{coef_network}{estimated coefficients matrix.}
  \item{lambda}{lambda}
  \item{v}{optimal network index}
  \item{opt_lambda}{optimal lambda}
}
\references{Network analysis for excess zero}
\author{
Park Beomjin, Choi hosik
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
