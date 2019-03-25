\name{mitface}
\alias{mitface}
\title{
cbcl face database
}
\description{
This is a database of faces and non-faces, that has been used
extensively at the Center for Biological and Computational Learning at
MIT.
}
\usage{
mitface
}
%- maybe also 'usage' for other objects documented here.
\format{
  A data frame with 6977 rows and 362 variables.
}

\details{
The database is freely available for research use.  If you use this
database in published work, you must reference
\itemize{
  \item face. 1 is human face, 0 is not face.
  \item pixel. This is 19 x 19 images containing 361 pixels(19 x 19)
}
}
\references{
MIT Center For Biological and Computation Learning
\url{http://cbcl.mit.edu/software-datasets/FaceData.html}
}

\examples{
data(mitface)

# sample of face == 1
face_ind = which(mitface$face == 1)
plot(c(0, 1), c(0, 1), type = "n", xlab = "", ylab = "", axes = FALSE)
rasterImage(matrix(as.numeric(mitface[face_ind[1], -1]/255), ncol = 19), 0, 0, 1, 1)

# sample of face == 0
face_ind = which(mitface$face == 0)
plot(c(0, 1), c(0, 1), type = "n", xlab = "", ylab = "", axes = FALSE)
rasterImage(matrix(as.numeric(mitface[face_ind[1], -1]/255), ncol = 19), 0, 0, 1, 1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
