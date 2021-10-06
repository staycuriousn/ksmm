\name{predict.ksmm}
\alias{predict.ksmm}
\docType{package}
\title{
 Predict Methods for Kernel Support Matrix Machines
}
\description{
 This function predicts values based upon a model trained by \code{ksmm}
}
\usage{
	predict.ksmm(object, new_x, nCores = 1)
}
\arguments{
 \item{object}{
 Object of class 'ksmm', created by \code{ksmm}
}
\item{new_x}{
 input data matrix to test
}
\item{nCores}{
the number of cores to use for parallel computing
}

}
\details{
  This code is built by KSMM paper. Detailed theory is included in the KSMM paper
}
\author{
 Kyuri Park
}
\references{
  Ye, Y. (2019). A nonlinear kernel support matrix machine. \emph{International Journal of Machine Learning and Cybernetics}.
}
\keyword{ package }
\seealso{
  \code{ksmm}
}
\examples{


require(ksmm)
data(nottingham)

X = as.matrix(nottingham[,-1])
y = ifelse(nottingham[,1] == 1, 1, -1)
train_x = X[c(1,100),]
train_y = y[c(1,100)]
test_x = X[2,]
test_y = y[2]

ksmm_fit = ksmm(train_x, train_y, c(200,200), 1, 1, 'rbf', 100, 1, 5e-2)
ksmm_pred_y = predict.ksmm(object = ksmm_fit, new_x = test_x, nCores = 1)
table(ksmm_pred_y, test_y)
}
