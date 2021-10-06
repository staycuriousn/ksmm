\name{kfold_ksmm}
\alias{kfold_ksmm}
\docType{package}
\title{
  Parameter Tuning of Functions Using Grid Search
}
\description{
 This generic function tunes hyperparameters of statistical methods using a grid search over supplied parameter ranges for ksmm.
}
\usage{
	kfold_ksmm(x, y, x_dim, cost_range, sigma_range, kernel, nFold = 5, nCores = 1, \dots) 
}
\arguments{
\item{x}{
 a input data matrix
}
\item{y}{
 a response vector with one label for each row/component of X. It must be 1 or -1
}
\item{x_dim}{
 the size of input matrix X
}
\item{cost_range}{
 input data matrix to test
}
\item{sigma_range}{
 input data matrix to test
}
\item{kernel}{
 the kernel used in training and predicting
}
\item{nFold}{
 the number of partitions for cross-valication
}
\item{nCores}{
the number of cores to use for parallel computing
}

}
\details{
 This function is used to tune the parameters of \code{ksmm}. Detailed theory is included in the KSMM paper.
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
cost_seq = 2^{-1:1}
sigma_seq = 10^{-1:1}
train_x = X[c(1,2,99,100),]
train_y = y[c(1,2,99,100)]

kfold_fit = kfold_ksmm(x= train_x, y = train_y, x_dim = c(200, 200), cost_range = cost_seq, 
sigma_range = sigma_seq, kernel='rbf', nCores = 1, maxit = 1e+3, epsilon = 5e-2, nFold = 1)
opt_cost = kfold_fit$opt_params[1]
opt_sigma = kfold_fit$opt_params[2]
    	
}
