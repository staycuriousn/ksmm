\name{ksmm}
\alias{ksmm}
\docType{package}
\title{
  Kernel Support Matrix Machine
}
\description{
 \code{ksmm} is used to train a kernel support matrix machine. It can be used to carry out classification.  
}
\usage{
	ksmm(x = NULL, y, K = NULL, x_dim, cost, kernel = c("linear", "rbf"), sigma = 1,
                maxit = 1e+4, epsilon = 5e-2, init_alpha = 1e-6, nCores = 1, \dots)
}
\arguments{
 \item{x}{
a input data matrix
}
 \item{y}{
a response vector with one label for each row / component of X. It must be 1 or -1
}
 \item{K}{
 the calculated kernel matrix. If you don't input anything, it works automatically.
}
 \item{x_dim}{
 the size of input matrix X
}
 \item{cost}{
 cost of constraints violation
}
 \item{kernel}{
 type of the kernel used in training.  
\itemize{
	\item{linear}{ : linear kernel.}
	\item{rbf}{ : Gaussian RBF kernel.}
}
}
 \item{sigma}{
 kernel width for the Radial Basis kernel function "rbf"
}
 \item{maxit}{
a maximum number of iteration for SMO algorithm
}
 \item{epsilon}{
 a minimum value for optimization 
}
 \item{init_alpha}{
initialization of alpha
}
 \item{nCores}{
the number of cores to use for parallel computing.
}
 \item{\dots}{
 currently not used.
}


}
\details{
  This function is built by KSMM paper. Detailed theory is included in the KSMM paper.
}
\author{
 Kyuri Park (based on Matlab code by Yunfei Ye)
}
\references{
  Ye, Y. (2019). A nonlinear kernel support matrix machine. \emph{International Journal of Machine Learning and Cybernetics}.
}
\keyword{ package }
\seealso{
  \code{make_ksmm_kernel}
}
\examples{
require(ksmm)
data(nottingham)

X = as.matrix(nottingham[,-1])
y = ifelse(nottingham[,1] == 1, 1, -1)
train_x = X[1,]
train_y = y[1]

ksmm_fit = ksmm(train_x, train_y, c(200,200), 1, 1, 'rbf', 100, 1, 5e-2)
    	
}
