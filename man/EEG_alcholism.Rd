\name{EEG_alcohol}
\alias{EEG_alcohol}
\title{
EEG alcholism database
}
\description{
Multiple electrode time series EEG recordings of control and alcoholic subjects.
}
\usage{
EEG_alcohol
}
%- maybe also 'usage' for other objects documented here.
\format{
  A data frame with 122 rows and 16384 variables.
}

\details{
This data arises from a large study to examine EEG correlates of genetic predisposition to alcoholism. 
It contains measurements from 64 electrodes placed on the scalp sampled at 256 Hz (3.9-msec epoch) for 1 second.

\itemize{
  \item group. 77 alcoholic and 45 control subjects.
  \item pixel. This is 256 x 64 images containing 16384 pixels(256 x 64)
}
}
\references{
Henri Begleiter
Neurodynamics Laboratory, 
State University of New York Health Center
Brooklyn, New York 

\url{http://kdd.ics.uci.edu/databases/eeg/eeg.data.html}
}

\examples{
data(EEG_alcohol)
colnames(eeg_alchol)[16385] = 'result'
}

