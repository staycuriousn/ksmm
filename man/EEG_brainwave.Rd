\name{brainwave}
\alias{brainwave}
\title{
EEG signal data from 10 college students while they watched MOOC video clips.
}
\description{

}
\usage{
brainwave
}
%- maybe also 'usage' for other objects documented here.
\format{
  A data frame with 100 rows and 1152 variables.
}

\details{

These data are collected from ten students, each watching ten videos. Therefore, it can be seen as only 100 data points for these 12000+ rows.
After each session, the student rated his/her confusion level on a scale of 1-7, 
These labels if further normalized into labels of whether the students are confused or not. 
Since the length of the EEG information collected for each image is different, the size of each image was unified with the length of the longest image.


\itemize{
  \item group. 50 confused group & 50 not confused group.
  \item pixel. This is 144 x 8 images containing 1152 pixels(144 x 8)
}
}
\references{
Wang, H., Li, Y., Hu, X., Yang, Y., Meng, Z., & Chang, K. M. (2013, June). Using EEG to Improve Massive Open Online Courses Feedback Interaction. In AIED Workshops. 

\url{https://www.kaggle.com/wanghaohan/confused-eeg}
}

\examples{
data(brainwave)
colnames(brainwave)[1] = 'result' 
}

