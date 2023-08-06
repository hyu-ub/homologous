\name{homologous_test}
\alias{homologous_test}
\title{Test of homologous hypothesis}
\arguments{
\item{x}{a a data.frame with four columns: x, the measurement such as tumor volume; time, the time of measurement; subject, the subject id; group, the subject's treatment group.}
\item{alternative}{alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".}
}
\value{
a data frame of the test statistics
}
\description{
Perform the test of homologous hypothesis for comparing growth curves}

\author{
Han Yu
}
