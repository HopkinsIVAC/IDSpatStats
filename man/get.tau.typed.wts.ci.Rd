% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatialfuncs_clusters.r
\name{get.tau.typed.wts.ci}
\alias{get.tau.typed.wts.ci}
\title{Bootstrap confidence interval for the \code{get.tau.typed.wts} values}
\usage{
get.tau.typed.wts.ci(
  posmat,
  typeA = -1,
  typeB = -1,
  r = 1,
  r.low = rep(0, length(r)),
  boot.iter = 1000,
  ci.low = 0.025,
  ci.high = 0.975
)
}
\arguments{
\item{posmat}{a matrix appropriate for input to \code{get.tau}}

\item{r}{the series of spatial distances wer are interested in}

\item{r.low}{the low end of each range....0  by default}

\item{boot.iter}{the number of bootstrap iterations}

\item{ci.low}{the low end of the ci...0.025 by default}

\item{ci.high}{the high end of the ci...0.975 by default}

\item{fun}{a function appropriate as input to \code{get.pi}}
}
\value{
tau values for all the distances examined
}
\description{
Wrapper to \code{get.tau.bootstrap} that takes care of calulating
the confidence intervals based on the bootstrapped values
}
\seealso{
Other get.tau: 
\code{\link{get.tau.bootstrap}()},
\code{\link{get.tau.ci}()},
\code{\link{get.tau.clustsurvey.bootstrap}()},
\code{\link{get.tau.clustsurvey.ci}()},
\code{\link{get.tau.clustsurvey.hh.wts.window}()},
\code{\link{get.tau.clustsurvey.hh.wts}()},
\code{\link{get.tau.clustsurvey.window}()},
\code{\link{get.tau.clustsurvey.wts.bootstrap}()},
\code{\link{get.tau.clustsurvey.wts.ci}()},
\code{\link{get.tau.clustsurvey.wts.window}()},
\code{\link{get.tau.clustsurvey.wts}()},
\code{\link{get.tau.clustsurvey}()},
\code{\link{get.tau.hh.typed.bootstrap}()},
\code{\link{get.tau.hh.typed}()},
\code{\link{get.tau.permute}()},
\code{\link{get.tau.typed.bootstrap.parallel}()},
\code{\link{get.tau.typed.bootstrap}()},
\code{\link{get.tau.typed.clustsurvey.hh.wts}()},
\code{\link{get.tau.typed.clustsurvey.wts}()},
\code{\link{get.tau.typed.permute}()},
\code{\link{get.tau.typed.survey}()},
\code{\link{get.tau.typed.wts.bootstrap.parallel}()},
\code{\link{get.tau.typed.wts.bootstrap}()},
\code{\link{get.tau.typed.wts.permute}()},
\code{\link{get.tau.typed.wts}()},
\code{\link{get.tau.typed}()},
\code{\link{get.tau}()}
}
\author{
Justin Lessler and Henrik Salje
}
\concept{get.tau}
