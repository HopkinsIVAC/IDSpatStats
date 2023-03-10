% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatialfuncs.r, R/spatialfuncs_housholds.r
\name{get.pi.typed.bootstrap}
\alias{get.pi.typed.bootstrap}
\title{runs bootstrapping on \code{get.pi.typed}}
\usage{
get.pi.typed.bootstrap(
  posmat,
  typeA = -1,
  typeB = -1,
  r = 1,
  r.low = rep(0, length(r)),
  boot.iter
)

get.pi.typed.bootstrap(
  posmat,
  typeA = -1,
  typeB = -1,
  r = 1,
  r.low = rep(0, length(r)),
  boot.iter
)
}
\arguments{
\item{posmat}{a matrix with columns hh, type, x and y}

\item{typeA}{the "from" type that we are interested in, -1 is wildcard}

\item{typeB}{the "to" type that we are interested i, -1 is wildcard}

\item{r}{the series of spatial distances we are interested in}

\item{r.low}{the low end of each range....0  by default}

\item{boot.iter}{the number of bootstrap iterations}
}
\value{
pi values for all the distances we looked at

pi values for all the distances we looked at
}
\description{
Bootstraps typed pi values. Makes sure distances between a sample and
another draw of itself are left out

Bootstraps typed pi values. Makes sure distances between a sample and
another draw of itself are left out
}
\examples{
\dontrun{
 R/examples/get_pi_typed_bootstrap.R
 }
 
}
\seealso{
Other get.pi: 
\code{\link{get.pi.a.clustsurvey.hh.wts}()},
\code{\link{get.pi.a.clustsurvey.wts.bootstrap}()},
\code{\link{get.pi.a.clustsurvey.wts}()},
\code{\link{get.pi.a.typed.clustsurvey.hh.wts}()},
\code{\link{get.pi.a.typed.clustsurvey.wts}()},
\code{\link{get.pi.b.clustsurvey.hh.wts}()},
\code{\link{get.pi.b.clustsurvey.wts.bootstrap}()},
\code{\link{get.pi.b.clustsurvey.wts}()},
\code{\link{get.pi.b.typed.clustsurvey.hh.wts}()},
\code{\link{get.pi.b.typed.clustsurvey.wts}()},
\code{\link{get.pi.bootstrap}()},
\code{\link{get.pi.ci}()},
\code{\link{get.pi.clustsurvey.bootstrap}()},
\code{\link{get.pi.clustsurvey.ci}()},
\code{\link{get.pi.clustsurvey.hh.wts.window}()},
\code{\link{get.pi.clustsurvey.window}()},
\code{\link{get.pi.clustsurvey.wts.window}()},
\code{\link{get.pi.clustsurvey}()},
\code{\link{get.pi.hh.typed}()},
\code{\link{get.pi.permute}()},
\code{\link{get.pi.typed.gridcells}()},
\code{\link{get.pi.typed.permute}()},
\code{\link{get.pi.typed.survey}()},
\code{\link{get.pi.typed.wts.bootstrap}()},
\code{\link{get.pi.typed.wts.ci}()},
\code{\link{get.pi.typed.wts.permute}()},
\code{\link{get.pi.typed.wts}()},
\code{\link{get.pi.typed}()},
\code{\link{get.pi}()}

Other get.pi: 
\code{\link{get.pi.a.clustsurvey.hh.wts}()},
\code{\link{get.pi.a.clustsurvey.wts.bootstrap}()},
\code{\link{get.pi.a.clustsurvey.wts}()},
\code{\link{get.pi.a.typed.clustsurvey.hh.wts}()},
\code{\link{get.pi.a.typed.clustsurvey.wts}()},
\code{\link{get.pi.b.clustsurvey.hh.wts}()},
\code{\link{get.pi.b.clustsurvey.wts.bootstrap}()},
\code{\link{get.pi.b.clustsurvey.wts}()},
\code{\link{get.pi.b.typed.clustsurvey.hh.wts}()},
\code{\link{get.pi.b.typed.clustsurvey.wts}()},
\code{\link{get.pi.bootstrap}()},
\code{\link{get.pi.ci}()},
\code{\link{get.pi.clustsurvey.bootstrap}()},
\code{\link{get.pi.clustsurvey.ci}()},
\code{\link{get.pi.clustsurvey.hh.wts.window}()},
\code{\link{get.pi.clustsurvey.window}()},
\code{\link{get.pi.clustsurvey.wts.window}()},
\code{\link{get.pi.clustsurvey}()},
\code{\link{get.pi.hh.typed}()},
\code{\link{get.pi.permute}()},
\code{\link{get.pi.typed.gridcells}()},
\code{\link{get.pi.typed.permute}()},
\code{\link{get.pi.typed.survey}()},
\code{\link{get.pi.typed.wts.bootstrap}()},
\code{\link{get.pi.typed.wts.ci}()},
\code{\link{get.pi.typed.wts.permute}()},
\code{\link{get.pi.typed.wts}()},
\code{\link{get.pi.typed}()},
\code{\link{get.pi}()}
}
\concept{get.pi}
