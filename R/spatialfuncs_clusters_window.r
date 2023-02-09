

##' Optimized version of \code{get.pi} for typed data.
##'
##' Version of the \code{get.pi} function that is optimized for statically typed data. That is
##' data where we are interested in the probability of points within some distance of points of
##' typeA are of typeB.
##'
##' @param posmat a matrix with columns type, x, y, and weight
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##'
##' @return pi values for all the distances we looked at
##'
##' @author Shaun Truelove
##'
##' @family get.pi
##' 
##' @export
##'
##'
get.pi.clustsurvey.window <- function(posmat,
                             r=1,
                             r.low=rep(0,length(r))) {
  
  return(.C("get_pi_clustsurvey_window",
            as.double(posmat[,"p"]),
            as.double(posmat[,"x"]),
            as.double(posmat[,"y"]),
            as.double(posmat[,"urban"]),
            as.integer(posmat[,"s"]),
            as.integer(nrow(posmat)),
            as.double(r.low),
            as.double(r),
            as.integer(length(r)),
            as.integer(1:nrow(posmat)),
            rc=double(length(r))
  )$rc)
}




##' Optimizewd version of \code{get.tau} for typed data
##'
##' Version of th e \code{get.tau} function that is optimized for
##' statically typed data. That is dat where we want the relationship between
##' points of type A and points of type B
##'
##' @param posmat a matrix with columns type, x, y, and weight
##' @param typeA the "from" type that we are interested in, -1 is wildcard
##' @param typeB the "to" type that we are interested i, -1 is wildcard
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##'
##' @return tau values for all the distances we looked at
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.tau
##' 
##' @export
##'
##'
get.tau.clustsurvey.window <- function(posmat,
                                    r=1,
                                    r.low=rep(0,length(r))) {
  return(.C("get_tau_clustsurvey_window",
            as.double(posmat[,"p"]),
            as.double(posmat[,"x"]),
            as.double(posmat[,"y"]),
            as.double(posmat[,"urban"]),
            as.integer(posmat[,"s"]),
            as.integer(nrow(posmat)),
            as.double(r.low),
            as.double(r),
            as.integer(length(r)),
            as.integer(1:nrow(posmat)),
            rc=double(length(r))
  )$rc)
}




##' Optimized version of \code{get.pi} for weighted cluster survey data.
##'
##' Version of the \code{get.pi} function that is optimized for cluster survey data. That is
##' data where individuals are grouped into clusters, as collected with the Demographics and Health Surveys.
##'
##' @param posmat a matrix with columns type, x, y, and weight
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##'
##' @return pi values for all the distances we looked at
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.pi
##' 
##' @export
##'
##'
get.pi.clustsurvey.wts.window <- function(posmat,
                                   r=1,
                                   r.low=rep(0,length(r))) {
    return(.C("get_pi_clustsurvey_wts_window",
              as.double(posmat[,"p"]),
              as.double(posmat[,"x"]),
              as.double(posmat[,"y"]),
              as.double(posmat[,"urban"]),
              as.integer(posmat[,"s"]),
              as.double(posmat[,"weight"]),
              as.double(posmat[,"alpha"]),
              as.double(posmat[,"delta"]),
              as.integer(nrow(posmat)),
              as.double(r.low),
              as.double(r),
              as.integer(length(r)),
              as.integer(1:nrow(posmat)),
              rc=double(length(r))
    )$rc)
}



##' Optimizewd version of \code{get.tau} for typed data
##'
##' Version of th e \code{get.tau} function that is optimized for
##' statically typed data. That is dat where we want the relationship between
##' points of type A and points of type B
##'
##' @param posmat a matrix with columns type, x, y, and weight
##' @param typeA the "from" type that we are interested in, -1 is wildcard
##' @param typeB the "to" type that we are interested i, -1 is wildcard
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##'
##' @return tau values for all the distances we looked at
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.tau
##'
##' @export
##'
get.tau.clustsurvey.wts.window <- function(posmat,
                                    r=1,
                                    r.low=rep(0,length(r))) {
    return(.C("get_tau_clustsurvey_wts_window",
              as.double(posmat[,"p"]),
              as.double(posmat[,"x"]),
              as.double(posmat[,"y"]),
              as.double(posmat[,"urban"]),
              as.integer(posmat[,"s"]),
              as.double(posmat[,"weight"]),
              as.double(posmat[,"alpha"]),
              as.double(posmat[,"delta"]),
              as.integer(nrow(posmat)),
              as.double(r.low),
              as.double(r),
              as.integer(length(r)),
              as.integer(1:nrow(posmat)),
              rc=double(length(r))
    )$rc)
}





##' Optimized version of \code{get.pi} for weighted cluster survey data.
##'
##' Version of the \code{get.pi} function that is optimized for cluster survey data. That is
##' data where individuals are grouped into clusters, as collected with the Demographics and Health Surveys.
##'
##' @param posmat a matrix with columns type, x, y, and weight
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##'
##' @return pi values for all the distances we looked at
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.pi
##'
##' @export
##'
get.pi.clustsurvey.hh.wts.window <- function(posmat,
                                      r=1,
                                      r.low=rep(0,length(r))) {
    return(.C("get_pi_clustsurvey_hh_wts_window",
              as.double(posmat[,"p"]),
              as.double(posmat[,"x"]),
              as.double(posmat[,"y"]),
              as.double(posmat[,"urban"]),
              as.integer(posmat[,"s"]),
              as.double(posmat[,"weight"]),
              as.double(posmat[,"alpha"]),
              as.double(posmat[,"delta"]),
              as.double(posmat[,"rho1"]),
              as.double(posmat[,"rho2"]),
              as.integer(nrow(posmat)),
              as.double(r.low),
              as.double(r),
              as.integer(length(r)),
              as.integer(1:nrow(posmat)),
              rc=double(length(r))
    )$rc)
}








##' Optimizewd version of \code{get.tau} for typed data
##'
##' Version of th e \code{get.tau} function that is optimized for
##' statically typed data. That is dat where we want the relationship between
##' points of type A and points of type B
##'
##' @param posmat a matrix with columns type, x, y, and weight
##' @param typeA the "from" type that we are interested in, -1 is wildcard
##' @param typeB the "to" type that we are interested i, -1 is wildcard
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##'
##' @return tau values for all the distances we looked at
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.tau
##' 
##' @export
##'
##'
get.tau.clustsurvey.hh.wts.window <- function(posmat,
                                       r=1,
                                       r.low=rep(0,length(r))) {
    return(.C("get_tau_clustsurvey_hh_wts_window",
              as.double(posmat[,"p"]),
              as.double(posmat[,"x"]),
              as.double(posmat[,"y"]),
              as.double(posmat[,"urban"]),
              as.integer(posmat[,"s"]),
              as.double(posmat[,"weight"]),
              as.double(posmat[,"alpha"]),
              as.double(posmat[,"delta"]),
              as.double(posmat[,"rho1"]),
              as.double(posmat[,"rho2"]),
              as.integer(nrow(posmat)),
              as.double(r.low),
              as.double(r),
              as.integer(length(r)),
              as.integer(1:nrow(posmat)),
              rc=double(length(r))
    )$rc)
}


