


##'
##' Optimized version of \code{get.pi} for typed data.
##'
##' Version of the \code{get.pi} function that is optimized for statically typed data. That is
##' data where we are interested in the probability of points within some distance of points of
##' typeA are of typeB.
##'
##' @param posmat a matrix with columns type, x and y
##' @param typeA the "from" type that we are interested in, -1 is wildcard
##' @param typeB the "to" type that we are interested i, -1 is wildcard
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
get.pi.typed.survey <- function(posmat,
                         typeA = -1,
                         typeB = -1,
                         r=1,
                         r.low=rep(0,length(r))) {
  
  return(.C("get_pi_typed_survey",
            as.integer(posmat[,"type"]),
            as.double(posmat[,"x"]),
            as.double(posmat[,"y"]),
            as.integer(nrow(posmat)),
            as.integer(typeA),
            as.integer(typeB),
            as.double(r.low),
            as.double(r),
            as.integer(length(r)),
            as.integer(1:nrow(posmat)),
            rc=double(length(r))
  )$rc)
}



##'
##' Optimized version of \code{get.theta} for typed data.
##'
##' Version of the \code{get.theta} function that is optimized for statically typed data. That is
##' data where we are interested in the odds that points within some distance of points of
##' typeA are of typeB.
##'
##' @param posmat a matrix with columns type, x and y
##' @param typeA the "from" type that we are interested in, -1 is wildcard
##' @param typeB the "to" type that we are interested i, -1 is wildcard
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##'
##' @return theta values for all the distances we looked at
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.theta
##' 
##' @export
##'
##'
get.theta.typed.survey <- function(posmat,
                            typeA = -1,
                            typeB = -1,
                            r=1,
                            r.low=rep(0,length(r))) {
  
  return(.C("get_theta_typed_survey",
            as.integer(posmat[,"type"]),
            as.double(posmat[,"x"]),
            as.double(posmat[,"y"]),
            as.integer(nrow(posmat)),
            as.integer(typeA),
            as.integer(typeB),
            as.double(r.low),
            as.double(r),
            as.integer(length(r)),
            as.integer(1:nrow(posmat)),
            rc=double(length(r))
  )$rc)
}





##'
##' Optimizewd version of \code{get.tau} for typed data
##'
##' Version of th e \code{get.tau} function that is optimized for
##' statically typed data. That is data where we want the relationship between
##' points of type A and points of type B
##'
##' @param posmat a matrix with columns type, x and y
##' @param typeA the "from" type that we are interested in, -1 is wildcard
##' @param typeB the "to" type that we are interested i, -1 is wildcard
##' @param r the series of spatial distances we are interested in
##' @param r.low the low end of each range....0  by default
##' @param comparison.type what type of points are included in the comparison set.
##' \itemize{
##'   \item "representative" if comparison set is representative of the underlying population
##'   \item "independent" if comparison set is cases/events coming from an indepedent process
##' }
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
get.tau.typed.survey <- function(posmat,
                          typeA = -1,
                          typeB = -1,
                          r=1,
                          r.low=rep(0,length(r)),
                          comparison.type = "representative") {
  
  if (comparison.type == "representative") {
    comp.type.int <- 0
  } else if (comparison.type == "independent") {
    comp.type.int <- 1
  } else {
    stop("unkown comparison type specified")
  }
  
  
  return(.C("get_tau_typed_survey",
            as.integer(posmat[,"type"]),
            as.double(posmat[,"x"]),
            as.double(posmat[,"y"]),
            as.integer(nrow(posmat)),
            as.integer(typeA),
            as.integer(typeB),
            as.double(r.low),
            as.double(r),
            as.integer(length(r)),
            as.integer(1:nrow(posmat)),
            as.integer(comp.type.int),
            rc=double(length(r))
  )$rc)
}

