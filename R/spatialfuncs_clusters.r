
##' Optimized version of \code{get.pi} for typed data with sampling weights.
##'
##' Version of the \code{get.pi} function that is optimized for statically typed data 
##' and sampling weights. That is data where we are interested in the probability of 
##' points within some distance of points of typeA are of typeB.
##'
##' @param posmat a matrix with columns type, x, y, and sample weight
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
##' 
##' @export
##'
get.pi.typed.wts <- function(posmat,
                             typeA = -1,
                             typeB = -1,
                             r=1,
                             r.low=rep(0,length(r))) {
    
    return(.C("get_pi_typed_wts",
              as.integer(posmat[,"type"]),
              as.double(posmat[,"x"]),
              as.double(posmat[,"y"]),
              as.double(posmat[,"weight"]),
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
get.pi.clustsurvey <- function(posmat,
                               r=1,
                               r.low=rep(0,length(r)),
                               coord_type_geo = FALSE,
                               dist_unit = "K") {
    
    dist_unit <- toupper(substr(dist_unit, 1, 1))
    
    return(.C("get_pi_clustsurvey",
              as.double(posmat[,"p"]),
              as.double(posmat[,"x"]),
              as.double(posmat[,"y"]),
              as.integer(posmat[,"s"]),
              as.integer(nrow(posmat)),
              as.double(r.low),
              as.double(r),
              as.integer(length(r)),
              as.integer(1:nrow(posmat)),
              as.double(coord_type_geo),
              as.character(dist_unit),
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
##' 
##' @export
##'
get.pi.a.clustsurvey.wts <- function(posmat,
                                     r=1,
                                     r.low=rep(0,length(r)),
                                     remove_self=FALSE) {
    return(.C("get_pi_a_clustsurvey_wts",
              as.double(posmat[,"p"]),
              as.double(posmat[,"x"]),
              as.double(posmat[,"y"]),
              as.integer(posmat[,"s"]),
              as.double(posmat[,"weight"]),
              as.double(posmat[,"delta"]),
              as.double(posmat[,"alpha"]),
              as.integer(nrow(posmat)),
              as.double(r.low),
              as.double(r),
              as.integer(length(r)),
              as.integer(1:nrow(posmat)),
              rc=double(length(r)),
              as.integer(remove_self)
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
##' 
##' @export
##'
get.pi.b.clustsurvey.wts <- function(posmat,
                                     r=1,
                                     r.low=rep(0,length(r)),
                                     remove_self=FALSE) {
    return(.C("get_pi_b_clustsurvey_wts",
              as.double(posmat[,"p"]),
              as.double(posmat[,"x"]),
              as.double(posmat[,"y"]),
              as.integer(posmat[,"s"]),
              as.double(posmat[,"weight"]),
              as.double(posmat[,"delta"]),
              as.double(posmat[,"alpha"]),
              as.integer(nrow(posmat)),
              as.double(r.low),
              as.double(r),
              as.integer(length(r)),
              as.integer(1:nrow(posmat)),
              rc=double(length(r)),
              as.integer(remove_self)
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
##' 
##' @export
##'
get.pi.a.typed.clustsurvey.wts <- function(posmat,
                                           typeA = -1,
                                           typeB = -1,
                                           r=1,
                                           r.low=rep(0,length(r)),
                                           remove_self=FALSE) {
    return(.C("get_pi_a_typed_clustsurvey_wts",
              as.double(posmat[,"p"]),
              as.integer(posmat[,"type_a"]),
              as.integer(posmat[,"type_b"]),
              as.double(posmat[,"x"]),
              as.double(posmat[,"y"]),
              as.integer(posmat[,"s"]),
              as.double(posmat[,"weight"]),
              as.double(posmat[,"delta"]),
              as.double(posmat[,"alpha"]),
              as.integer(nrow(posmat)),
              as.integer(typeA),
              as.integer(typeB),
              as.double(r.low),
              as.double(r),
              as.integer(length(r)),
              as.integer(1:nrow(posmat)),
              rc=double(length(r)),
              as.integer(remove_self)
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
##' 
##' @export
##'
get.pi.b.typed.clustsurvey.wts <- function(posmat,
                                           typeA = -1,
                                           typeB = -1,
                                           r=1,
                                           r.low=rep(0,length(r)),
                                           remove_self=FALSE) {
    return(.C("get_pi_b_typed_clustsurvey_wts",
              as.double(posmat[,"p"]),
              as.integer(posmat[,"type_a"]),
              as.integer(posmat[,"type_b"]),
              as.double(posmat[,"x"]),
              as.double(posmat[,"y"]),
              as.integer(posmat[,"s"]),
              as.double(posmat[,"weight"]),
              as.double(posmat[,"delta"]),
              as.double(posmat[,"alpha"]),
              as.integer(nrow(posmat)),
              as.integer(typeA),
              as.integer(typeB),
              as.double(r.low),
              as.double(r),
              as.integer(length(r)),
              as.integer(1:nrow(posmat)),
              rc=double(length(r)),
              as.integer(remove_self)
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
##' 
##' @export
##'
get.pi.a.typed.clustsurvey.hh.wts <- function(posmat,
                                              typeA = -1,
                                              typeB = -1,
                                              r=1,
                                              r.low=rep(0,length(r)),
                                              remove_self=FALSE) {
    return(.C("get_pi_a_typed_clustsurvey_hh_wts",
              as.double(posmat[,"p"]),
              as.integer(posmat[,"type_a"]),
              as.integer(posmat[,"type_b"]),
              as.double(posmat[,"x"]),
              as.double(posmat[,"y"]),
              as.integer(posmat[,"s"]),
              as.double(posmat[,"weight"]),
              as.double(posmat[,"delta"]),
              as.double(posmat[,"alpha"]),
              as.double(posmat[,"rho1"]),
              as.double(posmat[,"rho2"]),
              as.integer(nrow(posmat)),
              as.integer(typeA),
              as.integer(typeB),
              as.double(r.low),
              as.double(r),
              as.integer(length(r)),
              as.integer(1:nrow(posmat)),
              rc=double(length(r)),
              as.integer(remove_self)
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
##' 
##' @export
##'
get.pi.b.typed.clustsurvey.hh.wts <- function(posmat,
                                              typeA = -1,
                                              typeB = -1,
                                              r=1,
                                              r.low=rep(0,length(r)),
                                              remove_self=FALSE) {
    return(.C("get_pi_b_typed_clustsurvey_hh_wts",
              as.double(posmat[,"p"]),
              as.integer(posmat[,"type_a"]),
              as.integer(posmat[,"type_b"]),
              as.double(posmat[,"x"]),
              as.double(posmat[,"y"]),
              as.integer(posmat[,"s"]),
              as.double(posmat[,"weight"]),
              as.double(posmat[,"delta"]),
              as.double(posmat[,"alpha"]),
              as.double(posmat[,"rho1"]),
              as.double(posmat[,"rho2"]),
              as.integer(nrow(posmat)),
              as.integer(typeA),
              as.integer(typeB),
              as.double(r.low),
              as.double(r),
              as.integer(length(r)),
              as.integer(1:nrow(posmat)),
              rc=double(length(r)),
              as.integer(remove_self)
    )$rc)
}




##' Calculate bootstrapped confidence intervals for \code{get.pi.typed.wts} values.
##'
##' Wrapper to \code{get.pi.typed.wts.bootstrap} that takes care of calculateing the
##' confience interavals based on the bootstrapped values..
##'
##'
##' @param posmat a matrix with columns type, x and y
##' @param fun the function to decide relationships
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range. 0 by default
##' @param boot.iter the number of bootstrap iterations
##' @param ci.low the low end of the ci...0.025 by default
##' @param ci.high the high end of the ci...0.975 by default
##'
##' @return a matrix with a row for the high and low values and
##'     a column per distance
##'
##' @author Justin Lessler
##'
##' @family get.pi
##'
##' 
##' @export
##'
get.pi.typed.wts.ci <- function(posmat,
                                typeA = -1,
                                typeB = -1,
                                r=1,
                                r.low=rep(0,length(r)),
                                boot.iter = 1000,
                                ci.low=0.025,
                                ci.high=0.975) {
    boots <- get.pi.typed.wts.bootstrap(posmat, typeA = -1, typeB = -1, r, r.low, boot.iter)
    
    rc <- matrix(nrow=2,ncol=ncol(boots))
    
    rownames(rc) <- c(ci.low,ci.high)
    
    for (i in 1:ncol(rc)) {
        rc[,i] <- stats::quantile(boots[,i], probs=c(ci.low, ci.high))
    }
    
    return(rc)
}



##' Calculate bootstrapped confidence intervals for \code{get.pi.clustsurvey} values.
##'
##' Wrapper to \code{get.pi.clustsurvey.bootstrap} that takes care of calculateing the
##' confience interavals based on the bootstrapped values..
##'
##'
##' @param posmat a matrix with columns p, x, y, s, and weight
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range. 0 by default
##' @param boot.iter the number of bootstrap iterations
##' @param ci.low the low end of the ci...0.025 by default
##' @param ci.high the high end of the ci...0.975 by default
##'
##' @return a matrix with a row for the high and low values and
##'     a column per distance
##'
##' @author Justin Lessler
##'
##' @family get.pi
##'
##' 
##' @export
##'
get.pi.clustsurvey.ci <- function(posmat,
                                  r=1,
                                  r.low=rep(0,length(r)),
                                  boot.iter = 1000,
                                  ci.low=0.025,
                                  ci.high=0.975) {
    boots <- get.pi.clustsurvey.bootstrap(posmat, r, r.low, boot.iter)
    
    rc <- matrix(nrow=2,ncol=ncol(boots))
    
    rownames(rc) <- c(ci.low,ci.high)
    
    for (i in 1:ncol(rc)) {
        if (length(which(is.na(boots[,i])))==nrow(boots)) next
        rc[,i] <- stats::quantile(boots[,i], probs=c(ci.low, ci.high), na.rm=TRUE)
    }
    
    return(rc)
}





##' runs bootstrapping on \code{get.pi.typed.wts}
##'
##' Bootstraps typed pi values. Makes sure distances between a sample and
##' another draw of itself are left out
##'
##' @param boot.iter the number of bootstrap iterations
##' @param posmat a matrix with columns type, x and y
##' @param typeA the "from" type that we are interested in, -1 is wildcard
##' @param typeB the "to" type that we are interested i, -1 is wildcard
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##'
##' @return pi values for all the distances we looked at
##'
##' @family get.pi
##'
##' 
##' @export
##'
get.pi.typed.wts.bootstrap <- function(posmat,
                                       typeA = -1,
                                       typeB = -1,
                                       r=1,
                                       r.low=rep(0,length(r)),
                                       boot.iter) {
    
    rc <- matrix(nrow=boot.iter, ncol=length(r))
    for (i in 1:boot.iter) {
        inds <- sample(nrow(posmat), replace=T)
        rc[i,] <- .C("get_pi_typed_wts",
                     as.integer(posmat[inds,"type"]),
                     as.double(posmat[inds,"x"]),
                     as.double(posmat[inds,"y"]),
                     as.double(posmat[inds,"weight"]),
                     as.integer(nrow(posmat)),
                     as.integer(typeA),
                     as.integer(typeB),
                     as.double(r.low),
                     as.double(r),
                     as.integer(length(r)),
                     as.integer(inds),
                     rc=double(length(r))
        )$rc
    }
    return(rc)
}



##' runs bootstrapping on \code{get.pi.typed}
##'
##' Bootstraps typed pi values. Makes sure distances between a sample and
##' another draw of itself are left out
##'
##' @param boot.iter the number of bootstrap iterations
##' @param posmat a matrix with columns type, x and y
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##'
##' @return pi values for all the distances we looked at
##'
##' @family get.pi
##'
##' 
##' @export
##'
get.pi.clustsurvey.bootstrap <- function(posmat,
                                         r=1,
                                         r.low=rep(0,length(r)),
                                         boot.iter) {
    
    rc <- matrix(nrow=boot.iter, ncol=length(r))
    for (i in 1:boot.iter) {
        inds <- sample(nrow(posmat), replace=T)
        rc[i,] <- .C("get_pi_clustsurvey",
                     as.double(posmat[inds,"p"]),
                     as.double(posmat[inds,"x"]),
                     as.double(posmat[inds,"y"]),
                     as.integer(posmat[inds,"s"]),
                     as.integer(nrow(posmat)),
                     as.double(r.low),
                     as.double(r),
                     as.integer(length(r)),
                     as.integer(inds),
                     as.double(coord_type_geo),
                     as.character(dist_unit),
                     rc=double(length(r))
        )$rc
    }
    return(rc)
}




##' runs bootstrapping on \code{get.pi.clustsurvey.wts}
##'
##' Bootstraps typed pi values. Makes sure distances between a sample and
##' another draw of itself are left out
##'
##' @param boot.iter the number of bootstrap iterations
##' @param posmat a matrix with columns type, x and y
##' @param typeA the "from" type that we are interested in, -1 is wildcard
##' @param typeB the "to" type that we are interested i, -1 is wildcard
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##'
##' @return pi values for all the distances we looked at
##'
##' @family get.pi
##'
##' 
##' @export
##'
get.pi.a.clustsurvey.wts.bootstrap <- function(posmat,
                                               r=1,
                                               r.low=rep(0,length(r)),
                                               boot.iter,
                                               remove_self=FALSE) {
    
    rc <- matrix(nrow=boot.iter, ncol=length(r))
    for (i in 1:boot.iter) {
        inds <- sample(nrow(posmat), replace=T)
        rc[i,] <- .C("get_pi_a_clustsurvey_wts",
                     as.double(posmat[inds,"p"]),
                     as.double(posmat[inds,"x"]),
                     as.double(posmat[inds,"y"]),
                     as.integer(posmat[inds,"s"]),
                     as.double(posmat[inds,"weight"]),
                     as.double(posmat[,"delta"]),
                     as.double(posmat[,"alpha"]),
                     as.integer(nrow(posmat)),
                     as.double(r.low),
                     as.double(r),
                     as.integer(length(r)),
                     as.integer(inds),
                     rc=double(length(r)),
              as.integer(remove_self)
        )$rc
    }
    return(rc)
}


##' runs bootstrapping on \code{get.pi.clustsurvey.wts}
##'
##' Bootstraps typed pi values. Makes sure distances between a sample and
##' another draw of itself are left out
##'
##' @param boot.iter the number of bootstrap iterations
##' @param posmat a matrix with columns type, x and y
##' @param typeA the "from" type that we are interested in, -1 is wildcard
##' @param typeB the "to" type that we are interested i, -1 is wildcard
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##'
##' @return pi values for all the distances we looked at
##'
##' @family get.pi
##'
##' 
##' @export
##'
get.pi.b.clustsurvey.wts.bootstrap <- function(posmat,
                                               r=1,
                                               r.low=rep(0,length(r)),
                                               boot.iter,
                                               remove_self = FALSE) {
    
    rc <- matrix(nrow=boot.iter, ncol=length(r))
    for (i in 1:boot.iter) {
        inds <- sample(nrow(posmat), replace=T)
        rc[i,] <- .C("get_pi_b_clustsurvey_wts",
                     as.double(posmat[inds,"p"]),
                     as.double(posmat[inds,"x"]),
                     as.double(posmat[inds,"y"]),
                     as.integer(posmat[inds,"s"]),
                     as.double(posmat[inds,"weight"]),
                     as.double(posmat[,"delta"]),
                     as.double(posmat[,"alpha"]),
                     as.integer(nrow(posmat)),
                     as.double(r.low),
                     as.double(r),
                     as.integer(length(r)),
                     as.integer(inds),
                     rc=double(length(r)),
                     as.integer(remove_self)
        )$rc
    }
    return(rc)
}




##' get the null distribution of the \code{get.pi.typed.wts} function
##'
##' Does permutations to calculate the null distribution of get pi
##' if there were no spatial dependence. Randomly reassigns coordinates
##' to each observation permutations times
##'
##' @param posmat a matrix with columns type, x and y
##' @param typeA the "from" type that we are interested in, -1 is wildcard
##' @param typeB the "to" type that we are interested i, -1 is wildcard
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##' @param permutations the number of permute iterations
##'
##' @return pi values for all the distances we looked at
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.pi
##'
##' 
##' @export
##'
get.pi.typed.wts.permute <- function(posmat,
                                     typeA = -1,
                                     typeB = -1,
                                     r=1,
                                     r.low=rep(0,length(r)),
                                     permutations) {
    
    xcol <-  which(colnames(posmat)=="x")
    ycol <- which(colnames(posmat)=="y")
    
    #check that both columns exist
    if (length(xcol)!=1 & length(ycol)!=1) {
        stop("unique x and y columns must be defined")
    }
    
    rc <- matrix(nrow=permutations, ncol=length(r))
    for (i in 1:permutations) {
        inds <- sample(nrow(posmat))#, replace=T)
        rc[i,] <- .C("get_pi_typed_wts",
                     as.integer(posmat[,"type"]),
                     as.double(posmat[inds,"x"]),
                     as.double(posmat[inds,"y"]),
                     as.double(posmat[inds,"weight"]),
                     as.integer(nrow(posmat)),
                     as.integer(typeA),
                     as.integer(typeB),
                     as.double(r.low),
                     as.double(r),
                     as.integer(length(r)),
                     as.integer(1:nrow(posmat)),
                     rc=double(length(r))
        )$rc
    }
    return(rc)
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
##' 
##' @export
##'
get.tau.typed.wts <- function(posmat,
                              typeA = -1,
                              typeB = -1,
                              r=1,
                              r.low=rep(0,length(r))) {
    return(.C("get_tau_typed_wts",
              as.integer(posmat[,"type"]),
              as.double(posmat[,"x"]),
              as.double(posmat[,"y"]),
              as.double(posmat[,"weight"]),
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
##' 
##' @export
##'
get.tau.clustsurvey <- function(posmat,
                                r=1,
                                r.low=rep(0,length(r)),
                                coord_type_geo=FALSE,
                                dist_unit = "K") {
    
    dist_unit <- toupper(substr(dist_unit, 1, 1))
    
    return(.C("get_tau_clustsurvey",
              as.double(posmat[,"p"]),
              as.double(posmat[,"x"]),
              as.double(posmat[,"y"]),
              as.integer(posmat[,"s"]),
              as.integer(nrow(posmat)),
              as.double(r.low),
              as.double(r),
              as.integer(length(r)),
              as.integer(1:nrow(posmat)),
              as.double(coord_type_geo),
              as.character(dist_unit),
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
##' 
##' @export
##'
get.tau.clustsurvey.wts <- function(posmat,
                                    r=1,
                                    r.low=rep(0,length(r)),
                                    remove_self=FALSE) {
    return(.C("get_tau_clustsurvey_wts",
              as.double(posmat[,"p"]),
              as.double(posmat[,"x"]),
              as.double(posmat[,"y"]),
              as.integer(posmat[,"s"]),
              as.double(posmat[,"weight"]),
              as.double(posmat[,"delta"]),
              as.double(posmat[,"alpha"]),
              as.integer(nrow(posmat)),
              as.double(r.low),
              as.double(r),
              as.integer(length(r)),
              as.integer(1:nrow(posmat)),
              rc=double(length(r)),
              numerator=double(length(r)),
              denominator=double(length(r)),
              as.integer(remove_self)
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
##' 
##' @export
##'
get.tau.typed.clustsurvey.wts <- function(posmat,
                                          typeA = -1,
                                          typeB = -1,
                                          r=1,
                                          r.low=rep(0,length(r)),
                                          remove_self=FALSE) {
    return(.C("get_tau_typed_clustsurvey_wts",
              as.double(posmat[,"p"]),
              as.integer(posmat[,"type_a"]),
              as.integer(posmat[,"type_b"]),
              as.double(posmat[,"x"]),
              as.double(posmat[,"y"]),
              as.integer(posmat[,"s"]),
              as.double(posmat[,"weight"]),
              as.double(posmat[,"delta"]),
              as.double(posmat[,"alpha"]),
              as.integer(nrow(posmat)),
              as.integer(typeA),
              as.integer(typeB),
              as.double(r.low),
              as.double(r),
              as.integer(length(r)),
              as.integer(1:nrow(posmat)),
              rc=double(length(r)),
              numerator=double(length(r)),
              denominator=double(length(r)),
              as.integer(remove_self)
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
##' 
##' @export
##'
get.tau.typed.clustsurvey.hh.wts <- function(posmat,
                                             typeA = -1,
                                             typeB = -1,
                                             r=1,
                                             r.low=rep(0,length(r)),
                                             remove_self=FALSE) {
    return(.C("get_tau_typed_clustsurvey_hh_wts",
              as.double(posmat[,"p"]),
              as.integer(posmat[,"type_a"]),
              as.integer(posmat[,"type_b"]),
              as.double(posmat[,"x"]),
              as.double(posmat[,"y"]),
              as.integer(posmat[,"s"]),
              as.double(posmat[,"weight"]),
              as.double(posmat[,"delta"]),
              as.double(posmat[,"alpha"]),
              as.double(posmat[,"rho1"]),
              as.double(posmat[,"rho2"]),
              as.integer(nrow(posmat)),
              as.integer(typeA),
              as.integer(typeB),
              as.double(r.low),
              as.double(r),
              as.integer(length(r)),
              as.integer(1:nrow(posmat)),
              rc=double(length(r)),
              numerator=double(length(r)),
              denominator=double(length(r)),
              as.integer(remove_self)
    )$rc)
}







##' Bootstrap confidence interval for the \code{get.tau.typed.wts} values
##'
##' Wrapper to \code{get.tau.bootstrap} that takes care of calulating
##' the confidence intervals based on the bootstrapped values
##'
##' @param posmat a matrix appropriate for input to \code{get.tau}
##' @param fun a function appropriate as input to \code{get.pi}
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##' @param boot.iter the number of bootstrap iterations
##' @param ci.low the low end of the ci...0.025 by default
##' @param ci.high the high end of the ci...0.975 by default
##'
##' @return tau values for all the distances examined
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.tau
##'
##' 
##' @export
##'
get.tau.typed.wts.ci <- function(posmat,
                                 typeA = -1,
                                 typeB = -1,
                                 r=1,
                                 r.low=rep(0,length(r)),
                                 boot.iter = 1000,
                                 ci.low=0.025,
                                 ci.high=0.975) {
    boots <- get.tau.typed.wts.bootstrap(posmat, typeA = -1, typeB = -1, r, r.low, boot.iter)
    
    rc <- matrix(nrow=2, ncol=ncol(boots))
    
    rownames(rc) <- c(ci.low,ci.high)
    
    for (i in 1:ncol(rc)) {
        rc[,i] <- stats::quantile(boots[,i], probs=c(ci.low, ci.high), na.rm = T)
    }
    
    return(rc)
}


##' Bootstrap confidence interval for the \code{get.tau.clustsurvey} values
##'
##' Wrapper to \code{get.tau.bootstrap} that takes care of calulating
##' the confidence intervals based on the bootstrapped values
##'
##' @param posmat a matrix appropriate for input to \code{get.tau}
##' @param fun a function appropriate as input to \code{get.pi}
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##' @param boot.iter the number of bootstrap iterations
##' @param ci.low the low end of the ci...0.025 by default
##' @param ci.high the high end of the ci...0.975 by default
##'
##' @return tau values for all the distances examined
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.tau
##'
##' 
##' @export
##'
get.tau.clustsurvey.ci <- function(posmat,
                                   r=1,
                                   r.low=rep(0,length(r)),
                                   boot.iter = 1000,
                                   ci.low=0.025,
                                   ci.high=0.975) {
    boots <- get.tau.clustsurvey.bootstrap(posmat, r, r.low, boot.iter)
    
    rc <- matrix(nrow=2, ncol=ncol(boots))
    
    rownames(rc) <- c(ci.low,ci.high)
    
    for (i in 1:ncol(rc)) {
        if (length(which(is.na(boots[,i])))==nrow(boots)) next
        rc[,i] <- stats::quantile(boots[,i], probs=c(ci.low, ci.high), na.rm = T)
    }
    
    return(rc)
}



##' Bootstrap confidence interval for the \code{get.tau.clustsurvey.wts} values
##'
##' Wrapper to \code{get.tau.bootstrap} that takes care of calulating
##' the confidence intervals based on the bootstrapped values
##'
##' @param posmat a matrix appropriate for input to \code{get.tau}
##' @param fun a function appropriate as input to \code{get.pi}
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##' @param boot.iter the number of bootstrap iterations
##' @param ci.low the low end of the ci...0.025 by default
##' @param ci.high the high end of the ci...0.975 by default
##'
##' @return tau values for all the distances examined
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.tau
##'
##' 
##' @export
##'
get.tau.clustsurvey.wts.ci <- function(posmat,
                                       r=1,
                                       r.low=rep(0,length(r)),
                                       boot.iter = 1000,
                                       ci.low=0.025,
                                       ci.high=0.975) {
    boots <- get.tau.clustsurvey.wts.bootstrap(posmat, r, r.low, boot.iter)
    
    rc <- matrix(nrow=2, ncol=ncol(boots))
    
    rownames(rc) <- c(ci.low,ci.high)
    
    for (i in 1:ncol(rc)) {
        if (length(which(is.na(boots[,i])))==nrow(boots)) next
        rc[,i] <- stats::quantile(boots[,i], probs=c(ci.low, ci.high), na.rm = T)
    }
    
    return(rc)
}



##' runs bootstrapping for \code{get.tau.typed.wts}
##'
##' @param posmat a matrix with columns type, x and y
##' @param typeA the "from" type that we are interested in, -1 is wildcard
##' @param typeB the "to" type that we are interested i, -1 is wildcard
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##' @param boot.iter the number of bootstrap iterations
##'
##' @return tau values for all the distances we looked at
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.tau
##'
##'
##' @export
##'
get.tau.typed.wts.bootstrap <- function(posmat,
                                        typeA = -1,
                                        typeB = -1,
                                        r=1,
                                        r.low=rep(0,length(r)),
                                        boot.iter) {
    
    rc <- matrix(nrow=boot.iter, ncol=length(r))
    for (i in 1:boot.iter) {
        inds <- sample(nrow(posmat), replace=T)
        rc[i,] <- .C("get_tau_typed_wts",
                     as.integer(posmat[inds,"type"]),
                     as.double(posmat[inds,"x"]),
                     as.double(posmat[inds,"y"]),
                     as.double(posmat[inds,"weight"]),
                     as.integer(nrow(posmat)),
                     as.integer(typeA),
                     as.integer(typeB),
                     as.double(r.low),
                     as.double(r),
                     as.integer(length(r)),
                     as.integer(inds),
                     rc=double(length(r))
        )$rc
    }
    return(rc)
}





##' runs bootstrapping for \code{get.tau.typed} using parallel processing with package \code{doParallel}
##'
##' @param posmat a matrix with columns type, x and y
##' @param typeA the "from" type that we are interested in, -1 is wildcard
##' @param typeB the "to" type that we are interested i, -1 is wildcard
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##' @param boot.iter the number of bootstrap iterations
##'
##' @return tau values for all the distances we looked at
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.tau
##'
##' 
##' @export
##' @import doParallel foreach
##'
get.tau.typed.bootstrap.parallel <- function(cores, posmat,
                                             typeA = -1,
                                             typeB = -1,
                                             r=1,
                                             r.low=rep(0,length(r)),
                                             boot.iter,
                                             comparison.type = "representative") {
    
    if (comparison.type == "representative") {
        comp.type.int <- 0
    } else if (comparison.type == "independent") {
        comp.type.int <- 1
    } else {
        stop("unkown comparison type specified")
    }
    
    registerDoParallel(cores=cores)
    
    rc <- foreach(i=1:boot.iter, .errorhandling="pass", .combine=rbind, .multicombine=TRUE) %dopar% {
        inds <- sample(nrow(posmat), replace=T)
        .C("get_tau_typed",
           as.integer(posmat[inds,"type"]),
           as.double(posmat[inds,"x"]),
           as.double(posmat[inds,"y"]),
           as.integer(nrow(posmat)),
           as.integer(typeA),
           as.integer(typeB),
           as.double(r.low),
           as.double(r),
           as.integer(length(r)),
           as.integer(inds),
           as.integer(comp.type.int),
           rc=double(length(r))
        )$rc
    }
    return(rc)
}





##' runs bootstrapping for \code{get.tau.typed.wts} using parallel processing with package \code{doParallel}
##'
##' @param posmat a matrix with columns type, x and y
##' @param typeA the "from" type that we are interested in, -1 is wildcard
##' @param typeB the "to" type that we are interested i, -1 is wildcard
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##' @param boot.iter the number of bootstrap iterations
##'
##' @return tau values for all the distances we looked at
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.tau
##'
##' 
##' @export
##' @import doParallel
##'
get.tau.typed.wts.bootstrap.parallel <- function(cores, posmat,
                                                 typeA = -1,
                                                 typeB = -1,
                                                 r=1,
                                                 r.low=rep(0,length(r)),
                                                 boot.iter) {
    
    registerDoParallel(cores=cores)
    
    rc <- foreach(i=1:boot.iter, .errorhandling="pass", .combine=rbind, .multicombine=TRUE) %dopar% {
        
        inds <- sample(nrow(posmat), replace=T)
        .C("get_tau_typed_wts",
           as.integer(posmat[inds,"type"]),
           as.double(posmat[inds,"x"]),
           as.double(posmat[inds,"y"]),
           as.double(posmat[inds,"weight"]),
           as.integer(nrow(posmat)),
           as.integer(typeA),
           as.integer(typeB),
           as.double(r.low),
           as.double(r),
           as.integer(length(r)),
           as.integer(inds),
           rc=double(length(r))
        )$rc
    }
    return(rc)
}




##' runs bootstrapping for \code{get.tau.clustsurvey}
##'
##' @param posmat a matrix with columns type, x and y
##' @param typeA the "from" type that we are interested in, -1 is wildcard
##' @param typeB the "to" type that we are interested i, -1 is wildcard
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##' @param boot.iter the number of bootstrap iterations
##'
##' @return tau values for all the distances we looked at
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.tau
##'
##' 
##' @export
##'
get.tau.clustsurvey.bootstrap <- function(posmat,
                                          r=1,
                                          r.low=rep(0,length(r)),
                                          boot.iter) {
    
    rc <- matrix(nrow=boot.iter, ncol=length(r))
    for (i in 1:boot.iter) {
        inds <- sample(nrow(posmat), replace=T)
        rc[i,] <- .C("get_tau_clustsurvey",
                     as.double(posmat[inds,"p"]),
                     as.double(posmat[inds,"x"]),
                     as.double(posmat[inds,"y"]),
                     as.integer(posmat[inds,"s"]),
                     as.integer(nrow(posmat)),
                     as.double(r.low),
                     as.double(r),
                     as.integer(length(r)),
                     as.integer(inds),
                     as.double(coord_type_geo),
                     as.character(dist_unit),
                     rc=double(length(r))
        )$rc
    }
    return(rc)
}



##' runs bootstrapping for \code{get.tau.clustsurvey.wts}
##'
##' @param posmat a matrix with columns type, x and y
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##' @param boot.iter the number of bootstrap iterations
##'
##' @return tau values for all the distances we looked at
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.tau
##'
##' 
##' @export
##'
get.tau.clustsurvey.wts.bootstrap <- function(posmat,
                                              r=1,
                                              r.low=rep(0,length(r)),
                                              boot.iter) {
    
    rc <- matrix(nrow=boot.iter, ncol=length(r))
    for (i in 1:boot.iter) {
        inds <- sample(nrow(posmat), replace=T)
        rc[i,] <- .C("get_tau_clustsurvey_wts",
                     as.double(posmat[inds,"p"]),
                     as.double(posmat[inds,"x"]),
                     as.double(posmat[inds,"y"]),
                     as.integer(posmat[inds,"s"]),
                     as.double(posmat[inds,"weight"]),
                     as.double(posmat[,"delta"]),
                     as.double(posmat[,"alpha"]),
                     as.integer(nrow(posmat)),
                     as.double(r.low),
                     as.double(r),
                     as.integer(length(r)),
                     as.integer(inds),
                     rc=double(length(r)),
                     numerator=double(length(r)),
                     denominator=double(length(r)),
                     as.integer(remove_self)
        )$rc
    }
    return(rc)
}










##' get the null distribution for the \code{get.tau.typed} function
##'
##'
##' @param posmat a matrix with columns type, x and y
##' @param typeA the "from" type that we are interested in, -1 is wildcard
##' @param typeB the "to" type that we are interested i, -1 is wildcard
##' @param r the series of spatial distances wer are interested in
##' @param r.low the low end of each range....0  by default
##' @param permutations the number of permute iterations
##'
##' @return a matrix with permutation tau values for each distance specified
##'
##' @author Justin Lessler and Henrik Salje
##'
##' @family get.tau
##'
##'
##' @export
##' 
get.tau.typed.wts.permute <- function(posmat,
                                      typeA = -1,
                                      typeB = -1,
                                      r=1,
                                      r.low=rep(0,length(r)),
                                      permutations) {
    
    xcol <-  which(colnames(posmat)=="x")
    ycol <- which(colnames(posmat)=="y")
    
    #check that both columns exist
    if (length(xcol)!=1 & length(ycol)!=1) {
        stop("unique x and y columns must be defined")
    }
    
    rc <- matrix(nrow=permutations, ncol=length(r))
    for (i in 1:permutations) {
        inds <- sample(nrow(posmat))#, replace=T)
        rc[i,] <- .C("get_tau_typed_wts",
                     as.integer(posmat[,"type"]),
                     as.double(posmat[inds,"x"]),
                     as.double(posmat[inds,"y"]),
                     as.double(posmat[inds,"weight"]),
                     as.integer(nrow(posmat)),
                     as.integer(typeA),
                     as.integer(typeB),
                     as.double(r.low),
                     as.double(r),
                     as.integer(length(r)),
                     #as.integer(inds),
                     as.integer(1:nrow(posmat)),
                     rc=double(length(r))
        )$rc
    }
    return(rc)
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
##' 
##' @export
##'
get.pi.a.clustsurvey.hh.wts <- function(posmat,
                                        r=1,
                                        r.low=rep(0,length(r)),
                                        remove_self=FALSE) {
    return(.C("get_pi_a_clustsurvey_hh_wts",
              as.double(posmat[,"p"]),
              as.double(posmat[,"x"]),
              as.double(posmat[,"y"]),
              as.integer(posmat[,"s"]),
              as.double(posmat[,"weight"]),
              as.double(posmat[,"delta"]),
              as.double(posmat[,"alpha"]),
              as.double(posmat[,"rho1"]),
              as.double(posmat[,"rho2"]),
              as.integer(nrow(posmat)),
              as.double(r.low),
              as.double(r),
              as.integer(length(r)),
              as.integer(1:nrow(posmat)),
              rc=double(length(r)),
              as.integer(remove_self)
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
##' 
##' @export
##'
get.pi.b.clustsurvey.hh.wts <- function(posmat,
                                        r=1,
                                        r.low=rep(0,length(r)),
                                        remove_self=FALSE) {
    return(.C("get_pi_b_clustsurvey_hh_wts",
              as.double(posmat[,"p"]),
              as.double(posmat[,"x"]),
              as.double(posmat[,"y"]),
              as.integer(posmat[,"s"]),
              as.double(posmat[,"weight"]),
              as.double(posmat[,"delta"]),
              as.double(posmat[,"alpha"]),
              as.double(posmat[,"rho1"]),
              as.double(posmat[,"rho2"]),
              as.integer(nrow(posmat)),
              as.double(r.low),
              as.double(r),
              as.integer(length(r)),
              as.integer(1:nrow(posmat)),
              rc=double(length(r)),
              as.integer(remove_self)
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
##'
##' @export
##' 
get.tau.clustsurvey.hh.wts <- function(posmat,
                                       r=1,
                                       r.low=rep(0,length(r)),
                                       remove_self=FALSE) {
    return(.C("get_tau_clustsurvey_hh_wts",
              as.double(posmat[,"p"]),
              as.double(posmat[,"x"]),
              as.double(posmat[,"y"]),
              as.integer(posmat[,"s"]),
              as.double(posmat[,"weight"]),
              as.double(posmat[,"delta"]),
              as.double(posmat[,"alpha"]),
              as.double(posmat[,"rho1"]),
              as.double(posmat[,"rho2"]),
              as.integer(nrow(posmat)),
              as.double(r.low),
              as.double(r),
              as.integer(length(r)),
              as.integer(1:nrow(posmat)),
              rc=double(length(r)),
              numerator=double(length(r)),
              denominator=double(length(r)),
              as.integer(remove_self)
    )$rc)
}
