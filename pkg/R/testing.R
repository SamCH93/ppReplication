#' @title Bayes factor for testing effect size
#'
#' @description This function computes the Bayes factor contrasting \eqn{H_0 :
#'     \theta = 0}{H0: theta = 0} to \eqn{H_1 : \theta \sim f(\theta |
#'     \code{to}, \code{so}, \alpha)}{H1: theta ~ f(theta|original data, alpha)}
#'     for the replication data assuming a normal likelihood. The prior of the
#'     effect size \eqn{\theta}{theta} under \eqn{H_1}{H1} is the posterior of
#'     the effect size obtained from combining an initial flat prior with a
#'     normal likelihood of the original data raised to the power of
#'     \eqn{\alpha}{alpha}. Under \eqn{H_1}{H1}, the power parameter can either
#'     be fixed to some value between 0 and 1, or it can have a Beta
#'     distribution \eqn{\alpha | H_1 \sim \mbox{Beta}(\code{x},
#'     \code{y})}{alpha|H1 ~ Beta(\code{x}, \code{y})}.
#'
#'
#' @param tr Effect estimate of the replication study.
#' @param to Effect estimate of the original study.
#' @param sr Standard error of the replication effect estimate.
#' @param so Standard error of the replication effect estimate.
#' @param x Number of successes parameter \eqn{x}{x} for beta prior of power
#'     parameter under \eqn{H_1}{H1}. Defaults to 1. Is only taken into account
#'     when \code{alpha = NA}.
#' @param y Number of failures parameter \eqn{y}{y} for beta prior of power
#'     parameter under \eqn{H_1}{H1}. Defaults to 1. Is only taken into account
#'     when \code{alpha = NA}.
#' @param alpha Power parameter under \eqn{H_1}{H1}. Can be set to a number
#'     between 0 and 1. Defaults to \code{NA}.
#' @param ... Additional arguments for integration function.
#'
#' @return Bayes factor (BF > 1 indicates evidence for \eqn{H_0}{H0}, whereas BF
#'     < 1 indicates evidence for \eqn{H_1}{H1})
#'
#' @author Samuel Pawel
#'
#' @seealso \code{\link{bfPPalpha}}
#'
#' @examples
#' ## uniform prior on power parameter
#' bfPPtheta(tr = 0.09,  sr = 0.0518, to = 0.205, so = 0.0506)
#'
#' ## power parameter fixed to alpha = 1
#' bfPPtheta(tr = 0.090, sr = 0.0518, to = 0.205, so = 0.0506, alpha = 1)
#' @export
bfPPtheta <- function(tr, sr, to, so, x = 1, y = 1, alpha = NA, ...) {
    ## input checks
    stopifnot(
        length(tr) == 1,
        is.numeric(tr),
        is.finite(tr),

        length(to) == 1,
        is.numeric(to),
        is.finite(to),

        length(sr) == 1,
        is.numeric(sr),
        is.finite(sr),
        0 < sr,

        length(so) == 1,
        is.numeric(so),
        is.finite(so),
        0 < so,

        length(x) == 1,
        is.numeric(x),
        is.finite(x),
        0 <= x,

        length(y) == 1,
        is.numeric(y),
        is.finite(y),
        0 <= y,

        length(alpha) == 1,
        (is.na(alpha) |
         ((is.numeric(alpha)) &
          (is.finite(alpha)) &
          (0 <= alpha) &
          (alpha <= 1)))
    )

    ## marginal density under H0
    fH0 <- stats::dnorm(x = tr, mean = 0, sd = sr)

    ## marginal density under H1
    if (!is.na(alpha)) {
        fH1 <- stats::dnorm(x = tr, mean = to, sd = sqrt(sr^2 + so^2/alpha))
    } else {
        fH1 <- margLik(tr = tr, to = to, sr = sr, so = so, x = x, y = y,
                            ... = ...)
    }

    ## compute BF
    bf01 <- fH0/fH1
    return(bf01)
}


#' @title Bayes factor for testing power parameter
#'
#' @description This function computes the Bayes factor contrasting \eqn{H_0 :
#'     \alpha = 0}{H0: alpha = 0} to \eqn{H_0 : \alpha = 1}{H0: alpha = 1} for
#'     the replication data assuming a normal likelihood. The power parameter
#'     \eqn{\alpha}{alpha} indicates how much the normal likelihood of the
#'     original is raised and then incorporated in the power prior for the
#'     effect size \eqn{\theta}{theta} (e.g. for \eqn{\alpha = 1}{alpha = 1} the
#'     original data are completely discounted). An initial unit-information
#'     prior \eqn{\theta \sim \mathrm{N}(0, \mbox{uv})}{theta ~ N(0, uv)} is
#'     assumed for the effect size \eqn{\theta}{theta} under either hypothesis.
#'
#' @param tr Effect estimate of the replication study.
#' @param to Effect estimate of the original study.
#' @param sr Standard error of the replication effect estimate.
#' @param so Standard error of the replication effect estimate.
#' @param uv Variance of the unit-information prior.
#'
#' @return Bayes factor (BF > 1 indicates evidence for \eqn{H_0}{H0}, whereas BF
#'     < 1 indicates evidence for \eqn{H_1}{H1})
#'
#' @author Samuel Pawel
#'
#' @seealso \code{\link{bfPPtheta}}
#'
#' @examples
#' ## use unit variance of 2
#' bfPPalpha(tr = 0.09,  sr = 0.0518, to = 0.205, so = 0.0506, uv = 2)
#' @export
bfPPalpha <- function(tr, sr, to, so, uv) {
    ## input checks
    stopifnot(
        length(tr) == 1,
        is.numeric(tr),
        is.finite(tr),

        length(to) == 1,
        is.numeric(to),
        is.finite(to),

        length(sr) == 1,
        is.numeric(sr),
        is.finite(sr),
        0 < sr,

        length(so) == 1,
        is.numeric(so),
        is.finite(so),
        0 < so,

        length(uv) == 1,
        is.numeric(uv),
        is.finite(uv),
        0 < uv
    )

    ## marginal density under H0
    fH0 <- stats::dnorm(x = tr, mean = 0, sd = sqrt(sr^2 + uv))

    ## marginal density under H1:
    ## updating theta ~ N(0, uv) with to|theta ~ N(theta, so^2) leads to
    ## posterior theta | to ~ N(g/(1 + g)*to, g/(1 + g)*so^2) where
    g <- uv/so^2
    s <- g/(1 + g)
    fH1 <- stats::dnorm(x = tr, mean = s*to, sd = sqrt(sr^2 + s*so^2))

    ## compute BF
    bf01 <- fH0/fH1
    return(bf01)
}
