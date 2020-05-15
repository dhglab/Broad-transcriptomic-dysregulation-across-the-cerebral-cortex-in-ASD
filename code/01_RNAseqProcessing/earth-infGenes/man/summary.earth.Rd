\name{summary.earth}
\alias{summary.earth}
\alias{print.summary.earth}
\title{Summary method for earth objects}
\description{
Summary method for \code{\link{earth}} objects.
}
\usage{
\method{summary}{earth}(object = stop("no 'object' argument"),
        details = FALSE,  style = c("h", "pmax", "max", "C", "bf"),
        decomp = "anova", digits = getOption("digits"), fixed.point=TRUE,
        newdata = NULL, \dots)

\method{print}{summary.earth}(x = stop("no 'x' argument"),
        details = x$details,
        decomp = x$decomp, digits = x$digits, fixed.point = x$fixed.point,
        newdata = x$newdata, \dots)
}
\arguments{
  \item{object}{
    An \code{\link{earth}} object.
    This is the only required argument for \code{summary.earth}.
  }
  \item{x}{
    A \code{\link{summary.earth}} object.
    This is the only required argument for \code{print.summary.earth}.
  }
  \item{details}{
    Default is \code{FALSE}.
    Use \code{TRUE}
    to print more information about \code{\link{earth}}--\code{\link{glm}} models.
    But note that the displayed P-values of the GLM coefficients
    are meaningless because of the amount of preprocessing
    by earth to select the regression terms.
  }
  \item{style}{
    Formatting style. One of\cr
\code{"h"} (default) more compact\cr
\code{"pmax"} for those who prefer it and for compatibility with old versions of earth\cr
\code{"max"} is the same as \code{"pmax"} but prints \code{max} rather than \code{pmax}\cr
\code{"C"} C style expression with zero based indexing\cr
\code{"bf"} basis function format.
  }
  \item{decomp}{
    Specify how terms are ordered.
    Default is \code{"anova"}.
    Use \code{"none"} to order the terms as created by the forward.pass.
    See \code{\link{format.earth}} for a full description.
  }
  \item{digits}{
    The number of significant digits.\cr
    For \code{summary.earth}, the default is \code{getOption("digits")}.\cr
    For \code{print.summary.earth}, the default is the \code{$digits} component of \code{object}.
  }
  \item{fixed.point}{
    Method of printing numbers in matrices.
    Default is \code{TRUE} which prints like this
    (making it easier to compare coefficients):\preformatted{
        (Intercept)    15.029
        h(temp-58)      0.313
        h(234-ibt)     -0.046
        ...
} whereas \code{fixed.point=FALSE} prints like this (which is more usual in R):\preformatted{
        (Intercept)   1.5e+01
        h(temp-58)    3.1e-01
        h(234-ibt)   -4.6e-02
        ...
}
Matrices with two or fewer rows are never printed with a fixed point.
  }
  \item{newdata}{
Default \code{NULL}.\cr
Else print R-Squared for the new data
(and the returned object will have \code{newrsq} and \code{newdata} fields).
Additionally, if a variance model is present print the interval coverage table for the new data.\cr
  }
  \item{\dots}{
    Extra arguments are passed to \code{\link{format.earth}}.
  }
}
\value{
The value is the same as that returned by \code{\link{earth}}
but with the following extra components.\cr

  \item{strings}{
    String(s) created by \code{\link{format.earth}}.
    For multiple response models, a vector of strings.
  }
  \item{newrsq}{Only if \code{newdata} was passed to \code{summary.earth}.}
  \item{newdata}{Only if \code{newdata} was passed to \code{summary.earth}.}
  \item{digits}{}
  \item{details}{}
  \item{decomp}{}
  \item{fixed.point}{The corresponding arguments, passed on to \code{print.summary.earth}.}
}
\note{
The printed \code{Estimated importance} uses \code{\link{evimp}}
with the \code{nsubsets} criterion.
The most important predictor is printed first, and so on.
}
\seealso{
  \code{\link{earth}},
  \code{\link{evimp}},
  \code{\link{format.earth}}
}
\examples{
earth.mod <- earth(Volume~ ., data = trees)
summary(earth.mod, digits = 2)
}
\keyword{models}
