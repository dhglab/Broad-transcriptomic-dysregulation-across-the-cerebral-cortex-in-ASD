\name{residuals.earth}
\alias{residuals.earth}
\title{Residuals for an earth model}
\description{
Residuals of an \code{\link{earth}} model.
}
\usage{
\method{residuals}{earth}(object = stop("no 'object' argument"),
          type = NULL, warn = TRUE, \dots)
}
\arguments{
  \item{object}{
    An \code{\link{earth}} object.
    This is the only required argument.
  }
  \item{type}{
     One of\cr
        \code{"earth"} (default) Residuals (from the \code{\link{lm}} fit on \code{bx}).\cr
        \code{"standardize"} Residuals divided by \code{se * sqrt(1 - h_ii)}.
See the \code{standardize} argument of \code{\link{plot.earth}}.\cr
        \code{"delever"} Residuals divided by \code{sqrt(1 - h_ii)}.
See the \code{delever} argument of \code{\link{plot.earth}}.\cr
        \code{"deviance"}
            Residuals as above, unless the object has a \code{\link{glm}} component,
            in which case return the \code{\link{glm}} deviance residuals.\cr
        \code{"glm.pearson"}\cr
        \code{"glm.working"}\cr
        \code{"glm.response"}\cr
        \code{"glm.partial"}
            Return the corresponding \code{\link{glm}} residuals
            (from the \code{\link{glm}} fit on \code{bx}).
            Can be used only if the earth model has a \code{\link{glm}} component.
  }
  \item{warn}{
        This function gives warnings when the results are not what you may expect.
        Use \code{warn=FALSE} to turn of just these warnings.
  }
  \item{\dots}{
     Unused, but provided for generic/method consistency.
  }
}
\value{
The residual values (will be a matrix for multiple response models).
}
\seealso{
  \code{\link{earth}}\cr
  \code{\link{residuals}}\cr
  \code{\link{resid}} identical to \code{\link{residuals}}
}
\examples{
data(etitanic)
earth.mod <- earth(pclass ~ ., data=etitanic, glm=list(family=binomial))
head(resid(earth.mod, warn=FALSE))      # earth residuals, a column for each response
head(resid(earth.mod, type="earth"))    # same
head(resid(earth.mod, type="deviance")) # GLM deviance residuals, a column for each response
}
\keyword{models}
