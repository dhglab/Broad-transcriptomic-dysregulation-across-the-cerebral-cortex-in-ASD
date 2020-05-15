\name{predict.varmod}
\alias{predict.varmod}
\title{Predict with a varmod model}
\description{
You probably won't need to call this function directly.
It is called by \code{\link{predict.earth}}
when that function's \code{interval} argument is used.
}
\usage{
\method{predict}{varmod}(
    object  = stop("no 'object' argument"),
    newdata = NULL,
    type    = c("pint", "cint", "se", "abs.residual"),
    level   = .95,
    trace   = FALSE,
    \dots)
}
\arguments{
  \item{object}{
    A \code{varmod} object.
  }
  \item{newdata}{
    Make predictions using \code{newdata}.
    Default is NULL, meaning return values predicted from the training set.
  }
  \item{type}{
    Type of prediction. This is the \code{interval} argument of
\code{\link{predict.earth}}. One of
\cr\cr
\code{"pint"} Prediction intervals.
\cr\cr
\code{"cint"} Confidence intervals.
Cannot be used with \code{newdata}.
\cr\cr
\code{"se"} Standard error of the parent model residuals.
\cr\cr
\code{"abs.residual"} The absolute residuals of the parent model
on which the residual model regresses.
\cr\cr
  }
  \item{level}{
    Confidence level for the \code{interval} argument.
    Default is \code{.95}, meaning construct 95\% confidence bands
    (estimate the 2.5\% and 97.5\% levels).
  }
  \item{trace}{
     Currently unused.
 }
  \item{\dots}{
     Unused, but provided for generic/method consistency.
 }
}
\note{
\code{predict.varmod} is called by \code{predict.earth}
when its \code{interval} argument is used.
}
\examples{
data(ozone1)

set.seed(1) # optional, for cross validation reproducibility

# note: should really use ncross=30 below but for a quick demo we don't

earth.mod <- earth(O3~temp, data=ozone1, nfold=10, ncross=3, varmod.method="lm")

# call predict.earth, which calls predict.varmod

predict(earth.mod, newdata=ozone1[200:203,], interval="pint", level=.95)
}
\seealso{
  \code{\link{predict.earth}}
  \code{\link{varmod}}
}
