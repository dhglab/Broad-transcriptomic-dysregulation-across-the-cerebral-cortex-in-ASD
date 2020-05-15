\name{predict.earth}
\alias{predict.earth}
\title{Predict with an earth model}
\description{
Predict with an \code{\link{earth}} model.
}
\usage{
\method{predict}{earth}(object = stop("no 'object' argument"), newdata = NULL,
        type = c("link", "response", "earth", "class", "terms"),
        interval = "none", level = .95,
        thresh = .5, trace = FALSE, \dots)
}
\arguments{
  \item{object}{
    An \code{\link{earth}} object.
    This is the only required argument.
  }
  \item{newdata}{
    Make predictions using \code{newdata}, which
    can be a data frame, a matrix, or a vector with length equal to a
    multiple of the number of columns
    of the original input matrix \code{x}.\cr
    Default is NULL, meaning return values predicted from the training set.\cr
    NAs are allowed in \code{newdata} (and the predicted value will be NA
    unless the NAs are in variables that are unused in the earth model).
  }
  \item{type}{
     Type of prediction.
     One of \code{"link"} (default), \code{"response"}, \code{"earth"}, \code{"class"}, or \code{"terms"}.
     See the \bold{Note} below.
  }
  \item{interval}{
    Return prediction or confidence levels.
    Default is \code{"none"}.
    Use \code{interval="pint"} to get prediction intervals on new data.
\cr
    Requires that the earth model was built with \code{varmod.method}.
\cr
    This argument gets passed on as the \code{type} argument to \code{\link{predict.varmod}}.
    See its help page for details.
  }
  \item{level}{
    Confidence level for the \code{interval} argument.
    Default is \code{0.95}, meaning construct 95\% confidence bands
    (estimate the 2.5\% and 97.5\% levels).
  }
  \item{thresh}{
    Threshold, a value between 0 and 1 when predicting a probability.
    Only applies when \code{type="class"}.
    Default is 0.5.
    See the \bold{Note} below.
  }
  \item{trace}{
     Default \code{FALSE}. Set to \code{TRUE} to see which data, subset, etc.  \code{predict.earth} is using.
  }
  \item{\dots}{
     Unused, but provided for generic/method consistency.
  }
}
\value{
The predicted values (a matrix for multiple response models).

If \code{type="terms"}, a matrix with each column showing the contribution of a predictor.

If \code{interval="pint"} or \code{"cint"}, a matrix with three columns:\cr
\code{fit}: the predicted values\cr
\code{lwr}: the lower confidence or prediction limit\cr
\code{upr}: the upper confidence or prediction limit

If \code{interval="se"}, the standard errors.
}
\note{
\bold{Predicting with standard earth models}

    Use the default \code{type="link"}, or possibly \code{type="class"}.

    Actually, the \code{"link"}, \code{"response"}, and \code{"earth"}
    choices all return the same value unless the \code{glm} argument
    was used in the original call to \code{\link{earth}}.

\bold{Predicting with earth-GLM models}

    This section applies to earth models with a GLM component, i.e.,
    when the \code{glm} argument was used
    in the original call to \code{\link{earth}}.

    The \code{"link"} and \code{"response"} options:
    see \code{\link{predict.glm}} for a description of these.
    In brief: for logistic models
    use \code{type="link"} to get log-odds and
    \code{type="response"} to get probabilities.

    Use option \code{"earth"} to get the linear fit (this gives the prediction you would get
    if your original call to earth had no \code{glm} argument).

\bold{Predicting with "class"}

    Use option \code{"class"} to get the predicted class.
    With option \code{"class"}, this function first makes predictions with
    \code{type="response"} and then assigns the predicted values to classes as follows:

    (i) When the response is a \emph{logical}, predict \code{TRUE} if
    the predicted probability is greater than \code{thresh}.

    (ii) When the response is a \emph{numeric}, predict \code{TRUE} if
    the predicted value is greater than \code{thresh}.
    Actually, this is identical to the above case,
    although \code{thresh} here may legitimately be a value
    outside the 0...1 range.

    (iii) When the response is a \emph{two level factor},
    predict the second level if its probability is more than \code{thresh}.
    In other words, with the default \code{thresh=0.5} predict the most probable level.

    (iv) When the response is a \emph{three or more level factor},
    predict the most probable level (and \code{thresh} is ignored).

\bold{Predicting with "terms"}

    The \code{"terms"} option returns a \code{"link"} response suitable for \code{\link{termplot}}.
    Only the additive terms and the first response (for multi-response models) are returned.
    Also, \code{"terms"} always returns the earth terms, and ignores the GLM component
    of the model, if any.
}
\seealso{
  \code{\link{earth}},
  \code{\link{predict}}
}
\examples{
data(trees)
earth.mod <- earth(Volume ~ ., data = trees)
predict(earth.mod)           # same as earth.mod$fitted.values
predict(earth.mod, c(10,80)) # yields 16.8
}
\keyword{models}
