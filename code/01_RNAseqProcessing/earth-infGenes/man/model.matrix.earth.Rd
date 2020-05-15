\name{model.matrix.earth}
\alias{model.matrix.earth}
\title{Get the earth basis matrix}
\description{
Get the basis matrix of an \code{\link{earth}} model.
}
\usage{
\method{model.matrix}{earth}(object = stop("no 'object' argument"),
    x = NULL, subset = NULL, which.terms = NULL,
    trace = 0,
    \dots,
    Env = parent.frame(),
    Callers.name = "model.matrix.earth")
}
\arguments{
  \item{object}{
    An \code{\link{earth}} model.
    This is the only required argument.
  }
  \item{x}{
Default is NULL, meaning use the original data
used to build the \code{earth} model
(after taking the original \code{subset}, if any).

Else \code{x} can be a data frame, a matrix, or a vector with
length equal to a multiple of the number of columns
of the original input matrix \code{x}.
(There is some leniency here.  For example, column names aren't
necessary if \code{x} has the same number of predictors originally
used to build the \code{earth} model.)
  }
  \item{subset}{
    Which rows to use in \code{x}.
    Default is NULL, meaning use all of \code{x}.
  }
  \item{which.terms}{
    Which terms to use.
    Default is NULL, meaning all terms in the earth model
    (i.e. the terms in \code{object$selected.terms}).
  }
  \item{trace}{
     Default 0. Set to non-zero to see which data \code{model.matrix.earth} is using.
  }
  \item{\dots}{
     Unused, but provided for generic/method consistency.
  }
  \item{Env}{
     For internal use.
  }
  \item{Callers.name}{
     For internal use (used by earth in trace messages).
  }
}
\value{
A basis matrix \code{bx} of the same form returned by \code{\link{earth}}.
The format of \code{bx} is described in \code{\link{earth.object}}.

If \code{x}, \code{subset}, and \code{which.terms} are all NULL (the
default), this function returns the model's \code{bx}.  In this case,
it is perhaps easier to simply use \code{object$bx}.\cr

The matrix \code{bx} can be used
as the input matrix to \code{\link{lm}} or \code{\link{glm}},
as shown below in the example.
In fact, that is what earth does internally after the pruning pass ---
it calls \code{\link{lm.fit}},
and additionally \code{\link{glm}} if earth's \code{glm} argument is used.
}
\seealso{
  \code{\link{earth}}
}
\examples{
data(trees)
earth.mod <- earth(Volume ~ ., data = trees)
summary(earth.mod, decomp = "none") # "none" to print terms in same order as lm.mod below

bx <- model.matrix(earth.mod)         # equivalent to bx <- earth.mod$bx
lm.mod <- lm(trees$Volume ~ bx[,-1])  # -1 to drop intercept
summary(lm.mod)                       # yields same coeffs as above summary
                                      # displayed t values are not meaningful
}
\keyword{models}
