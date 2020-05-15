\name{varmod}
\alias{varmod}
\alias{summary.varmod}
\title{Variance models for estimating prediction intervals}
\description{
A \emph{variance model} estimates the variance of predicted values.
It can be used to estimate prediction intervals.
See the \code{interval} argument of \code{\link{predict.earth}}.

A variance model is built by \code{earth} if \code{earth}'s
\code{varmod.method} argument is specified.
Results are stored in the \code{$varmod} field of the \code{earth} model.
See the vignette \dQuote{\href{../doc/earth-varmod.pdf}{Variance models in earth}} for details.

You probably won't need to directly call
\code{print.varmod} or \code{summary.varmod}.
They get called internally by \code{\link{summary.earth}}.
}
\usage{
\method{summary}{varmod}(
    object  = stop("no 'object' argument"),
    level   = .95,
    style   = "standard",
    digits  = 2,
    newdata = NULL,
    \dots)
}
\arguments{
  \item{object}{
    A \code{varmod} object.  This is the only required argument.
  }
  \item{level}{
    Same as \code{predict.earth}'s \code{level} argument.
  }
  \item{style}{
    Determines how the coefficients of the \code{varmod} are printed
by \code{summary.varmod}:\cr
\code{"standard"} (default)\cr
\code{"unit"} for easy comparison normalize the coefficients
by dividing by the first coefficient.
  }
  \item{digits}{
    Number of digits to print.  Default is \code{2}.
  }
  \item{newdata}{
Default \code{NULL}.\cr
Else print the interval coverage table for the new data.
  }
  \item{\dots}{
     Dots are passed on.
 }
}
\note{
A \code{"varmod"} object has the following fields:

\itemize{
    \item{\code{call}}
    The call used internally in the parent model to build the \code{varmod} object.

    \item{\code{parent}} The parent \code{earth} model.

    \item{\code{method}} Copy of the \code{varmod.method} argument to the parent model.

    \item{\code{package}} NULL, unless \code{method="gam"}, in which case either \code{"gam"} or \code{"mgcv"}.

    \item{\code{exponent}} Copy of the \code{varmod.exponent} argument to the parent model.

    \item{\code{lambda}} Currently always 1, meaning use absolute residuals.

    \item{\code{rmethod}} Currently always "hc2", meaning correct the residuals with \code{1/(1-h_ii)}.

    \item{\code{converged}} Did the residual submodel IRLS converge?

    \item{\code{iters}} Number of residual model IRLS iterations (1 to 50).

    \item{\code{residmod}} The residual submodel.
    So for example, if \code{varmod.method="lm"}, this will be an \code{lm} object.

    \item{\code{min.sd}}
    The predicted residual standard deviation is clamped
    so it will always be at least this value.
    This prevents prediction of negative or absurdly small variances.
    See \code{earth}'s \code{varmod.clamp} argument.
    Clamping takes place in \code{predict.varmod}, which is called
    by \code{predict.earth} when estimating prediction intervals.

    \item{\code{model.var}}
    An n x 1 matrix.
    The \code{model.var} for an observation is the estimated model
    variance for that observation over all datasets, and is estimated with
    repeated cross validation.
    It is the variance of the mean out-of-fold prediction for that
    observation over \code{ncross} repetitions.

    \item{\code{abs.resids}}
    An n x 1 matrix.
    The absolute residuals used to build the residual model.

    \item{\code{parent.x}} {
    An n x p matrix. Parent earth model \code{x}.
    }
    \item{\code{parent.y}} {
    An n x 1 matrix. Parent earth model \code{y}.
    }
    \item{\code{iter.rsq}} {
    Weighted R-Squared of residual submodel \code{residmod},
    after IRLS iteration.
    }
    \item{\code{iter.stderr}} {
    Standard errors of the coefficients of the residual submodel \code{residmod},
    after IRLS iteration.
    }
  }
}
\examples{
data(ozone1)

set.seed(1) # optional, for cross validation reproducibility

# note: should really use ncross=30 below but for a quick demo we don't

earth.mod <- earth(O3~temp, data=ozone1, nfold=10, ncross=3, varmod.method="lm")

print(summary(earth.mod)) # note additional info on the variance model

old.mfrow <- par(mfrow=c(2,2), mar=c(3, 3, 3, 1), mgp=c(1.5, 0.5, 0))

plotmo(earth.mod, do.par=FALSE, response.col=1, level=.90, main="earth model: O3~temp")

plot(earth.mod, which=3, level=.90) # residual plot: note 90% pred and darker conf intervals

par(par=old.mfrow)
}
\seealso{
  \code{\link{plot.varmod}},
  \code{\link{predict.varmod}}
}
