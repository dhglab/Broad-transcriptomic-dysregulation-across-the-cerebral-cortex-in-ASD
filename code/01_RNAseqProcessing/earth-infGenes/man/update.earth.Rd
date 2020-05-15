\name{update.earth}
\alias{update.earth}
\title{Update an earth model}
\description{
Update an \code{\link{earth}} model.
}
\usage{
\method{update}{earth}(object = stop("no 'object' argument"),
       formula. = NULL, ponly = FALSE, \dots, evaluate = TRUE)
}
\arguments{
  \item{object}{ The earth object}
  \item{formula.}{
     The \code{formula.} argument is treated like earth's \code{formula} argument.\cr
  }
  \item{ponly}{
     Force pruning only, no forward pass.
     Default is \code{FALSE},
     meaning \code{update.earth} decides automatically if a forward pass is needed.
     See note below.
  }
  \item{\dots}{
     Arguments passed on to \code{\link{earth}}.
  }
  \item{evaluate}{
     If \code{TRUE} (default) evaluate the new call, else return the call.
     Mostly for compatibility with the generic \code{\link{update}}.
  }
}
\details{
If only the following arguments are used, a forward pass
is unnecessary, and \code{update.earth} will perform only the pruning pass.
This is usually much faster for large models.
\preformatted{
     object
     glm
     trace
     nprune
     pmethod
     Eval.model.subsets
     Print.pruning.pass
     Force.xtx.prune
     Use.beta.cache
     Endspan.penalty
     Get.leverages
}
This automatic determination to do a forward pass can be overridden
with the \code{ponly} argument.
If \code{ponly=TRUE} the forward pass will be skipped and only the pruning pass will be executed.
This is useful for doing a pruning pass with new data.
(Use earth's \code{data} argument to specify the new data.)
Typically in this scenario you would also specify \code{penalty=-1}.
This is because with sufficient new data, independent of the original
training data, the RSS not the GCV should be used for evaluating model
subsets
(The GCV approximates what the RSS would be on new data --- but here
we actually have new data, so why bother approximating.
This "use new data for pruning" approach is useful in situations where
you don't trust the GCV approximation for your data.)
By making \code{penalty=-1}, earth will calculate the RSS, not the GCV.
See also the description of \code{penalty} on the
\code{\link{earth}} help page.
\cr
\cr
Another (somewhat esoteric) use of \code{ponly=TRUE} is to do subset
selection with a different \code{penalty} from that used to build the
original model.
\cr
\cr
With \code{trace=1}, \code{update.earth} will tell you if earth's
forward pass was skipped.
\cr
\cr
If you used \code{keepxy=TRUE} in your original call to \code{earth}, then
\code{update.earth} will use the saved values of \code{x}, \code{y}, etc.,
unless you specify otherwise by arguments to \code{update.earth}.
It can be helpful to set \code{trace=1} to see which \code{x} and \code{y}
is used by \code{update.earth}.
}
\value{
   The value is the same as that returned by \code{\link{earth}}.
   If \code{object} is the only parameter then no changes are made
   --- the returned value will be the same as the original \code{object}.
}
\seealso{
  \code{\link{earth}}
}
\examples{
data(ozone1)

(earth.mod <- earth(O3 ~ ., data = ozone1, degree = 2))

update(earth.mod, formula = O3 ~ . - temp) # requires forward pass and pruning

update(earth.mod, nprune = 8)              # requires only pruning

update(earth.mod, penalty=1, ponly=TRUE)   # pruning pass only with a new penalty
}
\keyword{models}
