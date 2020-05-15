\name{mars.to.earth}
\alias{mars.to.earth}
\title{Convert a mars object from the mda package to an earth object}
\description{
Convert a \code{\link[mda]{mars}} object from the \code{mda} package to an \code{\link{earth}} object
}
\usage{
mars.to.earth(object, trace=TRUE)
}
\arguments{
  \item{object}{
    A \code{mars} object, created using \code{\link[mda]{mars}} in the \code{mda} package.
  }
  \item{trace}{
    If \code{TRUE} (default) print a summary of the conversion.
  }
}
\value{
  The value is the same format as that returned by \code{\link{earth}} but
  with skeletal versions of \code{rss.per.subset},
  \code{gcv.per.subset}, and \code{prune.terms}.\cr

  You can fully initialize these components by calling \code{\link{update.earth}}
  after \code{mars.to.earth}, but if you do this \code{selected.terms} may change.
  However with \code{pmethod="backward"} a change is unlikely ---
  \code{selected.terms} would change only if GCVs are so close that numerical errors
  have an effect.
}
\note{
\bold{Differences between mars and earth objects}

  Perhaps the most notable difference between
  \code{mars} and \code{earth} objects is that \code{mars} returns the
  MARS basis matrix in a field called "\code{x}"
  whereas \code{earth} returns "\code{bx}" with only the selected terms.
  Also, \code{earth} returns "\code{dirs}" rather than "\code{factors}",
  and in \code{earth} this matrix can have entries of value 2 for linear predictors.
  % The calculation of \code{minspan} in \code{earth} follows Friedman's paper more closely
  % and is slightly different from \code{mars}.

  For details of other differences between \code{mars} and \code{earth} objects,
  see the comments in the source code of \code{mars.to.earth}.

\bold{Weights}

  The \code{w} argument is silently ignored by \code{mars}.
  % The equivalent \code{earth} argument \code{weights} is also not yet supported,
  % and you will get a warning.

  \code{mars} normalizes \code{wp} to (euclidean) length 1;
  \code{earth} normalizes
  \code{wp} to length equal to the number of responses, i.e., the number
  of columns in \code{y}.  This change was made so an all ones \code{wp}
  (or in fact any all constant \code{wp}) is equivalent to using no \code{wp}.

  If the original call to \code{mars} used the \code{wp} argument,
  \code{mars.to.earth} will run \code{\link{update.earth}} to force consistency.
  This could modify the model, so a warning is issued.
}
\seealso{
  \code{\link{earth}},
  \code{\link[mda]{mars}}
}
\examples{
if(require(mda)) {
    mars.mod <- mars(trees[,-3], trees[,3])
    earth.mod <- mars.to.earth(mars.mod)
    # the standard earth functions can now be used
    # note the reconstructed call in the summary
    summary(earth.mod, digits = 2)
}
}
\keyword{models}
