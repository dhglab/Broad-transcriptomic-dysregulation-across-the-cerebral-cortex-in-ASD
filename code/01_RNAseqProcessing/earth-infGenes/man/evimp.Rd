\name{evimp}
\alias{evimp}
\title{Estimate variable importances in an earth object}
\description{
Estimate variable importances in an \code{\link{earth}} object
}
\usage{
evimp(object, trim=TRUE, sqrt.=TRUE)
}
\arguments{
  \item{object}{
    An \code{\link{earth}} object.
  }
  \item{trim}{
    If \code{TRUE} (default), delete rows in the returned matrix for
    variables that don't appear in any subsets.
  }
  \item{sqrt.}{
    Default is \code{TRUE},
    meaning take the \code{\link{sqrt}} of the GCV and RSS importances before
    normalizing to 0 to 100.
    Taking the square root gives a better indication of
    relative importances because the raw importances are calculated using
    a sum of squares.
    Use \code{FALSE} to not take the square root.
  }
}
\value{
This function returns a matrix showing the relative importances of the
variables in the model.  There is a row for each variable.  The row
name is the variable name, but with \code{-unused} appended if the
variable does not appear in the final model.

The columns of the matrix are (not all of these are printed by \code{print.evimp}):
\itemize{
\item \code{col}: Column index of the variable in the \code{x} argument to \code{earth}.
\item \code{used}: 1 if the variable is used in the final model, else 0.
Equivalently, 0 if the row name has an \code{-unused} suffix.
\item \code{nsubsets}: Variable importance using the "number of subsets" criterion.
Is the number of subsets that include the variable (see "Three Criteria" in the chapter
on \code{evimp} in the \code{earth} vignette
\dQuote{\href{../doc/earth-notes.pdf}{Notes on the earth package}}).
\item \code{gcv}: Variable importance using the GCV criterion (see "Three Criteria").
\item \code{gcv.match}: 1, except is
0 where the rank using the \code{gcv} criterion differs from
that using the \code{nsubsets} criterion.
In other words, there is a 0 for values that increase as you go
down the \code{gcv} column.
\item \code{rss}: Variable importance using the RSS criterion (see "Three Criteria").
\item \code{rss.match}: Like \code{gcv.match} but for the \code{rss}.
}
The rows are sorted on the \code{nsubsets} criterion.
This means that values in the \code{nsubsets} column decrease as you go down the column
(more accurately, they are non-increasing).
The values in the \code{gcv} and \code{rss} columns
are also non-increasing, except where the
\code{gcv} or \code{rss} rank differs from the \code{nsubsets} ranking.
}
\note{
There is a chapter on \code{evimp} in the \code{earth} package vignette
\dQuote{\href{../doc/earth-notes.pdf}{Notes on the earth package}}.

\bold{Acknowledgment}

Thanks to Max Kuhn for the original \code{evimp} code and for helpful discussions.
}
\seealso{
  \code{\link{earth}},
  \code{\link{plot.evimp}}
}
\examples{
data(ozone1)
earth.mod <- earth(O3 ~ ., data=ozone1, degree=2)
ev <- evimp(earth.mod, trim=FALSE)
plot(ev)
print(ev)
}
\keyword{models}
