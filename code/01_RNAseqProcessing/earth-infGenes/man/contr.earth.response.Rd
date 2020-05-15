\name{contr.earth.response}
\alias{contr.earth.response}
\title{Please ignore}
\description{
Contrasts function for factors in the \code{\link{earth}} response.
For internal use by earth.
}
\usage{
contr.earth.response(x, base, contrasts)
}
\arguments{
  \item{x}{
    a factor
  }
  \item{base}{
    unused
  }
  \item{contrasts}{
    unused
  }
}
\value{
Returns a diagonal matrix.
An example for a 3 level factor with levels \code{A}, \code{B}, and \code{C}:
\preformatted{
      A B C
    A 1 0 0
    B 0 1 0
    C 0 0 1
}
}
\note{
Earth uses this function internally.
You shouldn't need it.
It is made publicly available only because it seems that is necessary
for \code{model.matrix}.
}
\seealso{
  \code{\link{contrasts}}
}
\keyword{models}
