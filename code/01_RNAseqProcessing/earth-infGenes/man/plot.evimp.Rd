\name{plot.evimp}
\alias{plot.evimp}
\title{Plot an evimp object (created by the evimp function)}
\description{
    Plot an \code{\link{evimp}} object.
}
\usage{
\method{plot}{evimp}(x = stop("no 'x' argument"),
    cex.var = 1,
    type.nsubsets = "l", col.nsubsets = "black", lty.nsubsets = 1,
    type.gcv = "l", col.gcv = 2, lty.gcv = 1,
    type.rss = "l", col.rss = "gray60", lty.rss = 1,
    cex.legend = 1, x.legend = nrow(x), y.legend = x[1,"nsubsets"],
    rh.col = 1, do.par = TRUE, \dots)
}
\arguments{
  \item{x}{
     An \code{\link{evimp}} object.
  }
  \item{cex.var}{
     cex for variable names.
     Default is 1.
     Make smaller (say 0.8) if you have lots of variables.
  }
  \item{type.nsubsets}{
     Plot type for nsubsets graph.
     Default is "l".
     Use "n" for none, "b" looks good too.
  }
  \item{col.nsubsets}{
     Color of nsubsets line.
     Default is "black".
  }
  \item{lty.nsubsets}{
     Line type of nsubsets line.
     Default is 1.
  }
  \item{type.gcv,col.gcv,lty.gcv}{
     As above but for the gcv plot
  }
  \item{type.rss,col.rss,lty.rss}{
     As above but for the rss plot
  }
  \item{cex.legend}{
     cex for legend strings.
     Default is 1.
     Make smaller (say 0.8) if you want a smaller legend.
  }
  \item{x.legend}{
      x position of legend.
      Use 0 for no legend.
  }
  \item{y.legend}{
      y position of legend.
  }
  \item{rh.col}{
    Color of right hand axis label.  Use \code{rh.col=0} for no label,
    a workaround for when the label is mispositioned.
  }
  \item{do.par}{
    Call \code{par()} for global settings as appropriate.
    Default is \code{TRUE}, which sets\cr
    \code{oma=c(bottom.margin,0,0,3), cex=cex.var}.\cr
    Set to \code{FALSE} if you want to append figures to an existing plot.
  }
  \item{\dots}{
       Extra arguments passed to plotting functions.
  }
}
\seealso{
  \code{\link{earth}},
  \code{\link{evimp}},
  \code{\link{plot.earth.models}},
  \code{\link[plotmo]{plotmo}}
}
\examples{
data(ozone1)
earth.mod <- earth(O3 ~ ., data=ozone1, degree=2)
ev <- evimp(earth.mod)
plot(ev)
print(ev)
}
\keyword{models}
