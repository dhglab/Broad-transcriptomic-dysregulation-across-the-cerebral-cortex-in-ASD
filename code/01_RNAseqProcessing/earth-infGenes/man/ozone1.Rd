\name{ozone1}
\alias{ozone1}
\title{Ozone readings in Los Angeles with incomplete cases removed}
\description{
Ozone readings in Los Angeles, with incomplete cases removed.}
\format{
  A data frame with 330 observations on 10 variables.

  \tabular{ll}{
      \code{O3} \tab daily maximum of the hourly average ozone concentrations in Upland, CA\cr
      \code{vh} \tab 500 millibar pressure height, measured at the Vandenberg air force base\cr
      \code{wind} \tab wind speed in mph at LAX airport\cr
      \code{humidity} \tab humidity in percent at LAX\cr
      \code{temp} \tab Sandburg Air Force Base temperature in degrees Fahrenheit\cr
      \code{ibh} \tab temperature inversion base height in feet\cr
      \code{dpg} \tab pressure gradient from LAX to Daggert in mm Hg\cr
      \code{ibt} \tab inversion base temperature at LAX in degrees Fahrenheit\cr
      \code{vis} \tab visibility at LAX in miles\cr
      \code{doy} \tab day of the year\cr
  }
}
\source{
This dataset was copied from \code{library(faraway)} and the name
changed to \code{ozone1} to prevent a name clash.
The data were originally made available by Leo Breiman who was a consultant
on a project where the data were generated.
Example analyses using these data may be found in Faraway and in Hastie and Tibshirani.

\preformatted{    > ozone1
        O3   vh wind humidity temp  ibh dpg ibt vis doy
    1    3 5710    4       28   40 2693 -25  87 250  33
    2    5 5700    3       37   45  590 -24 128 100  34
    3    5 5760    3       51   54 1450  25 139  60  35
    ...
    330  1 5550    4       85   39 5000   8  44 100 390
}
}
\references{
Faraway (2005) \emph{Extending the Linear Model with R}
\url{http://www.maths.bath.ac.uk/~jjf23}

Hastie and Tibshirani (1990) \emph{Generalized Additive Models}
\url{http://web.stanford.edu/~hastie/pub.htm}
}
\seealso{
  \code{\link{earth}}

  \code{\link[datasets]{airquality}} a different set of ozone data
}
\keyword{datasets}
