\name{etitanic}
\alias{etitanic}
\title{Titanic data with incomplete cases removed}
\description{
Titanic data with incomplete cases,
passenger names, and other details removed.
}
\format{
  A data frame with 1046 observations on 6 variables.

  \tabular{ll}{
      \code{pclass}   \tab passenger class, unordered factor: 1st 2nd 3rd\cr
      \code{survived} \tab integer: 0 or 1\cr
      \code{sex}      \tab unordered factor: male female\cr
      \code{age}      \tab age in years, min 0.167 max 80.0\cr
      \code{sibsp}    \tab number of siblings or spouses aboard, integer: 0...8\cr
      \code{parch}    \tab number of parents or children aboard, integer: 0...6\cr
  }
}
\source{
This dataset is included in the earth package because it is a convenient
vehicle for illustrating earth's GLM and factor handling.

The dataset was compiled by Frank Harrell and Robert Dawson:
\url{http://biostat.mc.vanderbilt.edu/twiki/pub/Main/DataSets/titanic.html}\cr
See also:\cr
\url{http://biostat.mc.vanderbilt.edu/twiki/pub/Main/DataSets/titanic3info.txt}.

For this version of the Titanic data,
passenger details and incomplete cases were deleted
and the name changed to \code{etitanic} to minimize confusion with other
versions ("e" because it is part of the earth package).

Note that \code{survived} is an integer (it should arguably be a logical).

In this data the crew are conspicuous by their absence.

Contents of \code{etitanic}:
    \preformatted{
         pclass survived    sex    age sibsp parch
    1       1st        1 female 29.000     0     0
    2       1st        1   male  0.917     1     2
    3       1st        0 female  2.000     1     2
    4       1st        0   male 30.000     1     2
    5       1st        0 female 25.000     1     2
    ...
    1309    3rd        0   male 29.000     0     0
    }
How \code{etitanic} was built:
    \preformatted{
    load("titanic3") # from Harrell's web site
    # discard name, ticket, fare, cabin, embarked, body, home.dest
    etitanic <- titanic3[,c(1,2,4,5,6,7)]
    etitanic <- etitanic[!is.na(etitanic$age),]
    save(etitanic, file="etitanic.rda")
    }
}
\references{
Further details and analyses of the Titanic data may be found in:

F. Harrell (2001) \emph{
Regression Modeling Strategies with Applications to Linear Models, Logistic Regression, and Survival Analysis }
\url{http://biostat.mc.vanderbilt.edu/twiki/bin/view/Main/RmS}
}
\seealso{
  \code{\link{earth}}
}
\keyword{datasets}
