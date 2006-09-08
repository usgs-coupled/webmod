\name{phrGetLastErrorString}
\alias{phrGetLastErrorString}
\title{
    Return error string
}
\description{

    phrGetLastErrorString executes a phreeqc run
}
\synopsis{
    phrGetLastErrorString()
}
\usage{
    phrGetLastErrorString()
}
\arguments{
} %END arguments
\details{

    Any details about the operation of this function should go here.
} %END details
\value{
    This function returns the errors that occured during the previous
    phrReadDB, phrRun, phrRunFile call.
}
\references{
    Literature references and web URLs can go here.
}
\author{
    Scott Charlton <charlton@usgs.gov>

    Lopaka(Rob) Lee <rclee@usgs.gov>
} 
\seealso{
    \code{\link{phrReadString}},
    \code{\link{phrRun}},
    \code{\link{phrRunFile}},
    \code{\link{phrGetSelectedOutput}}
}

\examples{

    #phrGetLastErrorString()
}
\keyword{interface}
