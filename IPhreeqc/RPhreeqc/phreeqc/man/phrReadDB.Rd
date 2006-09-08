\name{phrReadDB}
\alias{phrReadDB}
\title{
    Read phreeqc database file
}
\description{

    phrReadDB reads a phreeqc database file into ...
}
\synopsis{
    phrReadDB(filename)
}
\usage{
    phrReadDB(filename)
}
\arguments{
\item{filename}{
    The file name of the database.
} 
} %END arguments
\details{

    Any details about the operation of this function should go here.

} %END details
\value{
    This function returns NULL.
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
    \code{\link{phrGetSelectedOutput}}
}

\examples{

    # This will only work if phreeqc.dat is in the current directory
    #phrReadDB("phreeqc.dat")
}
\keyword{interface}
