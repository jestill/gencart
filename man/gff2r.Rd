\name{gff2r}
\alias{gff2r}
\title{Import GFF File to R Format}
\description{
  Import GFF formatted text files into the base R format used by
  the GenCart package.
}
\usage{gff2r(file,header=TRUE)}
\arguments{
  \item{file}{The GFF formatted file to import}
  \item{header}{Boolean, TRUE if the import file contains a header line.}
}
\references{
  \url{http://www.sanger.ac.uk/Software/formats/GFF/}
  \url{http://www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml}
}
\seealso{
  nothing at the moment
}
\examples{
## Import a gff file
gff2r("my_gff_file.gff", header=FALSE)
}
\keyword{gff}