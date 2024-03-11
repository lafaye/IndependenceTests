\name{dna}
\alias{dna}
\docType{data}
\title{dna sequence}
\description{
This data from Whisenant et al. (1991) is a nucleotides sequence of 4156
base pairs (bp). The categorical variable represents the
nucleotide which is either one of the two purines (r), adenine (a) or
guanine (g), or one of the two pyrimidines (y), cytosine (c) or thymine (t).
}

\usage{
data(dna)
}

\format{
A character vector of length 70 representing 70 consecutive segments of
a dna strand of length 4156.
}

\references{
 Whisenant E.C., Rasheed B.K.A., Ostrer H., Bhatnagar 
  Y.M. (1991). Evolution and sequence analysis of a human Y-chromosomal
  DNA fragment, J. Mol. Evol., 33, 133-141.
}

\examples{
data(dna)
x2 <- dna[1]
for (i in 2:length(dna)) x2 <- paste(x2, dna[i], sep = "")
x <- unlist(strsplit(x2, ""))
}
\keyword{datasets}
