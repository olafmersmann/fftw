\name{planFFT}
\alias{planFFT}
\alias{planDCT}
\title{Create FFTW plan}
\description{see title}
\usage{
planFFT(n, effort=0)
planDCT(n, type=1, effort=0)
}
\arguments{
  \item{n}{size of transform}
  \item{type}{type of DCT}
  \item{effort}{how hard fftw tries to find an optimal plan (0 to 3)}
}
\author{Olaf Mersmann <olafm@statistik.uni-dortmund.de>}
\seealso{\code{\link{FFT}} and \code{\link{IFFT}}}
