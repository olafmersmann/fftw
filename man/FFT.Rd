\name{FFT}
\alias{FFT}
\alias{IFFT}
\alias{DCT}
\alias{IDCT}
\title{Calculate (inverse) DFT using the FFT method}
\description{see title}
\usage{
FFT(x, ..., plan)
IFFT(x, ..., plan, scale=TRUE)
DCT(x, ..., plan, type=1)
IDCT(x, ..., plan, type=1, scale=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{(complex) vector to process}
  \item{\dots}{ignored}
  \item{plan}{FFTW plan, can be missing}
  \item{scale}{scale results}
  \item{type}{type of DCT}
}
\author{Olaf Mersmann <olafm@statistik.uni-dortmund.de>}
\seealso{\code{\link{planFFT}}}
\examples{
n <- 2**16
x <- rnorm(n)
p <- planFFT(n)
y <- FFT(x, plan=p)
}