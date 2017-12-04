##
## fftw.R - glue between user and fftw.c
##
## Authors:
##  Sebastian Krey <skrey@statistik.tu-dortmund.de>
##  Olaf Mersmann  <olafm@statistik.tu-dortmund.de>
##

planFFT <- function(n, effort=0) {
  if (length(n) == 1)
    n <- as.integer(n)
  plan <- .Call("FFT_plan", n, as.integer(effort), PACKAGE="fftw")
  class(plan) <- "FFTplan"
  return (plan)
}

print.FFTplan <- function(x, ...)
  invisible(.Call("FFT_print_plan", x, PACKAGE="fftw"))

FFT <- function(x, ..., plan, inverse=FALSE) {
  if (missing(plan))
    plan <- planFFT(x)
  if (is.integer(x))
    x <- as.numeric(x)
  .Call("FFT_execute", plan, x, inverse, PACKAGE="fftw")
}

IFFT <- function(x, ..., plan, scale=TRUE) {
  if (missing(plan))
    plan <- planFFT(x)
  if (is.integer(x))
    x <- as.numeric(x)
  y <- .Call("FFT_execute", plan, x, TRUE, PACKAGE="fftw")
  if (scale)
    y <- y / length(y)
  return(y)
}

## DCT:
planDCT <- function(n, type=1, effort=0) {
  if (length(n) == 1)
    n <- as.integer(n)
  plan <- .Call("DCT_plan", n, as.integer(type), as.integer(effort), PACKAGE="fftw")
  class(plan) <- "DCTplan"
  return (plan)  
}

DCT <- function(x, ..., plan, type=1, inverse=FALSE) {
  if (missing(plan))
    plan <- planDCT(x, type)
  if (is.integer(x))
    x <- as.numeric(x)
  .Call("DCT_execute", plan, x, inverse, PACKAGE="fftw")
}

IDCT <- function(x, ..., plan, type=1, scale=TRUE) {
  if (missing(plan))
    plan <- planDCT(x, type)
  if (is.integer(x))
    x <- as.numeric(x)
  y <- .Call("DCT_execute", plan, x, TRUE, PACKAGE="fftw")  
  if (scale) {
    n <- length(y)
    N <- if (type == 1) 2*(n-1) else 2*n
    y <- y / N
  }
  return(y)
}
