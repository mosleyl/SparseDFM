linspace <- function(x1, x2, n=100) {
  stopifnot(is.numeric(x1), is.numeric(x2), length(x1)==1, length(x2)==1)
  n <- floor(n)
  if (n <= 1) x2
  else seq(x1, x2, length.out=n)
}

logspace <- function(x1, x2, n=50) {
  if (x2 == pi) x2 <- log10(x2)
  10^linspace(x1, x2, n)
}



