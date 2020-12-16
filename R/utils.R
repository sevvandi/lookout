nsphere <- function(n=3, r=1, npts=3145, sd=1){
  set.seed(sd)
  nrows <- npts
  x <- matrix(1, nrow = nrows, ncol=n)
  for(i in 1:(n-1)){
    j <- n - i
    theta <- seq(from = -1*pi, to = pi, length.out=npts)
    theta1 <- sample(theta, size=length(theta))
    x[ ,1:j] <- x[ ,1:j]*cos(theta1)
    x[ ,(j+1)] <- x[ ,(j+1)]*sin(theta1)
  }
  x <- x*r
  return(x)
}
