nYears <- 48
nMonths <- 12
error <- 0.05
foi <- runif(nYears, 0.1,0.4)
x <- seq(0,nMonths-1,by=1)/nMonths
knots <- c(0.33,0.66)
degree <- 2
nKnots <- length(knots) + degree + 1
theta <- runif(nKnots)
allDat <- NULL
for(i in 1:nYears){
  dat <- genSpline(x, knots, degree, theta)$dt
  dat$y.spline <- (dat$y.spline*foi[i])/sum(dat$y.spline)
  dat$x <- dat$x*12 + (i-1)*12
  allDat <- rbind(allDat, dat)
}
allDat$y.spline <- allDat$y.spline #+ rnorm(length(allDat$y),0,error)
plot(allDat,type='l')


#plot.spline(dat)

n_indiv <- 1000
nYears <- 48
infHist <- matrix(sample(c(0,1),n_indiv*nYears*12,prob=c(0.98,0.02),replace=TRUE),nrow=n_indiv)
ageMask <- rep(1,n_indiv)
library(microbenchmark)
Rprof(tmp <- tempfile())
## Run the MCMC using the inputs generated above
for(i in 1:100) calc_lambda_probs_monthly(foi, knots, theta, infHist, ageMask)

Rprof()
summaryRprof(tmp)
library(proftools)
plotProfileCallGraph(readProfileData(tmp),score = "total")


microbenchmark(bs(x,knots,3))
microbenchmark(bSpline(x,knots=knots,degree=3))
microbenchmark(calc_lambda_probs_monthly(foi, knots, theta, infHist, ageMask))
