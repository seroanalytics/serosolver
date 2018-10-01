library(splines)
library(data.table)
library(ggplot2)
library(broom)

genSpline <- function(x, knots, degree, theta, intercept=TRUE) {
  
  basis <- bs(x = x, knots = knots, degree = degree,
              Boundary.knots = c(0,1), intercept = intercept)
  
  y.spline <- basis %*% theta
  
  dt <- data.table(x, y.spline = as.vector(y.spline))
  
  return(list(dt = dt, basis = basis, knots = knots))
  
}
plot.basis <- function(basisdata) {
  
  dtbasis <- as.data.table(basisdata$basis)
  dtbasis[, x := seq(0, 1, length.out = .N)]
  dtmelt <- melt(data = dtbasis, id = "x", 
                 variable.name = "basis", variable.factor = TRUE)
  
  ggplot(data=dtmelt, aes(x=x, y=value, group = basis)) +
    geom_line(aes(color=basis), size = 1) +
    theme(legend.position = "none") +
    scale_x_continuous(limits = c(0, 1), 
                       breaks = c(0, basisdata$knots, 1)) +
    theme(panel.grid.minor = element_blank())
}

plot.spline <- function(basisdata, points = FALSE) {
  
  p <- ggplot(data = basisdata$dt)
  
  if (points) p <- p + geom_point(aes(x=x, y = y), color = "grey75")  
  
  p <- p + 
    geom_line(aes(x = x, y = y.spline), color = "red", size = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(limits = c(0, 1), breaks = knots) +
    theme(panel.grid.minor = element_blank())
  
  return(p)
  
}

x <- seq(0, 1, length.out = 1000)

knots <- c(0.25, 0.5, 0.75)
theta = c(0.6, 0.1, 0.3, 0.2, 0.9)

sdata <- genSpline(x, knots, 1, theta)
plot.spline(sdata)

theta = c(0.2, 0.3, 0.8, 0.2, 0.1)
sdata <- genSpline(x, knots, 1, theta)

plot.basis(sdata)
plot.spline(sdata)
knots <- c(0.25, 0.5, 0.75)
theta = c(0.6, 0.1, 0.5, 0.2, 0.8, 0.3)

sdata <- genSpline(x, knots, 2, theta)

plot.basis(sdata)
plot.spline(sdata)
knots <- c(0.333, 0.666)
theta = c(0.2, 0.4, 0.1, 0.9, 0.6)

sdata <- genSpline(x, knots, 2, theta)
plot.basis(sdata)
knots <- c(0.333, 0.666)
theta = c(0.1, 0.7, 0.1, 0.9, 0.2, 0.1)
sdata <- genSpline(x, knots, 3, theta)
plot.basis(sdata)
set.seed(5)
x <- runif(250)
sdata <- genSpline(x, knots, 3, theta)

sdata$dt[,  y := rnorm(.N, y.spline, 0.1)]

plot.spline(sdata, points = TRUE)


#knots <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
knots <- c(0.25,0.5,0.75)
degree <- 3
nKnots <- length(knots) + degree + 1
theta <- runif(nKnots)
x <- seq(0,1,by=0.01)
dat <- genSpline(x, knots, degree, theta)
plot.spline(dat)
