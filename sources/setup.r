###
# Given in the project statement
###

n=100

X = function(n=100) runif(n)
m = function(x) (sin(2*pi*x^3))^3
Y = function(X=runif(n)) m(X) + rnorm(length(X), 0, 0.5) # transformed into a function