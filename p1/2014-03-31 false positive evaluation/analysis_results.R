##
## Then run the commands with
## R --vanilla < analysis_results.R >! analysis_results.txt

#######analyze the result
result <- read.csv("allResults.csv", header=T, stringsAsFactors=F)

hist(result$p)
mean(result$p < 0.05)
hist(result$t, freq=F)
curve(dnorm, add=TRUE, col=2)
mean(result$t)

