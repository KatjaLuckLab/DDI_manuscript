# Set the working directory to the specified path
setwd('~/Documents/IMB/Johanna/')

# Source a custom R script named 'funs.R' from the 'model/scripts/' directory
source('model/scripts/funs.R')

# Define a function 'f' that takes a variable 'avar' and a vector 'x' as arguments
f <- function(avar, x) {
  # Apply a function to each element 'y' in vector 'x'
  sapply(x, function(y) mean(avar > y))
}

# Define a function 'g' that takes a vector 'x' as an argument
g <- function(x) {
  # Call the function 'f' with 'whole$X3did.z.score' as 'avar' and vector 'x'
  f(whole$X3did.z.score, x)
}

# Vectorize the function 'g' to make it applicable to vectors
h <- Vectorize(g)

{
# Create a PDF file for the plot
pdf('model/figs/ecdfish.pdf', width = 2, height = 2)
  
# Set plotting parameters for the graph
#par(cex = 1.5)
par(mar = c(3, 3, 1, 0), mgp = c(2, .5, 0), cex = .9)
  
# Create a curve by applying the vectorized function 'h' over the range from -3 to 12

curve(h, -3, 12, bty = 'l', xlab = '3did Z score (mean)', ylab = 'Proportion of DDIs', las=1)

dev.off()
}

# Remaining proportion of DDI types after applying thresholds
round(1-h(thrs),2)

