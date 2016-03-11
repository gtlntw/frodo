## Not run:
## paste_mainSim_results.R
## Put the following code in paste_mainSim_results.R (or just
## type them in), and run
## R --vanilla < paste_mainSim_results.R
##
## Pastes all of the csv files created in 'run_mainSim_parallel'
## together.
mergeCsv <- function (every = 1, outfile = "allResults.csv", quote = FALSE) 
{
  csvFiles <- setdiff(sort(dir(pattern = ".csv")), outfile)
  if (length(csvFiles) == 0) 
    stop("pasteCsv: No csv file in the current directory.")
  print(csvFiles)
  res <- read.csv(csvFiles[1], stringsAsFactors = FALSE)
  for (i in 2:length(csvFiles)) {
    temp <- read.csv(csvFiles[i])
    if (ncol(temp) != ncol(res) || !all(names(temp) == names(res))) {
      cat(paste("Names of csv file (", paste(names(temp), 
                                             collapse = ","), ") do not match names of current data (", 
                paste(names(res)), ").", sep = ""), "\n")
      cat("Current dataset:\n")
      print(res)
      cat("Dataset trying to merge in:\n")
      print(temp)
      if (ncol(temp) != ncol(res)) {
        stop("mergeCsv Error, dimensions not compatible")
      }
      else {
        names(temp) <- names(res)
        cat("warning, names of new dataset do not match, but were overwritten")
      }
    }
    res <- rbind(res, temp)
  }
  wh <- which(names(res) == "seed")
  if (length(wh) == 1) 
    res <- res[order(res$seed), ]
  if (every > 1) {
    i <- 1
    while (i < nrow(res)) {
      rows <- i + 1:every - 1
      for (col in 1:ncol(res)) {
        if (is.numeric(res[rows, col])) {
          res[i, col] <- mean(res[rows, col])
        }
      }
      res <- res[-rows[-1], ]
      i <- i + 1
    }
  }
  if (nchar(outfile) > 0) 
    write.csv(res, outfile, quote = quote, row.names = FALSE)
  return(res)
}
mergeCsv()
## End(Not run