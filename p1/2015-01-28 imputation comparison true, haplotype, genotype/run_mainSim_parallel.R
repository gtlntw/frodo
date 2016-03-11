## Selects a variety of parameter combinations to run
## mainSim.R in parallel on a cluster.
##
## Then run the commands with
## R --vanilla < run_mainSim_parallel.R >! run_mainSim_parallel.txt

parseCommandArgs <- function (evaluate = TRUE) 
{
  arglist <- list()
  args <- commandArgs()
  i <- which(args == "--args")
  if (length(i) == 0 || length(args) < 1) 
    return(invisible())
  args <- args[(i + 1):length(args)]
  for (i in seq(1, length(args), by = 2)) {
    value <- NA
    tryCatch(value <- as.double(args[i + 1]), warning = function(e) {
    })
    if (is.na(value)) {
      value <- args[i + 1]
      if (substr(value, 1, 2) == "c(") 
        value <- eval(parse(text = args[i + 1]))
    }
    if (evaluate) 
    assign(args[i], value, inherits = TRUE)
    arglist[[length(arglist) + 1]] <- value
    names(arglist)[length(arglist)] <- args[i]
  }
  return(arglist)
}

#need to take in parameter
## passed in from the command prompt.
parseCommandArgs()
## Will disply the default values on the first run,
##initialize
seed=1000
#number of replications
n_rep=200
#relative risk
#r

print(getwd())

##command line
parallel <- function(...) {
  names.args <- names(list(...))
  value.args <- as.vector(list(...))
  ##run command
  for(i in 1:5) {
    rfile="mainSim.R"
    dir=paste("/net/frodo", substr(getwd(), 0, 500), sep="")
    cmd <- paste("nohup mosbatch -E", dir, " -j2,4-8 R --vanilla --args seed ", seed, " ", paste(names.args, value.args, collapse=" "),
                 " < ", rfile, " > mainSim_", paste(value.args, collapse="_"), ".Rout", seed, " &", sep="")
    print(cmd)
    system(cmd)
    
    #add seed number
    seed<<-seed+1
    Sys.sleep(5)
  }
}
system("rm *.csv")
system("rm *.Rout*")
for(i in seq(1,2, by=0.2)) {
  parallel(n_rep=n_rep, r=i)
}

# parallel(n_rep=n_rep, r=1)


  
## Choose a high enough seed for later for pasting the results together
## (1,...,9,10) sorts not the way you want, for example.
## R.lellarap_miSniam_nur
## End(Not run)
## Not run: