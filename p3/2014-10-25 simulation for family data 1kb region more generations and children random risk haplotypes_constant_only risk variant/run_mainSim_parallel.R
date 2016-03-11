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
n_rep=400
#number of family
print(n_family) #=500
#null or alternative hypothesis
print(null) #null=FALSE
#family structure
print(family_strct) #family_strct='family_strct.2g2c'

print(getwd())

##command line
parallel <- function(seed, n_rep, f, n_family, p_dis, null, family_strct) {
  ##run command
  for(i in 1:5) {
    rfile="mainSim.R"
    dir=paste("/net/frodo", substr(getwd(), 0, 500), sep="")
    cmd <- paste("nohup mosbatch -E", dir, " -b -q R --vanilla --args seed ", seed, " n_rep ", n_rep," f ", f , " n_family ", n_family, " p_dis ", p_dis, " null ", null, " family_strct ", family_strct,
                 " < ", rfile, " > mainSim_", f, "_", n_family, "_", p_dis, ".Rout", seed, " &", sep="")
    print(cmd)
    system(cmd)
    
    #add seed number
    seed=seed+1
    Sys.sleep(10)
  }
}
parallel(seed=seed, n_rep=n_rep, f=0.01, n_family=n_family, p_dis=0.05, null=null, family_strct=family_strct)


  
## Choose a high enough seed for later for pasting the results together
## (1,...,9,10) sorts not the way you want, for example.
## R.lellarap_miSniam_nur
## End(Not run)
## Not run: