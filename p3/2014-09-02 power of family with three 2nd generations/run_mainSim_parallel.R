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


##initialize
seed=1000
#number of replications
n_rep=250
#allele frquency
#f=0.20
#relative risk
# r
#variance 
v=0
#number of family
n_family=400
#delta to check Lyapunov condition
delta=1

print(getwd())

##command line
parallel <- function(seed, n_rep, f, r, v, n_family, p_dis, delta) {
  run <- function(run=1, rfile="mainSim.R", dir=paste("/net/frodo", substr(getwd(), 0, 500), sep=""), seed, n_rep, f, r, v, n_family, p_dis, delta) {
    cmd <- paste("nohup mosbatch -E", dir, " -b -q R --vanilla --args seed ", seed, " n_rep ", n_rep," f ", f , " r ", r, " v ", v, " n_family ", n_family, " delta ", delta, " p_dis ", p_dis, 
                 " < ", rfile, " > mainSim_", f, "_", r, "_", v, "_", n_family, "_", p_dis, "_", delta, ".Rout", seed, " &", sep="")
    print(cmd)
    if(run==1) system(cmd)
    else print("not run")
  }
  
  ##run command
  for(i in 1:4) {
    Sys.sleep(10)
    run(run=1, rfile="mainSim.R", seed=seed, n_rep=n_rep, f=f, r=r, v=v, n_family=n_family, p_dis=p_dis, delta=delta)
    #add seed number
    seed=seed+1  
  }
}


for(i in seq(1,2, length=11)) {
  parallel(seed=seed, n_rep=n_rep, f=0.01, r=i, v=v, n_family=n_family, p_dis=0.1, delta=delta)
}


  
## Choose a high enough seed for later for pasting the results together
## (1,...,9,10) sorts not the way you want, for example.
## R.lellarap_miSniam_nur
## End(Not run)
## Not run: