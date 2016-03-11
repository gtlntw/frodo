##load project related functions
source('function.r')

##setup function
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
print(seed)
print(n_rep)
print(f)
print(r) #relative risk
print(v) #variance of relative risk
print(n_family) #number of family
print(p_dis) #prevelance
print(delta) #delta to check Lyapunov condition
## ... your simualtion code
sim_result <- sim_family_new(m=r, var=v, f=f, SRR=5, n_family=n_family, p_dis=p_dis, rep=n_rep, delta=delta, structure="three2nd")

## Write out your results to a csv file
write.csv(data.frame(seed=seed, f=f, r=r, v=v, p_dis=p_dis, sum=sim_result[1,], sd=sim_result[2,], T=sim_result[3,], 
                     p_value=sim_result[4,], Lyapunov=sim_result[5,]),
          paste("res_",f,"_",r,"_",v,"_",n_family,"_",p_dis,"_",delta,"_",seed,".csv",sep=""), row.names=FALSE)
## R.miSniam
## End(Not run)
## Not run: