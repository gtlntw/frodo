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

#########################################################################
##passed in from the command prompt.
parseCommandArgs()
## Will disply the default values on the first run,
##initialize
seed=1000
#number of replications
n_rep=250
#number of family
n_family=1000
#prevalence
p_dis=0.1
# print(null) #null=FALSE
#family structure
# print(family_strct) #family_strct='family_strct.2g2c2a'

print(getwd())

##command line
parallel <- function(...) {
  names.args <- names(list(...))
  value.args <- as.vector(list(...))
  ##commands to run
  for(i in 1:6) {
    rfile="mainSim.R"
    cmd <- paste("R --vanilla --args seed ", seed, " ", paste(names.args, value.args, collapse=" "),
                 " < ", rfile, " > ", "mainSim_", paste(value.args, collapse="_"), ".Rout", seed, " 2>&1", sep="")
    print(cmd)
    cat(cmd, file=fileConn, sep="\n") #write jobs to text
    #     writeLines(cmd, fileConn) #write jobs to text
    
    #add seed number
    seed<<-seed+1
  }
}
##clean-up & initialization
system("rm *.csv")
system("rm *.Rout*")
system("rm *.out*")
system("rm *.sh")
system("rm result*.txt")
system("rm mendelian*.txt")
system("rm data*.ped")
fileConn<-file("jobs.txt", "w")

# ##create your jobs here
# for(i in seq(1,1.5, length.out = 11)) {
#   parallel(n_rep=n_rep, r=i, n_family=n_family, p_dis=p_dis, risk.variant="\"c(4,16,42)\"") #common
# }
# for(i in seq(1,2.1, length.out = 11)) {
#   parallel(n_rep=n_rep, r=i, n_family=n_family, p_dis=p_dis, risk.variant="\"c(4,15,32)\"") #rare
# }
for(i in seq(1,2.8, length.out = 11)) {
  parallel(n_rep=n_rep, r=i, n_family=n_family, p_dis=p_dis, risk.variant=39) #super rare
}

##write jobs.txt
close(fileConn)
##submit the jobs.txt using runslurm.pl
system("runslurm.pl -replace -logfile jobs.log -copts \"--time=24:00:00 --mem=1024\" jobs.txt")
##submit jobs.txt via runbatch.pl i.e. mosic mode
# dir=paste("/net/frodo", substr(getwd(), 0, 500), sep="")
# cmd <- paste("runbatch.pl -replace -logfile jobs.log -concurrent 40 -copts \"-E", dir," -j2,4-8\" jobs.txt &", sep="")
# print(cmd)
# system(cmd)
