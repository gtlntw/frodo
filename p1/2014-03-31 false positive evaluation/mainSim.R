source('/net/frodo/home/khlin/p1/function.r')

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
## ... your simualtion code
sim <- function(rep=1000, var=0, f=0.01, n_sibpair=1000, SRR=5, r=1) {
  sim.result <- NULL
  r <<- r
  for(i in r) {
    print(paste("Generating Data...m=", i, "var=", var, "f=", f))
    mark <<- 1 #reset the current iteration
    print(mark)  #start iteration
    sim.result <- replicate(rep, {
      data.sample <- gene.sibpair(m=i, var=var, f=f, SRR=SRR, n_sibpair=n_sibpair)
      ##sibpair internal control
      attach(data.sample)
      n_chr_s <- sum(S1==1) + 2*sum(S2==1)
      n_chr_ns <- 4*sum(S0==1) + 2*sum(S1==1)
      
      case.count.mis <- sum(HS.mis)
      control.count.mis <- sum(HN.mis)
      
      #chi-square test on misclassification data
      temp <- matrix(c(case.count.mis, n_chr_s - case.count.mis, control.count.mis, n_chr_ns - control.count.mis), ncol = 2)
      test.mis <- chisq.test(temp, correct=F) ##prevent anti-conservative
      
      
      #       #chi-square test on corrected data by EM
      corr_fac <- EM(data.sample)
      #       
      count.HS1_HN0_S1 <- sum(S1==TRUE & HS.mis==1 & HN.mis==0)
      #       
      case.count.mis.cor <- corr_fac*count.HS1_HN0_S1 + (case.count.mis - count.HS1_HN0_S1)
      control.count.mis.cor <- control.count.mis + 2*(1-corr_fac)*count.HS1_HN0_S1
      #       
      #       temp <- matrix(c(case.count.mis.cor, n_chr_s - case.count.mis.cor, control.count.mis.cor, n_chr_ns - control.count.mis.cor), ncol = 2)
      #       test.mis.cor <- chisq.test(temp, correct=F) ##prevent anti-conservative
      #             
      #chi-square test on corrected data by Multiple Imputation
      test.mis.cor.MI <- MI(corr_fac, count.HS1_HN0_S1, case.count.mis, control.count.mis, n_chr_s, n_chr_ns)
      
      #power under true model
      case.count.true <- sum(HS)
      control.count.true <- sum(HN)
      
      temp <- matrix(c(case.count.true, n_chr_s - case.count.true, control.count.true, n_chr_ns - control.count.true), ncol = 2)
      test.true <- chisq.test(temp, correct=F) ##prevent anti-conservative
      
      detach(data.sample)
      if(mark%%10 == 0) print(mark) #print every 1000 iterations
      mark <<- mark + 1
      
      #print(c(case.count.true, case.count.mis.cor, control.count.true, control.count.mis.cor))
      c(true=test.true$p.value, naive=test.mis$p.value, mi=test.mis.cor.MI)
    }
    )
  }
  sim.result
}

#sim.result.1.1 <- sim(rep=1000000, f=0.01, var=1, sibpair=1000, SRR=5, r=c(1.1, 1.3, 1.5))
#sim.result.5.1 <- sim(rep=1000000, f=0.05, var=1, sibpair=1000, SRR=5, r=c(1.1, 1.3, 1.5))

sim.result <- sim(rep=n_rep, f=f, var=0, n_sibpair=1000, SRR=5, r=r)

## Write out your results to a csv file
write.csv(data.frame(seed=seed, f=f, r=r, true=sim.result[1,], naive=sim.result[2,],
                     mi=sim.result[3,]), paste("res_",f,"_",r,"_",seed,".csv",sep=""), row.names=FALSE)
## R.miSniam
## End(Not run)
## Not run: