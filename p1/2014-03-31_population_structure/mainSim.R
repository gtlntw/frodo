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
print(prop) #relative risk
## ... your simualtion code
#####################Simulation on population structure for 1. traditional case-control 2. internal control 3. selected cases design
sim_ps <- function(n_sb=1000, rep=100, prop=0) {
  pop1.f <- 0.01
  pop1.p_dis <- 0.01
  pop2.f <- 0.05
  pop2.p_dis <- 0.01
  alpha <- 0.05
  print(paste("proportion =", prop))
  mark <<- 1 #reset the current iteration
  print(mark)  #start iteration
  simresult <- replicate(rep, {
    #generate two populations
    pop1.case <- gene.pop(m=1, var=0, f=pop1.f, SRR=5, p_dis=pop1.p_dis, case=T, n_sample=n_sb*prop)
    pop1.control <-  gene.pop(m=1, var=0, f=pop1.f, SRR=5, p_dis=pop1.p_dis, case=F, n_sample=n_sb*0.5)
    family.sample.pop1 <- gene.sibpair(m=1, var=0, f=pop1.f, SRR=5, p_dis=pop1.p_dis, n_sibpair=n_sb*prop)
    pop2.case <- gene.pop(m=1, var=0, f=pop2.f, SRR=5, p_dis=pop2.p_dis, case=T, n_sample=n_sb*(1-prop))
    pop2.control <-  gene.pop(m=1, var=0, f=pop2.f, SRR=5, p_dis=pop2.p_dis, case=F, n_sample=n_sb*0.5)
    family.sample.pop2 <- gene.sibpair(m=1, var=0, f=pop2.f, SRR=5, p_dis=pop2.p_dis, n_sibpair=n_sb*(1-prop))
    
    #1. traditional case-control
    sample.control <- rbind(pop1.control, pop2.control)
    sample.case <- rbind(pop1.case, pop2.case)
    
    sample.case.count <- sum(sample.case$H1, sample.case$H2)
    sample.control.count <- sum(sample.control$H1, sample.control$H2)
    temp <- matrix(c(sample.case.count, n_sb*2 - sample.case.count, sample.control.count, n_sb*2 - sample.control.count), ncol = 2)
    test.case_control <- chisq.test(temp, correct=F)
    
    #2. internal control
    data.sample <- rbind(family.sample.pop1, family.sample.pop2)
    attach(data.sample)
    n_chr_s <- sum(S1==1) + 2*sum(S2==1)
    n_chr_ns <- 4*sum(S0==1) + 2*sum(S1==1)
    case.count.mis <- sum(HS.mis)
    control.count.mis <- sum(HN.mis)
    
    #chi-square test on corrected data by MI
    corr_fac <- EM(data.sample)
    count.HS1_HN0_S1 <- sum(S1==TRUE & HS.mis==1 & HN.mis==0)
    test.mis.cor.MI <- MI(corr_fac, count.HS1_HN0_S1, case.count.mis, control.count.mis, n_chr_s, n_chr_ns)
    detach(data.sample)
    
    
    #3. selected case design
    sample.selected.case <- rbind(family.sample.pop1, family.sample.pop2)
    sample.selected.case.count <- sum(sample.selected.case$H1, sample.selected.case$H2)
    sample.selected.control <- rbind(pop1.control, pop2.control)
    sample.selected.control.count <- sum(sample.selected.control$H1, sample.selected.control$H2)
    temp <- matrix(c(sample.selected.case.count, n_sb*2 - sample.selected.case.count, sample.selected.control.count, n_sb*2 - sample.selected.control.count), ncol = 2)
    test.selected_case <- chisq.test(temp, correct=F)
    
    
    if(mark%%10 == 0) print(mark) #print every 1000 iterations
    mark <<- mark + 1
    
    c(IC=test.mis.cor.MI, SC=test.selected_case$p.value, CC=test.case_control$p.value)
  })
  simresult
}

sim.result <- sim_ps(n_sb=1000, rep=n_rep, prop=prop)

write.csv(data.frame(seed=seed, prop=prop, IC=sim.result[1,], SC=sim.result[2,],
                     CC=sim.result[3,]), paste("res_",prop,"_",seed,".csv",sep=""), row.names=FALSE)
## R.miSniam
## End(Not run)
## Not run: