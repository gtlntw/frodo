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
# print(r) #relative risk
## ... your simualtion code
##read in template haplotype
haplotype <- read.table("out_100k_10k_1kb.hap-1", header=F)
colnames(haplotype) <- c("HAP", "CHROM", paste("SNP", 1:(ncol(haplotype)-2), sep=""))
snp <-read.table("out_100k_10k_1kb.pos-1", header=T)
#make allele 1 is the minor allele, 2 is the common allele
temp.idx <- snp$FREQ1 > snp$FREQ2
temp.freq <- snp$FREQ2
snp$FREQ2[temp.idx] <- snp$FREQ1[temp.idx]
snp$FREQ1[temp.idx] <- temp.freq[temp.idx]
#also change the genotype file
haplotype[,which(temp.idx==T)+2] <- 3 - haplotype[,which(temp.idx==T)+2]

#allele frequency
nrow(snp) #total number of snp
sum(snp$FREQ1 < 0.05) # number of snp with f < 0.05
sum(snp$FREQ1 <= 0.01) # number of snp with f < 0.01
sum(snp$FREQ1 == 0.0001) # number of singletons

##assign risk variants and the corresponding effect size (proportional to allele frequency)
null <- FALSE
n_haplo <- 10000
n_snp <- ncol(haplotype)-2
prevalence <- 0.05
b0_sqrt <- sqrt(prevalence)   #baseline

# set.seed(0110)
##generate risk haplotypes
# variant.id <- sort(sample(which(snp$FREQ1 < 0.02), 20))
variant.id <- c(4,5,6,8,10,13,15,23,25,26,27,32,34,36,37,39,43,44,45,46)
print(variant.id)
risk.haplo.id <- which(apply(2-haplotype[, variant.id+2], 1, sum)>0)
(risk.haplo.f <- mean(apply(2-haplotype[, variant.id+2], 1, sum)>0)) #carrier haplotype frequency
# f= 0.0503
# risk.variant.id <- sort(sample(variant.id, 10))
risk.variant.id <- c(5,8,25,26,32,37,39,44,45,46)
print(risk.variant.id) #print risk variants
b <- rep(1, length=(ncol(haplotype)-2)) #initialize every haplotype is neutral
if(null==FALSE) b[risk.variant.id] <- abs(log10(snp$FREQ1[risk.variant.id])) #effect size is log10 of its allele frequency
#calculate the haplotype variants p(A|h)
haplotype.risk <<- apply(2-haplotype[, -c(1:2)], 1, function(x) prod(b^x)*b0_sqrt)
#calculate the probability of drawing given the affected status of founder p(h|A)
# haplotype.draw_prob <- (haplotype.risk %o% haplotype.risk)/n_haplo^2/prevalence
(mu <<- mean(haplotype.risk[risk.haplo.id]/b0_sqrt)) #mean relative risk
(sigma2 <<- var(haplotype.risk[risk.haplo.id]/b0_sqrt)) #variance relative risk

##gene drop simulation for three generations
family_strct.2g2c <- data.frame(family=c(1,1,1,1), person=c(1,2,3,4), father=c(0,0,1,1),
                                mother=c(0,0,2,2), sex=c(1,2,1,1), affect=c(0,0,2,2)) #1=male, 2=female, 1=unaffected, 2=affected

mark=0
sim.result <- replicate(n_rep, {
  if(mark %% 10 == 0) print(mark)
  mark <<- mark+1
  ##remove father and mother
  family_generated <- gene_family(family_strct=family_strct.2g2c, n=1000, haplotype.risk=haplotype.risk)
  temp.idx <- which(family_generated$data_family$person==1:2) #take out founders
  family_generated$data_family <- family_generated$data_family[-temp.idx, ]
  family_generated$tran_vec <- family_generated$tran_vec[-temp.idx, ]
  
  #convert sibpair data from Merlin format to each haplotype a row
  # family_generated$data_family
  sibpair_haplotype <- as.data.frame(matrix(-999, nrow=4000, ncol=6+nrow(snp), dimnames=list(NULL, colnames(family_generated$data_family)[1:56])))
  sibpair_haplotype[seq(1,4000, by=2),1:6] <- family_generated$data_family[,1:6]
  sibpair_haplotype[seq(2,4000, by=2),1:6] <- family_generated$data_family[,1:6]
  sibpair_haplotype[seq(1,4000, by=2),7:56] <- family_generated$data_family[,7:56]
  sibpair_haplotype[seq(2,4000, by=2),7:56] <- family_generated$data_family[,57:106]
  
  sibpair_tran_vec <- as.data.frame(matrix(-999, nrow=4000, ncol=1, dimnames=list(NULL, "Chrom")))
  sibpair_tran_vec[seq(1,4000, by=2),1] <- family_generated$tran_vec [,2]
  sibpair_tran_vec[seq(2,4000, by=2),1] <- family_generated$tran_vec [,3]
  
  #allele frequency on the shared and non-shared chromosome
  S_hap <- rep(NA, nrow(sibpair_tran_vec)) #indicate the Shared (only count once) and non-shared chromosome
  for(i in seq(1,nrow(sibpair_tran_vec), by=4)) { #S:shared  NS:non-shared  D: to dump shared
    S_hap[i] <- ifelse(sibpair_tran_vec[i,1]==sibpair_tran_vec[i+2,1], "S", "NS")
    S_hap[i+1] <- ifelse(sibpair_tran_vec[i+1,1]==sibpair_tran_vec[i+1+2,1], "S", "NS")
    S_hap[i+2] <- ifelse(sibpair_tran_vec[i+2,1]==sibpair_tran_vec[i,1], "D", "NS")
    S_hap[i+3] <- ifelse(sibpair_tran_vec[i+3,1]==sibpair_tran_vec[i+1,1], "D", "NS")
  }
  table(S_hap) #no. of shared and non-shared chromosome
  # tapply(apply(2-sibpair_haplotype[, 6+risk.variant.id], 1, sum), S_hap, sum) #allele on shared and non-shared
  # tapply(apply(2-sibpair_haplotype[, 6+risk.variant.id], 1, sum)>0, S_hap, sum) #number of risk haplotype based on collaping framework
  
  #test statistics
  true.result <- prop.test(table(S_hap[which(!S_hap=="D")], apply(2-sibpair_haplotype[, 6+variant.id], 1, sum)[which(!S_hap=="D")]>0), correct=F)
  
  #create S vector
  S_sibpair <- (sibpair_tran_vec[seq(1,length(S_hap), by=4),] == sibpair_tran_vec[seq(3,length(S_hap), by=4),]) +
    (sibpair_tran_vec[seq(2,length(S_hap), by=4),] == sibpair_tran_vec[seq(4,length(S_hap), by=4),])
  table(S_sibpair)
  
  #need EM for double-het in S=1
  #at a position what is the allele freq. on share and non-shared chromosome
  ##EM algorithm for imputation
  n_sample=1000
  data_sibpair <- 2-sibpair_haplotype[, 6+variant.id] #only causal snp
  cn <- 4*sum(S_sibpair==0) + 2*sum(S_sibpair==1)  #no. of non-shared chromosomes
  cs <- sum(S_sibpair==1) + 2*sum(S_sibpair==2) #no. of shared chromosome
  
  ##count u, cs, cn
  para <- array(NA, c(3, ncol(data_sibpair)), list(c("u", "kn", "ks"), colnames(data_sibpair)))
  amb_sibpair <- array(FALSE, c(1000,ncol(data_sibpair)))
  for(j in 1:ncol(data_sibpair)) {
    u <- kn <- ks <- 0
    for(i in 1:n_sample) {
      idx <- (i-1)*4+1
      if(S_sibpair[i]==0) {
        kn <- kn + sum(data_sibpair[c(idx+0:3), j])
      }
      
      if(S_sibpair[i]==1) {
        sib1 <- sum(data_sibpair[c(idx+0:1), j])
        sib2 <- sum(data_sibpair[c(idx+2:3), j])
        sib_sum <- sib1+sib2
        if(sib1==1 & sib2==1) {
          u <- u + 1
          #         print(data_sibpair[c(idx+0:3), j])
          amb_sibpair[i,j] <- T
        }
        else {
          if(sib_sum==1) kn <- kn + 1
          if(sib_sum==3) {
            kn <- kn + 1
            ks <- ks + 1
          }
          if(sib_sum==4) {
            kn <- kn + 2
            ks <- ks + 1
          }
        }
      }
      
      if(S_sibpair[i]==2) {
        ks <- ks + sum(data_sibpair[c(idx+0:1), j])
      }
      #   u
      # kn
      # ks
    }
    para[,j] <- c(u, kn, ks)
  }
  para
  # sum(apply(amb_sibpair, 1, sum)>0) # no. of S=1 sibpair with ambiguous genotype
  # sum(apply(amb_sibpair, 2, sum)>0) # no. of SNP with S=1 sibpair with ambiguous genotype
  
  #estimate the probaility of having a shared variant
  EM <- function(para, cn, cs) {
    factor <- rep(NA, ncol(para))
    for(i in 1:ncol(para)) {#i <- 1 iterate over the positions
      u <- para[1,i] #number of unknown configuration (Double hets in IBD 1)
      if(u==0) {
        factor[i] <- NA
        next
      }
      
      #initialization
      kn <- para[2,i] #known non-shared variants (On ibd 0 or single variants on ibd 1)
      ks <- para[3,i] #known shared variants  (On ibd 2 or more than two variants on ibd 1)
      cn <- cn #total number of non-shared chromosomes
      cs <- cs # total number of shared chromosomes
      
      pn.init <- kn/(cn-u*2) #probability of rare variant on non-shared chromosome
      pn.cur <- ifelse(pn.init==0, runif(1), pn.init)
      ps.init <- ks/(cs-u) #probability of rare variant on shared chromosome
      ps.cur <- ifelse(ps.init==0, runif(1), ps.init)
      delta <- Inf
      iter <- 1
      
      while(delta > 10^-6) {
        #E step
        #us <- u*ps.cur/(pn.cur+ps.cur)
        #un <- u*pn.cur/(pn.cur+ps.cur)
        us <- u* ps.cur*(1-pn.cur)^2 / (ps.cur*(1-pn.cur)^2 + (1-ps.cur)*pn.cur^2)
        un <- u* (1-ps.cur)*pn.cur^2 / (ps.cur*(1-pn.cur)^2 + (1-ps.cur)*pn.cur^2)
        #M-step
        pn.new <- (kn + 2*un)/cn
        ps.new <- (ks+us)/cs
        #print(c(mu.new, sigma2.new, f.new, cor.factor, optim.result$value))
        
        #check convergence
        delta <- max(abs(pn.cur - pn.new), abs(ps.cur - ps.new))
        pn.cur <- pn.new
        ps.cur <- ps.new
        
        #print(c(pn.cur, ps.cur, iter))
        #iter <- iter + 1
      }
      #c(pn.init, ps.init)
      factor[i] <- result <- c(ps.cur*(1-pn.cur)^2 / (ps.cur*(1-pn.cur)^2 + (1-ps.cur)*pn.cur^2))
    }
    #output a correction factor for each position
    factor
  }
  
  #assume the phase is know but still need to solve ambiguity, assuming the phase is know for
  prob_shared <- EM(para=para, cn=cn, cs=cs) #the probability of being a shared variant
  amb_sibpair_idx <- which(amb_sibpair==T, TRUE) #index of which sibpair and snp is ambiguous
  impute <- function() {
    data_sibpair_tmp <- data_sibpair
    if(nrow(amb_sibpair_idx)>0) {
      for(a in 1:nrow(amb_sibpair_idx)) {
        sibpair_idx <- amb_sibpair_idx[a, 1]  #which sibpair
        haplotype_idx <- (sibpair_idx-1)*4+1:4 #which haplotypes
        snp_idx <- amb_sibpair_idx[a, 2] #which snp
        s_ns_status <- rbinom(1, 1, prob_shared[snp_idx]) #impute if the variant is shared or not
        s_chrm_idx <- which(S_hap[haplotype_idx]=="S")
        data_sibpair_tmp[haplotype_idx, snp_idx] <- if(s_ns_status==T) {#if shared
          if(s_chrm_idx==1) { #update the genotype data according to where the shared chromosome is
            c(1,0,1,0)
          }else{
            c(0,1,0,1)
          }
        } else{ # the variant is not shared
          if(s_chrm_idx==1) { #update the genotype data according to where the shared chromosome is
            c(0,1,0,1)
          }else{
            c(1,0,1,0)
          }
        }
      }
    }
    return(data_sibpair_tmp)
  }
  # data_sibpair_tmp <- impute()
  # tapply(apply(data_sibpair[,], 1, sum), S_hap, sum) #true
  # tapply(apply(data_sibpair_tmp[,], 1, sum), S_hap, sum) #imputed
  
  #apply test with multiple imputation
  n_chr_s <- cs
  n_chr_ns <- cn
  MI <- function() {
    diff <- NULL
    var <- NULL
    p1_D <- NULL
    p2_D <- NULL
    D <- 10
    
    for(i in 1:D) {
      data_sibpair_tmp <- impute()
      impute_result <- tapply(apply(data_sibpair_tmp, 1, sum)>0, S_hap, sum) #count no. of risk haplotype
      xns <- impute_result[2]
      xs <- impute_result[3]
      
      p1 <- xs/n_chr_s
      p2 <- xns/n_chr_ns
      p <- (xs+xns)/(n_chr_s+n_chr_ns)
      
      p1_D <- cbind(p1_D,p1)
      p2_D <- cbind(p2_D,p2)
      diff <- cbind(diff, p1-p2)
      var <- cbind(var, p*(1-p)*(1/n_chr_s+1/n_chr_ns))
    }
    
    TD <- mean(diff)
    VARD <- mean(var) + (1+1/D)*sum((diff-TD)^2)/(D-1)
    c(mean(p1_D), mean(p2_D),pchisq(TD^2/VARD, df=1, lower=F))
  }
  haplotype.result <- MI()
  # t.test(apply(2-sibpair_haplotype[, 6+risk.variant.id], 1, sum)[which(!S_hap=="D")]>0 ~ S_hap[which(!S_hap=="D")], var.equal=TRUE)
  
  
  #use genotype form and do the test in a pair for chromosome
  prob_shared <- EM(para=para, cn=cn, cs=cs) #the probability of being a shared variant
  amb_sibpair_idx <- which(amb_sibpair==T, TRUE) #index of which sibpair and snp is ambiguous
  # count number of allele b{y a sibpair
  impute_geno <- function() {
    allele_sibpair <- array(NA, c(1000,4), list(NULL, c("s", "ns1", "ns2", "ambiguous")))
    for(i in 1:1000) {
      allele_count_sib1 <- apply(data_sibpair[(i-1)*4+1:2,], 2, sum)
      allele_count_sib2 <- apply(data_sibpair[(i-1)*4+3:4,], 2, sum)
      allele_count <- allele_count_sib1 + allele_count_sib2
      allele_sibpair[i, ] <- if(S_sibpair[i]==0) {
        c(NA,sum(allele_count_sib1),sum(allele_count_sib2),NA)
      }else{if(S_sibpair[i]==2){
        c(sum(allele_count)/2, NA, NA, NA)
      }else{
        c(sum(allele_count %in% c(3,4)), sum(allele_count %in% c(1,3)) +
            2*sum(allele_count ==4), NA, sum(allele_count ==2)) #special treatment to count for S=1 sibpair
      }
      }
    }
    cbind(allele_sibpair, S_sibpair)
    allele_sibpair_impute <- allele_sibpair
    if(nrow(amb_sibpair_idx)>0) {
      for(i in 1:nrow(amb_sibpair_idx)){
        sibpair_idx <- amb_sibpair_idx[i, 1]  #which sibpair
        snp_idx <- amb_sibpair_idx[i, 2] #which snp
        s_ns_status <- rbinom(1, 1, prob_shared[snp_idx]) #impute if the variant is shared or not
        allele_sibpair_impute[sibpair_idx, ] <- if(s_ns_status==1) {
          allele_sibpair_impute[sibpair_idx, ] + c(1,0, NA, -1)#update no. of shared variant
        }else{
          allele_sibpair_impute[sibpair_idx, ] + c(0,2, NA, -1)
        }
      }
    }
    allele_sibpair_impute
  }
  allele_sibpair_impute <- impute_geno()
  cbind(allele_sibpair_impute, S_sibpair)
  apply(allele_sibpair_impute, 2, sum, na.rm=T)
  
  #apply test with multiple imputation using random pairing for S=1 sibpairs with ambiguity
  MI_geno <- function() {
    diff <- NULL
    var <- NULL
    p1_D <- NULL
    p2_D <- NULL
    D <- 10
    S1_idx <- which(S_sibpair==1) #which sibpair is S=1 for ramdom pairing
    no_S1 <- length(which(S_sibpair==1))
    n_case <- length(which(S_sibpair==2)) + (no_S1-(no_S1 %% 2))/2
    n_control <- 2*length(which(S_sibpair==0)) + (no_S1-(no_S1 %% 2))
    
    for(i in 1:D) {
      allele_sibpair_impute <- impute_geno()
      
      xns <- sum(allele_sibpair_impute[which(S_sibpair==0), c("ns1", "ns2")]>0) + sum(allele_sibpair_impute[which(S_sibpair==1), "ns1"]>0)
      xs <- sum(allele_sibpair_impute[which(S_sibpair==2), "s"]>0) + sum((allele_sibpair_impute[S1_idx[seq(1, no_S1-(no_S1 %% 2), by=2)], "s"] + allele_sibpair_impute[S1_idx[seq(2, no_S1-(no_S1 %% 2), by=2)], "s"])>0) ##S=1 discards the last sibpair when there are even number of sibpairs
      
      p1 <- xs/n_case
      p2 <- xns/n_control
      p <- (xs+xns)/(n_case+n_control)
      
      p1_D <- cbind(p1_D,p1)
      p2_D <- cbind(p2_D,p2)
      diff <- cbind(diff, p1-p2)
      var <- cbind(var, p*(1-p)*(1/n_case+1/n_control))
    }
    
    TD <- mean(diff)
    VARD <- mean(var) + (1+1/D)*sum((diff-TD)^2)/(D-1)
    c(mean(p1_D), mean(p2_D),pchisq(TD^2/VARD, df=1, lower=F))
  }
  genotype.result <- MI_geno()
  
  
  ##CC simulation
  generated_case_control <- gene_case_control(n_case_control_pair=1000, haplotype.risk)
  ##test CC
  ##hypothesis testing count the number of carrier haplotypes transmitted to the offsprings
  n_case_control_pair <-1000
  #check if founder's haplotype carries any variant's with f < 0.1
  snp2look.idx <- variant.id # snp to look for
  
  #check if carrier haplotype
  carrier <- apply(generated_case_control, 1, function(x) {
    h1 <- x[3:(2+n_snp)]
    h2 <- x[-(1:(2+n_snp))]
    #     sum(h1[snp2look.idx]==1)>0
    #     sum(h2[snp2look.idx]==1)>0
    c(h1.carrier=sum(h1[snp2look.idx]==1)>0, h2.carrier=sum(h2[snp2look.idx]==1)>0)
  })
  
  carrier.case <- sum(carrier[,1:n_case_control_pair])
  carrier.control <- sum(carrier[,(n_case_control_pair+1):(2*n_case_control_pair)])
  cc.haplotype.result <- prop.test(c(carrier.case, carrier.control), c(2*n_case_control_pair, 2*n_case_control_pair), correct=F) #turn off correct to avoid a conservative test 
  
  carrier.genotype.case <- sum(apply(carrier[,1:n_case_control_pair], 2, sum)>0)
  carrier.genotype.control <- sum(apply(carrier[,(n_case_control_pair+1):(2*n_case_control_pair)], 2, sum)>0)
  cc.genotype.result <- prop.test(c(carrier.genotype.case, carrier.genotype.control), c(n_case_control_pair, n_case_control_pair), correct=F) #turn off correct to avoid a conservative test 
  
  c(trafic.true=true.result$p.value, trafic.haplotype=haplotype.result[3], trafic.genotype=genotype.result[3], cc.haplotype=cc.haplotype.result$p.value, cc.genotype=cc.genotype.result$p.value)
})

sim.result

write.csv(data.frame(seed=seed, mu=mu, sigma2=sigma2, risk.haplo.f=risk.haplo.f,
                     trafic.true=sim.result[1,], trafic.haplotype=sim.result[2,], trafic.genotype=sim.result[3,],
                     cc.haplotype=sim.result[4,], cc.genotype=sim.result[5,])
          , paste("res_",round(mu,3),"_",round(sigma2,3),"_",seed,".csv",sep=""), row.names=FALSE)
## R.miSniam
## End(Not run)
## Not run: