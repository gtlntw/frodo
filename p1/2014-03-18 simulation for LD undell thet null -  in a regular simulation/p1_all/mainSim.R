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
print(bbeta)
print(n_rep)
print(n_sib)
## ... your simualtion code
#read in haplotype
hap <- 2 - read.table("out_10k_10k.hap-1", stringsAsFactors=F, header=F)[, -(1:2)] #ignore the first two columns and covert to 0/1
snp <- read.table("out_10k_10k.pos-1", stringsAsFactors=F, header=T)
colnames(hap) <- snp$CHROM_POS #assign colnames

#calculate effect size
m=1; var=0; f=0.01; SRR=5; p_dis=0.01
adj <- ifelse(var==0 & m==1, 1, ifelse(var==0, m, 1))
c = (m - 1) #shift to accommondate gamma distribution
beta <- ifelse(var==0 & m==1, 0, var/c) #prevent not a number error when var==0 & m==1
alpha = c/beta
KL <- (1+f*(m-1))^2 #contribution from locus
KLKLR <- 0.25*(1+f*(m-1))^4 + 0.5*(1+f*(m-1))^2*(1+f*(m^2+var-1)) + 0.25*(1+f*(m^2+var-1))^2
KG <- p_dis/KL #contribution from other genome
SRR <- SRR #sibling relaive risk
KGKGR <- SRR*p_dis*p_dis/KLKLR #implement the heriatbility from other locus given SRR


sim <- function(freq="2", rep=100, n_sib=1000) {
  i <- 0
  sim_result <- replicate(rep, {
    i <<- i+1
    if(i%%100==0) print(i)
    
    #generate sibpair data
    n_snp <- ncol(hap)
    n_sib <- n_sib
    hap_sib <- array(NA, c(n_sib*4, (3+n_snp)), dimnames=list(x=NULL, y=c("sibpair", "S", "inh", snp$CHROM_POS)))
    ##generate sibpairs
    counter=1 
    while(counter<=n_sib) {
      
      h_founder <- sample(1:n_snp, 4)
      H1_f <- h_founder[1]
      H2_f <- h_founder[2]
      H1_m <- h_founder[3]
      H2_m <- h_founder[4]
      H1_rr_f <- 1
      H2_rr_f <- 1
      H1_rr_m <- 1
      H2_rr_m <- 1
      
      H1_inh_s <- sample(c(1,2),1) #sibling 1
      H1_s <- ifelse(H1_inh_s==1, H1_f, H2_f)
      H1_rr_s <- ifelse(H1_inh_s==1, H1_rr_f, H2_rr_f)
      
      H2_inh_s <- sample(c(1,2),1) 
      H2_s <- ifelse(H2_inh_s==1, H1_m, H2_m)
      H2_rr_s <- ifelse(H2_inh_s==1, H1_rr_m, H2_rr_m)
      
      H1_inh <- sample(c(1,2), 1) #sibling 2
      H1 <- ifelse(H1_inh==1, H1_f, H2_f)
      H1_rr <- ifelse(H1_inh==1, H1_rr_f, H2_rr_f)
      
      H2_inh <- sample(c(1,2), 1) 
      H2 <- ifelse(H2_inh==1, H1_m, H2_m)
      H2_rr <- ifelse(H2_inh==1, H1_rr_m, H2_rr_m)
      
      S0 <- (H1_inh_s != H1_inh) & (H2_inh_s != H2_inh)
      S2 <- (H1_inh_s == H1_inh) & (H2_inh_s == H2_inh)
      S1 <- ((H1_inh_s == H1_inh) & (H2_inh_s != H2_inh)) | ((H1_inh_s != H1_inh) & (H2_inh_s == H2_inh))
      S11 <- ((H1_inh_s == H1_inh) & (H2_inh_s != H2_inh))
      S12 <- ((H1_inh_s != H1_inh) & (H2_inh_s == H2_inh))
      
      S <- (S0==TRUE)*0 + (S1==TRUE)*1 + (S2==TRUE)*2
      
      penetrance <- H1_rr*H2_rr*H1_rr_s*H2_rr_s*KGKGR #penetrance of disease given haplotypes of both siblings
      penetrance <- ifelse(penetrance>1, 1, penetrance)
      dis <- rbinom(length(penetrance),1,penetrance) #disease status of both affected
      
      if(dis==1) {
        #         print(counter)
        idx <- 4*(counter-1)+1
        hap_sib[idx,] <- unlist(c(counter, S, H1_s, hap[H1_s,]))
        hap_sib[idx+1,] <- unlist(c(counter, S, H2_s, hap[H2_s,]))##this is wrong
        hap_sib[idx+2,] <- unlist(c(counter, S, H1, hap[H1,]))
        hap_sib[idx+3,] <- unlist(c(counter, S, H2, hap[H2,]))
        counter <- counter + 1
      }
    }
    
    #allele frequency and Sharing status
    snp_f1 <- match(snp$CHROM_POS[which(snp$FREQ1 <0.01)], colnames(hap)) #only f<0.01 snps
    snp_f5 <- match(snp$CHROM_POS[which(snp$FREQ1 <0.05)], colnames(hap)) #only f<0.05 snps
    
    
    #subset the data by allele frequency
    if(freq=="2") data_sample <- hap_sib[, snp_f1] #only f<0.01 snps
    if(freq=="3") data_sample <- hap_sib[, snp_f5] #only f<0.05 snps
    if(freq=="4") data_sample <- hap_sib[, -c(1:3)] #all snp
    
    #determine sharing status
    S_hap <- rep(NA, nrow(hap_sib)) #indicate the Shared (only count once) and non-shared chromosome
    for(i in seq(1,nrow(hap_sib), by=4)) { #S:shared  NS:non-shared  D: to dump shared
      S_hap[i] <- ifelse(hap_sib[i,3]==hap_sib[i+2,3], "S", "NS") 
      S_hap[i+1] <- ifelse(hap_sib[i+1,3]==hap_sib[i+1+2,3], "S", "NS")
      S_hap[i+2] <- ifelse(hap_sib[i+2,3]==hap_sib[i,3], "D", "NS")
      S_hap[i+3] <- ifelse(hap_sib[i+3,3]==hap_sib[i+1,3], "D", "NS")
    }
    
    test_idx <- which(!S_hap=="D") #chromsome to keep
    #1. count the number of variant allele
    result <- t.test(apply(data_sample[,], 1, sum)[test_idx] ~ S_hap[test_idx])
    #2. count the number of risk haplotypes
    #   result <- t.test(apply(data_sample[,], 1, sum)[test_idx]>0 ~ S_hap_sample[test_idx])
    c(p=result$p.value,   result$statistic)
  }
  )
  sim_result
}
#f<all
sim_result_2 <- sim("4", rep=n_rep, n_sib=n_sib)
# hist(sim_result_2[1,])
mean(sim_result_2[1,] < 0.05)
# hist(sim_result_2[2,], freq=F)
# curve(dnorm, add=TRUE, col=2)
mean(sim_result_2[2,])

## Write out your results to a csv file
write.csv(data.frame(seed=seed, bbeta=paste(bbeta,collapse="~"), p=sim_result_2[1,], t=sim_result_2[2,]),
          paste("res",seed,".csv",sep=""), row.names=FALSE)
## R.miSniam
## End(Not run)
## Not run:
