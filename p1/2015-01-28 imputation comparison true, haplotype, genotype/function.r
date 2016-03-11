##Theoratical allele frequency: E(HS|AAr, S=1,2) and E(HN|AAr, S=0,1)on a chromosome
EH <- function(mu=1, sigma2=0, f=0.01) {
  mu <- mu
  sigma2 <- sigma2
  p <- f
  P_AA <- 0.25*(1+p*(mu-1))^4 + 0.5*(1+p*(mu-1))^2*(1+p*(mu^2+sigma2-1)) + 0.25*(1+p*(mu^2+sigma2-1))^2
  
  E_S_C <- (0.5*(1+p*(mu-1))^2*(1+p*(mu^2+sigma2-1))) / P_AA + 2*(0.25*(1+p*(mu^2+sigma2-1))^2) / P_AA
  E_NS_C <- 4*(0.25*(1+p*(mu-1))^4)/P_AA + 2*0.5*(1+p*(mu-1))^2*(1+p*(mu^2+sigma2-1))/P_AA
  
  P_HS1 <- ((mu^2+sigma2)*((1-p)+mu*p)^2*p*0.5+(mu^2+sigma2)*2*p*(1-p)*0.25)/P_AA
  P_HS2 <- (mu^2+sigma2)^2*p^2*0.25/P_AA
  E_HS <- P_HS1 + 2*P_HS2
  
  P_HN1 <- (mu*4*p*(1-p)^3*0.25 + (mu*(1-p)+mu*(mu^2+sigma2)*p)*2*p*(1-p)*0.5)/P_AA
  P_HN2 <- (mu^2*6*p^2*(1-p)^2*0.25 + (mu^2*(1-p)+mu^2*(mu^2+sigma2)*p)*p^2*0.5)/P_AA
  P_HN3 <- (mu^3*4*p^3*(1-p)*0.25)/P_AA
  P_HN4 <- (mu^4*p^4*0.25)/P_AA
  E_HN <- P_HN1 + 2*P_HN2 + 3*P_HN3 + 4*P_HN4
  
  c(E_HS/E_S_C, E_HN/E_NS_C)
}
EH(mu=5, sigma2=0, f=0.02)

## E(N_C|AAr) E(N_S_C|AAr) E(N_NS_C|AAr) no. of chromosomes
E_C <- function(mu=1, sigma2=0, f=0.01) {
  mu <- mu
  sigma2 <- sigma2
  p <- f
  P_AA <- 0.25*(1+p*(mu-1))^4 + 0.5*(1+p*(mu-1))^2*(1+p*(mu^2+sigma2-1)) + 0.25*(1+p*(mu^2+sigma2-1))^2
  
  E_S_C <- (1*0.5*(1+p*(mu-1))^2*(1+p*(mu^2+sigma2-1)) + 2*0.25*(1+p*(mu^2+sigma2-1))^2)/P_AA
  E_NS_C <- (4*(0.25*(1+p*(mu-1))^4) + 2*0.5*(1+p*(mu-1))^2*(1+p*(mu^2+sigma2-1)))/P_AA
  
  E_C <- (4*(0.25*(1+p*(mu-1))^4) + 3*0.5*(1+p*(mu-1))^2*(1+p*(mu^2+sigma2-1)) + 2*0.25*(1+p*(mu^2+sigma2-1))^2)/P_AA
  
  c(E_S_C, E_NS_C, E_C)
}

power.sas <- function(mu=1.2, sigma2=3, f=0.01, n_sb=50, alpha=10^-6) {
  e_c <- E_C(mu=mu, sigma2=sigma2, f=f) ## expected no. of chromosomes per sibpair
  N <- n_sb*e_c[3] #no. of  total chromosomes
  n1 <- n_sb*e_c[1]  ##no. of shared chromosomes
  n2 <- n_sb*e_c[2]  #no. of non-shared chromosomes
  
  P_A <- EH(mu=mu, sigma2=sigma2, f=f) ##proportion under HA
  p1 <- P_A[1]
  p2 <- P_A[2]
  p0 <- 0
  w1 <- n1/N
  w2 <- n2/N
  
  pnorm( ((p2-p1-p0)*(N*w1*w2)^.5 -  qnorm(1-alpha/2)*((w1*p1+w2*p2)*(1-w1*p1-w2*p2))^.5)/(w2*p1*(1-p1)+w1*p2*(1-p2))^.5 ) +
    pnorm( (-(p2-p1-p0)*(N*w1*w2)^.5 -  qnorm(1-alpha/2)*((w1*p1+w2*p2)*(1-w1*p1-w2*p2))^.5)/(w2*p1*(1-p1)+w1*p2*(1-p2))^.5) 
}
power.sas(mu=3.0, sigma2=0, f=0.01, n_sb=1000, alpha=10^-6)

##use chisquare distribution with noncentrality parameter to calculate power
power.sas.ncp <- function(mu=1.2, sigma2=3, f=0.01, n_sb=50, alpha=10^-6) {
  e_c <- E_C(mu=mu, sigma2=sigma2, f=f) ## expected no. of chromosomes per sibpair
  N <- n_sb*e_c[3] #no. of  total chromosomes
  n1 <- n_sb*e_c[1]  ##no. of shared chromosomes
  n2 <- n_sb*e_c[2]  #no. of non-shared chromosomes
  
  P_A <- EH(mu=mu, sigma2=sigma2, f=f) ##proportion under HA
  p1 <- P_A[1]
  p2 <- P_A[2]
  p0 <- 0
  w1 <- n1/N
  w2 <- n2/N
  
  pchisq((-qnorm(1-alpha/2)*((w1*p1+w2*p2)*(1-w1*p1-w2*p2))^.5/(w2*p1*(1-p1)+w1*p2*(1-p2))^.5)^2, df=1, ncp=(-(p2-p1-p0)*(N*w1*w2)^.5/(w2*p1*(1-p1)+w1*p2*(1-p2))^.5)^2, lower=F)
}
power.sas.ncp(mu=3.0, sigma2=0, f=0.01, n_sb=1000, alpha=10^-6)

##sample size calculation -- doesn't really use
samplesize.sas <- function(mu=1.2, sigma2=3, f=0.01, power=0.8, alpha=10^-6) {
  e_c <- E_C(mu=mu, sigma2=sigma2, f=f) ## expected no. of chromosomes per sibpair
  P_A <- EH(mu=mu, sigma2=sigma2, f=f) ##proportion under HA
  p1 <- P_A[1]
  p2 <- P_A[2]
  p0 <- 0
  
  obj_fn <- function(n_sb) {
    N <- n_sb*e_c[3] #no. of  total chromosomes
    n1 <- n_sb*e_c[1]  ##no. of shared chromosomes
    n2 <- n_sb*e_c[2]  #no. of non-shared chromosomes
    w1 <- n1/N
    w2 <- n2/N
    
    (pnorm( ((p2-p1-p0)*(N*w1*w2)^.5 -  qnorm(1-alpha/2)*((w1*p1+w2*p2)*(1-w1*p1-w2*p2))^.5)/(w2*p1*(1-p1)+w1*p2*(1-p2))^.5 ) +
       pnorm( (-(p2-p1-p0)*(N*w1*w2)^.5 -  qnorm(1-alpha/2)*((w1*p1+w2*p2)*(1-w1*p1-w2*p2))^.5)/(w2*p1*(1-p1)+w1*p2*(1-p2))^.5) 
     - power)^2
  }
  
  init_n <- (qnorm(1-alpha)*((0.33*p1+0.67*p2)*(1-0.33*p1-0.67*p2))^0.5 + qnorm(power)*(0.67*p1*(1-p1)+0.33*p2*(1-p2))^0.5)^2 / (0.33*0.67*(p2-p1-p0)^2)
  ceiling(optimize(obj_fn, c(1,init_n))$minimum)
  #ceiling(optim(1000, obj_fn)$par)
}
samplesize.sas(mu=5, sigma2=0, f=0.01, power=0.8, alpha=10^-6)


##Theoratical allele frequency: Unaffected and Affected, E(H|A) and E(H|AAr) on a chromosome for comparison
EH_comp <- function(mu=1, sigma2=0, f=0.01, out="A") {
  mu <- mu
  sigma2 <- sigma2
  p <- f
  P_A <- (1+p*(mu-1))^2
  P_AA <- 0.25*(1+p*(mu-1))^4 + 0.5*(1+p*(mu-1))^2*(1+p*(mu^2+sigma2-1)) + 0.25*(1+p*(mu^2+sigma2-1))^2
  if(out=="A") {
    return(c(EH_A <- (mu*2*p*(1-p)+2*mu^2*p^2)/P_A/2, p))
  }
  else {
    return(c(EH_AA <- (0.25*2*p*(1-p)*(mu*(1+p*(mu-1))^2) + 0.5*2*p*(1-p)*(0.5*((1-p)*mu+p*mu^2)+0.5*((1-p)*(mu^2+sigma2)+p*mu*(mu^2+sigma2))) 
                       + 0.25*2*p*(1-p)*(mu^2+sigma2) + 2*(0.25*p^2*mu^2*(1+p*(mu-1))^2+0.5*p^2*((1-p)*mu*(mu^2+sigma2)+p*mu^2*(mu^2+sigma2))+0.25*p^2*(mu^2+sigma2)^2) )/P_AA/2,p))
  }  		   
}

EH_comp(mu=5, sigma2=0, f=0.02, out="A")


power_comp.sas <- function(mu=1.2, sigma2=3, f=0.01, n_pair=50, alpha=10^-6, out="A") {
  N <- 4*n_pair
  P_A <- EH_comp(mu=mu, sigma2=sigma2, f=f, out=out) ##proportion under HA
  p1 <- P_A[1]
  p2 <- P_A[2]
  p0 <- 0
  w1 <- 0.5
  w2 <- 0.5
  
  pnorm( ((p2-p1-p0)*(N*w1*w2)^.5 -  qnorm(1-alpha/2)*((w1*p1+w2*p2)*(1-w1*p1-w2*p2))^.5)/(w2*p1*(1-p1)+w1*p2*(1-p2))^.5 ) +
    pnorm( (-(p2-p1-p0)*(N*w1*w2)^.5 -  qnorm(1-alpha/2)*((w1*p1+w2*p2)*(1-w1*p1-w2*p2))^.5)/(w2*p1*(1-p1)+w1*p2*(1-p2))^.5) 
}
power_comp.sas(mu=3.0, sigma2=0, f=0.01, n_pair=1000, alpha=10^-6, out="A")

samplesize_comp.sas <- function(mu=1.2, sigma2=3, f=0.01, power=0.8, alpha=10^-6, out="A") {
  P_A <- EH_comp(mu=mu, sigma2=sigma2, f=f, out=out) ##proportion under HA
  p1 <- P_A[1]
  p2 <- P_A[2]
  p0 <- 0
  w1 <- 0.5
  w2 <- 0.5
  
  obj_fn <- function(n_pair) {
    
    N <- 4*n_pair    
    (pnorm( ((p2-p1-p0)*(N*w1*w2)^.5 -  qnorm(1-alpha/2)*((w1*p1+w2*p2)*(1-w1*p1-w2*p2))^.5)/(w2*p1*(1-p1)+w1*p2*(1-p2))^.5 ) +
       pnorm( (-(p2-p1-p0)*(N*w1*w2)^.5 -  qnorm(1-alpha/2)*((w1*p1+w2*p2)*(1-w1*p1-w2*p2))^.5)/(w2*p1*(1-p1)+w1*p2*(1-p2))^.5) 
     - power)^2
  }
  
  init_n <- (qnorm(1-alpha)*((0.5*p1+0.5*p2)*(1-0.5*p1-0.5*p2))^0.5 + qnorm(power)*(0.5*p1*(1-p1)+0.5*p2*(1-p2))^0.5)^2 / (0.5*0.5*(p2-p1-p0)^2)
  ceiling(optimize(obj_fn, c(1,init_n))$minimum)
  #ceiling(optim(1000, obj_fn)$par)
  #c(init_n, optimize(obj_fn, c(1,init_n)))
}
samplesize_comp.sas(mu=5, sigma2=0, f=0.01, power=0.8, alpha=10^-6)



### Using Gamma distribution to generate the efffect of haplotype in a batch
##simulate a case/control
gene.data <- function(m=1, var=0, f=0.01, SRR=5, p_dis=0.01, pop=FALSE) {
  adj <- ifelse(var==0 & m==1, 1, ifelse(var==0, m, 1))
  c = (m - 1) #shift to accommondate gamma distribution
  beta <- ifelse(var==0 & m==1, 0, var/c) #prevent not a number error when var==0 & m==1
  alpha = c/beta
  KL <- (1+f*(m-1))^2 #contribution from locus
  KLKLR <- 0.25*(1+f*(m-1))^4 + 0.5*(1+f*(m-1))^2*(1+f*(m^2+var-1)) + 0.25*(1+f*(m^2+var-1))^2
  KG <- p_dis/KL #contribution from other genome
  SRR <- SRR #sibling relaive risk
  KGKGR <- SRR*p_dis*p_dis/KLKLR #implement the heriatbility from other locus given SRR
  
  ##generate case and control individuals from population
  if(pop){ #skip individual part
    n_pop <- 5000000
    H1 <- rbinom(n_pop,1,f) #if the first haplotype carries risk raviant
    H2 <- rbinom(n_pop,1,f) #if the second haplotype carries risk raviant
    H1_rr <- ifelse(H1 == 1, rgamma(length(H1), alpha, scale=beta) + adj, 1) #RR of the first haplotype
    H2_rr <- ifelse(H2 == 1, rgamma(length(H2), alpha, scale=beta) + adj, 1) #RR of the second haplotype
    penetrance <- H1_rr*H2_rr*KG #penetrance of disease given haplotypes
    penetrance <- ifelse(penetrance>1, 1, penetrance)
    dis <- rbinom(length(penetrance),1,penetrance) #disease status
    data.pop <<- data.frame(H1=H1, H2=H2, H1_rr=H1_rr, H2_rr=H2_rr, penetrance=penetrance, dis=dis)
  }
  
  ##generate sibpairs
  n_family <- 50000000
  H1_f <- rbinom(n_family,1,f) #if the first haplotype carries risk raviant
  H2_f <- rbinom(n_family,1,f) #if the second haplotype carries risk raviant
  H1_rr_f <- ifelse(H1_f == 1, rgamma(length(H1_f), alpha, scale=beta) + adj, 1) #RR of the first haplotype
  H2_rr_f <- ifelse(H2_f == 1, rgamma(length(H2_f), alpha, scale=beta) + adj, 1) #RR of the second haplotype
  H1_m <- rbinom(n_family,1,f) #if the first haplotype carries risk raviant
  H2_m <- rbinom(n_family,1,f) #if the second haplotype carries risk raviant
  H1_rr_m <- ifelse(H1_m == 1, rgamma(length(H1_m), alpha, scale=beta) + adj, 1) #RR of the first haplotype
  H2_rr_m <- ifelse(H2_m == 1, rgamma(length(H2_m), alpha, scale=beta) + adj, 1) #RR of the second haplotype
  
  H1_inh_s <- sample(c(1,2),n_family, replace=TRUE) 
  H1_s <- ifelse(H1_inh_s==1, H1_f, H2_f)
  H1_rr_s <- ifelse(H1_inh_s==1, H1_rr_f, H2_rr_f)
  
  H2_inh_s <- sample(c(1,2),n_family, replace=TRUE) 
  H2_s <- ifelse(H2_inh_s==1, H1_m, H2_m)
  H2_rr_s <- ifelse(H2_inh_s==1, H1_rr_m, H2_rr_m)
  
  #penetrance_s <- H1_rr_s*H2_rr_s #penetrance of disease given haplotypes
  #penetrance_s <- ifelse(penetrance_s>1, 1, penetrance_s)
  #dis_s <- rbinom(length(penetrance_s),1,penetrance_s) #disease status
  
  H1_inh <- sample(c(1,2), n_family, replace=TRUE) 
  H1 <- ifelse(H1_inh==1, H1_f, H2_f)
  H1_rr <- ifelse(H1_inh==1, H1_rr_f, H2_rr_f)
  
  H2_inh <- sample(c(1,2), n_family, replace=TRUE) 
  H2 <- ifelse(H2_inh==1, H1_m, H2_m)
  H2_rr <- ifelse(H2_inh==1, H1_rr_m, H2_rr_m)
  
  rm(H1_f, H2_f, H1_rr_f, H2_rr_f, H1_m, H2_m, H1_rr_m, H2_rr_m) #save memory
  
  S0 <- (H1_inh_s != H1_inh) & (H2_inh_s != H2_inh)
  S2 <- (H1_inh_s == H1_inh) & (H2_inh_s == H2_inh)
  S1 <- ((H1_inh_s == H1_inh) & (H2_inh_s != H2_inh)) | ((H1_inh_s != H1_inh) & (H2_inh_s == H2_inh))
  S11 <- ((H1_inh_s == H1_inh) & (H2_inh_s != H2_inh))
  S12 <- ((H1_inh_s != H1_inh) & (H2_inh_s == H2_inh))
  
  S <- (S0==TRUE)*0 + (S1==TRUE)*1 + (S2==TRUE)*2
  HS.mis <- HS <- (S1==TRUE&S11)*(H1) + (S1==TRUE&S12)*(H2) + (S2==TRUE)*(H1+H2)
  HN.mis <- HN <- (S0==TRUE)*(H1+H2+H1_s+H2_s) + (S1==TRUE&S11)*(H2+H2_s) + (S1==TRUE&S12)*(H1+H1_s)
  HS.mis <- ifelse(S1==TRUE & HS==0 & HN==2, 1, HS.mis) 
  HN.mis <- ifelse(S1==TRUE & HS==0 & HN==2, 0, HN.mis) 
  
  penetrance <- H1_rr*H2_rr*H1_rr_s*H2_rr_s*KGKGR #penetrance of disease given haplotypes of both siblings
  penetrance <- ifelse(penetrance>1, 1, penetrance)
  dis <- rbinom(length(penetrance),1,penetrance) #disease status of both affected
  
  #data.family <- data.frame(H1_f, H2_f, H1_rr_f, H2_rr_f,   #for varification
  #                          H1_m, H2_m, H1_rr_m, H2_rr_m,
  #                          H1_inh, H1, H1_rr, H2_inh, H2, H2_rr, penetrance, dis,
  #                          H1_inh_s, H1_s, H1_rr_s, H2_inh_s, H2_s, H2_rr_s, penetrance_s, dis_s)
  
  idx <- which(dis==1)
  
  data.family <<- data.frame(H1=H1[idx], H1_rr=H1_rr[idx], H2=H2[idx], H2_rr=H2_rr[idx], dis=dis[idx],
                             H1_s=H1_s[idx], H1_rr_s=H1_rr_s[idx], H2_s=H2_s[idx], H2_rr_s=H2_rr_s[idx],
                             S0=S0[idx], S1=S1[idx], S2=S2[idx], S11=S11[idx], S12=S12[idx],
                             HS=HS[idx],HN=HN[idx],S=S[idx], HS.mis=HS.mis[idx],HN.mis=HN.mis[idx])
}

#generate one by one
gene.pop <- function(m=1, var=0, f=0.01, SRR=5, p_dis=0.01, case=T, n_sample=1000, g_gamma=T) {
  adj <- ifelse(var==0 & m==1, 1, ifelse(var==0, m, 1))
  c = (m - 1) #shift to accommondate gamma distribution
  beta <- ifelse(var==0 & m==1, 0, var/c) #prevent not a number error when var==0 & m==1
  alpha = c/beta
  KL <- (1+f*(m-1))^2 #contribution from locus
  KLKLR <- 0.25*(1+f*(m-1))^4 + 0.5*(1+f*(m-1))^2*(1+f*(m^2+var-1)) + 0.25*(1+f*(m^2+var-1))^2
  KG <- p_dis/KL #contribution from other genome
  SRR <- SRR #sibling relaive risk
  KGKGR <- SRR*p_dis*p_dis/KLKLR #implement the heriatbility from other locus given SRR
  
  ##generate case and control individuals from population
  data.pop <- as.data.frame(array(NA, c(n_sample, 5), dimnames=list(x=NULL, y=c("H1", "H1_rr", "H2", "H2_rr", "dis"))))
  ##generate sibpairs
  n_counter <- 1
  while(n_counter <= n_sample){
    H1 <- rbinom(1,1,f) #if the first haplotype carries risk raviant
    H2 <- rbinom(1,1,f) #if the second haplotype carries risk raviant
    if(g_gamma==T){
      H1_rr <- ifelse(H1 == 1, rgamma(1, alpha, scale=beta) + adj, 1) #RR of the first haplotype
      H2_rr <- ifelse(H2 == 1, rgamma(1, alpha, scale=beta) + adj, 1) #RR of the second haplotype
    } else{
  
    }
    penetrance <- H1_rr*H2_rr*KG #penetrance of disease given haplotypes
    penetrance <- ifelse(penetrance>1, 1, penetrance)
    dis <- rbinom(1,1,penetrance) #disease status
    
    if(dis==(case==T)) {
      #         print(counter)
      data.pop[n_counter,] <- c(H1, H1_rr, H2, H2_rr, dis)
      
      n_counter <- n_counter + 1
    }
  }
  data.pop
}
gene.pop(m=1, var=0, f=0.01, SRR=5, p_dis=0.01, case=T, n_sample=1000)

#generate one by one
gene.sibpair <- function(m=1, var=0, f=0.01, SRR=5, p_dis=0.01, n_sibpair=1000, g_gamma=T) {
  adj <- ifelse(var==0 & m==1, 1, ifelse(var==0, m, 1))
  c = (m - 1) #shift to accommondate gamma distribution
  beta <- ifelse(var==0 & m==1, 0, var/c) #prevent not a number error when var==0 & m==1
  alpha = c/beta
  KL <- (1+f*(m-1))^2 #contribution from locus
  KLKLR <- 0.25*(1+f*(m-1))^4 + 0.5*(1+f*(m-1))^2*(1+f*(m^2+var-1)) + 0.25*(1+f*(m^2+var-1))^2
  KG <- p_dis/KL #contribution from other genome
  SRR <- SRR #sibling relaive risk
  KGKGR <- SRR*p_dis*p_dis/KLKLR #implement the heriatbility from other locus given SRR
    
  data.sibpair <- as.data.frame(array(NA, c(n_sibpair, 19), dimnames=list(x=NULL, y=c("H1", "H1_rr", "H2", "H2_rr", "dis", "H1_s", "H1_rr_s", "H2_s", "H2_rr_s", "S0", "S1", "S2", "S11", "S12",
                                                                        "HS","HN","S", "HS.mis","HN.mis"))))
  ##generate sibpairs
  n_counter <- 1
  while(n_counter <= n_sibpair){
    H1_f <- rbinom(1,1,f) #if the first haplotype carries risk raviant
    H2_f <- rbinom(1,1,f) #if the second haplotype carries risk raviant
    H1_m <- rbinom(1,1,f) #if the first haplotype carries risk raviant
    H2_m <- rbinom(1,1,f) #if the second haplotype carries risk raviant
    if(g_gamma==T) {
      H1_rr_f <- ifelse(H1_f == 1, rgamma(1, alpha, scale=beta) + adj, 1) #RR of the first haplotype
      H2_rr_f <- ifelse(H2_f == 1, rgamma(1, alpha, scale=beta) + adj, 1) #RR of the second haplotype
      H1_rr_m <- ifelse(H1_m == 1, rgamma(1, alpha, scale=beta) + adj, 1) #RR of the first haplotype
      H2_rr_m <- ifelse(H2_m == 1, rgamma(1, alpha, scale=beta) + adj, 1) #RR of the second haplotype
    } else{
      H1_rr_f <- ifelse(H1_f == 1, rnorm(1, mean=m, sd=sqrt(var)), 1) #RR of the first haplotype
      H2_rr_f <- ifelse(H2_f == 1, rnorm(1, mean=m, sd=sqrt(var)), 1) #RR of the second haplotype
      H1_rr_m <- ifelse(H1_m == 1, rnorm(1, mean=m, sd=sqrt(var)), 1) #RR of the first haplotype
      H2_rr_m <- ifelse(H2_m == 1, rnorm(1, mean=m, sd=sqrt(var)), 1) #RR of the second haplotype 
    }    
    H1_inh_s <- sample(c(1,2),1, replace=TRUE) 
    H1_s <- ifelse(H1_inh_s==1, H1_f, H2_f)
    H1_rr_s <- ifelse(H1_inh_s==1, H1_rr_f, H2_rr_f)
    
    H2_inh_s <- sample(c(1,2),1, replace=TRUE) 
    H2_s <- ifelse(H2_inh_s==1, H1_m, H2_m)
    H2_rr_s <- ifelse(H2_inh_s==1, H1_rr_m, H2_rr_m)
    
    #penetrance_s <- H1_rr_s*H2_rr_s #penetrance of disease given haplotypes
    #penetrance_s <- ifelse(penetrance_s>1, 1, penetrance_s)
    #dis_s <- rbinom(length(penetrance_s),1,penetrance_s) #disease status
    
    H1_inh <- sample(c(1,2), 1, replace=TRUE) 
    H1 <- ifelse(H1_inh==1, H1_f, H2_f)
    H1_rr <- ifelse(H1_inh==1, H1_rr_f, H2_rr_f)
    
    H2_inh <- sample(c(1,2), 1, replace=TRUE) 
    H2 <- ifelse(H2_inh==1, H1_m, H2_m)
    H2_rr <- ifelse(H2_inh==1, H1_rr_m, H2_rr_m)
    
    S0 <- (H1_inh_s != H1_inh) & (H2_inh_s != H2_inh)
    S2 <- (H1_inh_s == H1_inh) & (H2_inh_s == H2_inh)
    S1 <- ((H1_inh_s == H1_inh) & (H2_inh_s != H2_inh)) | ((H1_inh_s != H1_inh) & (H2_inh_s == H2_inh))
    S11 <- ((H1_inh_s == H1_inh) & (H2_inh_s != H2_inh))
    S12 <- ((H1_inh_s != H1_inh) & (H2_inh_s == H2_inh))
    
    S <- (S0==TRUE)*0 + (S1==TRUE)*1 + (S2==TRUE)*2
    HS.mis <- HS <- (S1==TRUE&S11)*(H1) + (S1==TRUE&S12)*(H2) + (S2==TRUE)*(H1+H2)
    HN.mis <- HN <- (S0==TRUE)*(H1+H2+H1_s+H2_s) + (S1==TRUE&S11)*(H2+H2_s) + (S1==TRUE&S12)*(H1+H1_s)
    HS.mis <- ifelse(S1==TRUE & HS==0 & HN==2, 1, HS.mis) 
    HN.mis <- ifelse(S1==TRUE & HS==0 & HN==2, 0, HN.mis) 
    
    penetrance <- H1_rr*H2_rr*H1_rr_s*H2_rr_s*KGKGR #penetrance of disease given haplotypes of both siblings
    penetrance <- ifelse(penetrance>1, 1, penetrance)
    dis <- rbinom(1,1,penetrance) #disease status of both affected
    
    if(dis==1) {
      #         print(counter)
      data.sibpair[n_counter,] <- c(H1, H1_rr, H2, H2_rr, dis, H1_s, H1_rr_s, H2_s, H2_rr_s, S0, S1, S2, S11, S12,
                                    HS,HN,S, HS.mis,HN.mis)
      n_counter <- n_counter + 1
    }
  }
  data.sibpair
}
gene.sibpair(n_=10)

#generate one by one
gene.sibpair.linkage <- function(m=1, var=0, f=0.01, lambda_mz=11, lambda_o=4, p_dis=0.01, n_sibpair=1000, g_gamma=T) {
  lambda_s <- (lambda_mz+2*lambda_o+1)/4
  print(paste("lambda_s =", (lambda_mz+2*lambda_o+1)/4))
  z0 <- 0.25*1/lambda_s
  z1 <- 0.5*lambda_o/lambda_s
  z2 <- 0.25*lambda_mz/lambda_s
  theta <- (n_sibpair*z0)*log10(z0)+(n_sibpair*z1)*log10(z1)+(n_sibpair*z2)*log10(z2) - ((n_sibpair*z0)*log10(1/4)+(n_sibpair*z1)*log10(1/2)+(n_sibpair*z2)*log10(1/4))
  print(paste("z0 =", z0, "z1 =", z1, "z2 =", z2))
  print(paste("LOD =", theta))
  adj <- ifelse(var==0 & m==1, 1, ifelse(var==0, m, 1))
  c = (m - 1) #shift to accommondate gamma distribution
  beta <- ifelse(var==0 & m==1, 0, var/c) #prevent not a number error when var==0 & m==1
  alpha = c/beta
  K <- p_dis
  
  data.sibpair <- as.data.frame(array(NA, c(n_sibpair, 19), dimnames=list(x=NULL, y=c("H1", "H1_rr", "H2", "H2_rr", "dis", "H1_s", "H1_rr_s", "H2_s", "H2_rr_s", "S0", "S1", "S2", "S11", "S12",
                                                                                      "HS","HN","S", "HS.mis","HN.mis"))))
  ##generate sibpairs
  n_counter <- 1
  while(n_counter <= n_sibpair){
    H1_f <- rbinom(1,1,f) #if the first haplotype carries risk raviant
    H2_f <- rbinom(1,1,f) #if the second haplotype carries risk raviant
    H1_m <- rbinom(1,1,f) #if the first haplotype carries risk raviant
    H2_m <- rbinom(1,1,f) #if the second haplotype carries risk raviant
    if(g_gamma==T) {
      H1_rr_f <- ifelse(H1_f == 1, rgamma(1, alpha, scale=beta) + adj, 1) #RR of the first haplotype
      H2_rr_f <- ifelse(H2_f == 1, rgamma(1, alpha, scale=beta) + adj, 1) #RR of the second haplotype
      H1_rr_m <- ifelse(H1_m == 1, rgamma(1, alpha, scale=beta) + adj, 1) #RR of the first haplotype
      H2_rr_m <- ifelse(H2_m == 1, rgamma(1, alpha, scale=beta) + adj, 1) #RR of the second haplotype
    } else{
      H1_rr_f <- ifelse(H1_f == 1, rnorm(1, mean=m, sd=sqrt(var)), 1) #RR of the first haplotype
      H2_rr_f <- ifelse(H2_f == 1, rnorm(1, mean=m, sd=sqrt(var)), 1) #RR of the second haplotype
      H1_rr_m <- ifelse(H1_m == 1, rnorm(1, mean=m, sd=sqrt(var)), 1) #RR of the first haplotype
      H2_rr_m <- ifelse(H2_m == 1, rnorm(1, mean=m, sd=sqrt(var)), 1) #RR of the second haplotype 
    }    
    H1_inh_s <- sample(c(1,2),1, replace=TRUE) 
    H1_s <- ifelse(H1_inh_s==1, H1_f, H2_f)
    H1_rr_s <- ifelse(H1_inh_s==1, H1_rr_f, H2_rr_f)
    
    H2_inh_s <- sample(c(1,2),1, replace=TRUE) 
    H2_s <- ifelse(H2_inh_s==1, H1_m, H2_m)
    H2_rr_s <- ifelse(H2_inh_s==1, H1_rr_m, H2_rr_m)
    
    #penetrance_s <- H1_rr_s*H2_rr_s #penetrance of disease given haplotypes
    #penetrance_s <- ifelse(penetrance_s>1, 1, penetrance_s)
    #dis_s <- rbinom(length(penetrance_s),1,penetrance_s) #disease status
    
    H1_inh <- sample(c(1,2), 1, replace=TRUE) 
    H1 <- ifelse(H1_inh==1, H1_f, H2_f)
    H1_rr <- ifelse(H1_inh==1, H1_rr_f, H2_rr_f)
    
    H2_inh <- sample(c(1,2), 1, replace=TRUE) 
    H2 <- ifelse(H2_inh==1, H1_m, H2_m)
    H2_rr <- ifelse(H2_inh==1, H1_rr_m, H2_rr_m)
    
    S0 <- (H1_inh_s != H1_inh) & (H2_inh_s != H2_inh)
    S2 <- (H1_inh_s == H1_inh) & (H2_inh_s == H2_inh)
    S1 <- ((H1_inh_s == H1_inh) & (H2_inh_s != H2_inh)) | ((H1_inh_s != H1_inh) & (H2_inh_s == H2_inh))
    S11 <- ((H1_inh_s == H1_inh) & (H2_inh_s != H2_inh))
    S12 <- ((H1_inh_s != H1_inh) & (H2_inh_s == H2_inh))
    
    S <- (S0==TRUE)*0 + (S1==TRUE)*1 + (S2==TRUE)*2
    HS.mis <- HS <- (S1==TRUE&S11)*(H1) + (S1==TRUE&S12)*(H2) + (S2==TRUE)*(H1+H2)
    HN.mis <- HN <- (S0==TRUE)*(H1+H2+H1_s+H2_s) + (S1==TRUE&S11)*(H2+H2_s) + (S1==TRUE&S12)*(H1+H1_s)
    HS.mis <- ifelse(S1==TRUE & HS==0 & HN==2, 1, HS.mis) 
    HN.mis <- ifelse(S1==TRUE & HS==0 & HN==2, 0, HN.mis) 
    
    penetrance <- H1_rr*H2_rr*H1_rr_s*H2_rr_s*K*K*ifelse(S==2, lambda_mz, ifelse(S==1, lambda_o, 1)) #penetrance of disease given haplotypes of both siblings
    penetrance <- ifelse(penetrance>1, 1, penetrance)
    dis <- rbinom(1,1,penetrance) #disease status of both affected
    
    if(dis==1) {
      #         print(counter)
      data.sibpair[n_counter,] <- c(H1, H1_rr, H2, H2_rr, dis, H1_s, H1_rr_s, H2_s, H2_rr_s, S0, S1, S2, S11, S12,
                                    HS,HN,S, HS.mis,HN.mis)
      n_counter <- n_counter + 1
    }
  }
  data.sibpair
}
gene.sibpair.linkage(n_sibpair=10)


##EM algorithm for imputation
EM <- function(data) {
  #initialization
  pn.init <- sum(subset(data, S==0)$HN.mis)/(4*sum(data$S==0)) #probability of rare variant on shared chromosome
  pn.cur <- ifelse(pn.init==0, runif(1), pn.init)
  ps.init <- sum(subset(data, S==2)$HS.mis)/(2*sum(data$S==2)) #probability of rare variant on non-shared chromosome
  ps.cur <- ifelse(ps.init==0, runif(1), ps.init)
  kn <- sum(subset(data, S==0 | (S==1 & !(HS.mis==1 &HN.mis==0)))$HN.mis) #known non-shared variants (On ibd 0 or single variants on ibd 1)
  ks <- sum(subset(data, (S==1 & !(HS.mis==1 &HN.mis==0)) | S==2)$HS.mis) #known shared variants  (On ibd 2 or more than two variants on ibd 1)
  cn <- 4*sum(data$S==0) + 2*sum(data$S==1) # total number of shared chromosomes
  cs <- sum(data$S==1) + 2*sum(data$S==2) #total number of non-shared chromosomes
  u <- sum(data$S1==TRUE & data$HS.mis==1 & data$HN.mis==0) #number of unknown configuration (Double hets in IBD 1)
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
  c(ps.cur*(1-pn.cur)^2 / (ps.cur*(1-pn.cur)^2 + (1-ps.cur)*pn.cur^2))
}

EM2 <- function(data) {
  #initialization
  p.init <- (sum(subset(data, S==0)$HN.mis)+sum(subset(data, S==2)$HS.mis))/((2*sum(data$S==2))+(4*sum(data$S==0))) #probability of rare variant on shared chromosome
  p.cur <- ifelse(p.init==0, runif(1), p.init)
  k <- sum(subset(data, S==0 | (S==1 & !(HS.mis==1 &HN.mis==0)))$HN.mis) + sum(subset(data, (S==1 & !(HS.mis==1 &HN.mis==0)) | S==2)$HS.mis)
  c <- 4*sum(data$S==0) + 2*sum(data$S==1) + sum(data$S==1) + 2*sum(data$S==2) 
  u <- sum(data$S1==TRUE & data$HS.mis==1 & data$HN.mis==0) #number of unknown configuration (Double hets in IBD 1)
  delta <- Inf
  iter <- 1
  
  while(delta > 10^-6) {
    us <- 1-p.cur
    un <- p.cur
    #M-step
    p.new <- (k + us + 2*un)/c
    #print(c(mu.new, sigma2.new, f.new, cor.factor, optim.result$value))
    
    #check convergence
    delta <- max(abs(p.cur - p.new))
    p.cur <- p.new
    
    #print(c(pn.cur, ps.cur, iter))
    #iter <- iter + 1
  }
  #c(pn.init, ps.init)
  c(1-p.cur)
}

MI <- function(corr_fac, count.HS1_HN0_S1, case.count.mis, control.count.mis, n_chr_s, n_chr_ns) {
  diff <- NULL
  var <- NULL
  u <- count.HS1_HN0_S1
  N <- 10
  
  for(i in 1:N) {
    us <- rbinom(1, u, corr_fac)
    uns <- u - us
    
    xs <- case.count.mis - u + us
    xns <- control.count.mis + 2*uns
    
    p1 <- xs/n_chr_s
    p2 <- xns/n_chr_ns
    p <- (xs+xns)/(n_chr_s+n_chr_ns)
    
    diff <- cbind(diff, p1-p2)
    var <- cbind(var, p*(1-p)*(1/n_chr_s+1/n_chr_ns))
  }
  
  TD <- mean(diff)
  VARD <- mean(var) + (1+1/N)*sum((diff-TD)^2)/(N-1)
  pchisq(TD^2/VARD, df=1, lower=F)
}

#Power curve simulation
#count allele on shared and non-shared chromosomes
#sibpair
sim <- function(rep=100000, var=0, f=0.01, sibpair=1000, SRR=5, r=seq(1, 3.5, by=0.1)) {
  sim.result <- NULL
  r <<- r
  for(i in r) {
    print(paste("Generating Data...m=", i, "var=", var, "f=", f))
    gene.data(m=i, var=var, f=f, SRR=SRR)
    mark <<- 1 #reset the current iteration
    print(mark)  #start iteration
    sim <- replicate(rep, {
      ##sibpair internal control
      data.sample <- data.family[sample(1:nrow(data.family), sibpair), ]
      attach(data.sample)
      n_chr_s <- sum(S1==1) + 2*sum(S2==1)
      n_chr_ns <- 4*sum(S0==1) + 2*sum(S1==1)
      
      case.count.mis <- sum(HS.mis)
      control.count.mis <- sum(HN.mis)
      
      #chi-square test on misclassification data
      temp <- matrix(c(case.count.mis, n_chr_s - case.count.mis, control.count.mis, n_chr_ns - control.count.mis), ncol = 2)
      test.mis <- chisq.test(temp, correct=F) ##prevent anti-conservative
      
      
      #chi-square test on corrected data by EM
      corr_fac <- EM(data.sample)
      
      count.HS1_HN0_S1 <- sum(S1==TRUE & HS.mis==1 & HN.mis==0)
      
      case.count.mis.cor <- corr_fac*count.HS1_HN0_S1 + (case.count.mis - count.HS1_HN0_S1)
      control.count.mis.cor <- control.count.mis + 2*(1-corr_fac)*count.HS1_HN0_S1
      
      temp <- matrix(c(case.count.mis.cor, n_chr_s - case.count.mis.cor, control.count.mis.cor, n_chr_ns - control.count.mis.cor), ncol = 2)
      test.mis.cor <- chisq.test(temp, correct=F) ##prevent anti-conservative
      
      #chi-square test on corrected data by Multiple Imputation
      test.mis.cor.MI <- MI(corr_fac, count.HS1_HN0_S1, case.count.mis, control.count.mis, n_chr_s, n_chr_ns)
      
      #power under true model
      case.count.true <- sum(HS)
      control.count.true <- sum(HN)
      
      temp <- matrix(c(case.count.true, n_chr_s - case.count.true, control.count.true, n_chr_ns - control.count.true), ncol = 2)
      test.true <- chisq.test(temp, correct=F) ##prevent anti-conservative
      
      detach(data.sample)
      if(mark%%1000 == 0) print(mark) #print every 1000 iterations
      mark <<- mark + 1
      c(test.mis$p.value, test.mis.cor.MI, test.true$p.value)
    }
    )
    sim.result <- cbind(sim.result, rowMeans(sim<10^-6))
    #sim.result
  }
  sim.result
}

sim_var <- function(rep=100000, mu=1, f=0.01, sibpair=1000, SRR=5, r=seq(1, 3.5, by=0.1)) {
  sim.result <- NULL
  r <<- r
  for(i in r) {
    print(paste("Generating Data...var=", i, "mu=", mu, "f=", f))
    gene.data(m=mu, var=i, f=f, SRR=SRR)
    mark <<- 1 #reset the current iteration
    print(mark)  #start iteration
    sim <- replicate(rep, {
      ##sibpair internal control
      data.sample <- data.family[sample(1:nrow(data.family), sibpair), ]
      attach(data.sample)
      n_chr_s <- sum(S1==1) + 2*sum(S2==1)
      n_chr_ns <- 4*sum(S0==1) + 2*sum(S1==1)
      
      case.count.mis <- sum(HS.mis)
      control.count.mis <- sum(HN.mis)
      
      #chi-square test on misclassification data
      temp <- matrix(c(case.count.mis, n_chr_s - case.count.mis, control.count.mis, n_chr_ns - control.count.mis), ncol = 2)
      test.mis <- chisq.test(temp, correct=F) ##prevent anti-conservative
      
      
      #chi-square test on corrected data by EM
      corr_fac <- EM(data.sample)
      
      count.HS1_HN0_S1 <- sum(S1==TRUE & HS.mis==1 & HN.mis==0)
      
      case.count.mis.cor <- corr_fac*count.HS1_HN0_S1 + (case.count.mis - count.HS1_HN0_S1)
      control.count.mis.cor <- control.count.mis + 2*(1-corr_fac)*count.HS1_HN0_S1
      
      temp <- matrix(c(case.count.mis.cor, n_chr_s - case.count.mis.cor, control.count.mis.cor, n_chr_ns - control.count.mis.cor), ncol = 2)
      test.mis.cor <- chisq.test(temp, correct=F) ##prevent anti-conservative
      
      #chi-square test on corrected data by Multiple Imputation
      test.mis.cor.MI <- MI(corr_fac, count.HS1_HN0_S1, case.count.mis, control.count.mis, n_chr_s, n_chr_ns)
      
      #power under true model
      case.count.true <- sum(HS)
      control.count.true <- sum(HN)
      
      temp <- matrix(c(case.count.true, n_chr_s - case.count.true, control.count.true, n_chr_ns - control.count.true), ncol = 2)
      test.true <- chisq.test(temp, correct=F) ##prevent anti-conservative
      
      detach(data.sample)
      if(mark%%1000 == 0) print(mark) #print every 1000 iterations
      mark <<- mark + 1
      c(test.mis$p.value, test.mis.cor.MI, test.true$p.value)
    }
    )
    sim.result <- cbind(sim.result, rowMeans(sim<10^-6))
    #sim.result
  }
  sim.result
}




###############Interaction model########################################
inter_model <- function(p=c(0.01, 0.05), trr_alpha=2, trr_beta=6.2, gamma=1, bl=0.0000001) {
  p <- p #Allele frequencies
  
  phi <- function(x) ifelse(x==1, p[1], 1-p[1])
  phj <- function(x) ifelse(x==1, p[1], 1-p[1])
  pgs <- function(x) ifelse(x==1, p[2], 1-p[2])
  pgt <- function(x) ifelse(x==1, p[2], 1-p[2])
  
  trr_alpha <- trr_alpha #target relative risk
  trr_beta <- trr_beta
  gamma <- gamma
  bl <- bl
  PofS <- c(0.25, 0.5, 0.25)
  
  ##new obj function
  obj <- function(x, trr_alpha, trr_beta, gamma) {
    beta_L <- x[1] 
    beta_G <- x[2]
    
    i=0;j=0;k=0;l=0;
    for(hj in 0:1) {
      for(gs in 0:1) {
        for(gt in 0:1){
          i=i+(beta_G*gamma)^(gs+gt)*gamma^(hj*(gs+gt))*pgs(gs)*pgt(gt)
          j=j+(beta_G)^(gs+gt)*gamma^(hj*(gs+gt))*pgs(gs)*pgt(gt)
        }
      }
      k=k+beta_L^(hj)*phj(hj)*i
      l=l+beta_L^(hj)*phj(hj)*j
    }
    rr_alpha=(beta_L*k)/l
    
    i=0;j=0;k=0;l=0;
    for(gt in 0:1) {
      for(hi in 0:1) {
        for(hj in 0:1){
          i=i+(beta_L*gamma)^(hi+hj)*gamma^(gt*(hi+hj))*phi(hi)*phj(hj)
          j=j+(beta_L)^(hi+hj)*gamma^(gt*(hi+hj))*phi(hi)*phj(hj)
        }
      }
      k=k+beta_G^(gt)*pgt(gt)*i
      l=l+beta_G^(gt)*pgt(gt)*j
    }
    rr_beta=(beta_G*k)/l
    
    return(abs(rr_alpha-trr_alpha)+abs(rr_beta - trr_beta))
  }
  #obj(c(1.671843,1.93), 2, 2, 2)
  
  result <- optim(c(1,1), obj, trr_alpha=trr_alpha, trr_beta=trr_beta, gamma=gamma)
  result$par
  
  #solution of beta_L and beta_G
  alpha <- result$par[1]
  beta <- result$par[2]
  
  #penetrance of one individual P(D|A,B)
  Pgt <- array(NA, c(3,3))
  for(a in 0:2) {
    for (b in 0:2) {
      Pgt[a+1,b+1] <- bl*alpha^a*beta^b*gamma^(a*b)
      if(Pgt[a+1,b+1]>1) Pgt[a+1,b+1] <- 1
      #print(Pgt[a+1,b+1])
    }
  }
  
  # P(i, j|A,B,S,loc)  given the first sibling and sharing status what is the genotype of the second sibling
  Ptr <- array(NA, c(3,3,3,2)) 
  for (loc in 0:1) {
    for (a in 0:2) {
      for (i in 0:2) {
        if(a==i) {
          Ptr[a+1, i+1, 2+1, loc+1]=1;
        }
        else
          Ptr[a+1, i+1, 2+1, loc+1]=0;
        Ptr[a+1, i+1, 0+1, loc+1]=p[loc+1]^i*(1-p[loc+1])^(2-i)
      }
      Ptr[a+1, 1+1, 0+1, loc+1]=2*Ptr[a+1, 1+1, 0+1, loc+1]
    }
    Ptr[2+1, 0+1, 1+1, loc+1]=0
    Ptr[0+1, 2+1, 1+1, loc+1]=0
    Ptr[2+1, 1+1, 1+1, loc+1]=1-p[loc+1]
    Ptr[2+1, 2+1, 1+1, loc+1]=p[loc+1]
    Ptr[1+1, 2+1, 1+1, loc+1]=0.5*p[loc+1]
    Ptr[1+1, 1+1, 1+1, loc+1]=0.5
    Ptr[1+1, 0+1, 1+1, loc+1]=0.5*(1-p[loc+1])
    Ptr[0+1, 1+1, 1+1, loc+1]=p[loc+1]
    Ptr[0+1, 0+1, 1+1, loc+1]=1-p[loc+1]
  }
  
  #penetrance of P(DD|A,B)
  Pab <- array(0, c(3,3)) 
  for(a in 0:2){
    for(b in 0:2){
      for(i in 0:2){
        for(j in 0:2){
          Pab[a+1, b+1]=Pab[a+1, b+1]+Pgt[i+1, j+1]*(PofS[0+1]*Ptr[a+1, i+1, 0+1, 0+1]+PofS[1+1]*Ptr[a+1, i+1, 1+1, 0+1]+PofS[2+1]*Ptr[a+1, i+1, 2+1, 0+1])*(PofS[0+1]*Ptr[b+1, j+1, 0+1, 1+1]+PofS[1+1]*Ptr[b+1, j+1, 1+1, 1+1]+PofS[2+1]*Ptr[b+1, j+1, 2+1, 1+1])
        }
      }
      Pab[a+1, b+1]=Pab[a+1, b+1]*Pgt[a+1, b+1]
    }
  }
  
  #probability of observed genotype P(A,B)
  Pf <- array(0, c(3,3))
  norm <- 0
  for(a in 0:2){
    for(b in 0:2) {
      Pf[a+1, b+1]=p[0+1]^a*(1-p[0+1])^(2-a)*p[1+1]^b*(1-p[1+1])^(2-b)
      if(a==1) Pf[a+1, b+1]=2*Pf[a+1, b+1]
      if(b==1) Pf[a+1, b+1]=2*Pf[a+1, b+1]
      norm = norm + Pf[a+1, b+1]*Pab[a+1, b+1]
    }
  }
  
  #prob. of P(A|DD)
  prob <- array(0, 3)
  for(a in 0:2){
    for(b in 0:2){
      prob[a+1] = prob[a+1]+Pf[a+1, b+1]*Pab[a+1, b+1]/norm
    }
  }
  
  #prevalence
  prev <- 0
  for(a in 0:2){
    for(b in 0:2){
      prev = prev+Pf[a+1, b+1]*Pgt[a+1, b+1]
    }
  }
  
  #sibling relative risk
  srr <- 0
  for(a in 0:2){
    for(b in 0:2){
      srr = srr+Pf[a+1, b+1]*Pab[a+1, b+1]
    }
  }
  prev_dd <- srr ##P(DD)
  srr=srr/prev^2
  
  #print(c(alpha, beta, srr))
  
  #baseline
  prob_b<-array(0, 3)
  norm_b <- 0
  for(a in 0:2){
    for(b in 0:2){
      norm_b = norm_b + Pf[a+1, b+1]*Pgt[a+1, b+1] ##P(D)
      prob_b[a+1]=prob_b[a+1] + Pf[a+1, b+1]*Pgt[a+1, b+1]; ##P(D,A)
    }
  }
  
  for(a in 0:2){
    prob_b[a+1]=prob_b[a+1]/norm_b;  ##P(A|D)
    #print "$a $prob[$a] $prob_b[$a]\n";
  }
  
  Ex_sib=prob[1+1]+2*prob[2+1]; Ex_unrelated=prob_b[1+1]+2*prob_b[2+1];
  #print "Expectations sib-pair $Ex unrelated $Ex_b\t";
  
  #Calculate P(AA|vns1, vns2,vs) summing over P(b,b')
  PAAs <- array(0, c(3,3,3));
  for(vns1 in 0:2){
    for(vns2 in 0:2){
      for(vs in 0:2){
        if(vns1+vs<3 & vns2+vs<3)
        {
          for(b in 0:2){
            for(i in 0:2){
              
              ah=vns1+vs;
              ah2=vns2+vs;
              p_sib=0;
              for(S in 0:2){
                p_sib=p_sib+PofS[S+1]*Ptr[b+1, i+1, S+1, 1+1]*dbinom(b,2,p[1+1]); ##contribution from loci B
                #print "tt $S $PofS[$S] $Ptr[$b, $i, $S, 1]\n";
              }
              #print "$b $i $p_sib\n";
              PAAs[vns1+1, vns2+1, vs+1]=PAAs[vns1+1, vns2+1, vs+1]+Pgt[ah+1, b+1]*Pgt[ah2+1, i+1]*p_sib;
            }}
          #print "test $vns1 $vns2 $vs $PAAs[vns1+1, vns2+1, vs+1]\n";
        }
        else {PAAs[vns1+1, vns2+1, vs+1]=-1;}
      }}}
  #Calculate the prob of shared/nonshared variants based on S
  Psns <- array(0, c(3,3,3))
  norm <- 0;
  for(vns1 in 0:2){
    for(vns2 in 0:2){
      for(vs in 0:2){
        if(vns1+vs<3 & vns2+vs<3)
        {
          for(S in 0:2){
            Psns[vns1+1, vns2+1, vs+1]=Psns[vns1+1, vns2+1, vs+1]+dbinom(vs,S,p[0+1])*dbinom(vns1,2-S,p[0+1])*dbinom(vns2,2-S,p[0+1])*PofS[S+1];
          }
          Psns[vns1+1, vns2+1, vs+1]=Psns[vns1+1, vns2+1, vs+1]*PAAs[vns1+1, vns2+1, vs+1]; ##P(DD, ns1, ns2, ns)
          norm=norm+Psns[vns1+1, vns2+1, vs+1]  ##P(DD)
        }else  {Psns[vns1+1, vns2+1, vs+1]=-1;}
      }}}
  #Expectation of shared and non-shared variants
  vncount <- 0;vscount <-0;
  for(vns1 in 0:2){
    for(vns2 in 0:2){
      for(vs in 0:2){
        if(vns1+vs<3 & vns2+vs<3)
        {
          Psns[vns1+1, vns2+1, vs+1]=Psns[vns1+1, vns2+1, vs+1]/norm;  ##P(vns1, vns2, ns|DD)
          vncount=vncount+(vns1+vns2)*Psns[vns1+1, vns2+1, vs+1]; ##E(vns1+vns2|DD)
          vscount=vscount+vs*Psns[vns1+1, vns2+1, vs+1]; ##E(vs|DD)
          
          #print "$vns1 $vns2 $vs $Psns[vns1+1, vns2+1, vs+1]\n"
        }
      }}}
  
  
  ##calculate expectation of S given DD
  #penetrance of P(DD|A,B,S)
  Pab_S <- array(0, c(3,3,3))
  for(S in 0:2) {
    for(a in 0:2){
      for(b in 0:2){
        for(i in 0:2){
          for(j in 0:2){
            Pab_S[a+1, b+1, S+1]=Pab_S[a+1, b+1, S+1]+Pgt[i+1, j+1]*Ptr[a+1, i+1, S+1, 0+1]*(PofS[0+1]*Ptr[b+1, j+1, 0+1, 1+1]+PofS[1+1]*Ptr[b+1, j+1, 1+1, 1+1]+PofS[2+1]*Ptr[b+1, j+1, 2+1, 1+1])
          }
        }
        Pab_S[a+1, b+1, S+1]=Pab_S[a+1, b+1, S+1]*Pgt[a+1, b+1]
      }
    }
  }
  
  
  #Calculate P(AA|vns1, vns2,vs, S) summing over P(b,b')
  PAAs_S <- array(0, c(3,3,3,3));
  for(S in 0:2) {
    for(vns1 in 0:2){
      for(vns2 in 0:2){
        for(vs in 0:2){
          if(vns1+vs<3 & vns2+vs<3)
          {
            for(b in 0:2){
              for(i in 0:2){
                
                ah=vns1+vs;
                ah2=vns2+vs;
                p_sib=0;
                for(S2 in 0:2){
                  p_sib=p_sib+PofS[S2+1]*Ptr[b+1, i+1, S2+1, 1+1]*dbinom(b,2,p[1+1]); ##contribution from loci B
                  #print "tt $S $PofS[$S] $Ptr[$b, $i, $S, 1]\n";
                }
                #print "$b $i $p_sib\n";
                PAAs_S[vns1+1, vns2+1, vs+1, S+1]=PAAs_S[vns1+1, vns2+1, vs+1, S+1]+Pgt[ah+1, b+1]*Pgt[ah2+1, i+1]*p_sib;
              }}
            #print "test $vns1 $vns2 $vs $PAAs_S[vns1+1, vns2+1, vs+1]\n";
          }
          else {PAAs_S[vns1+1, vns2+1, vs+1, S+1]=-1;}
        }}}
  }
  
  #P(S|AA)
  Ps <- array(0, 3)
  norm <- 0;
  for(S in 0:2){
    for(vns1 in 0:2){
      for(vns2 in 0:2){
        for(vs in 0:2){
          if(vns1+vs<3 & vns2+vs<3)
          {
            Ps[S+1]=Ps[S+1]+PAAs_S[vns1+1, vns2+1, vs+1, S+1]*dbinom(vs,S,p[0+1])*dbinom(vns1,2-S,p[0+1])*dbinom(vns2,2-S,p[0+1])*PofS[S+1];##P(DD, ns1, ns2, ns)
          }
        }}}
    norm=norm+Ps[S+1]  ##P(DD)
  }
  #Expectation of shared and non-shared chromosome
  ES <- 0;ENS <-0;
  for(S in 0:2){
    Ps[S+1]=Ps[S+1]/norm;  ##P(S|DD)
    ES=ES+S*Ps[S+1]; ##E(vns1+vns2|DD)
    ENS=ENS+(2-S)*2*Ps[S+1]; ##E(vs|DD)
  }
  return(list(EH=c(vscount/ES, vncount/ENS), E_c=c(ES, ENS, ES+ENS), alpha=alpha, beta=beta, prevalence=prev, srr=srr, Ex_sib=Ex_sib, Ex_unrelated=Ex_unrelated))
}

power_inter.sas <- function(p=c(0.01, 0.05), trr_alpha=2, trr_beta=6.2, gamma=1, n_sb=50, alpha=10^-6) {
  result <- inter_model(p, trr_alpha, trr_beta, gamma)
  e_c <- result$E_c  ## expected no. of chromosomes per sibpair
  N <- n_sb*e_c[3] #no. of  total chromosomes
  n1 <- n_sb*e_c[1]  ##no. of shared chromosomes
  n2 <- n_sb*e_c[2]  #no. of non-shared chromosomes
  
  P_A <- result$EH ##proportion under HA
  p1 <- P_A[1]
  p2 <- P_A[2]
  p0 <- 0
  w1 <- n1/N
  w2 <- n2/N
  
  pnorm( ((p2-p1-p0)*(N*w1*w2)^.5 -  qnorm(1-alpha/2)*((w1*p1+w2*p2)*(1-w1*p1-w2*p2))^.5)/(w2*p1*(1-p1)+w1*p2*(1-p2))^.5 ) +
    pnorm( (-(p2-p1-p0)*(N*w1*w2)^.5 -  qnorm(1-alpha/2)*((w1*p1+w2*p2)*(1-w1*p1-w2*p2))^.5)/(w2*p1*(1-p1)+w1*p2*(1-p2))^.5) 
}
power_inter.sas(p=c(0.01, 0.05), trr_alpha=2, trr_beta=6.2, gamma=1, n_sb=1000, alpha=10^-6)

power_inter_comp.sas <- function(p=c(0.01, 0.05), trr_alpha=2, trr_beta=6.2, gamma=1, n_sb=50, alpha=10^-6, out="A") {
  result <- inter_model(p, trr_alpha, trr_beta, gamma)
  N <- n_sb*4 #no. of  total chromosomes
  n1 <- n_sb*2  ##no. of case chromosomes
  n2 <- n_sb*2  #no. of control chromosomes
  
  P_A <- ifelse(out=="A", result$Ex_unrelated / 2,  result$Ex_sib / 2)##proportion under HA
  p1 <- P_A
  p2 <- p[1]
  p0 <- 0
  w1 <- n1/N
  w2 <- n2/N
  
  pnorm( ((p2-p1-p0)*(N*w1*w2)^.5 -  qnorm(1-alpha/2)*((w1*p1+w2*p2)*(1-w1*p1-w2*p2))^.5)/(w2*p1*(1-p1)+w1*p2*(1-p2))^.5 ) +
    pnorm( (-(p2-p1-p0)*(N*w1*w2)^.5 -  qnorm(1-alpha/2)*((w1*p1+w2*p2)*(1-w1*p1-w2*p2))^.5)/(w2*p1*(1-p1)+w1*p2*(1-p2))^.5) 
}
power_inter_comp.sas(p=c(0.01, 0.05), trr_alpha=2, trr_beta=6.2, gamma=1, n_sb=1000, alpha=10^-6)


##realistic example based on Cosi
#simulate the tranmission vector by generation and determine the affected status
gene_family <- function(family_strct=family_strct.2g2c, n_family=100, haplotype.risk) {
  n_family_member <- length(family_strct$person)
  data_family <- matrix(NA, nrow=n_family*n_family_member, ncol=(6+n_snp*2))
  tran_vec <- matrix(NA, nrow=n_family*n_family_member, ncol=3)
  family.haplo <- matrix(NA, nrow=n_family_member, ncol=2)
  #basic strategy is generate each individual one by one
  data_family.idx <- 1 #index of the current family member being generated
  n_family.idx <- 1 #index of how many families have been generated
  while(n_family.idx <= n_family) {
    for(i in 1:n_family_member) {
      disease <<- 0
      while(disease!=family_strct$affect[i]) {
        if(family_strct$father[i]==0 & family_strct$mother[i]==0) { #if founder
          #if(family_strct$affect[i]==1) { #if unaffected directly draw from population haplotype
          haplo.id <- sample(1:n_haplo, 2, replace=T)
          disease_prob <- prod(haplotype.risk[haplo.id])
          #           print(haplo.id)
          disease_prob <- ifelse(disease_prob <= 1, disease_prob, 1)
          #           debug.1 <<- haplo.id
          disease <<- rbinom(1,1, prob=disease_prob) + 1
          
        }
        else{ #if not founder
          haplo.id <- c(sample(family.haplo[family_strct$father[i],], 1), sample(family.haplo[family_strct$mother[i],], 1))
          disease_prob <- prod(haplotype.risk[haplo.id])
          #           print(haplo.id)
          disease_prob <- ifelse(disease_prob <= 1, disease_prob, 1)
          #           debug.2 <<- haplo.id
          disease <<- rbinom(1,1, prob=disease_prob) + 1
        }
      }
      #store haplotype's id
      family.haplo[i,] <- haplo.id
    }
    #save the haplotype file
    letter.idx <- 1 #indicator used in the transmission vector
    for(i in 1:n_family_member) {
      #store transmission vector
      data_family[data_family.idx, ] <- unlist(c(n_family.idx, family_strct[i,2:6], haplotype[family.haplo[i,1],-c(1:2)], haplotype[family.haplo[i,2],-c(1:2)]))
      if(family_strct$father[i]==0 & family_strct$mother[i]==0) { #if founder
        tran_vec[data_family.idx,] <- c(n_family.idx, LETTERS[letter.idx], LETTERS[letter.idx+1])
        letter.idx <- letter.idx + 2
      }
      else{ #if not founder then compare with his/her parents
        current_row <- (n_family.idx-1)*n_family_member
        #print(current_row)
        tran_vec[data_family.idx,] <- c(n_family.idx, ifelse(family.haplo[i,1] == family.haplo[family_strct$father[i],1], tran_vec[family_strct$father[i]+current_row, 2], tran_vec[family_strct$father[i]+current_row, 3]) 
                                        ,ifelse(family.haplo[i,2] == family.haplo[family_strct$mother[i],1], tran_vec[family_strct$mother[i]+current_row, 2], tran_vec[family_strct$mother[i]+current_row, 3])) 
      }
      data_family.idx <- data_family.idx + 1
    }
    #     print(n_family.idx)
    n_family.idx <- n_family.idx + 1
  }
  colnames(data_family) <- c("family","person","father","mother","sex","affect",rep(paste("SNP", 1:n_snp, sep=""),2))
  colnames(tran_vec) <- c("family","h1","h2")
  return(list(data_family=data.frame(data_family, stringsAsFactors=F), tran_vec=data.frame(tran_vec, stringsAsFactors=F)))
}
# family_generated <- gene_family()    
# family_generated 
