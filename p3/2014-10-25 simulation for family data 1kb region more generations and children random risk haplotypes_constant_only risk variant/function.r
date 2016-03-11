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



### Using Gamma distribution to generate the efffect of haplotype
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

### Using Gamma distribution to generate the efffect of haplotype
##simulate family
gene.data.family <- function(m=1, var=0, f=0.01, SRR=5, p_dis=0.01, n_family=1000, structure="") {
  if(!(structure %in% c("2g2c", "3g", "2g3c"))) stop("need to specify the family structure")
  adj <- ifelse(var==0 & m==1, 1, ifelse(var==0, m, 1))
  c = (m - 1) #shift to accommondate gamma distribution
  beta <- ifelse(var==0 & m==1, 0, var/c) #prevent not a number error when var==0 & m==1
  alpha = c/beta
  KL <- (1+f*(m-1))^2 #contribution from locus
  KLKLR <- 0.25*(1+f*(m-1))^4 + 0.5*(1+f*(m-1))^2*(1+f*(m^2+var-1)) + 0.25*(1+f*(m^2+var-1))^2
  KG <- p_dis/KL #contribution from other genome
  SRR <- SRR #sibling relaive risk
  KGKGR <- SRR*p_dis*p_dis/KLKLR #implement the heriatbility from other locus given SRR
  f_father <- EH_comp(mu=m, sigma2=var, f=f)[2]
  f_mother <- EH_comp(mu=m, sigma2=var, f=f)[1]
   
  #print(c(f_father, f_mother))
  
  ##generate sibpairs
  if(structure=="2g2c") {
    data_family <- array(NA, c(n_family, 29), dimnames=list(NULL,c("H1_f", "H1_rr_f", "H2_f", "H2_rr_f", 
                                                              "H1_m", "H1_rr_m", "H2_m", "H2_rr_m",
                                                              "H1", "H1_rr", "H2", "H2_rr",
                                                              "H1_s", "H1_rr_s", "H2_s", "H2_rr_s",
                                                              "S","S0", "S1", "S2", "S11", "S12",
                                                              "penetrance", "dis", "H1_f_count", "H2_f_count", 
                                                              "H1_m_count", "H2_m_count","T_stat")))
    n_family_count <- 1
      while(n_family_count <= n_family) {
      H1_f <- rbinom(1,1,f_father) #if the first haplotype carries risk raviant
      H2_f <- rbinom(1,1,f_father) #if the second haplotype carries risk raviant
      H1_rr_f <- ifelse(H1_f == 1, rgamma(1, alpha, scale=beta) + adj, 1) #RR of the first haplotype
      H2_rr_f <- ifelse(H2_f == 1, rgamma(1, alpha, scale=beta) + adj, 1) #RR of the second haplotype
      H1_m <- rbinom(1,1,f_mother) #if the first haplotype carries risk raviant
      H2_m <- rbinom(1,1,f_mother) #if the second haplotype carries risk raviant
      H1_rr_m <- ifelse(H1_m == 1, rgamma(1, alpha, scale=beta) + adj, 1) #RR of the first haplotype
      H2_rr_m <- ifelse(H2_m == 1, rgamma(1, alpha, scale=beta) + adj, 1) #RR of the second haplotype
      
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
      
      #rm(H1_f, H2_f, H1_rr_f, H2_rr_f, H1_m, H2_m, H1_rr_m, H2_rr_m) #save memory
      
      S0 <- (H1_inh_s != H1_inh) & (H2_inh_s != H2_inh)
      S2 <- (H1_inh_s == H1_inh) & (H2_inh_s == H2_inh)
      S1 <- ((H1_inh_s == H1_inh) & (H2_inh_s != H2_inh)) | ((H1_inh_s != H1_inh) & (H2_inh_s == H2_inh))
      S11 <- ((H1_inh_s == H1_inh) & (H2_inh_s != H2_inh))
      S12 <- ((H1_inh_s != H1_inh) & (H2_inh_s == H2_inh))
      
      S <- (S0==TRUE)*0 + (S1==TRUE)*1 + (S2==TRUE)*2
      
      #modified to unaffected father, affected mother and both affected children
      penetrance <- H1_rr*H2_rr*H1_rr_s*H2_rr_s*KGKGR #penetrance of disease given haplotypes of both siblings
      penetrance <- ifelse(penetrance>1, 1, penetrance)
      dis <- rbinom(1,1,penetrance) #disease status of both affected
      
      if(dis==1) {
        #print(n_family_count)
        #count the T statistics
        T_stat <- (H1_f==1)*((H1_inh_s==1)+(H1_inh==1)) + (H2_f==1)*((H1_inh_s==2)+(H1_inh==2)) + (H1_m==1)*((H2_inh_s==1)+(H2_inh==1)) + (H2_m==1)*((H2_inh_s==2)+(H2_inh==2)) + H1_m + H2_m
        H1_f_count <- ((H1_inh_s==1)+(H1_inh==1))
        H2_f_count <- ((H1_inh_s==2)+(H1_inh==2))
        H1_m_count <- ((H2_inh_s==1)+(H2_inh==1)) + 1
        H2_m_count <- ((H2_inh_s==2)+(H2_inh==2)) + 1
        
        data_family[n_family_count,] <- c(H1_f, H1_rr_f, H2_f, H2_rr_f, 
                                   H1_m, H1_rr_m, H2_m, H2_rr_m,
                                   H1, H1_rr, H2, H2_rr,
                                   H1_s, H1_rr_s, H2_s, H2_rr_s,
                                   S,S0, S1, S2, S11, S12,
                                   penetrance, dis, H1_f_count, H2_f_count, H1_m_count, H2_m_count, T_stat)
        n_family_count <- n_family_count + 1
      }
    }
  }
  
  ##generate families with three generations
  if(structure=="3g") {
    data_family <- array(NA, c(n_family, 34), dimnames=list(NULL,c("H1_f", "H1_rr_f", "H2_f", "H2_rr_f", 
                                                                   "H1_m", "H1_rr_m", "H2_m", "H2_rr_m",
                                                                   "H1", "H1_rr", "H2", "H2_rr",
                                                                   "H1_s", "H1_rr_s", "H2_s", "H2_rr_s",
                                                                   "H1_3rd_inh", "H1_3rd", "H1_rr_3rd", "H2_3rd", "H2_rr_3rd",
                                                                   "S","S0", "S1", "S2", "S11", "S12",
                                                                   "penetrance", "dis", "H1_f_count", "H2_f_count", 
                                                                   "H1_m_count", "H2_m_count","T_stat")))
    n_family_count <- 1
    while(n_family_count <= n_family) {
      H1_f <- rbinom(1,1,f_father) #if the first haplotype carries risk raviant
      H2_f <- rbinom(1,1,f_father) #if the second haplotype carries risk raviant
      H1_rr_f <- ifelse(H1_f == 1, rgamma(1, alpha, scale=beta) + adj, 1) #RR of the first haplotype
      H2_rr_f <- ifelse(H2_f == 1, rgamma(1, alpha, scale=beta) + adj, 1) #RR of the second haplotype
      H1_m <- rbinom(1,1,f_mother) #if the first haplotype carries risk raviant
      H2_m <- rbinom(1,1,f_mother) #if the second haplotype carries risk raviant
      H1_rr_m <- ifelse(H1_m == 1, rgamma(1, alpha, scale=beta) + adj, 1) #RR of the first haplotype
      H2_rr_m <- ifelse(H2_m == 1, rgamma(1, alpha, scale=beta) + adj, 1) #RR of the second haplotype
      
      #the second generation
      H1_inh_s <- sample(c(1,2),1, replace=TRUE) 
      H1_s <- ifelse(H1_inh_s==1, H1_f, H2_f)
      H1_rr_s <- ifelse(H1_inh_s==1, H1_rr_f, H2_rr_f)
      
      H2_inh_s <- sample(c(1,2),1, replace=TRUE) 
      H2_s <- ifelse(H2_inh_s==1, H1_m, H2_m)
      H2_rr_s <- ifelse(H2_inh_s==1, H1_rr_m, H2_rr_m)
      
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
      
      H1_3rd_inh <- sample(c(1,2), 1, replace=TRUE)
      H1_3rd <- ifelse(H1_3rd_inh==1, H1_s, H2_s)
      H1_rr_3rd <- ifelse(H1_3rd_inh==1, H1_rr_s, H2_rr_s)
      H2_3rd <- rbinom(1,1, EH_comp(mu=m, sigma2=var, f=f)[1]) #assume affected
      H2_rr_3rd <- ifelse(H2_3rd==1, rgamma(1, alpha, scale=beta) + adj, 1)
      S <- (S0==TRUE)*0 + (S1==TRUE)*1 + (S2==TRUE)*2
      
      #modified to unaffected father, affected mother, two affected second generation children and one affected third generation
      penetrance <- H1_rr*H2_rr*H1_rr_s*H2_rr_s*KGKGR*H1_rr_3rd #penetrance of disease given haplotypes of both siblings
      penetrance <- ifelse(penetrance>1, 1, penetrance)
#       penetrance_3rd <- H1_rr_3rd*H2_rr_3rd*KG
#       penetrance_3rd <- ifelse(penetrance_3rd>1, 1, penetrance_3rd)
#       print(c(penetrance, penetrance_3rd, penetrance*penetrance_3rd))
      dis <- rbinom(1,1,penetrance) #disease status of both affected second-generation siblings and an affected third-generation offspring
      
      if(dis==1) {
        #print(n_family_count)
        #count the inherited haplotypes of affected
        H1_f_count <- (H1_inh_s==1)+(H1_inh==1)+(H1_inh_s==1)*(H1_3rd_inh==1)
        H2_f_count <- (H1_inh_s==2)+(H1_inh==2)+(H1_inh_s==2)*(H1_3rd_inh==1)
        H1_m_count <- (H2_inh_s==1)+(H2_inh==1)+(H2_inh_s==1)*(H1_3rd_inh==2) + 1
        H2_m_count <- (H2_inh_s==2)+(H2_inh==2)+(H2_inh_s==2)*(H1_3rd_inh==2) + 1
        #count the T statistics
        T_stat <- (H1_f==1)*H1_f_count + (H2_f==1)*H2_f_count + (H1_m==1)*H1_m_count + (H2_m==1)*H2_m_count
        
        data_family[n_family_count,] <- c(H1_f, H1_rr_f, H2_f, H2_rr_f, 
                                          H1_m, H1_rr_m, H2_m, H2_rr_m,
                                          H1, H1_rr, H2, H2_rr,
                                          H1_s, H1_rr_s, H2_s, H2_rr_s,
                                          H1_3rd_inh, H1_3rd, H1_rr_3rd, H2_3rd, H2_rr_3rd,
                                          S,S0, S1, S2, S11, S12,
                                          penetrance, dis, H1_f_count, H2_f_count, H1_m_count, H2_m_count, T_stat)
        n_family_count <- n_family_count + 1
      }
    }
  }
  
  ##generate families with three generations
  if(structure=="2g3c") {
    data_family <- array(NA, c(n_family, 34), dimnames=list(NULL,c("H1_f", "H1_rr_f", "H2_f", "H2_rr_f", 
                                                                   "H1_m", "H1_rr_m", "H2_m", "H2_rr_m",
                                                                   "H1", "H1_rr", "H2", "H2_rr",
                                                                   "H1_s", "H1_rr_s", "H2_s", "H2_rr_s",
                                                                   "H1_3rdsib_inh", "H1_3rdsib", "H1_rr_3rdsib", "H2_3rdsib", "H2_rr_3rdsib",
                                                                   "S","S0", "S1", "S2", "S11", "S12",
                                                                   "penetrance", "dis", "H1_f_count", "H2_f_count", 
                                                                   "H1_m_count", "H2_m_count","T_stat")))
    n_family_count <- 1
    while(n_family_count <= n_family) {
      H1_f <- rbinom(1,1,f_father) #if the first haplotype carries risk raviant
      H2_f <- rbinom(1,1,f_father) #if the second haplotype carries risk raviant
      H1_rr_f <- ifelse(H1_f == 1, rgamma(1, alpha, scale=beta) + adj, 1) #RR of the first haplotype
      H2_rr_f <- ifelse(H2_f == 1, rgamma(1, alpha, scale=beta) + adj, 1) #RR of the second haplotype
      H1_m <- rbinom(1,1,f_mother) #if the first haplotype carries risk raviant
      H2_m <- rbinom(1,1,f_mother) #if the second haplotype carries risk raviant
      H1_rr_m <- ifelse(H1_m == 1, rgamma(1, alpha, scale=beta) + adj, 1) #RR of the first haplotype
      H2_rr_m <- ifelse(H2_m == 1, rgamma(1, alpha, scale=beta) + adj, 1) #RR of the second haplotype
      
      #the second generation
      H1_inh_s <- sample(c(1,2),1, replace=TRUE) 
      H1_s <- ifelse(H1_inh_s==1, H1_f, H2_f)
      H1_rr_s <- ifelse(H1_inh_s==1, H1_rr_f, H2_rr_f)
      
      H2_inh_s <- sample(c(1,2),1, replace=TRUE) 
      H2_s <- ifelse(H2_inh_s==1, H1_m, H2_m)
      H2_rr_s <- ifelse(H2_inh_s==1, H1_rr_m, H2_rr_m)
      
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
      
      H1_3rdsib_inh <- sample(c(1,2), 1, replace=TRUE)
      H1_3rdsib <- ifelse(H1_3rdsib_inh==1, H1_f, H2_f)
      H1_rr_3rdsib <- ifelse(H1_3rdsib_inh==1, H1_rr_f, H2_rr_f)
      H2_3rdsib_inh <- sample(c(1,2), 1, replace=TRUE)
      H2_3rdsib <- ifelse(H2_3rdsib_inh==1, H1_m, H2_m)
      H2_rr_3rdsib <- ifelse(H2_3rdsib_inh==1, H1_rr_m, H2_rr_m)
      S <- (S0==TRUE)*0 + (S1==TRUE)*1 + (S2==TRUE)*2
      
      #modified to unaffected father, affected mother, three affected second generation children
      penetrance <- H1_rr*H2_rr*H1_rr_s*H2_rr_s*H1_rr_3rdsib*H2_rr_3rdsib*p_dis^3*SRR #penetrance of disease given haplotypes of both siblings
      penetrance <- ifelse(penetrance>1, 1, penetrance)
      #       penetrance_3rd <- H1_rr_3rd*H2_rr_3rd*KG
      #       penetrance_3rd <- ifelse(penetrance_3rd>1, 1, penetrance_3rd)
      #       print(c(penetrance, penetrance_3rd, penetrance*penetrance_3rd))
      dis <- rbinom(1,1,penetrance) #disease status of three affected second-generation siblings
      
      if(dis==1) {
        #print(n_family_count)
        #count the inherited haplotypes of affected
        H1_f_count <- (H1_inh_s==1)+(H1_inh==1)+(H1_3rdsib_inh==1)
        H2_f_count <- (H1_inh_s==2)+(H1_inh==2)+(H1_3rdsib_inh==2)
        H1_m_count <- (H2_inh_s==1)+(H2_inh==1)+(H2_3rdsib_inh==1) + 1
        H2_m_count <- (H2_inh_s==2)+(H2_inh==2)+(H2_3rdsib_inh==2) + 1
        #count the T statistics
        T_stat <- (H1_f==1)*H1_f_count + (H2_f==1)*H2_f_count + (H1_m==1)*H1_m_count + (H2_m==1)*H2_m_count
        
        data_family[n_family_count,] <- c(H1_f, H1_rr_f, H2_f, H2_rr_f, 
                                          H1_m, H1_rr_m, H2_m, H2_rr_m,
                                          H1, H1_rr, H2, H2_rr,
                                          H1_s, H1_rr_s, H2_s, H2_rr_s,
                                          H1_3rdsib_inh, H1_3rdsib, H1_rr_3rdsib, H2_3rdsib, H2_rr_3rdsib,
                                          S,S0, S1, S2, S11, S12,
                                          penetrance, dis, H1_f_count, H2_f_count, H1_m_count, H2_m_count, T_stat)
        n_family_count <- n_family_count + 1
      }
    }
  }  

  as.data.frame(data_family)
}
#test <- gene.data.family(m=1, var=0, f=0.01, SRR=5, n_family=1, structure="sibpair")
#test <- gene.data.family(m=1, var=0, f=0.01, SRR=5, n_family=1, structure="3rd")
#test <- gene.data.family(m=1, var=0, f=0.01, SRR=5, n_family=1, structure="three2nd")

exp_var <- function(data, delta=1) {
  e <- 0
  v <- 0
  attach(data)
  if(H1_f+H2_f+H1_m+H2_m==1) {
    e <- (H1_f_count+H2_f_count+H1_m_count+H2_m_count)/4
    v <- ((H1_f_count-e)^2+(H2_f_count-e)^2+(H1_m_count-e)^2+(H2_m_count-e)^2)/4
    moment <- (abs(H1_f_count-e)^(2+delta)+abs(H2_f_count-e)^(2+delta)+abs(H1_m_count-e)^(2+delta)+abs(H2_m_count-e)^(2+delta))/4
  }
  
  if(H1_f+H2_f+H1_m+H2_m==2) {
    e <- ((H1_f_count+H2_f_count)+(H1_f_count+H1_m_count)+(H1_f_count+H2_m_count)+(H2_f_count+H1_m_count)+
          (H2_f_count+H2_m_count)+(H1_m_count+H2_m_count))/6
    v <- (((H1_f_count+H2_f_count)-e)^2+((H1_f_count+H1_m_count)-e)^2+((H1_f_count+H2_m_count)-e)^2+
            ((H2_f_count+H1_m_count)-e)^2+((H2_f_count+H2_m_count)-e)^2+((H1_m_count+H2_m_count)-e)^2)/6
    moment <- (abs((H1_f_count+H2_f_count)-e)^(2+delta)+abs((H1_f_count+H1_m_count)-e)^(2+delta)+abs((H1_f_count+H2_m_count)-e)^(2+delta)+
                 abs((H2_f_count+H1_m_count)-e)^(2+delta)+abs((H2_f_count+H2_m_count)-e)^(2+delta)+abs((H1_m_count+H2_m_count)-e)^(2+delta))/6
  }
  if(H1_f+H2_f+H1_m+H2_m==3) {
    e <- ((H1_f_count+H2_f_count+H1_m_count)+(H1_f_count+H2_f_count+H2_m_count)+(H1_f_count+H1_m_count+H2_m_count)+(H2_f_count+H1_m_count+H2_m_count))/4
    v <- (((H1_f_count+H2_f_count+H1_m_count)-e)^2+((H1_f_count+H2_f_count+H2_m_count)-e)^2+
            ((H1_f_count+H1_m_count+H2_m_count)-e)^2+((H2_f_count+H1_m_count+H2_m_count)-e)^2)/4
    moment <- (abs((H1_f_count+H2_f_count+H1_m_count)-e)^(2+delta)+abs((H1_f_count+H2_f_count+H2_m_count)-e)^(2+delta)+
                 abs((H1_f_count+H1_m_count+H2_m_count)-e)^(2+delta)+abs((H2_f_count+H1_m_count+H2_m_count)-e)^(2+delta))/4
  }
  detach(data)
  return(c(e=e, v=v, moment=moment))
}

sim_family_new <- function(m=1, var=0, f=0.01, SRR=5, n_family=500, p_dis=0.01, rep=10, delta=1, structure) {
  print(paste("Generating Data...m=", m, "var=", var, "f=", f))
  mark <<- 1 #reset the current iteration
  print(mark)  #start iteration
  result <- replicate(rep, {
    family.data <- gene.data.family(m=m, var=var, f=f, SRR=SRR, n_family=n_family, p_dis=p_dis, structure=structure)
    attach(family.data)
    idx <- which(!((H1_f==0 & H2_f==0 & H1_m==0 & H2_m==0) | (H1_f==1 & H2_f==1 & H1_m==1 & H2_m==1))) #keep only informative families
    detach()
    if(mark%%10 == 0) print(mark) #print every 1000 iterations
    mark <<- mark + 1
    
    if(length(idx)>=1) {
      test.data <- family.data[idx,]
      n_test.data <- nrow(test.data)
#             for(i in 1:n_test.data) {
#               print((exp_var(test.data[i,])["e"]))
#               e_sum <- e_sum + (exp_var(test.data[i,])["e"]) #sum of expectation for every family
#             }
      temp <- sapply(1:n_test.data, function(x) exp_var(test.data[x,], delta))
      e <- temp["e", ]
      v <- temp["v", ]
      moment <- temp["moment", ]
      v_sq <- sum(v)
      #e_avg <- mean(e) #overall expectation
      T <- sum(test.data[, "T_stat"] - e)/sqrt(v_sq)
      Lyapunov <- sum(moment)/v_sq^((2+delta)/2)
      
      #c(t, n_test.data)#T and number of informative families
      p.value <- 2*pnorm(abs(T), lower=F) #p-value
      c(final.test.stat=T, sum_e=sum(e), se=sqrt(v_sq), mean_observed=mean(test.data[, "T_stat"]), mean_mean=mean(e), mean_var=mean(v), p.value=p.value, n_info_family=n_test.data)
#       c(sum(e), sqrt(v_sq), T, p.value, Lyapunov)  old output
    }
  }
  )
  result
}
#sim_family_new(m=1, var=0, f=0.01, SRR=5, n_family=500, p_dis=0.01, rep=2, structure="sibpair")
#sim_family_new(m=1, var=0, f=0.01, SRR=5, n_family=500, p_dis=0.01, rep=2, structure="3rd")
#sim_family_new(m=1, var=0, f=0.01, SRR=5, n_family=500, p_dis=0.01, rep=2, structure="three2nd")

#calculate only 0,1,2,3,4 variant in founders and dump 0 and 4, fix the sibling affected status
lyapunov <- function(affect=c(c(0,1),c(1,1)), n_family=500, f=0.01, mu=1, sigma2=0, p_dis=0.01, delta=1, alpha=10^-6, stat=F) {
  #I go over the possible parents' genotypes and the combination of children's genotypes
  #but not the affected status, this effect will be included in the combination of children genotyps
  
  #calculate KG
  KL <- (1+f*(mu-1))^2
  KG=p_dis/KL
  
  P_A_F <- function(F, A, mu) { #prob of being affected for founder
    if(is.na(A)) return(1)
    if(A==1) return(mu^F)
    if(A==0) return((1-mu^F*KG))
  }
  
  P_AA_H <- function(HS1, HS2, AA, mu, sharing_IBD, H_sib) { #prob of being double affected siblings
    if(sum(AA)==2) {
      if(sharing_IBD==0)  result <- mu^HS1
      if(sharing_IBD==1)  result <- mu^HS1*(mu^2+sigma2)^HS2
      if(sharing_IBD==2)  result <- (mu^2+sigma2)^HS2
    }
    if(sum(AA)==1) {
      if(sharing_IBD==0)  {
        #         if(HS1==0) result <- 1
        #         if(HS1==1) result <- 0.5*mu+0.5*(1-mu*KG)
        #         if(HS1==2) result <- 1/6*mu^2 + 4/6*mu*(1-mu*KG) + 1/6*(1-mu*KG)^2
        #         if(HS1==3) result <- 2/4*mu^2*(1-mu*KG) + 2/4*mu*(1-mu*KG)^2
        #         if(HS1==4) result <- mu^2*(1-mu*KG)^2
        result <- prod(ifelse(rep(AA, each=2)==1, mu, (1-mu^H_sib*KG))) #multiple each haplotype's effect which determined by the affected status
      }
      if(sharing_IBD==1)  {
        if(HS1==0) result <- (mu^HS2*(1-mu^HS2*KG))
        #         if(HS1==1) result <- 0.5*mu*(mu*(1-mu*KG))^HS2 + 0.5*(1-mu*KG)*(mu*(1-mu*KG))^HS2
        if(HS1==1) result <- ifelse(AA[1]==1, mu^(1+HS2)*(1-mu^(HS2)*KG), mu^(HS2)*(1-mu^(1+HS2)*KG)) #see if the affected carries the non-shared variant
        if(HS1==2) result <- mu^(1+HS2)*(1-mu^(1+HS2)*KG)
      }
      if(sharing_IBD==2)  {
        result <- (mu^HS2*(1-mu^HS2*KG))
      }
    }
    if(sum(AA)==0) {
      if(sharing_IBD==0)  result <- (1-mu^HS1*KG*KG)
      if(sharing_IBD==1)  result <- (1-mu^HS1*(mu^2)^HS2*KG*KG)
      if(sharing_IBD==2)  result <- (1-(mu^2)^HS2*KG*KG)
    }
    result
  }
  
  #calculate expectation and variance given sharing status of chromosomes and no. of carrier haplotype in founders
  exp_var <- function(sce=c(0,0,0,1), IBD_Str=c(1,3,2,4), affect=c(c(0,1),c(1,1)), delta=1) {
    n_variant = sum(sce)
    if(n_variant==0 | n_variant==4) return(c(0,0))
    if(n_variant==1) {
      founder <- matrix(c(1,0,0,0,
                          0,1,0,0,
                          0,0,1,0,
                          0,0,0,1), byrow=T, ncol=4)
    }
    if(n_variant==2) {
      founder <- matrix(c(1,1,0,0,
                          1,0,1,0,
                          1,0,0,1,
                          0,1,1,0,
                          0,1,0,1,
                          0,0,1,1), byrow=T, ncol=4)
    }
    if(n_variant==3) {
      founder <- matrix(c(1,1,1,0,
                          1,1,0,1,
                          1,0,1,1,
                          0,1,1,1), byrow=T, ncol=4)
    }
    S <- apply(founder, 1, function(x) {
      sce <- as.vector(x) #founder's haplotype
      h1 <- sce[1]
      h2 <- sce[2]
      h3 <- sce[3]
      h4 <- sce[4]        
      H_sib <- sce[IBD_Str] #variant or not on each haplotype
      S <- h1*affect[1] + h2*affect[1] + h3*affect[2] + h4*affect[2] + sum(H_sib*rep(affect[3:4], each=2)) + ifelse(is.na(affect[5]), 0, H_sib[4]*affect[5])
      }
    )
    #print(S)
    #print(c(mean(S), sum((S-mean(S))^2)/nrow(founder), sum(abs(S-mean(S))^(2+delta))/nrow(founder)))
    return(c(mean(S), sum((S-mean(S))^2)/nrow(founder), sum(abs(S-mean(S))^(2+delta))/nrow(founder)))
  }
  
  #affected status
  affect <- affect  
  
  #initialization
  p_A_sum=0 #sum of prob. under the alternative
  p_N_sum=0 #sum of prob. under the null
  p_A_n0_sum=0 #sum of prob. under the alternative non 0 in founders
  p_N_n0_sum=0 #sum of prob. under the null non 0 in founders
  Sn_sq_A=0 #pooled variance
  Sn_sq_N=0 #pooled variance
  numerator_A=0 #numerator of Lyapunov condition
  numerator_N=0 #numerator of Lyapunov condition
  T_A=0 #numerator of test statistics
  T_N=0 #numerator of test statistics
  S_A=0 #numerator of observed count
  S_N=0 #numerator of observed count
  sum_A=0 # expectation
  sum_N=0 # expectation
  debug=0 #
  p_IBD = c(0.25, 0.5, 0.25) #PROB.of IBD

  for(IBD in 0:2) {
    #create founders' haplotypes
    n_founder=4
    result <- sapply(0:n_founder, function(y) {
      combn(n_founder, y, function(x) {
        temp <- rep(0,n_founder)
        sapply(x, function(x) {
          temp[x] <<- 1
        } )
        temp
      })
    })  
    founder <- matrix(unlist(result), ncol=n_founder, byrow=T)
    for(h in 1:nrow(founder)) {
      #initialization
      sce <- as.vector(founder[h,1:4]) #founder's haplotype
      h1 <- sce[1]
      h2 <- sce[2]
      h3 <- sce[3]
      h4 <- sce[4]
      p_A <- 0
      p_N <- 0
      S <- NA
      
      HS1 <- HS2 <- 0 #shared variant
      F1 <- sum(sce[1:2]) #founder: father
      F2 <- sum(sce[3:4]) #founder: mothrt
      FB <- sum(sce[1:4]) #founder's haplotype
      p_f <- f^(FB)*(1-f)^(4-FB) #allele frequency of founders
      
      if(IBD==0) {
        IBD_Str=matrix(c(1,2,3,4,
                         1,2,4,3,
                         1,3,2,4,
                         1,3,4,2,
                         1,4,2,3,
                         1,4,3,2,
                         2,1,3,4,
                         2,1,4,3,
                         2,3,1,4,
                         2,3,4,1,
                         2,4,1,3,
                         2,4,3,1,
                         3,1,2,4,
                         3,1,4,2,
                         3,2,1,4,
                         3,2,4,1,
                         3,4,1,2,
                         3,4,2,1,
                         4,1,2,3,
                         4,1,3,2,
                         4,2,1,3,
                         4,2,3,1,
                         4,3,1,2,
                         4,3,2,1), byrow=T, ncol=4)
        for(i in 1:nrow(IBD_Str)) {
          HS1 <- FB
          HS2 <- 0
          H_sib <- sce[IBD_Str[i,]] #variant or not on each haplotype
          
          p_A <- P_A_F(F1, affect[1], mu)*P_A_F(F2, affect[2], mu)*P_AA_H(HS1,HS2,AA=affect[3:4], mu,IBD,H_sib)*p_f*p_IBD[IBD+1]*1/nrow(IBD_Str)*P_A_F(H_sib[4], affect[5], mu) #(1/2)^HS1
          p_N <- P_A_F(F1, affect[1], 1)*P_A_F(F2, affect[2], 1)*P_AA_H(HS1,HS2,AA=affect[3:4], 1,IBD,H_sib)*p_f*p_IBD[IBD+1]*1/nrow(IBD_Str)*P_A_F(H_sib[4], affect[5], 1)
          
          result <- exp_var(sce=sce, IBD_Str=IBD_Str[i,], affect=affect, delta=delta)
          S <- h1*affect[1] + h2*affect[1] + h3*affect[2] + h4*affect[2] + sum(H_sib*rep(affect[3:4], each=2)) + ifelse(is.na(affect[5]), 0, H_sib[4]*affect[5])
          #p_A <- P_A_F(F2, mu)*P_AA_H(HS1,HS2, mu) #(1/2)^HS1
          #p_N <- P_A_F(F2, 1)*P_AA_H(HS1,HS2, 1)
          p_A_sum <- p_A_sum+p_A
          p_N_sum <- p_N_sum+p_N
          if(h!=1 & h!=16) {
            Sn_sq_A <- Sn_sq_A + result[2]*p_A*n_family
            Sn_sq_N <- Sn_sq_N + result[2]*p_N*n_family

            numerator_A <- numerator_A + result[3]*p_A*n_family
            numerator_N <- numerator_N + result[3]*p_N*n_family
            
            T_A <- T_A + (S-result[1])*p_A*n_family
            T_N <- T_N + (S-result[1])*p_N*n_family

            S_A <- S_A + S*p_A
            S_N <- S_N + S*p_N
            
            sum_A <- sum_A + result[1]*p_A*n_family
            sum_N <- sum_N + result[1]*p_N*n_family
            
#             debug <- debug + (S-result[1])*p_N*n_family
#             print(c(S,result[1], S-result[1], debug, p_N, P_A_F(F1, affect[1], 1),P_A_F(F2, affect[2], 1),P_AA_H(HS1,HS2,AA=affect[3:4], 1,IBD,H_sib)))
            
            p_A_n0_sum <- p_A_n0_sum + p_A
            p_N_n0_sum <- p_N_n0_sum + p_N
          }
        }
      }
      #print(IBD)
      if(IBD==1) {
        IBD_Str=matrix(c(1,3,1,4,
                         1,4,1,3,
                         3,1,4,1,
                         4,1,3,1,
                         2,3,2,4,
                         2,4,2,3,
                         3,2,4,2,
                         4,2,3,2,
                         3,1,3,2,
                         3,2,3,1,
                         1,3,2,3,
                         2,3,1,3,
                         4,1,4,2,
                         4,2,4,1,
                         1,4,2,4,
                         2,4,1,4), byrow=T, ncol=4)
        for(i in 1:nrow(IBD_Str)) {
          H_sib <- sce[IBD_Str[i,]] #variant or not on each haplotype
          
          HS1 <- sum(table(match(IBD_Str[i,], which(sce==1)))==1) #how many variant shared once
          HS2 <- sum(table(match(IBD_Str[i,], which(sce==1)))==2) #how many variant shared twice
          p_A <- P_A_F(F1, affect[1], mu)*P_A_F(F2, affect[2], mu)*P_AA_H(HS1,HS2,AA=affect[3:4], mu,IBD,H_sib)*p_f*p_IBD[IBD+1]*1/nrow(IBD_Str)*P_A_F(H_sib[4], affect[5], mu) #(1/2)^HS1
          p_N <- P_A_F(F1, affect[1], 1)*P_A_F(F2, affect[2], 1)*P_AA_H(HS1,HS2,AA=affect[3:4], 1,IBD,H_sib)*p_f*p_IBD[IBD+1]*1/nrow(IBD_Str)*P_A_F(H_sib[4], affect[5], 1)
          
          S <- h1*affect[1] + h2*affect[1] + h3*affect[2] + h4*affect[2] + sum(H_sib*rep(affect[3:4], each=2)) + ifelse(is.na(affect[5]), 0, H_sib[4]*affect[5])
          result <- exp_var(sce=sce, IBD_Str=IBD_Str[i,], affect=affect, delta=delta)
          p_A_sum <- p_A_sum+p_A
          p_N_sum <- p_N_sum+p_N
          if(h!=1 & h!=16) {
            Sn_sq_A <- Sn_sq_A + result[2]*p_A*n_family
            Sn_sq_N <- Sn_sq_N + result[2]*p_N*n_family
            
            numerator_A <- numerator_A + result[3]*p_A*n_family
            numerator_N <- numerator_N + result[3]*p_N*n_family
            
            T_A <- T_A + (S-result[1])*p_A*n_family
            T_N <- T_N + (S-result[1])*p_N*n_family
            
            S_A <- S_A + S*p_A
            S_N <- S_N + S*p_N
            
            sum_A <- sum_A + result[1]*p_A*n_family
            sum_N <- sum_N + result[1]*p_N*n_family

            p_A_n0_sum <- p_A_n0_sum + p_A
            p_N_n0_sum <- p_N_n0_sum + p_N
          }
        }
      }
      if(IBD==2) {
        IBD_Str=matrix(c(1,3,1,3,
                         1,4,1,4,
                         3,1,3,1,
                         4,1,4,1,
                         2,3,2,3,
                         2,4,2,4,
                         3,2,3,2,
                         4,2,4,2), byrow=T, ncol=4)
        for(i in 1:nrow(IBD_Str)) {
          H_sib <- sce[IBD_Str[i,]] #variant or not on each haplotype
          HS1 <- sum(table(match(IBD_Str[i,], which(sce==1)))==1) #how many variant shared once
          HS2 <- sum(table(match(IBD_Str[i,], which(sce==1)))==2) #how many variant shared twice
          p_A <- P_A_F(F1, affect[1], mu)*P_A_F(F2, affect[2], mu)*P_AA_H(HS1,HS2,AA=affect[3:4], mu,IBD)*p_f*p_IBD[IBD+1]*1/nrow(IBD_Str)*P_A_F(H_sib[4], affect[5], mu) #(1/2)^HS1
          p_N <- P_A_F(F1, affect[1], 1)*P_A_F(F2, affect[2], 1)*P_AA_H(HS1,HS2,AA=affect[3:4], 1,IBD)*p_f*p_IBD[IBD+1]*1/nrow(IBD_Str)*P_A_F(H_sib[4], affect[5], 1)
          
          S <- h1*affect[1] + h2*affect[1] + h3*affect[2] + h4*affect[2] + sum(H_sib*rep(affect[3:4], each=2)) + ifelse(is.na(affect[5]), 0, H_sib[4]*affect[5])
          result <- exp_var(sce=sce, IBD_Str=IBD_Str[i,], affect=affect, delta=delta)
          p_A_sum <- p_A_sum+p_A
          p_N_sum <- p_N_sum+p_N
          if(h!=1 & h!=16) {
            Sn_sq_A <- Sn_sq_A + result[2]*p_A*n_family
            Sn_sq_N <- Sn_sq_N + result[2]*p_N*n_family
            
            numerator_A <- numerator_A + result[3]*p_A*n_family
            numerator_N <- numerator_N + result[3]*p_N*n_family
            
            T_A <- T_A + (S-result[1])*p_A*n_family
            T_N <- T_N + (S-result[1])*p_N*n_family
            
            S_A <- S_A + S*p_A
            S_N <- S_N + S*p_N
            
            sum_A <- sum_A + result[1]*p_A*n_family
            sum_N <- sum_N + result[1]*p_N*n_family
                        
            p_A_n0_sum <- p_A_n0_sum + p_A
            p_N_n0_sum <- p_N_n0_sum + p_N
          }
        }
      }
    }
  }
  
  #
#   print(c(p_A_n0_sum, p_A_sum))
#   print(c(p_N_n0_sum, p_N_sum))
  Sn_sq_A <- Sn_sq_A/p_A_sum #normalized by sum
  Sn_sq_N <- Sn_sq_N/p_N_sum #normalized by sum
  numerator_A <- numerator_A/p_A_sum #normalized by sum
  numerator_N <- numerator_N/p_N_sum #normalized by sum
  
  T_A <- T_A/p_A_sum #normalized by sum
  T_N <- T_N/p_N_sum #normalized by sum

  S_A <- S_A/p_A_sum #normalized by sum
  S_N <- S_N/p_N_sum #normalized by sum

  sum_A <- sum_A/p_A_sum
  sum_N <- sum_N/p_N_sum

  p_A_n0 <- p_A_n0_sum/p_A_sum
  p_N_n0 <- p_N_n0_sum/p_N_sum
  
  if(stat==T) {
    print("sum_N, sqrt(Sn_sq_N), p_N_n0, T_N, S_N/p_N_n0, T_N/sqrt(Sn_sq_N), sum_A, sqrt(Sn_sq_A), p_A_n0, T_A, S_A/p_A_n0, T_A/sqrt(Sn_sq_A)")
    print(c(sum_N, sqrt(Sn_sq_N), p_N_n0, T_N, S_N/p_N_n0, T_N/sqrt(Sn_sq_N), sum_A, sqrt(Sn_sq_A), p_A_n0, T_A, S_A/p_A_n0, T_A/sqrt(Sn_sq_A)))
  } 

  #test statistics
     
  #c(t, n_test.data)#T and number of informative families
#   correct p.value calculation
  crit_L <- qnorm(alpha/2)
  crit_H <- qnorm(1-alpha/2)
  p.value_A <- pnorm(crit_L, T_A/sqrt(Sn_sq_A),  lower=T) + pnorm(crit_H, T_A/sqrt(Sn_sq_A),  lower=F) #p-value
  #wrong p.value calculation
#   crit_L <- qnorm(alpha/2)*sqrt(Sn_sq_N)+T_N
#   crit_H <- qnorm(1-alpha/2)*sqrt(Sn_sq_N)+T_N
#   p.value_A <- pnorm(crit_L, T_A, sqrt(Sn_sq_A),  lower=T) + pnorm(crit_H, T_A, sqrt(Sn_sq_A),  lower=F) #p-value
# #   
  p.value_N <- pnorm(crit_L, T_N/sqrt(Sn_sq_N),  lower=T) + pnorm(crit_H, T_N/sqrt(Sn_sq_N),  lower=F) #p-value

  return(c(Lya_A=1/Sn_sq_A^((2+delta)/2)*numerator_A, Lya_N=1/Sn_sq_N^((2+delta)/2)*numerator_N, p.value_A=p.value_A, p.value_N=p.value_N))
}
lyapunov(affect=c(c(0,1),c(1,1)), n_family=500, f=0.01, mu=1, p_dis=0.01, delta=1, alpha=10^-6, stat=F)


#calculate only 0,1,2,3,4 variant in founders and dump 0 and 4, fix the sibling affected status
lyapunov.conditional <- function(affect=c(c(0,1),c(1,1)), n_family=500, f=0.01, mu=1, sigma2=0, p_dis=0.01, delta=1, alpha=10^-6, stat=F) {
  #I go over the possible parents' genotypes and the combination of children's genotypes
  #but not the affected status, this effect will be included in the combination of children genotyps
  
  #calculate KG
  KL <- (1+f*(mu-1))^2
  KG=p_dis/KL
  
  P_A_F <- function(F, A, mu) { #prob of being affected for founder
    if(is.na(A)) return(1)
    if(A==1) return(mu^F)
    if(A==0) return((1-mu^F*KG))
  }
  
  P_AA_H <- function(HS1, HS2, AA, mu, sharing_IBD, H_sib) { #prob of being double affected siblings
    if(sum(AA)==2) {
      if(sharing_IBD==0)  result <- mu^HS1
      if(sharing_IBD==1)  result <- mu^HS1*(mu^2+sigma2)^HS2
      if(sharing_IBD==2)  result <- (mu^2+sigma2)^HS2
    }
    if(sum(AA)==1) {
      if(sharing_IBD==0)  {
        #         if(HS1==0) result <- 1
        #         if(HS1==1) result <- 0.5*mu+0.5*(1-mu*KG)
        #         if(HS1==2) result <- 1/6*mu^2 + 4/6*mu*(1-mu*KG) + 1/6*(1-mu*KG)^2
        #         if(HS1==3) result <- 2/4*mu^2*(1-mu*KG) + 2/4*mu*(1-mu*KG)^2
        #         if(HS1==4) result <- mu^2*(1-mu*KG)^2
        result <- prod(ifelse(rep(AA, each=2)==1, mu, (1-mu^H_sib*KG))) #multiple each haplotype's effect which determined by the affected status
      }
      if(sharing_IBD==1)  {
        if(HS1==0) result <- (mu^HS2*(1-mu^HS2*KG))
        #         if(HS1==1) result <- 0.5*mu*(mu*(1-mu*KG))^HS2 + 0.5*(1-mu*KG)*(mu*(1-mu*KG))^HS2
        if(HS1==1) result <- ifelse(AA[1]==1, mu^(1+HS2)*(1-mu^(HS2)), mu^(HS2)*(1-mu^(1+HS2))) #see if the affected carries the non-shared variant
        if(HS1==2) result <- mu^(1+HS2)*(1-mu^(1+HS2)*KG)
      }
      if(sharing_IBD==2)  {
        result <- (mu^HS2*(1-mu^HS2*KG))
      }
    }
    if(sum(AA)==0) {
      if(sharing_IBD==0)  result <- (1-mu^HS1*KG*KG)
      if(sharing_IBD==1)  result <- (1-mu^HS1*(mu^2)^HS2*KG*KG)
      if(sharing_IBD==2)  result <- (1-(mu^2)^HS2*KG*KG)
    }
    result
  }
  
  #calculate expectation and variance given sharing status of chromosomes and no. of carrier haplotype in founders
  exp_var <- function(sce=c(0,0,0,1), IBD_Str=c(1,3,2,4), affect=c(c(0,1),c(1,1)), delta=1) {
    n_variant = sum(sce)
    if(n_variant==0 | n_variant==4) return(c(0,0))
    if(n_variant==1) {
      founder <- matrix(c(1,0,0,0,
                          0,1,0,0,
                          0,0,1,0,
                          0,0,0,1), byrow=T, ncol=4)
    }
    if(n_variant==2) {
      founder <- matrix(c(1,1,0,0,
                          1,0,1,0,
                          1,0,0,1,
                          0,1,1,0,
                          0,1,0,1,
                          0,0,1,1), byrow=T, ncol=4)
    }
    if(n_variant==3) {
      founder <- matrix(c(1,1,1,0,
                          1,1,0,1,
                          1,0,1,1,
                          0,1,1,1), byrow=T, ncol=4)
    }
    S <- apply(founder, 1, function(x) {
      sce <- as.vector(x) #founder's haplotype
      h1 <- sce[1]
      h2 <- sce[2]
      h3 <- sce[3]
      h4 <- sce[4]        
      H_sib <- sce[IBD_Str] #variant or not on each haplotype
      S <- h1*affect[1] + h2*affect[1] + h3*affect[2] + h4*affect[2] + sum(H_sib*rep(affect[3:4], each=2)) + ifelse(is.na(affect[5]), 0, H_sib[4]*affect[5])
    }
    )
    #print(S)
    #print(c(mean(S), sum((S-mean(S))^2)/nrow(founder), sum(abs(S-mean(S))^(2+delta))/nrow(founder)))
    return(c(mean(S), sum((S-mean(S))^2)/nrow(founder), sum(abs(S-mean(S))^(2+delta))/nrow(founder)))
  }
  
  #affected status
  affect <- affect  
  
  #initialization
  p_A_sum=0 #sum of prob. under the alternative
  p_A_n0_sum=0 #sum of prob. under the alternative non 0 in founders
  Sn_sq_A=0 #pooled variance
  numerator_A=0 #numerator of Lyapunov condition
  T_A=0 #numerator of test statistics
  S_A=0 #numerator of observed count
  sum_A=0 # expectation
  debug=0 #
  p_IBD = c(0.25, 0.5, 0.25) #PROB.of IBD
  
  #create founders' haplotypes
  n_founder=4
  result <- sapply(0:n_founder, function(y) {
    combn(n_founder, y, function(x) {
      temp <- rep(0,n_founder)
      sapply(x, function(x) {
        temp[x] <<- 1
      } )
      temp
    })
  })  
  founder <- matrix(unlist(result), ncol=n_founder, byrow=T)
  for(h in 1:nrow(founder)) {
    #initialization
    sce <- as.vector(founder[h,1:4]) #founder's haplotype
    h1 <- sce[1]
    h2 <- sce[2]
    h3 <- sce[3]
    h4 <- sce[4]
    p_A <- 0
    p_N <- 0
    S <- NA
    
    HS1 <- HS2 <- 0 #shared variant
    F1 <- sum(sce[1:2]) #founder: father
    F2 <- sum(sce[3:4]) #founder: mothrt
    FB <- sum(sce[1:4]) #founder's haplotype
    #       p_f <- 0.05 #allele frequency of founders
    p_f <- P_A_F(F1, affect[1], mu)*P_A_F(F2, affect[2], mu)*f^(FB)*(1-f)^(4-FB) #allele frequency of founders
    #       (EH_comp(mu=mu, sigma2=sigma2, f=f)[2])^F1*(1-EH_comp(mu=mu, sigma2=sigma2, f=f)[2])^(2-F1)*EH_comp(mu=mu, sigma2=sigma2, f=f)[1]^F2*(1-EH_comp(mu=mu, sigma2=sigma2, f=f)[1])^(2-F2)
    
    p_A_sib_temp=0 #sum of prob. under the alternative
    p_A_n0_sum_temp=0 #sum of prob. under the alternative non 0 in founders
    Sn_sq_A_temp=0 #pooled variance
    numerator_A_temp=0 #numerator of Lyapunov condition
    T_A_temp=0 #numerator of test statistics
    S_A_temp=0 #numerator of observed count
    sum_A_temp=0 # expectation
    for(IBD in 0:2) {
      if(IBD==0) {
        IBD_Str=matrix(c(1,2,3,4,
                         1,2,4,3,
                         1,3,2,4,
                         1,3,4,2,
                         1,4,2,3,
                         1,4,3,2,
                         2,1,3,4,
                         2,1,4,3,
                         2,3,1,4,
                         2,3,4,1,
                         2,4,1,3,
                         2,4,3,1,
                         3,1,2,4,
                         3,1,4,2,
                         3,2,1,4,
                         3,2,4,1,
                         3,4,1,2,
                         3,4,2,1,
                         4,1,2,3,
                         4,1,3,2,
                         4,2,1,3,
                         4,2,3,1,
                         4,3,1,2,
                         4,3,2,1), byrow=T, ncol=4)
      }
      if(IBD==1) {
        IBD_Str=matrix(c(1,3,1,4,
                         1,4,1,3,
                         3,1,4,1,
                         4,1,3,1,
                         2,3,2,4,
                         2,4,2,3,
                         3,2,4,2,
                         4,2,3,2,
                         3,1,3,2,
                         3,2,3,1,
                         1,3,2,3,
                         2,3,1,3,
                         4,1,4,2,
                         4,2,4,1,
                         1,4,2,4,
                         2,4,1,4), byrow=T, ncol=4)
      }
      if(IBD==2) {
        IBD_Str=matrix(c(1,3,1,3,
                         1,4,1,4,
                         3,1,3,1,
                         4,1,4,1,
                         2,3,2,3,
                         2,4,2,4,
                         3,2,3,2,
                         4,2,4,2), byrow=T, ncol=4)
      }
      for(i in 1:nrow(IBD_Str)) {
        H_sib <- sce[IBD_Str[i,]] #variant or not on each haplotype
        HS1 <- sum(table(match(IBD_Str[i,], which(sce==1)))==1) #how many variant shared once
        HS2 <- sum(table(match(IBD_Str[i,], which(sce==1)))==2) #how many variant shared twice
        
        p_A <- p_f*p_IBD[IBD+1]*1/nrow(IBD_Str) #(1/2)^HS1
        
        result <- exp_var(sce=sce, IBD_Str=IBD_Str[i,], affect=affect, delta=delta)
        S <- h1*affect[1] + h2*affect[1] + h3*affect[2] + h4*affect[2] + sum(H_sib*rep(affect[3:4], each=2)) + ifelse(is.na(affect[5]), 0, H_sib[4]*affect[5])
        #p_A <- P_A_F(F2, mu)*P_AA_H(HS1,HS2, mu) #(1/2)^HS1
        #p_N <- P_A_F(F2, 1)*P_AA_H(HS1,HS2, 1)
        p_A_sib_temp <- p_A_sib_temp+P_AA_H(HS1,HS2,AA=affect[3:4], mu,IBD,H_sib)*p_IBD[IBD+1]*1/nrow(IBD_Str)
        p_A_sum <- p_A_sum+p_A
        
        if(h!=1 & h!=16) {
          Sn_sq_A_temp <- Sn_sq_A_temp + result[2]*p_A*P_AA_H(HS1,HS2,AA=affect[3:4], mu,IBD,H_sib)*n_family
          
          numerator_A_temp <- numerator_A_temp + result[3]*p_A*P_AA_H(HS1,HS2,AA=affect[3:4], mu,IBD,H_sib)*n_family
          
          T_A_temp <- T_A_temp + (S-result[1])*p_A*P_AA_H(HS1,HS2,AA=affect[3:4], mu,IBD,H_sib)*n_family
          
          S_A_temp <- S_A_temp + S*p_A*P_AA_H(HS1,HS2,AA=affect[3:4], mu,IBD,H_sib)*n_family
          
          sum_A_temp <- sum_A_temp + result[1]*p_A*P_AA_H(HS1,HS2,AA=affect[3:4], mu,IBD,H_sib)*n_family
          
          p_A_n0_sum <- p_A_n0_sum + p_A
        }
      }
    }
    Sn_sq_A <- Sn_sq_A + Sn_sq_A_temp/p_A_sib_temp
    numerator_A <- numerator_A + numerator_A_temp/p_A_sib_temp
    T_A <- T_A + T_A_temp/p_A_sib_temp
    S_A <- S_A + S_A_temp/p_A_sib_temp
    sum_A <- sum_A + sum_A_temp/p_A_sib_temp
  }
  
  
  #
  #   print(c(p_A_n0_sum, p_A_sum))
  #   print(c(p_N_n0_sum, p_N_sum))
  Sn_sq_A <- Sn_sq_A/p_A_sum #normalized by sum
  numerator_A <- numerator_A/p_A_sum #normalized by sum
  
  T_A <- T_A/p_A_sum #normalized by sum
  
  S_A <- S_A/p_A_sum #normalized by sum
  
  sum_A <- sum_A/p_A_sum
  
  p_A_n0 <- p_A_n0_sum/p_A_sum
  
  if(stat==T) {
    print("sum_A, sqrt(Sn_sq_A), p_A_n0, T_A, S_A/p_A_n0, T_A/sqrt(Sn_sq_A)")
    print(c(sum_A, sqrt(Sn_sq_A), p_A_n0, T_A, S_A/p_A_n0, T_A/sqrt(Sn_sq_A)))
  } 
  
  #test statistics
  
  #c(t, n_test.data)#T and number of informative families
  #   correct p.value calculation
  crit_L <- qnorm(alpha/2)
  crit_H <- qnorm(1-alpha/2)
  p.value_A <- pnorm(crit_L, T_A/sqrt(Sn_sq_A),  lower=T) + pnorm(crit_H, T_A/sqrt(Sn_sq_A),  lower=F) #p-value
  #wrong p.value calculation
  #   crit_L <- qnorm(alpha/2)*sqrt(Sn_sq_N)+T_N
  #   crit_H <- qnorm(1-alpha/2)*sqrt(Sn_sq_N)+T_N
  #   p.value_A <- pnorm(crit_L, T_A, sqrt(Sn_sq_A),  lower=T) + pnorm(crit_H, T_A, sqrt(Sn_sq_A),  lower=F) #p-value
  # #   
  
  return(c(Lya_A=1/Sn_sq_A^((2+delta)/2)*numerator_A, p.value_A=p.value_A))
}
lyapunov.conditional(affect=c(c(0,1),c(1,1)), n_family=500, f=0.01, mu=1, p_dis=0.01, delta=1, alpha=10^-6, stat=F)



##gene drop simulation given the number of childern in every generation
family_strct.2g2c <- data.frame(family=c(1,1,1,1), person=c(1,2,3,4), father=c(0,0,1,1), 
                                   mother=c(0,0,2,2), sex=c(1,2,1,1), affect=c(1,2,2,2)) #1=male, 2=female, 1=unaffected, 2=affected
#use the file format as in Merlin
#simulate the tranmission vector by generation and determine the affected status
gene_family <- function(family_strct=family_strct.2g2c, n_family=100) {
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

family.test <- function(data=family_generated, f=0.01, method=1) {
  ##hypothesis testing count the number of carrier haplotypes transmitted to the offsprings
  n_family <- max(data$data_family$family)
  n_family_member <- table(data$data_family$family)
  #check if founder's haplotype carries any variant's with f < 0.1
  if(length(f)==1 & f[1] <1) { #support allele frequency or list of snps
    snp2look.idx <-  which(snp$FREQ1 < f) # snp to look for
  } else(snp2look.idx <-  f)

  
  #start looking at each family
  test.stat <- sapply(1:n_family, function(x) {
    family.idx=x 
    current_row=sum(n_family_member[1:family.idx]) - n_family_member[family.idx]
    h1 <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),7:(6+n_snp)]
    h2 <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),-c(1:(6+n_snp))]
    #adjust here for more founders
    affect <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),"affect"]
    tran_vec <- data$tran_vec[(current_row+1):(current_row+n_family_member[family.idx]),]
    family_strct <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),1:6]
    person <- c(0,family_strct[,2]) #adding 0 to check if the chromosome if from missing founder
    
    #define who is a founder
    founder <- rep(0, length=n_family_member[family.idx])
    for(i in 1:n_family_member[family.idx]) {
      founder[i] <- ifelse((family_strct$father[i]==0 & family_strct$mother[i]==0), 1,  #full founder
                           ifelse(!(family_strct$father[i] %in% person), 0.5, #half founder from father
                                   ifelse(!(family_strct$mother[i] %in% person), -0.5, 0))) #half founder from mother
    }
    founder_idx <- which(founder!=0)
    n_founder <- sum(abs(founder))
    carrier <- as.vector(sapply(founder_idx, function(y){
      if(founder[y]==1) {
        c(ifelse(sum(h1[y, snp2look.idx]==1)>0, tran_vec[y, "h1"], NA),
          ifelse(sum(h2[y, snp2look.idx]==1)>0, tran_vec[y, "h2"], NA))
      }else{
        if(founder[y]==0.5) { #from father
          ifelse(sum(h1[y, snp2look.idx]==1)>0, tran_vec[y, "h1"], NA)
        }else{
          if(founder[y]==-0.5) { #from mother
            ifelse(sum(h2[y, snp2look.idx]==1)>0, tran_vec[y, "h2"], NA)
          }
         }
       }  
    }))
#     print(carrier)
    
#     carrier <- c( #indicator of which haplotypes is carrier
#       ifelse(sum(h1[1, snp2look.idx]==1)>0, tran_vec[1, "h1"], NA) #check first founder's h1
#       ,ifelse(sum(h2[1, snp2look.idx]==1)>0, tran_vec[1, "h2"], NA) #check first founder's h2
#       ,ifelse(sum(h1[2, snp2look.idx]==1)>0, tran_vec[2, "h1"], NA) #check second founder's h1
#       ,ifelse(sum(h2[2, snp2look.idx]==1)>0, tran_vec[2, "h2"], NA) #check second founder's h2
#     )
    n_carrier <- sum(!is.na(carrier)) #no. of carrier haplotypes
    
    if(!(n_carrier==(2*n_founder) | n_carrier==0)) { #skip families with 0 or 4 carrier haplotypes in founders
      IBD_haplotype_observed = 0
      if(method==1) { #original idea of Seb
        for(i in 1:n_family_member[family.idx]) {
          #print(c(sum(tran_vec[current_row+i, c("h1","h2")] %in% carrier),(affect[current_row+i]-1)))
            IBD_haplotype_observed <- IBD_haplotype_observed + sum(tran_vec[i, c("h1","h2")] %in% carrier)*(affect[i]-1)
        }
      }
      if(method==2) { #G's idea count carrier chromosome shared between at least two affected members
        match.result <- table(match(as.matrix(tran_vec[as.logical(affect-1),c("h1","h2")]), carrier))
        IBD_haplotype_observed <- sum(match.result[which(match.result>1)])
      }
      observed <- IBD_haplotype_observed
      
      #calculate expectation and variance
      founder <- t(combn(LETTERS[1:(2*n_founder)], n_carrier))
      S <- apply(founder, 1, function(x) {
        carrier <- x #founder's haplotype
        #       print(carrier)
        IBD_haplotype_observed = 0
        if(method==1) { #original idea of Seb
          for(i in 1:n_family_member[family.idx]) {
            #print(c(sum(tran_vec[current_row+i, c("h1","h2")] %in% carrier),(affect[current_row+i]-1)))
            IBD_haplotype_observed <- IBD_haplotype_observed + sum(tran_vec[i, c("h1","h2")] %in% carrier)*(affect[i]-1)
          }
        }
        if(method==2) { #G's idea count carrier chromosome shared between at least two affected members
          match.result <- table(match(as.matrix(tran_vec[as.logical(affect-1),c("h1","h2")]), carrier))
          IBD_haplotype_observed <- sum(match.result[which(match.result>1)])
        }
        IBD_haplotype_observed
      }
      )
      
      c(observed=observed, mean=mean(S), var=sum((S-mean(S))^2)/nrow(founder), n_carrier=n_carrier)
    }
  }
  )
  test.stat <- data.frame(do.call(rbind, test.stat))
  
  v <- test.stat$var
  se <- sqrt(sum(v))
  #e_avg <- mean(e) #overall expectation
  final.test.stat <- sum(test.stat$observed - test.stat$mean)/se
  
  #c(t, n_test.data)#T and number of informative families
  p.value <- 2*pnorm(abs(final.test.stat), lower=F) #p-value
  c(final.test.stat=final.test.stat, sum_e=sum(test.stat$mean), se=se, mean_observed=mean(test.stat$observed), mean_mean=mean(test.stat$mean), mean_var=mean(test.stat$var), p.value=p.value, n_info_family=length(test.stat$n_carrier))
}
# family.test()

gene_case_control <- function(n_case_control_pair=1000) {
  data_case_control <- matrix(NA, nrow=2*n_case_control_pair, ncol=(2+n_snp*2))
  #generate cases first
  current_row <- 1 #which row is generating
  n_case.idx <- 1 #index of how many families have been generated
  while(n_case.idx <= n_case_control_pair) {
    disease <<- 0
    while(disease!=2) {
      haplo.id <- sample(1:n_haplo, 2, replace=T)
      disease_prob <- prod(haplotype.risk[haplo.id])
      disease_prob <- ifelse(disease_prob <= 1, disease_prob, 1)
      disease <<- rbinom(1,1, prob=disease_prob) + 1
    }
    #save the haplotype file
    data_case_control[current_row, ] <- unlist(c(n_case.idx, disease, haplotype[haplo.id[1],-c(1:2)], haplotype[haplo.id[2],-c(1:2)]))
    current_row <- current_row + 1
    n_case.idx <- n_case.idx + 1
  }
  #generate controls first
  n_control.idx <- 1 #index of how many families have been generated
  while(n_control.idx <= n_case_control_pair) {
    disease <<- 0
    while(disease!=1) {
      haplo.id <- sample(1:n_haplo, 2, replace=T)
      disease_prob <- prod(haplotype.risk[haplo.id])
      disease_prob <- ifelse(disease_prob <= 1, disease_prob, 1)
      disease <<- rbinom(1,1, prob=disease_prob) + 1
    }
    #save the haplotype file
    data_case_control[current_row, ] <- unlist(c(n_control.idx, disease, haplotype[haplo.id[1],-c(1:2)], haplotype[haplo.id[2],-c(1:2)]))
    current_row <- current_row + 1
    n_control.idx <- n_control.idx + 1
  }
  
  colnames(data_case_control) <- c("id","affect", rep(paste("SNP", 1:n_snp, sep=""),2))
  return(data_case_control=data.frame(data_case_control, stringsAsFactors=F))
}
# generated_case_control <- gene_case_control()

case_control.test <- function(data=generated_case_control, f=0.01) {
  ##hypothesis testing count the number of carrier haplotypes transmitted to the offsprings
  n_case_control_pair <- nrow(data)/2
  #check if founder's haplotype carries any variant's with f < 0.1
  if(length(f)==1 & f[1] <1) { #support allele frequency or list of snps
    snp2look.idx <-  which(snp$FREQ1 < f) # snp to look for
  } else(snp2look.idx <-  f)
  
  #check if carrier haplotype
  carrier <- apply(data, 1, function(x) {
    h1 <- x[3:(2+n_snp)]
    h2 <- x[-(1:(2+n_snp))]
    #     sum(h1[snp2look.idx]==1)>0
    #     sum(h2[snp2look.idx]==1)>0
    c(h1.carrier=sum(h1[snp2look.idx]==1)>0, h2.carrier=sum(h2[snp2look.idx]==1)>0)
  })
  
  carrier.case <- sum(carrier[,1:n_case_control_pair])
  carrier.control <- sum(carrier[,(n_case_control_pair+1):(2*n_case_control_pair)])
  test.result <- prop.test(c(carrier.case, carrier.control), c(2*n_case_control_pair, 2*n_case_control_pair), correct=F) #turn off correct to avoid a conservative test 
  c(test.result$estimate[1],test.result$estimate[2], p.value=test.result$p.value)
}
# case_control.test()