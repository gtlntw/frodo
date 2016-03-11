##Theoratical allele frequency: E(HS|AAr, S=1,2) and E(HN|AAr, S=0,1) on a chromosome
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

EH_comp_unaffected <- function(mu=1, sigma2=0, f=0.01, out="A", p_dis=0.01) {
  
  #calculate KG
  KL <- (1+f*(mu-1))^2
  KG=p_dis/KL
  
  EH_U <- ((1-mu*KG)*2*(1-f)*f + 2*(1-mu^2*KG)*f^2)/(1-p_dis)/2 #allele frequency on controls
  
  mu <- mu
  sigma2 <- sigma2
  p <- f
  P_A <- (1+p*(mu-1))^2
  P_AA <- 0.25*(1+p*(mu-1))^4 + 0.5*(1+p*(mu-1))^2*(1+p*(mu^2+sigma2-1)) + 0.25*(1+p*(mu^2+sigma2-1))^2
  if(out=="A") {
    return(c(EH_A <- (mu*2*p*(1-p)+2*mu^2*p^2)/P_A/2, EH_U))
  }
  else {
    return(c(EH_AA <- (0.25*2*p*(1-p)*(mu*(1+p*(mu-1))^2) + 0.5*2*p*(1-p)*(0.5*((1-p)*mu+p*mu^2)+0.5*((1-p)*(mu^2+sigma2)+p*mu*(mu^2+sigma2))) 
                       + 0.25*2*p*(1-p)*(mu^2+sigma2) + 2*(0.25*p^2*mu^2*(1+p*(mu-1))^2+0.5*p^2*((1-p)*mu*(mu^2+sigma2)+p*mu^2*(mu^2+sigma2))+0.25*p^2*(mu^2+sigma2)^2) )/P_AA/2,EH_U))
  }         
}
EH_comp_unaffected(mu=5, sigma2=0, f=0.02, out="A", p_dis=0.01)

power_comp.sas <- function(mu=1.2, sigma2=3, f=0.01, n_pair=50, alpha=10^-6, out="A", unaffected=F, p_dis=0.01) {
  N <- 4*n_pair
  ifelse(unaffected==F, P_A <-EH_comp(mu=mu, sigma2=sigma2, f=f, out=out),
         P_A <- EH_comp_unaffected(mu=mu, sigma2=sigma2, f=f, out=out, p_dis=p_dis)) ##proportion under HA
  p1 <- P_A[1]
  p2 <- P_A[2]
  p0 <- 0
  w1 <- 0.5
  w2 <- 0.5
  
  pnorm( ((p2-p1-p0)*(N*w1*w2)^.5 -  qnorm(1-alpha/2)*((w1*p1+w2*p2)*(1-w1*p1-w2*p2))^.5)/(w2*p1*(1-p1)+w1*p2*(1-p2))^.5 ) +
    pnorm( (-(p2-p1-p0)*(N*w1*w2)^.5 -  qnorm(1-alpha/2)*((w1*p1+w2*p2)*(1-w1*p1-w2*p2))^.5)/(w2*p1*(1-p1)+w1*p2*(1-p2))^.5) 
}
power_comp.sas(mu=3.0, sigma2=0, f=0.01, n_pair=1000, alpha=10^-6, out="A")

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
#this calculation is to check why simulated power is lower than the theoretical power
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
# family_strct.2g2c <- data.frame(family=c(1,1,1,1), person=c(1,2,3,4), father=c(0,0,1,1), 
#                                    mother=c(0,0,2,2), sex=c(1,2,1,1), affect=c(1,2,2,2)) #1=male, 2=female, 1=unaffected, 2=affected
#use the file format as in Merlin
#simulate the tranmission vector by generation and determine the affected status
#this generate the whole family first and then keep only those matched input
gene_family <- function(family_strct=family_strct.2g2c, n_family=100, haplotype.risk) {
  n_family_member <- length(family_strct$person)
  data_family <- matrix(NA, nrow=n_family*n_family_member, ncol=(6+n_snp*2))
  tran_vec <- matrix(NA, nrow=n_family*n_family_member, ncol=3)
  #the strategy is generate each individual one by one and check affected status if fail then restart from the first individual
  data_family.idx <- 1 #index of the current family member being generated
  n_family.idx <- 1 #index of how many families have been generated
#   affect_idx <- which(family_strct$affect!=0) #who is affected status is not missing
  while(n_family.idx <= n_family) {
    disease_vec <- matrix(NA, nrow=n_family_member, ncol=1) #store generated affected status
    family.haplo <- matrix(NA, nrow=n_family_member, ncol=2) #store haplotype info.
    ind.idx <- 1 # which individual is generating
    while(ind.idx <= n_family_member) { #until the generated matches the input
      #generate the current individual
      if(family_strct$father[ind.idx]==0 & family_strct$mother[ind.idx]==0) { #if founder
        haplo.id <- sample(1:n_haplo, 2, replace=T)
        disease_prob <- prod(haplotype.risk[haplo.id])
        disease_prob <- ifelse(disease_prob <= 1, disease_prob, 1)
        disease_vec[ind.idx] <- rbinom(1,1, prob=disease_prob) + 1
      }
      else{ #if not founder
        haplo.id <- c(sample(family.haplo[family_strct$father[ind.idx],], 1), sample(family.haplo[family_strct$mother[ind.idx],], 1))
        disease_prob <- prod(haplotype.risk[haplo.id])
        disease_prob <- ifelse(disease_prob <= 1, disease_prob, 1)
        disease_vec[ind.idx] <- rbinom(1,1, prob=disease_prob) + 1
      }
      #check if match disease status o.w. restart 
      if(disease_vec[ind.idx] == family_strct$affect[ind.idx]) {
        #store haplotype's id
        family.haplo[ind.idx,] <- haplo.id
        
        ind.idx <- ind.idx + 1 #move to the next
      } else{
        ind.idx <- 1 #restart from the first
      }
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

##regular family test trap
family.test <- function(data=family_generated, f=risk.variant.id, method=1) {
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
#     print(family.idx)
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
    carrier <- as.vector(unlist(sapply(founder_idx, function(y){
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
    })))
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
      
      c(observed=observed, mean=mean(S), var=sum((S-mean(S))^2)/nrow(founder), n_carrier=n_carrier, family.idx=family.idx)
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

##family test in TRAFIC spirit
family.test.trafic.ext <- function(data=family_generated, f=risk.variant.id, nofounderpheno=F, test_sq=F) {
  ##hypothesis testing count the number of carrier haplotypes transmitted to the offsprings
  n_family <- max(data$data_family$family)
  n_family_member <- table(data$data_family$family)
  #check if founder's haplotype carries any variant's with f < 0.1
  if(length(f)==1 & f[1] <1) { #support allele frequency or list of snps
    snp2look.idx <-  which(snp$FREQ1 < f) # snp to look for
  } else(snp2look.idx <-  f)
  
  
  #start looking at each family
  data.family <- lapply(1:n_family, function(x) {
    family.idx=x 
    #     print(family.idx)
    current_row=sum(n_family_member[1:family.idx]) - n_family_member[family.idx]
    h1 <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),7:(6+n_snp)]
    h2 <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),-c(1:(6+n_snp))]
    #adjust here for more founders
    affect <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),"affect"]
    tran_vec <- data$tran_vec[(current_row+1):(current_row+n_family_member[family.idx]),]
    family_strct <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),1:6]
    person <- c(0,family_strct[,2]) #adding 0 to check if the chromosome if from missing founder
    
    #number of unique haplotype among affected family members
    if(nofounderpheno==F) {
      affect_id <- which(affect==2)
    }else {
      affect_id <- which(affect[-c(1:2)]==2)+2
    }
    unique_haplotype <-sort(unique(as.vector(as.matrix(tran_vec[affect_id,c("h1", "h2")]))))
    
    
    #calculate for each unique haplotype the number of affeted member and report if it is a carrier haplotype
    stat <- lapply(unique_haplotype, function(x) {
      idx <- which(tran_vec[,c("h1", "h2")]==x, arr.ind=T)[1,]
      haplotype <- switch(idx[2], "1"=h1[idx[1],], "2"=h2[idx[1],])
      carrier <- sum(haplotype[, snp2look.idx]==1)>0
      haplotype_on_affect <- sum(tran_vec[affect_id,c("h1", "h2")] == x)
      data.frame(family.idx=family.idx, x=x, haplotype_on_affect=haplotype_on_affect, carrier=carrier)
    })
    
    
    return(do.call(rbind,stat))
  })
  data.family <- do.call(rbind,data.family) #convert a list of dataframes to a dataframe  
  data.family <- transform(data.family, haplotype_on_affect_sq=haplotype_on_affect^2)
#   table(data.family$carrier, data.family$haplotype_on_affect)

  #fit a logistic regression
  if(test_sq==F) {
    glm.result <- summary(glm(carrier ~ haplotype_on_affect, family=binomial(link = "logit"), data=data.family))  
  }else{
    glm.result <- summary(glm(carrier ~ haplotype_on_affect_sq, family=binomial(link = "logit"), data=data.family))
  }
  
  
  c(p.value=glm.result$coefficients[2, "Pr(>|z|)"])
}
# family.test.trafic.ext()


##look at the allele frequency on shared(2,3,4) and nonshared chromosomes
family.maf <- function(data=family_generated, f=risk.variant.id) {
  ##hypothesis testing count the number of carrier haplotypes transmitted to the offsprings
  n_family <- max(data$data_family$family)
  n_family_member <- table(data$data_family$family)
  #check if founder's haplotype carries any variant's with f < 0.1
  if(length(f)==1 & f[1] <1) { #support allele frequency or list of snps
    snp2look.idx <-  which(snp$FREQ1 < f) # snp to look for
  } else(snp2look.idx <-  f)
  
  
  #start looking at each family
  data.family <- lapply(1:n_family, function(x) {
    family.idx=x 
    #     print(family.idx)
    current_row=sum(n_family_member[1:family.idx]) - n_family_member[family.idx]
    h1 <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),7:(6+n_snp)]
    h2 <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),-c(1:(6+n_snp))]
    #adjust here for more founders
    affect <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),"affect"]
    tran_vec <- data$tran_vec[(current_row+1):(current_row+n_family_member[family.idx]),]
    family_strct <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),1:6]
    person <- c(0,family_strct[,2]) #adding 0 to check if the chromosome if from missing founder
    
    #number of unique haplotype among affected family members
    affect_id <- which(affect==2)
    unique_haplotype <-sort(unique(as.vector(as.matrix(tran_vec[affect_id,c("h1", "h2")]))))
    
    
    #calculate for each unique haplotype the number of affeted member and report if it is a carrier haplotype
    stat <- lapply(unique_haplotype, function(x) {
      idx <- which(tran_vec[,c("h1", "h2")]==x, arr.ind=T)[1,]
      haplotype <- switch(idx[2], "1"=h1[idx[1],], "2"=h2[idx[1],])
      carrier <- sum(haplotype[, snp2look.idx]==1)>0
      haplotype_on_affect <- sum(tran_vec[affect_id,c("h1", "h2")] == x)
      data.frame(family.idx=family.idx, x=x, haplotype_on_affect=haplotype_on_affect, carrier=carrier)
    })
    
    
    return(do.call(rbind,stat))
  })
  data.family <- do.call(rbind,data.family) #convert a list of dataframes to a dataframe
  
  
  #fit a logistic regression
  result.table <- with(data.family, table(haplotype_on_affect, carrier))
  result.prop <- prop.table(result.table,1)
  row.sums <- rowSums(result.table)
  cbind(haplotype_on_affect=as.numeric(rownames(result.table)), result.table, total=row.sums, carrier.prop=result.prop[,"TRUE"])
}
# family.maf()


##test for TRAP 2g2c but do not consider founder's phenotype
##need to be fixed for three siblings cases
family.test.nofounderphenotype <- function(data=family_generated_2g2c, f=risk.variant.id, weight=F) {
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
    #     print(family.idx)
    current_row=sum(n_family_member[1:family.idx]) - n_family_member[family.idx]
    h1 <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),7:(6+n_snp)] #the first haplotype
    h2 <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),-c(1:(6+n_snp))] #the second haplotype
    #adjust here for more founders
    affect <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),"affect"]
    tran_vec <- data$tran_vec[(current_row+1):(current_row+n_family_member[family.idx]),]
    family_strct <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),1:6]
    person <- c(0,family_strct[,2]) #adding 0 to check if the chromosome if from missing founder
    
    #define who is a founder
    founder <- rep(0, length=n_family_member[family.idx])
    for(i in 1:n_family_member[family.idx]) {
      founder[i] <- ifelse((family_strct$father[i]==0 & family_strct$mother[i]==0) | 
                             (!(family_strct$father[i] %in% person) & !(family_strct$mother[i] %in% person)), 1,  #full founder
                           ifelse(!(family_strct$father[i] %in% person), 0.5, #half founder from father
                                  ifelse(!(family_strct$mother[i] %in% person), -0.5, 0))) #half founder from mother
    }
    founder_idx <- which(founder!=0)
    n_founder <- sum(abs(founder))
    carrier <- as.vector(unlist(sapply(founder_idx, function(y){
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
    })))
    n_carrier <- sum(!is.na(carrier)) #no. of carrier haplotypes
    
    if(!(n_carrier==0)) { #skip families with 0 carrier haplotypes in founder
      IBD_haplotype_observed = 0
      for(i in 3:n_family_member[family.idx]) { #do not include founders' phenotype
        #print(c(sum(tran_vec[current_row+i, c("h1","h2")] %in% carrier),(affect[current_row+i]-1)))
        IBD_haplotype_observed <- IBD_haplotype_observed + sum(tran_vec[i, c("h1","h2")] %in% carrier)*(affect[i]-1)
      }
      observed <- IBD_haplotype_observed
      
      #calculate expectation and variance
      n_unique_carrier <- sum(unique(carrier) %in% LETTERS)
      founder <- t(combn(LETTERS[1:(2*n_founder)], n_unique_carrier)) #use the number of founder's carrier haplotypes
      S <- apply(founder, 1, function(x) {
        carrier <- x #founder's haplotype
        #       print(carrier)
        IBD_haplotype_observed = 0
        for(i in 3:n_family_member[family.idx]) { #do not include founders' phenotype
          #print(c(sum(tran_vec[current_row+i, c("h1","h2")] %in% carrier),(affect[current_row+i]-1)))
          IBD_haplotype_observed <- IBD_haplotype_observed + sum(tran_vec[i, c("h1","h2")] %in% carrier)*(affect[i]-1)
        }
        IBD_haplotype_observed
      }
      )
      
      #       print(S)
#       S <- S[which(S!=0)] #only count the configuration when there is at least one variant observed in the sibpair
      c(family.idx=family.idx, observed=observed, mean=mean(S), var=sum((S-mean(S))^2)/length(S), n_carrier=n_carrier, n_unique_carrier=n_unique_carrier)
    }
  }
  )
  test.stat <- data.frame(do.call(rbind, test.stat))
  
  #use variance from each family as weight
  if(weight==F){ #regular TRAP
    v <- test.stat$var
    se <- sqrt(sum(v))
    #e_avg <- mean(e) #overall expectation
    final.test.stat <- sum(test.stat$observed - test.stat$mean)/se
    #c(t, n_test.data)#T and number of informative families
    p.value <- 2*pnorm(abs(final.test.stat), lower=F) #p-value
    #   print(c(final.test.stat=final.test.stat, sum_e=sum(test.stat$mean), se=se, mean_observed=mean(test.stat$observed), mean_mean=mean(test.stat$mean), mean_var=mean(test.stat$var), p.value=p.value, n_info_family=sum(test.stat$var!=0)))
    list(final.test.stat=final.test.stat, sum_e=sum(test.stat$mean), se=se, mean_observed=mean(test.stat$observed), mean_mean=mean(test.stat$mean), mean_var=mean(test.stat$var), p.value=p.value, n_info_family=sum(test.stat$var!=0))
  }else{
    #weighted test
    informative.idx <- which(test.stat$var!=0)
    n_informative <- length(informative.idx)
    final.test.stat <- sum(test.stat$observed[informative.idx] - test.stat$mean[informative.idx]/sqrt(test.stat$var[informative.idx]))/n_informative
    p.value <- 2*pnorm(abs(final.test.stat), sd=sqrt(1/n_informative), lower=F)
    list(final.test.stat=final.test.stat, p.value=p.value, n_info_family=n_informative)
  }

}

##test for TRAP 2c
family.test.nofounder <- function(data=family_generated_2c, f=risk.variant.id) {
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
    #     print(family.idx)
    current_row=sum(n_family_member[1:family.idx]) - n_family_member[family.idx]
    h1 <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),7:(6+n_snp)] #the first haplotype
    h2 <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),-c(1:(6+n_snp))] #the second haplotype
    #adjust here for more founders
    affect <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),"affect"]
    tran_vec <- data$tran_vec[(current_row+1):(current_row+n_family_member[family.idx]),]
    family_strct <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),1:6]
    person <- c(0,family_strct[,2]) #adding 0 to check if the chromosome if from missing founder
    
    #define who is a founder
    founder <- rep(0, length=n_family_member[family.idx])
    #!!!!!!!!!!!!!this part and the realated need to be fixed
    for(i in 1:n_family_member[family.idx]) {
      founder[i] <- ifelse((family_strct$father[i]==0 & family_strct$mother[i]==0) | 
                             (!(family_strct$father[i] %in% person) & !(family_strct$mother[i] %in% person)), 1,  #full founder
                           ifelse(!(family_strct$father[i] %in% person), 0.5, #half founder from father
                                  ifelse(!(family_strct$mother[i] %in% person), -0.5, 0))) #half founder from mother
    }
    founder_idx <- which(founder!=0)
    n_founder <- 2 #!!!!!!!!set two unknown founders for ease of testing run
    carrier <- as.vector(unlist(sapply(founder_idx, function(y){
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
    })))
    #     print(carrier)
    
    #     carrier <- c( #indicator of which haplotypes is carrier
    #       ifelse(sum(h1[1, snp2look.idx]==1)>0, tran_vec[1, "h1"], NA) #check first founder's h1
    #       ,ifelse(sum(h2[1, snp2look.idx]==1)>0, tran_vec[1, "h2"], NA) #check first founder's h2
    #       ,ifelse(sum(h1[2, snp2look.idx]==1)>0, tran_vec[2, "h1"], NA) #check second founder's h1
    #       ,ifelse(sum(h2[2, snp2look.idx]==1)>0, tran_vec[2, "h2"], NA) #check second founder's h2
    #     )
    n_carrier <- sum(!is.na(carrier)) #no. of carrier haplotypes
    
    if(!(n_carrier==0)) { #skip families with 0 carrier haplotypes in sibpair
      IBD_haplotype_observed = 0
        for(i in 1:n_family_member[family.idx]) {
          #print(c(sum(tran_vec[current_row+i, c("h1","h2")] %in% carrier),(affect[current_row+i]-1)))
          IBD_haplotype_observed <- IBD_haplotype_observed + sum(tran_vec[i, c("h1","h2")] %in% carrier)*(affect[i]-1)
        }
      observed <- IBD_haplotype_observed
      
      #calculate expectation and variance
      n_unique_carrier <- sum(unique(carrier) %in% LETTERS)
      founder <- t(combn(LETTERS[1:(2*n_founder)], n_unique_carrier)) #use minimum carrier haplotype in sibpair as an estimate for the number of founder's carrier haplotypes
      S <- apply(founder, 1, function(x) {
        carrier <- x #founder's haplotype
        #       print(carrier)
        IBD_haplotype_observed = 0
        for(i in 1:n_family_member[family.idx]) {
          #print(c(sum(tran_vec[current_row+i, c("h1","h2")] %in% carrier),(affect[current_row+i]-1)))
          IBD_haplotype_observed <- IBD_haplotype_observed + sum(tran_vec[i, c("h1","h2")] %in% carrier)*(affect[i]-1)
        }
        IBD_haplotype_observed
      }
      )
      
#       print(S)
      S <- S[which(S!=0)] #only count the configuration when there is at least one variant observed in the sibpair
      c(observed=observed, mean=mean(S), var=sum((S-mean(S))^2)/length(S), n_carrier=n_carrier, n_unique_carrier=n_unique_carrier)
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
  c(final.test.stat=final.test.stat, sum_e=sum(test.stat$mean), se=se, mean_observed=mean(test.stat$observed), mean_mean=mean(test.stat$mean), mean_var=mean(test.stat$var), p.value=p.value, n_info_family=sum(test.stat$var!=0))
}

##test for TRAP based on imputaation of the haplotype carrier status
##assuming the IBD status is known
family.test.nofounder.impute <- function(data=family_generated_3c, f=risk.variant.id, method, null.maf=NULL) {
  #check
  stopifnot(method %in% c("trap", "trafic"))
  
  ##hypothesis testing count the number of carrier haplotypes transmitted to the offsprings
  n_family <- max(data$data_family$family)
  n_family_member <- table(data$data_family$family)
  #check if founder's haplotype carries any variant's with f < 0.1
  if(length(f)==1 & f[1] <1) { #support allele frequency or list of snps
    snp2look.idx <-  which(snp$FREQ1 < f) # snp to look for
  } else(snp2look.idx <-  f)
  
  match.pattern.list <- list() #to store if the haplotype carrying a risk variant on multiple affected
  haplo.unique.count <- list()
  for(family.idx in 1:n_family) {
    current_row=sum(n_family_member[1:family.idx]) - n_family_member[family.idx]
    ##assign transmission vector based on IBD
    #create a IBD matrix n_ind*n_ind
    IBD.matrix <- matrix(NA, nrow=n_family_member[family.idx], ncol=n_family_member[family.idx])
    for(a in 1:n_family_member[family.idx]) {
      for(b in 1:n_family_member[family.idx]) {
        IBD.matrix[a,b] <- sum(data$tran_vec[a, c("h1","h2")] %in% data$tran_vec[b, c("h1","h2")])
      }
    }
    #assign a suitable transmission, not trivial to code up the algorithm assume is known for now
    tran_vec <- data$tran_vec[(current_row+1):(current_row+n_family_member[family.idx]),]
    ##tally the (un)ambiguous carrier haplotype
    h1 <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),7:(6+n_snp)] #the first haplotype
    h2 <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),-c(1:(6+n_snp))] #the second haplotype
    #observed allele count for each individual
    obersed.allele <- vector("numeric", n_family_member[family.idx])
    for(a in 1:n_family_member[family.idx]) {
      obersed.allele[a] <- (h1[a, snp2look.idx]==1) + (h2[a, snp2look.idx]==1)
    }
#     (tran_vec[, c("h1","h2")])
    #go over every possible haplotype carrier configuration with early terminatation rule
    #need to determine the number of founder haplotypes, assume know for now
    #can use dynamic programming, skip for now
    match.pattern <- NULL
    #the number of each unique haplotype occured on affecteds
    haplo.count <- table(factor(as.matrix(tran_vec[ ,c("h1","h2")]), LETTERS[1:4]))
    #save the number of each observed unique haplotype occured on affecteds for later use, i.e. can be A,B,D
    haplo.unique.count[[family.idx]] <- table(factor(as.matrix(tran_vec[ ,c("h1","h2")])))
    for(a in 0:as.numeric(haplo.count[1]>0)) { #only loop through observed haplotypes
      for(b in 0:as.numeric(haplo.count[2]>0)) {
        for(c in 0:as.numeric(haplo.count[3]>0)) {
          for(d in 0:as.numeric(haplo.count[4]>0)) {
            allele.count <- as.numeric(apply(tran_vec[ ,c("h1","h2")], 1, function(x) {
                temp <- c(0,0,0,0)
                temp[match(x, c(LETTERS[1:4]))] <- 1
                sum(c(a,b,c,d)*temp)
                }
              ))
            
#             print(c(allele.count,obersed.allele,identical(allele.count,obersed.allele)))
            result <- haplo.count * c(a,b,c,d)
            result[which(haplo.count==0)] <- -9 #replace non-obsered haplotype with -9
            if(identical(allele.count,obersed.allele)) match.pattern <- rbind(match.pattern, result)
          }
        }
      }
    }
    colnames(match.pattern) <- LETTERS[1:4]
    match.pattern.list[[family.idx]] <-  match.pattern
  }

  ##Use EM to impute the carrier haplotype
  #get ks, kn for unambuiguous families
  count.kskn <- lapply(match.pattern.list, function(x) if(nrow(x)==1) {
    return(c(ks=sum(x>1), kn=sum(x==1))) #count of kc, knc for each family
  })
  count.kskn.sum <- colSums(do.call(rbind, count.kskn))
  #total number of cn, cs
  count.cscn <- sapply(1:n_family, function(family.idx) {
    c(cs=sum(haplo.unique.count[[family.idx]]>1), cn=sum(haplo.unique.count[[family.idx]]==1))
  })
  count.cscn.sum <- rowSums(count.cscn)
  #count u and the number of kn ks cn cs in different configurations
  u.list <- list()
  for(a in 1:n_family) {
    n.row.config <- nrow(match.pattern.list[[a]])
    if(n.row.config > 1) {
      result <- NULL
      for(b in 1:n.row.config) {
        family.match.pattern <- match.pattern.list[[a]][b,] #current match pattern
        ks=sum(family.match.pattern>1); kn=sum(family.match.pattern==1) #count ks kn
        result <- rbind(result, as.vector(c(family=a, config=b, family.match.pattern, ks, kn, t(count.cscn[,a]))))
      }
      colnames(result) <- c("family", "config", "A", "B", "C", "D", "ks", "kn", "cs", "cn")
      u.list <- c(u.list ,list(as.data.frame(result)))
    }
  }
#   colnames(u.list) <- c("family", "config", "A", "B", "C", "D", "ks", "kn", "cs", "cn")
# match.pattern.list
  #assuming only one locus only for now
  EM <- function(count.kskn.sum=count.kskn.sum, count.cscn.sum=count.cscn.sum) {
    #initialization
    ks <- count.kskn.sum[1] #known shared variants  (variants occured On ibd 2 or more)
    kn <- count.kskn.sum[2] #known non-shared variants (variants occured on ibd 0)
    cs <- count.cscn.sum[1] # total number of shared chromosomes
    cn <- count.cscn.sum[2] #total number of non-shared chromosomes
          
    pn.init <- kn/(cn) #probability of rare variant on non-shared chromosome
    pn.cur <- ifelse(pn.init==0, runif(1), pn.init)
    ps.init <- ks/(cs) #probability of rare variant on shared chromosome
    ps.cur <- ifelse(ps.init==0, runif(1), ps.init)
    delta <- Inf
    iter <- 1
    
    while(delta > 10^-6) {
      #E step
      #us <- u*ps.cur/(pn.cur+ps.cur)
      #un <- u*pn.cur/(pn.cur+ps.cur)
      if(!length(u.list)==0) { #prevent the error when there's no ambiguous famlies
        for(a in 1:length(u.list)) {
          temp.prob <- with(u.list[[a]], ps.cur^ks*pn.cur^kn*(1-ps.cur)^(cs-ks)*(1-pn.cur)^(cn-kn)) #go through u.list and calculate the prob. for each config
          sum.prob <- sum(temp.prob)
          u.list[[a]]$prob <- temp.prob/sum.prob
        }
        u <- do.call(rbind, u.list)
        us <- sum(with(u, ks*prob)) #update us
        un <- sum(with(u, kn*prob)) #update un
      }else {
        u <- us <- un <- 0
      }
      #M-step
      pn.new <- (kn+un)/cn
      ps.new <- (ks+us)/cs
      #print(c(mu.new, sigma2.new, f.new, cor.factor, optim.result$value))
      
      #check convergence
      delta <- max(abs(pn.cur - pn.new), abs(ps.cur - ps.new))
      pn.cur <- pn.new
      ps.cur <- ps.new
      
#       print(c(pn.cur, ps.cur, iter))
      iter <- iter + 1
#       print(iter)
    }
    #c(pn.init, ps.init)
    c(ps.cur=ps.cur, pn.cur=pn.cur)
  }

  ##family test in TRAFIC spirit take input of imputed dataset
  family.test.trafic.ext.impute <- function(ps, pn) {
    ##hypothesis testing count the number of carrier haplotypes transmitted to the offsprings

    #extract info from unambiguous family
    data.family.list <- list()
    for(family.idx in 1:length(match.pattern.list)) {
      if(nrow(match.pattern.list[[family.idx]])==1) {
        match.pattern <- match.pattern.list[[family.idx]][which(match.pattern.list[[family.idx]]!=-9)]
        data.family.list <- c(data.family.list, list(t(rbind(haplotype_on_affect=haplo.unique.count[[family.idx]], carrier=match.pattern>=1))))
      }
    }
    
    #extract info from ambiguous family
    #update the configuration probability
    if(length(u.list)!=0) { #prevent the error when there's no ambiguous families
      for(a in 1:length(u.list)) {
        temp.prob <- with(u.list[[a]], ps^ks*pn^kn*(1-ps)^(cs-ks)*(1-pn)^(cn-kn)) #go through u.list and calculate the prob. for each config
        sum.prob <- sum(temp.prob)
        prob <- temp.prob/sum.prob #config probability
        sample.config <- sample(1:nrow(u.list[[a]]), 1, prob=prob) #which config
        match.pattern <- u.list[[a]][sample.config,c("A","B","C","D")]
        match.pattern <- match.pattern[which(match.pattern!=-9)]
        
        family.idx <- u.list[[a]]$family[1]
      
        data.family.list <- c(data.family.list, list(t(rbind(haplotype_on_affect=haplo.unique.count[[family.idx]], carrier=match.pattern>=1))))
      }
    }

    data.family <- as.data.frame(do.call(rbind, data.family.list))
    #fit a logistic regression
    glm.result <- summary(glm(carrier ~ haplotype_on_affect, family=binomial(link = "logit"), data=data.family))
    
    c(p.value=glm.result$coefficients["haplotype_on_affect", "Pr(>|z|)"], 
      est=glm.result$coefficients["haplotype_on_affect", "Estimate"],
      var=glm.result$coefficients["haplotype_on_affect", "Std. Error"]^2)
  }
  # family.test.trafic.ext.impute()

  ##test for 3c only use infor from children
  family.test.trap.impute <- function(ps, pn) {
    ##hypothesis testing count the number of carrier haplotypes transmitted to the offsprings
    #extract info from unambiguous family
    data.family.list <- list()
    for(family.idx in 1:length(match.pattern.list)) {
      if(nrow(match.pattern.list[[family.idx]])==1) {
        match.pattern <- match.pattern.list[[family.idx]][which(match.pattern.list[[family.idx]]!=-9)]
        data.family.list <- c(data.family.list, list(t(rbind(haplotype_on_affect=haplo.unique.count[[family.idx]], carrier=match.pattern>=1))))
      }
    }
    
    #extract info from ambiguous family
    #update the configuration probability
    if(length(u.list)!=0) { #prevent the error when there's no ambiguous families
      for(a in 1:length(u.list)) {
        temp.prob <- with(u.list[[a]], ps^ks*pn^kn*(1-ps)^(cs-ks)*(1-pn)^(cn-kn)) #go through u.list and calculate the prob. for each config
        sum.prob <- sum(temp.prob)
        prob <- temp.prob/sum.prob #config probability
        sample.config <- sample(1:nrow(u.list[[a]]), 1, prob=prob) #which config
        match.pattern <- u.list[[a]][sample.config,c("A","B","C","D")]
        match.pattern <- match.pattern[which(match.pattern!=-9)]
        
        family.idx <- u.list[[a]]$family[1]
        
        data.family.list <- c(data.family.list, list(t(rbind(haplotype_on_affect=haplo.unique.count[[family.idx]], carrier=match.pattern>=1))))
      }
    }
    
    #start looking at each family
    test.stat <- sapply(data.family.list, function(x) {
      n_founder <- 2
      n_carrier <- sum(x[,2]) #no. of carrier haplotypes
      if(!(n_carrier==0)) { #skip families with 0 carrier haplotypes in sibpair
        observed <- sum(x[,1] * x[, 2])
        
        #calculate expectation and variance
        n_unique_carrier <- sum(x[,2]==1)
        founder <- t(combn(LETTERS[1:(2*n_founder)], n_unique_carrier)) #use minimum carrier haplotype in sibpair as an estimate for the number of founder's carrier haplotypes
        S <- apply(founder, 1, function(y) {
          carrier <- as.numeric(rownames(x) %in% y) #founder's haplotype
          #       print(carrier)

          IBD_haplotype_observed <- sum(x[,1] * carrier)
        }
        )
        
        #       print(S)
        S <- S[which(S!=0)] #only count the configuration when there is at least one variant observed in the sibpair
        c(observed=observed, mean=mean(S), var=sum((S-mean(S))^2)/length(S), n_carrier=n_carrier, n_unique_carrier=n_unique_carrier)
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
    c(p.value=p.value, est=sum(test.stat$observed - test.stat$mean), var=se^2)
  }#family.test.trap.impute


  #choose trap or trafic to use
  impute.method <- switch(method, trap=family.test.trap.impute, trafic=family.test.trafic.ext.impute)
  #apply test with multiple imputation
  MI <- function() {
    #initialization
    diff <- NULL
    var <- NULL
    D <- 10
    #run EM to estimate the allele frequency
    if(!is.null(null.maf)) {
      pspn <- c(null.maf, null.maf)
    }else {
      pspn <- EM(count.kskn.sum=count.kskn.sum, count.cscn.sum=count.cscn.sum)  
    }
    for(i in 1:D) {
      mi.result <- impute.method(ps=pspn[1], pn=pspn[2])
      diff <- cbind(diff, mi.result[2])
      var <- cbind(var, mi.result[3])
    }
    
    TD <- mean(diff)
    VARD <- mean(var) + (1+1/D)*sum((diff-TD)^2)/(D-1)
    c(pchisq(TD^2/VARD, df=1, lower=F))
  }
  result <- MI()
  result
}


##test for TRAP based on imputaation of the haplotype carrier status
##assuming the IBD status is known
##do not aggregate the info of shared chromosomes
family.test.nofounder.impute.noaggregate <- function(data=family_generated_3c, f=risk.variant.id, method, null.maf=NULL) {
  #check
  stopifnot(method %in% c("trap", "trafic"))
  
  ##hypothesis testing count the number of carrier haplotypes transmitted to the offsprings
  n_family <- max(data$data_family$family)
  n_family_member <- table(data$data_family$family)
  #check if founder's haplotype carries any variant's with f < 0.1
  if(length(f)==1 & f[1] <1) { #support allele frequency or list of snps
    snp2look.idx <-  which(snp$FREQ1 < f) # snp to look for
  } else(snp2look.idx <-  f)
  
  match.pattern.list <- list() #to store if the haplotype carrying a risk variant on multiple affected
  haplo.unique.count <- list()
  for(family.idx in 1:n_family) {
    current_row=sum(n_family_member[1:family.idx]) - n_family_member[family.idx]
    ##assign transmission vector based on IBD
    #create a IBD matrix n_ind*n_ind
    IBD.matrix <- matrix(NA, nrow=n_family_member[family.idx], ncol=n_family_member[family.idx])
    for(a in 1:n_family_member[family.idx]) {
      for(b in 1:n_family_member[family.idx]) {
        IBD.matrix[a,b] <- sum(data$tran_vec[a, c("h1","h2")] %in% data$tran_vec[b, c("h1","h2")])
      }
    }
    #assign a suitable transmission, not trivial to code up the algorithm assume is known for now
    tran_vec <- data$tran_vec[(current_row+1):(current_row+n_family_member[family.idx]),]
    ##tally the (un)ambiguous carrier haplotype
    h1 <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),7:(6+n_snp)] #the first haplotype
    h2 <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),-c(1:(6+n_snp))] #the second haplotype
    #observed allele count for each individual
    obersed.allele <- vector("numeric", n_family_member[family.idx])
    for(a in 1:n_family_member[family.idx]) {
      obersed.allele[a] <- (h1[a, snp2look.idx]==1) + (h2[a, snp2look.idx]==1)
    }
    #     (tran_vec[, c("h1","h2")])
    #go over every possible haplotype carrier configuration with early terminatation rule
    #need to determine the number of founder haplotypes, assume know for now
    #can use dynamic programming, skip for now
    match.pattern <- NULL
    #the number of each unique haplotype occured on affecteds
    haplo.count <- table(factor(as.matrix(tran_vec[ ,c("h1","h2")]), LETTERS[1:4]))
    #save the number of each observed unique haplotype occured on affecteds for later use, i.e. can be A,B,D
    haplo.unique.count[[family.idx]] <- table(factor(as.matrix(tran_vec[ ,c("h1","h2")])))
    for(a in 0:as.numeric(haplo.count[1]>0)) { #only loop through observed haplotypes
      for(b in 0:as.numeric(haplo.count[2]>0)) {
        for(c in 0:as.numeric(haplo.count[3]>0)) {
          for(d in 0:as.numeric(haplo.count[4]>0)) {
            allele.count <- as.numeric(apply(tran_vec[ ,c("h1","h2")], 1, function(x) {
              temp <- c(0,0,0,0)
              temp[match(x, c(LETTERS[1:4]))] <- 1
              sum(c(a,b,c,d)*temp)
            }
            ))
            
            #             print(c(allele.count,obersed.allele,identical(allele.count,obersed.allele)))
            result <- haplo.count * c(a,b,c,d)
            result[which(haplo.count==0)] <- -9 #replace non-obsered haplotype with -9
            if(identical(allele.count,obersed.allele)) match.pattern <- rbind(match.pattern, result)
          }
        }
      }
    }
    colnames(match.pattern) <- LETTERS[1:4]
    match.pattern.list[[family.idx]] <-  match.pattern
  }
  
  ##Use EM to impute the carrier haplotype
  #get ks, kn for unambuiguous families
  count.kskn <- lapply(match.pattern.list, function(x) if(nrow(x)==1) {
    return(c(ks3=sum(x==3), ks2=sum(x==2), kn=sum(x==1))) #count of kc, knc for each family
  })
  count.kskn.sum <- colSums(do.call(rbind, count.kskn))
  #total number of cn, cs
  count.cscn <- sapply(1:n_family, function(family.idx) {
    c(cs3=sum(haplo.unique.count[[family.idx]]==3),
      cs2=sum(haplo.unique.count[[family.idx]]==2), cn=sum(haplo.unique.count[[family.idx]]==1))
  })
  count.cscn.sum <- rowSums(count.cscn)
  #count u and the number of kn ks cn cs in different configurations
  u.list <- list()
  for(a in 1:n_family) {
    n.row.config <- nrow(match.pattern.list[[a]])
    if(n.row.config > 1) {
      result <- NULL
      for(b in 1:n.row.config) {
        family.match.pattern <- match.pattern.list[[a]][b,] #current match pattern
        ks3=sum(family.match.pattern==3)
        ks2=sum(family.match.pattern==2)
        kn=sum(family.match.pattern==1) #count ks kn
        result <- rbind(result, as.vector(c(family=a, config=b, family.match.pattern, ks3, ks2, kn, t(count.cscn[,a]))))
      }
      colnames(result) <- c("family", "config", "A", "B", "C", "D", "ks3", "ks2", "kn", "cs3", "cs2", "cn")
      u.list <- c(u.list ,list(as.data.frame(result)))
    }
  }
  #   colnames(u.list) <- c("family", "config", "A", "B", "C", "D", "ks", "kn", "cs", "cn")
  # match.pattern.list
  #assuming only one locus only for now
  EM <- function(count.kskn.sum=count.kskn.sum, count.cscn.sum=count.cscn.sum) {
    #initialization
    ks3 <- count.kskn.sum[1] #known shared variants  (variants occured On ibd 2 or more)
    ks2 <- count.kskn.sum[2] #known shared variants  (variants occured On ibd 2 or more)
    kn <- count.kskn.sum[3] #known non-shared variants (variants occured on ibd 0)
    cs3 <- count.cscn.sum[1] # total number of shared chromosomes
    cs2 <- count.cscn.sum[2] # total number of shared chromosomes
    cn <- count.cscn.sum[3] #total number of non-shared chromosomes
    
    pn.init <- kn/(cn) #probability of rare variant on non-shared chromosome
    pn.cur <- ifelse(pn.init==0, runif(1), pn.init)
    ps3.init <- ks3/(cs3) #probability of rare variant on shared chromosome
    ps3.cur <- ifelse(ps3.init==0, runif(1), ps3.init)
    ps2.init <- ks2/(cs2) #probability of rare variant on shared chromosome
    ps2.cur <- ifelse(ps2.init==0, runif(1), ps2.init)
    delta <- Inf
    iter <- 1
    
    while(delta > 10^-6) {
      #E step
      #us <- u*ps.cur/(pn.cur+ps.cur)
      #un <- u*pn.cur/(pn.cur+ps.cur)
      if(length(u.list)!=0) { #prevent the error when there's no ambiguous families
        for(a in 1:length(u.list)) {
          temp.prob <- with(u.list[[a]], ps3.cur^ks3*ps2.cur^ks2*pn.cur^kn*(1-ps3.cur)^(cs3-ks3)*(1-ps2.cur)^(cs2-ks2)*(1-pn.cur)^(cn-kn)) #go through u.list and calculate the prob. for each config
          sum.prob <- sum(temp.prob)
          u.list[[a]]$prob <- temp.prob/sum.prob
        }
        u <- do.call(rbind, u.list)
        us3 <- sum(with(u, ks3*prob)) #update us
        us2 <- sum(with(u, ks2*prob)) #update us
        un <- sum(with(u, kn*prob)) #update un
      }else{
        u<-us3<-us2<-un<-0
      }
      #M-step
      pn.new <- (kn+un)/cn
      ps2.new <- (ks2+us2)/cs2
      ps3.new <- (ks3+us3)/cs3
      #print(c(mu.new, sigma2.new, f.new, cor.factor, optim.result$value))
      
      #check convergence
      delta <- max(abs(pn.cur - pn.new), abs(ps2.cur - ps2.new), abs(ps3.cur - ps3.new))
      pn.cur <- pn.new
      ps2.cur <- ps2.new
      ps3.cur <- ps3.new
      
      #       print(c(pn.cur, ps.cur, iter))
      iter <- iter + 1
      #       print(iter)
    }
    #c(pn.init, ps.init)
    c(ps3.cur=ps3.cur, ps2.cur=ps2.cur, pn.cur=pn.cur)
  }
  
  ##family test in TRAFIC spirit take input of imputed dataset
  family.test.trafic.ext.impute <- function(ps3, ps2, pn) {
    ##hypothesis testing count the number of carrier haplotypes transmitted to the offsprings
    
    #extract info from unambiguous family
    data.family.list <- list()
    for(family.idx in 1:length(match.pattern.list)) {
      if(nrow(match.pattern.list[[family.idx]])==1) {
        match.pattern <- match.pattern.list[[family.idx]][which(match.pattern.list[[family.idx]]!=-9)]
        data.family.list <- c(data.family.list, list(t(rbind(haplotype_on_affect=haplo.unique.count[[family.idx]], carrier=match.pattern>=1))))
      }
    }
    
    #extract info from ambiguous family
    #update the configuration probability
    if(length(u.list)!=0) { #prevent the error when there's no ambiguous families
      for(a in 1:length(u.list)) {
        temp.prob <- with(u.list[[a]], ps3^ks3*ps2^ks2*pn^kn*(1-ps3)^(cs3-ks3)*(1-ps2)^(cs2-ks2)*(1-pn)^(cn-kn)) #go through u.list and calculate the prob. for each config
        sum.prob <- sum(temp.prob)
        prob <- temp.prob/sum.prob #config probability
        sample.config <- sample(1:nrow(u.list[[a]]), 1, prob=prob) #which config
        match.pattern <- u.list[[a]][sample.config,c("A","B","C","D")]
        match.pattern <- match.pattern[which(match.pattern!=-9)]
        
        family.idx <- u.list[[a]]$family[1]
        
        data.family.list <- c(data.family.list, list(t(rbind(haplotype_on_affect=haplo.unique.count[[family.idx]], carrier=match.pattern>=1))))
      }
    }
    
    data.family <- as.data.frame(do.call(rbind, data.family.list))
    #fit a logistic regression
    glm.result <- summary(glm(carrier ~ haplotype_on_affect, family=binomial(link = "logit"), data=data.family))
    
    c(p.value=glm.result$coefficients["haplotype_on_affect", "Pr(>|z|)"], 
      est=glm.result$coefficients["haplotype_on_affect", "Estimate"],
      var=glm.result$coefficients["haplotype_on_affect", "Std. Error"]^2)
  }
  # family.test.trafic.ext.impute()
  
  ##test for 3c only use infor from children
  family.test.trap.impute <- function(ps3, ps2, pn) {
    ##hypothesis testing count the number of carrier haplotypes transmitted to the offsprings
    #extract info from unambiguous family
    data.family.list <- list()
    for(family.idx in 1:length(match.pattern.list)) {
      if(nrow(match.pattern.list[[family.idx]])==1) {
        match.pattern <- match.pattern.list[[family.idx]][which(match.pattern.list[[family.idx]]!=-9)]
        data.family.list <- c(data.family.list, list(t(rbind(haplotype_on_affect=haplo.unique.count[[family.idx]], carrier=match.pattern>=1))))
      }
    }
    
    #extract info from ambiguous family
    #update the configuration probability
    if(length(u.list)!=0) { #prevent the error when there's no ambiguous families
      for(a in 1:length(u.list)) {
        temp.prob <- with(u.list[[a]], ps3^ks3*ps2^ks2*pn^kn*(1-ps3)^(cs3-ks3)*(1-ps2)^(cs2-ks2)*(1-pn)^(cn-kn)) #go through u.list and calculate the prob. for each config
        sum.prob <- sum(temp.prob)
        prob <- temp.prob/sum.prob #config probability
        sample.config <- sample(1:nrow(u.list[[a]]), 1, prob=prob) #which config
        match.pattern <- u.list[[a]][sample.config,c("A","B","C","D")]
        match.pattern <- match.pattern[which(match.pattern!=-9)]
        
        family.idx <- u.list[[a]]$family[1]
        
        data.family.list <- c(data.family.list, list(t(rbind(haplotype_on_affect=haplo.unique.count[[family.idx]], carrier=match.pattern>=1))))
      }
    }
    
    #start looking at each family
    test.stat <- sapply(data.family.list, function(x) {
      n_founder <- 2
      n_carrier <- sum(x[,2]) #no. of carrier haplotypes
      if(!(n_carrier==0)) { #skip families with 0 carrier haplotypes in sibpair
        observed <- sum(x[,1] * x[, 2])
        
        #calculate expectation and variance
        n_unique_carrier <- sum(x[,2]==1)
        founder <- t(combn(LETTERS[1:(2*n_founder)], n_unique_carrier)) #use minimum carrier haplotype in sibpair as an estimate for the number of founder's carrier haplotypes
        S <- apply(founder, 1, function(y) {
          carrier <- as.numeric(rownames(x) %in% y) #founder's haplotype
          #       print(carrier)
          
          IBD_haplotype_observed <- sum(x[,1] * carrier)
        }
        )
        
        #       print(S)
        S <- S[which(S!=0)] #only count the configuration when there is at least one variant observed in the sibpair
        c(observed=observed, mean=mean(S), var=sum((S-mean(S))^2)/length(S), n_carrier=n_carrier, n_unique_carrier=n_unique_carrier)
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
    c(p.value=p.value, est=sum(test.stat$observed - test.stat$mean), var=se^2)
  }#family.test.trap.impute
  
  
  #choose trap or trafic to use
  impute.method <- switch(method, trap=family.test.trap.impute, trafic=family.test.trafic.ext.impute)
  #apply test with multiple imputation
  MI <- function() {
    #initialization
    diff <- NULL
    var <- NULL
    D <- 10
    #run EM to estimate the allele frequency
    if(!is.null(null.maf)) {
      pspn <- c(null.maf, null.maf, null.maf)
    }else {
      pspn <- EM(count.kskn.sum=count.kskn.sum, count.cscn.sum=count.cscn.sum)  
    }
    for(i in 1:D) {
      mi.result <- impute.method(ps3=pspn[1], ps2=pspn[2], pn=pspn[3])
      diff <- cbind(diff, mi.result[2])
      var <- cbind(var, mi.result[3])
    }
    
    TD <- mean(diff)
    VARD <- mean(var) + (1+1/D)*sum((diff-TD)^2)/(D-1)
    c(pchisq(TD^2/VARD, df=1, lower=F))
  }
  result <- MI()
  result
}#family.test.nofounder.impute.noaggregate


##test for TRAP based on imputaation of the haplotype carrier status
##assuming the IBD status is known
##do not aggregate the info of shared chromosomes
##the first half is family of three affected, and the second half is family of two affected
family.test.nofounder.impute.noaggregate.seperate <- function(data=family_generated, f=risk.variant.id, method, null.maf=NULL) {
  #check
  stopifnot(method %in% c("trap", "trafic"))
  
  ##hypothesis testing count the number of carrier haplotypes transmitted to the offsprings
  n_family <- max(data$data_family$family)
  n_family_member <- table(data$data_family$family)
  family.strct.idx <- ifelse(n_family_member==3, 1, 2)
  #check if founder's haplotype carries any variant's with f < 0.1
  if(length(f)==1 & f[1] <1) { #support allele frequency or list of snps
    snp2look.idx <-  which(snp$FREQ1 < f) # snp to look for
  } else(snp2look.idx <-  f)
  
  match.pattern.list <- list() #to store if the haplotype carrying a risk variant on multiple affected
  haplo.unique.count <- list()
  for(family.idx in 1:n_family) {
    current_row=sum(n_family_member[1:family.idx]) - n_family_member[family.idx]
    ##assign transmission vector based on IBD
    #create a IBD matrix n_ind*n_ind
    IBD.matrix <- matrix(NA, nrow=n_family_member[family.idx], ncol=n_family_member[family.idx])
    for(a in 1:n_family_member[family.idx]) {
      for(b in 1:n_family_member[family.idx]) {
        IBD.matrix[a,b] <- sum(data$tran_vec[a, c("h1","h2")] %in% data$tran_vec[b, c("h1","h2")])
      }
    }
    #assign a suitable transmission, not trivial to code up the algorithm assume is known for now
    tran_vec <- data$tran_vec[(current_row+1):(current_row+n_family_member[family.idx]),]
    ##tally the (un)ambiguous carrier haplotype
    h1 <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),7:(6+n_snp)] #the first haplotype
    h2 <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),-c(1:(6+n_snp))] #the second haplotype
    #observed allele count for each individual
    obersed.allele <- vector("numeric", n_family_member[family.idx])
    for(a in 1:n_family_member[family.idx]) {
      obersed.allele[a] <- (h1[a, snp2look.idx]==1) + (h2[a, snp2look.idx]==1)
    }
    #     (tran_vec[, c("h1","h2")])
    #go over every possible haplotype carrier configuration with early terminatation rule
    #need to determine the number of founder haplotypes, assume know for now
    #can use dynamic programming, skip for now
    match.pattern <- NULL
    #the number of each unique haplotype occured on affecteds
    haplo.count <- table(factor(as.matrix(tran_vec[ ,c("h1","h2")]), LETTERS[1:4]))
    #save the number of each observed unique haplotype occured on affecteds for later use, i.e. can be A,B,D
    haplo.unique.count[[family.idx]] <- table(factor(as.matrix(tran_vec[ ,c("h1","h2")])))
    for(a in 0:as.numeric(haplo.count[1]>0)) { #only loop through observed haplotypes
      for(b in 0:as.numeric(haplo.count[2]>0)) {
        for(c in 0:as.numeric(haplo.count[3]>0)) {
          for(d in 0:as.numeric(haplo.count[4]>0)) {
            allele.count <- as.numeric(apply(tran_vec[ ,c("h1","h2")], 1, function(x) {
              temp <- c(0,0,0,0)
              temp[match(x, c(LETTERS[1:4]))] <- 1
              sum(c(a,b,c,d)*temp)
            }
            ))
            
            #             print(c(allele.count,obersed.allele,identical(allele.count,obersed.allele)))
            result <- haplo.count * c(a,b,c,d)
            result[which(haplo.count==0)] <- -9 #replace non-obsered haplotype with -9
            if(identical(allele.count,obersed.allele)) match.pattern <- rbind(match.pattern, result)
          }
        }
      }
    }
    colnames(match.pattern) <- LETTERS[1:4]
    match.pattern.list[[family.idx]] <-  match.pattern
  }
  
  ##Use EM to impute the carrier haplotype
  #get ks, kn for unambuiguous families
  count.kskn <- lapply(match.pattern.list, function(x) if(nrow(x)==1) {
    return(c(ks3=sum(x==3), ks2=sum(x==2), kn=sum(x==1))) #count of kc, knc for each family
  })
  count.kskn.sum.1 <- colSums(do.call(rbind, count.kskn[which(family.strct.idx==1)]))
  count.kskn.sum.2 <- colSums(do.call(rbind, count.kskn[which(family.strct.idx==2)]))
  #total number of cn, cs
  count.cscn <- sapply(1:n_family, function(family.idx) {
    c(cs3=sum(haplo.unique.count[[family.idx]]==3),
      cs2=sum(haplo.unique.count[[family.idx]]==2), cn=sum(haplo.unique.count[[family.idx]]==1))
  })
  count.cscn.sum.1 <- rowSums(count.cscn[,which(family.strct.idx==1)])
  count.cscn.sum.2 <- rowSums(count.cscn[,which(family.strct.idx==2)])
  #count u and the number of kn ks cn cs in different configurations
  u.list <- list()
  for(a in 1:n_family) {
    n.row.config <- nrow(match.pattern.list[[a]])
    if(n.row.config > 1) {
      result <- NULL
      for(b in 1:n.row.config) {
        family.match.pattern <- match.pattern.list[[a]][b,] #current match pattern
        ks3=sum(family.match.pattern==3)
        ks2=sum(family.match.pattern==2)
        kn=sum(family.match.pattern==1) #count ks kn
        result <- rbind(result, as.vector(c(family=a, config=b, family.match.pattern, ks3, ks2, kn, t(count.cscn[,a]))))
      }
      colnames(result) <- c("family", "config", "A", "B", "C", "D", "ks3", "ks2", "kn", "cs3", "cs2", "cn")
      u.list <- c(u.list ,list(as.data.frame(result)))
    }
  }
  #   colnames(u.list) <- c("family", "config", "A", "B", "C", "D", "ks", "kn", "cs", "cn")
  # match.pattern.list
  #assuming only one locus only for now
  EM <- function(count.kskn.sum.1=count.kskn.sum.1, count.kskn.sum.2=count.kskn.sum.2,
                 count.cscn.sum.1=count.cscn.sum.1, count.cscn.sum.2=count.cscn.sum.2) {
    #initialization
    ks3.1 <- count.kskn.sum.1[1] #known shared variants  (variants occured On ibd 2 or more)
    ks2.1 <- count.kskn.sum.1[2] #known shared variants  (variants occured On ibd 2 or more)
    ks2.2 <- count.kskn.sum.2[2] #known shared variants  (variants occured On ibd 2 or more)
    kn.1 <- count.kskn.sum.1[3] #known non-shared variants (variants occured on ibd 0)
    kn.2 <- count.kskn.sum.2[3] #known non-shared variants (variants occured on ibd 0)
    cs3.1 <- count.cscn.sum.1[1] # total number of shared chromosomes
    cs2.1 <- count.cscn.sum.1[2] # total number of shared chromosomes
    cs2.2 <- count.cscn.sum.2[2] # total number of shared chromosomes
    cn.1 <- count.cscn.sum.1[3] #total number of non-shared chromosomes
    cn.2 <- count.cscn.sum.2[3] #total number of non-shared chromosomes
    
    pn.init.1 <- kn.1/(cn.1) #probability of rare variant on non-shared chromosome
    pn.init.2 <- kn.2/(cn.2) #probability of rare variant on non-shared chromosome
    pn.cur.1 <- ifelse(pn.init.1==0, runif(1), pn.init.1)
    pn.cur.2 <- ifelse(pn.init.2==0, runif(1), pn.init.2)
    ps3.init.1 <- ks3.1/(cs3.1) #probability of rare variant on shared chromosome
    ps3.cur.1 <- ifelse(ps3.init.1==0, runif(1), ps3.init.1)
    ps2.init.1 <- ks2.1/(cs2.1) #probability of rare variant on shared chromosome
    ps2.init.2 <- ks2.2/(cs2.2) #probability of rare variant on shared chromosome
    ps2.cur.1 <- ifelse(ps2.init.1==0, runif(1), ps2.init.1)
    ps2.cur.2 <- ifelse(ps2.init.2==0, runif(1), ps2.init.2)
    delta <- Inf
    iter <- 1
    
    u.idx <- unlist(sapply(u.list, function(x) family.strct.idx[x$family]))
    
    while(delta > 10^-6) {
      print("x")
      print(c(ps3.cur.1, ps2.cur.1, ps2.cur.2, 
      pn.cur.1, pn.cur.2))
      #E step
      #us <- u*ps.cur/(pn.cur+ps.cur)
      #un <- u*pn.cur/(pn.cur+ps.cur)
      if(length(u.list)!=0) { #prevent the error when there's no ambiguous families
        for(a in 1:length(u.list)) {
          temp.prob <- with(u.list[[a]], ifelse(family.strct.idx[family]==1,
                                                ps3.cur.1^ks3*ps2.cur.1^ks2*pn.cur.1^kn*(1-ps3.cur.1)^(cs3-ks3)*(1-ps2.cur.1)^(cs2-ks2)*(1-pn.cur.1)^(cn-kn),
                                                ps2.cur.2^ks2*pn.cur.2^kn*(1-ps2.cur.2)^(cs2-ks2)*(1-pn.cur.2)^(cn-kn)) 
                            ) #go through u.list and calculate the prob. for each config
          sum.prob <- sum(temp.prob)
          u.list[[a]]$prob <- temp.prob/sum.prob
        }
        u <- do.call(rbind, u.list)
        us3.1 <- sum(with(u[which(u.idx==1),], ks3*prob)) #update us
        us2.1 <- sum(with(u[which(u.idx==1),], ks2*prob)) #update us
        us2.2 <- sum(with(u[which(u.idx==2),], ks2*prob)) #update us
        un.1 <- sum(with(u[which(u.idx==1),], kn*prob)) #update un
        un.2 <- sum(with(u[which(u.idx==2),], kn*prob)) #update un
      }else{
        u<-us3.1<-us2.1<-us2.2<-un.1<-un.2<-0
      }
      #M-step
      pn.new.1 <- (kn.1+un.1)/cn.1
      pn.new.2 <- (kn.2+un.2)/cn.2
      ps2.new.1 <- (ks2.1+us2.1)/cs2.1
      ps2.new.2 <- (ks2.2+us2.2)/cs2.2
      ps3.new.1 <- (ks3.1+us3.1)/cs3.1
      #print(c(mu.new, sigma2.new, f.new, cor.factor, optim.result$value))
      
      #check convergence
      delta <- max(abs(pn.cur.1 - pn.new.1), abs(pn.cur.2 - pn.new.2),
                   abs(ps2.cur.1 - ps2.new.1), abs(ps2.cur.2 - ps2.new.2),
                   abs(ps3.cur.1 - ps3.new.1))
      pn.cur.1 <- pn.new.1
      pn.cur.2 <- pn.new.2
      ps2.cur.1 <- ps2.new.1
      ps2.cur.2 <- ps2.new.2
      ps3.cur.1 <- ps3.new.1
      
      #       print(c(pn.cur, ps.cur, iter))
      iter <- iter + 1
      #       print(iter)
    }
    #c(pn.init, ps.init)
    result <- c(ps3.cur.1, ps2.cur.1, ps2.cur.2, 
                pn.cur.1, pn.cur.2)
    names(result) <- c("ps3.1", "ps2.1", "ps2.2", "pn.1", "pn.2")
    return(result)
  }
  
  ##family test in TRAFIC spirit take input of imputed dataset
  family.test.trafic.ext.impute <- function(ps3.1, ps2.1, ps2.2, pn.1, pn.2) {
    ##hypothesis testing count the number of carrier haplotypes transmitted to the offsprings
    
    #extract info from unambiguous family
    data.family.list <- list()
    for(family.idx in 1:length(match.pattern.list)) {
      if(nrow(match.pattern.list[[family.idx]])==1) {
        match.pattern <- match.pattern.list[[family.idx]][which(match.pattern.list[[family.idx]]!=-9)]
        data.family.list <- c(data.family.list, list(t(rbind(haplotype_on_affect=haplo.unique.count[[family.idx]], carrier=match.pattern>=1))))
      }
    }
    
    #extract info from ambiguous family
    #update the configuration probability
    if(length(u.list)!=0) { #prevent the error when there's no ambiguous families
      for(a in 1:length(u.list)) {
        temp.prob <- with(u.list[[a]], ifelse(family.strct.idx[family]==1,
                                              ps3.cur.1^ks3*ps2.cur.1^ks2*pn.cur.1^kn*(1-ps3.cur.1)^(cs3-ks3)*(1-ps2.cur.1)^(cs2-ks2)*(1-pn.cur.1)^(cn-kn),
                                              ps2.cur.2^ks2*pn.cur.2^kn*(1-ps2.cur.2)^(cs2-ks2)*(1-pn.cur.2)^(cn-kn)) 
        )
        sum.prob <- sum(temp.prob)
        prob <- temp.prob/sum.prob #config probability
        sample.config <- sample(1:nrow(u.list[[a]]), 1, prob=prob) #which config
        match.pattern <- u.list[[a]][sample.config,c("A","B","C","D")]
        match.pattern <- match.pattern[which(match.pattern!=-9)]
        
        family.idx <- u.list[[a]]$family[1] #extract family id
        
        data.family.list <- c(data.family.list, list(t(rbind(haplotype_on_affect=haplo.unique.count[[family.idx]], carrier=match.pattern>=1))))
      }
    }
    
    data.family <- as.data.frame(do.call(rbind, data.family.list))
    #fit a logistic regression
    glm.result <- summary(glm(carrier ~ haplotype_on_affect, family=binomial(link = "logit"), data=data.family))
    
    c(p.value=glm.result$coefficients["haplotype_on_affect", "Pr(>|z|)"], 
      est=glm.result$coefficients["haplotype_on_affect", "Estimate"],
      var=glm.result$coefficients["haplotype_on_affect", "Std. Error"]^2)
  }
  # family.test.trafic.ext.impute()
  
  ##test for 3c only use infor from children
  family.test.trap.impute <- function(ps3.1, ps2.1, ps2.2, pn.1, pn.2) {
    ##hypothesis testing count the number of carrier haplotypes transmitted to the offsprings
    #extract info from unambiguous family
    data.family.list <- list()
    for(family.idx in 1:length(match.pattern.list)) {
      if(nrow(match.pattern.list[[family.idx]])==1) {
        match.pattern <- match.pattern.list[[family.idx]][which(match.pattern.list[[family.idx]]!=-9)]
        data.family.list <- c(data.family.list, list(t(rbind(haplotype_on_affect=haplo.unique.count[[family.idx]], carrier=match.pattern>=1))))
      }
    }
    
    #extract info from ambiguous family
    #update the configuration probability
    if(length(u.list)!=0) { #prevent the error when there's no ambiguous families
      for(a in 1:length(u.list)) {
        temp.prob <- with(u.list[[a]], ifelse(family.strct.idx[family]==1,
                                              ps3.cur.1^ks3*ps2.cur.1^ks2*pn.cur.1^kn*(1-ps3.cur.1)^(cs3-ks3)*(1-ps2.cur.1)^(cs2-ks2)*(1-pn.cur.1)^(cn-kn),
                                              ps2.cur.2^ks2*pn.cur.2^kn*(1-ps2.cur.2)^(cs2-ks2)*(1-pn.cur.2)^(cn-kn)) 
        )
        sum.prob <- sum(temp.prob)
        prob <- temp.prob/sum.prob #config probability
        sample.config <- sample(1:nrow(u.list[[a]]), 1, prob=prob) #which config
        match.pattern <- u.list[[a]][sample.config,c("A","B","C","D")]
        match.pattern <- match.pattern[which(match.pattern!=-9)]
        
        family.idx <- u.list[[a]]$family[1]
        
        data.family.list <- c(data.family.list, list(t(rbind(haplotype_on_affect=haplo.unique.count[[family.idx]], carrier=match.pattern>=1))))
      }
    }
    
    #start looking at each family
    test.stat <- sapply(data.family.list, function(x) {
      n_founder <- 2
      n_carrier <- sum(x[,2]) #no. of carrier haplotypes
      if(!(n_carrier==0)) { #skip families with 0 carrier haplotypes in sibpair
        observed <- sum(x[,1] * x[, 2])
        
        #calculate expectation and variance
        n_unique_carrier <- sum(x[,2]==1)
        founder <- t(combn(LETTERS[1:(2*n_founder)], n_unique_carrier)) #use minimum carrier haplotype in sibpair as an estimate for the number of founder's carrier haplotypes
        S <- apply(founder, 1, function(y) {
          carrier <- as.numeric(rownames(x) %in% y) #founder's haplotype
          #       print(carrier)
          
          IBD_haplotype_observed <- sum(x[,1] * carrier)
        }
        )
        
        #       print(S)
        S <- S[which(S!=0)] #only count the configuration when there is at least one variant observed in the sibpair
        c(observed=observed, mean=mean(S), var=sum((S-mean(S))^2)/length(S), n_carrier=n_carrier, n_unique_carrier=n_unique_carrier)
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
    c(p.value=p.value, est=sum(test.stat$observed - test.stat$mean), var=se^2)
  }#family.test.trap.impute
  
  
  #choose trap or trafic to use
  impute.method <- switch(method, trap=family.test.trap.impute, trafic=family.test.trafic.ext.impute)
  #apply test with multiple imputation
  MI <- function() {
    #initialization
    diff <- NULL
    var <- NULL
    D <- 10
    #run EM to estimate the allele frequency
    if(!is.null(null.maf)) {
      pspn <- c(null.maf, null.maf, null.maf)
    }else {
      pspn <- EM(count.kskn.sum.1=count.kskn.sum.1, count.kskn.sum.2=count.kskn.sum.2,
                 count.cscn.sum.1=count.cscn.sum.1, count.cscn.sum.2=count.cscn.sum.2)  
    }
    for(i in 1:D) {
      mi.result <- impute.method(ps3.1=pspn[1], ps2.1=pspn[2], ps2.2=pspn[3],
                                 pn.1=pspn[4], pn.2=pspn[5])
      diff <- cbind(diff, mi.result[2])
      var <- cbind(var, mi.result[3])
    }
    
    TD <- mean(diff)
    VARD <- mean(var) + (1+1/D)*sum((diff-TD)^2)/(D-1)
    c(pchisq(TD^2/VARD, df=1, lower=F))
  }
  result <- MI()
  result
}#family.test.nofounder.impute.noaggregate.seperate


##imputation test for TRAP using the sibpair who shares zero IBD chromosome
##assuming the IBD status is known
family.test.trap.impute.founder <- function(data=family_generated_3c, f=risk.variant.id, D=10, true.f=F, nontran.est.f=F, nontran.true.f=F, conditional.test=F) {
  ##hypothesis testing count the number of carrier haplotypes transmitted to the offsprings
  n_family <- max(data$data_family$family)
  n_family_member <- table(data$data_family$family)
  #check if founder's haplotype carries any variant's with f < 0.1
  if(length(f)==1 & f[1] <1) { #support allele frequency or list of snps
    snp2look.idx <-  which(snp$FREQ1 < f) # snp to look for
  } else(snp2look.idx <-  f)
  
  #disguitish those siblings with four founder's haplotype
  if(true.f==T | nontran.true.f==T | nontran.est.f==T) {
    data.founder <- family_generated_2g3c
    n_family_member.w.founder <- n_family_member + 2
    #calculate the allele frequency
    founder.list <- list() #to store if the haplotype carrying a risk variant on multiple affected
    carrier.count.2 <- 0 
    chromosome.count.2 <- 0
    carrier.count.2.mis <- 0 
    chromosome.count.2.mis <- 0
    carrier.count.2.obs <- 0 
    chromosome.count.2.obs <- 0
    carrier.count.3 <- 0 
    chromosome.count.3 <- 0
    carrier.count.3.mis <- 0 
    chromosome.count.3.mis <- 0
    carrier.count.3.obs <- 0 
    chromosome.count.3.obs <- 0
    carrier.count.4 <- 0 
    chromosome.count.4 <- 0
    haplo.unique.list <- list()
    for(family.idx in 1:n_family) {
      current_row=sum(n_family_member.w.founder[1:family.idx]) - n_family_member.w.founder[family.idx]
      #assign a suitable transmission, not trivial to code up the algorithm assume is known for now
      tran_vec <- data.founder$tran_vec[(current_row+1):(current_row+n_family_member.w.founder[family.idx]),]
      ##tally the (un)ambiguous carrier haplotype
      h1 <- data.founder$data_family[(current_row+1):(current_row+n_family_member.w.founder[family.idx]),7:(6+n_snp)] #the first haplotype
      h2 <- data.founder$data_family[(current_row+1):(current_row+n_family_member.w.founder[family.idx]),-c(1:(6+n_snp))] #the second haplotype
      #observed allele count for each individual
      carrier.founder <- c(0,0,0,0)
      carrier.founder.mis <- c(0,0,0,0)
      carrier.founder.obs <- c(0,0,0,0)
      carrier.offspring <- c(0,0,0,0)
      for(a in 3:n_family_member.w.founder[family.idx]) { #offsprings
        idx.h1 <- match(tran_vec[a,"h1"], c("A", "B", "C", "D"))
        if(h1[a, snp2look.idx]==1) carrier.offspring[idx.h1] <- carrier.offspring[idx.h1] + 1
        idx.h2 <- match(tran_vec[a,"h2"], c("A", "B", "C", "D"))
        if(h2[a, snp2look.idx]==1) carrier.offspring[idx.h2] <- carrier.offspring[idx.h2] + 1
      }
      #the number of each unique haplotype occured on affected offsprings
      haplo.unique <- unique(as.vector(as.matrix(tran_vec[3:n_family_member.w.founder[family.idx] ,c("h1","h2")])))
      haplo.unique.count <- length(haplo.unique)
      haplo.unique.list[[family.idx]] <- haplo.unique.count
      haplo.mis <- which(is.na(match(c("A", "B", "C", "D"), haplo.unique)))
      haplo.obs <- which(!is.na(match(c("A", "B", "C", "D"), haplo.unique)))
      
      for(a in 1:2) { #founder's 
        idx.h1 <- match(tran_vec[a,"h1"], c("A", "B", "C", "D"))
        if(h1[a, snp2look.idx]==1) carrier.founder[idx.h1] <- carrier.founder[idx.h1] + 1
        idx.h2 <- match(tran_vec[a,"h2"], c("A", "B", "C", "D"))
        if(h2[a, snp2look.idx]==1) carrier.founder[idx.h2] <- carrier.founder[idx.h2] + 1
        
        if(idx.h1 %in% haplo.mis & h1[a, snp2look.idx]==1) carrier.founder.mis[idx.h1] <- carrier.founder.mis[idx.h1] + 1
        if(idx.h2 %in% haplo.mis & h2[a, snp2look.idx]==1) carrier.founder.mis[idx.h2] <- carrier.founder.mis[idx.h2] + 1
        if(idx.h1 %in% haplo.obs & h1[a, snp2look.idx]==1) carrier.founder.obs[idx.h1] <- carrier.founder.obs[idx.h1] + 1
        if(idx.h2 %in% haplo.obs & h2[a, snp2look.idx]==1) carrier.founder.obs[idx.h2] <- carrier.founder.obs[idx.h2] + 1
      }
      
      #observed founder carrier based on offsprings
      founder.list[[family.idx]] <- c(observed=sum(carrier.offspring), carrier=sum(carrier.offspring>0), haplo.unique.count=haplo.unique.count)
      #steps to calculate true founder allele ferquency by the number of observed founder chromosome
      if(haplo.unique.count==2) {
        carrier.count.2 <- carrier.count.2 + sum(carrier.founder>=1)
        chromosome.count.2 <- chromosome.count.2 + 4
        carrier.count.2.mis <- carrier.count.2.mis + sum(carrier.founder.mis>=1)
        chromosome.count.2.mis <- chromosome.count.2.mis + 2
        carrier.count.2.obs <- carrier.count.2.obs + sum(carrier.founder.obs>=1)
        chromosome.count.2.obs <- chromosome.count.2.obs + 2
      }
      if(haplo.unique.count==3) {
        carrier.count.3 <- carrier.count.3 + sum(carrier.founder>=1)
        chromosome.count.3 <- chromosome.count.3 + 4
        carrier.count.3.mis <- carrier.count.3.mis + sum(carrier.founder.mis>=1)
        chromosome.count.3.mis <- chromosome.count.3.mis + 1
        carrier.count.3.obs <- carrier.count.3.obs + sum(carrier.founder.obs>=1)
        chromosome.count.3.obs <- chromosome.count.3.obs + 3      }
      if(haplo.unique.count==4) {
        carrier.count.4 <- carrier.count.4 + sum(carrier.founder>=1)
        chromosome.count.4 <- chromosome.count.4 + 4
      }
    }
    #the allele frequency
    (founder.freq.offspring.2 <- carrier.count.2/chromosome.count.2)
    (founder.freq.offspring.2.mis <- carrier.count.2.mis/chromosome.count.2.mis)
    (founder.freq.offspring.2.obs <- carrier.count.2.obs/chromosome.count.2.obs)
    (founder.freq.offspring.3 <- carrier.count.3/chromosome.count.3)
    (founder.freq.offspring.3.mis <- carrier.count.3.mis/chromosome.count.3.mis)
    (founder.freq.offspring.3.obs <- carrier.count.3.obs/chromosome.count.3.obs)
    (founder.freq.offspring.4 <- carrier.count.4/chromosome.count.4)
    
    #   print(founder.freq.offspring.4 <- carrier.count.4/chromosome.count.4)
    #   print(round(dbinom(0:4, 4, founder.freq.offspring.4), 3))
    #   print((founder.freq.offspring.all <- carrier.count.all/chromosome.count.all))
    #   print(round(dbinom(0:4, 4, founder.freq.offspring.all), 3))
    founder.list <- as.data.frame(do.call(rbind, founder.list))
    if(true.f==T) {
      impute.prob <- sapply(unlist(haplo.unique.list), function(x) {
        switch(as.character(x), 
               "4" = founder.freq.offspring.4,
               "3" = founder.freq.offspring.3,
               "2" = founder.freq.offspring.2)
      })
    }
    if(nontran.true.f==T) {
      impute.prob <- sapply(unlist(haplo.unique.list), function(x) {
        switch(as.character(x), 
               "4" = founder.freq.offspring.4,
               "3" = founder.freq.offspring.3.mis,
               "2" = founder.freq.offspring.2.mis)
      })
    }
  }
  
  
  #use family which inherits all 4 founder haplotypes
  if(!(true.f==T | nontran.true.f==T | nontran.est.f==T)) {
    #calculate the allele frequency
    founder.list <- list() #to store if the haplotype carrying a risk variant on multiple affected
    carrier.count.4 <- 0 
    chromosome.count.4 <- 0
    carrier.count.all <- 0 
    chromosome.count.all <- 0
    for(family.idx in 1:n_family) {
      current_row=sum(n_family_member[1:family.idx]) - n_family_member[family.idx]
      #assign a suitable transmission, not trivial to code up the algorithm assume is known for now
      tran_vec <- data$tran_vec[(current_row+1):(current_row+n_family_member[family.idx]),]
      ##tally the (un)ambiguous carrier haplotype
      h1 <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),7:(6+n_snp)] #the first haplotype
      h2 <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),-c(1:(6+n_snp))] #the second haplotype
      #observed allele count for each individual
      carrier <- c(0,0,0,0)
      for(a in 1:n_family_member[family.idx]) {
        idx.h1 <- match(tran_vec[a,"h1"], c("A", "B", "C", "D"))
        if(h1[a, snp2look.idx]==1) carrier[idx.h1] <- carrier[idx.h1] + 1
        idx.h2 <- match(tran_vec[a,"h2"], c("A", "B", "C", "D"))
        if(h2[a, snp2look.idx]==1) carrier[idx.h2] <- carrier[idx.h2] + 1
      }
      #the number of each unique haplotype occured on affecteds
      haplo.unique <- length(unique(as.vector(as.matrix(tran_vec[ ,c("h1","h2")]))))
      
      founder.list[[family.idx]] <- c(observed=sum(carrier), carrier=sum(carrier>0), haplo.unique=haplo.unique)
      if(haplo.unique==4) {
        carrier.count.4 <- carrier.count.4 + sum(carrier>=1)
        chromosome.count.4 <- chromosome.count.4 + haplo.unique
      }
      carrier.count.all <- carrier.count.all + sum(carrier>=1)
      chromosome.count.all <- chromosome.count.all + haplo.unique
    }
    #the allele frequency
    (founder.freq.offspring.4 <- carrier.count.4/chromosome.count.4)
    #   print(founder.freq.offspring.4 <- carrier.count.4/chromosome.count.4)
    #   print(round(dbinom(0:4, 4, founder.freq.offspring.4), 3))
    (founder.freq.offspring.all <- carrier.count.all/chromosome.count.all)
    #   print((founder.freq.offspring.all <- carrier.count.all/chromosome.count.all))
    #   print(round(dbinom(0:4, 4, founder.freq.offspring.all), 3))
    founder.list <- as.data.frame(do.call(rbind, founder.list))
    impute.prob <- rep(founder.freq.offspring.4, length.out = n_family)
  }

#   print(table(founder.list$carrier))
#   print(table(4-founder.list$haplo.unique))
  
  trap.impute <- function() {
    
    #redraw the estimation of allele frequency in each imputation iteration
    if(nontran.est.f==T) {
      #do drawing of allele frequenc to take into account the estimation uncertainty
      carrier.count.4.draw <- rbinom(1, chromosome.count.4, founder.freq.offspring.4)
      founder.freq.offspring.4.draw <- carrier.count.4.draw/chromosome.count.4
      carrier.count.2.obs.draw <- rbinom(1, chromosome.count.2.obs, founder.freq.offspring.2.obs)
      founder.freq.offspring.2.obs.draw <- carrier.count.2.obs.draw/chromosome.count.2.obs
      carrier.count.3.obs.draw <- rbinom(1, chromosome.count.3.obs, founder.freq.offspring.3.obs)
      founder.freq.offspring.3.obs.draw <- carrier.count.3.obs.draw/chromosome.count.3.obs
      
      (founder.freq.offspring.2.mis.est <- 2*(founder.freq.offspring.4.draw - 0.5*founder.freq.offspring.2.obs.draw))
      founder.freq.offspring.2.mis.est <- ifelse(founder.freq.offspring.2.mis.est > 0, founder.freq.offspring.2.mis.est, 0)
      (founder.freq.offspring.3.mis.est <- 4*(founder.freq.offspring.4.draw - 3/4*founder.freq.offspring.3.obs.draw))
      founder.freq.offspring.3.mis.est <- ifelse(founder.freq.offspring.3.mis.est > 0, founder.freq.offspring.3.mis.est, 0)
      
      impute.prob <- sapply(unlist(haplo.unique.list), function(x) {
        switch(as.character(x), 
               "4" = founder.freq.offspring.4,
               "3" = founder.freq.offspring.3.mis.est,
               "2" = founder.freq.offspring.2.mis.est)
      })
      print(c(founder.freq.offspring.2, founder.freq.offspring.3, founder.freq.offspring.4))
      print(c(founder.freq.offspring.2.obs, founder.freq.offspring.3.obs))
      print(c(founder.freq.offspring.2.mis, founder.freq.offspring.2.mis.est, 
              founder.freq.offspring.3.mis, founder.freq.offspring.3.mis.est))
    }
    
#     founder.carrier.impute <- mapply(sample, mapply(seq, 0, 4-founder.list$haplo.unique), as.list(rep(1, length.out = n_family)), as.list(rep(FALSE, length.out = n_family)), mapply(impute.prob, as.list(founder.list$carrier), as.list(4-founder.list$haplo.unique), as.list(rep(founder.freq.offspring.4, length.out = n_family)))) 
    founder.carrier.impute <- rbinom(n = rep(1, length.out = n_family), size = (4-founder.list$haplo.unique), prob = impute.prob)
    founder.carrier.impute <- founder.carrier.impute + founder.list$carrier 
    
#     table(founder.list$carrier)
#     table(founder.carrier.impute)
#     print(prop.table(table(founder.carrier.impute)))
#     founder.carrier.impute <- founder.list$carrier
#     
#     founder.carrier.impute <- rbinom(n = rep(1, length.out = n_family), size = rep(4, length.out = n_family), prob = founder.freq.offspring.4)
    
    
    #test statistic procedure
    #start looking at each family
    test.stat <- lapply(1:n_family, function(x) {
      family.idx=x 
      n_founder <- 2 #!!!!!!!!set two unknown founders for ease of testing run
      n_carrier <- founder.list$observed[family.idx] #no. of carrier haplotypes
      
      #!!!!!!!!using all families not just >0 carrier founder haplotype families
      if(conditional.test==F) {
        criteria <- !(founder.carrier.impute[family.idx]==0 | founder.carrier.impute[family.idx]==4)
      }else{
        criteria <- n_carrier>0
      }
        
      if(criteria) { #skip families with 0 carrier haplotypes in sibpair
        current_row=sum(n_family_member[1:family.idx]) - n_family_member[family.idx]
        #assign a suitable transmission, not trivial to code up the algorithm assume is known for now
        affect <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),"affect"]
        tran_vec <- data$tran_vec[(current_row+1):(current_row+n_family_member[family.idx]),]
        observed <- n_carrier
        
        #calculate expectation and variance
        n_unique_carrier <- founder.carrier.impute[family.idx]
        founder <- t(combn(LETTERS[1:(2*n_founder)], n_unique_carrier)) #use minimum carrier haplotype in sibpair as an estimate for the number of founder's carrier haplotypes
        S <- apply(founder, 1, function(x) {
          carrier <- x #founder's haplotype
          #       print(carrier)
          IBD_haplotype_observed = 0
          for(i in 1:n_family_member[family.idx]) {
            #print(c(sum(tran_vec[current_row+i, c("h1","h2")] %in% carrier),(affect[current_row+i]-1)))
            IBD_haplotype_observed <- IBD_haplotype_observed + sum(tran_vec[i, c("h1","h2")] %in% carrier)*(affect[i]-1)
          }
          IBD_haplotype_observed
        })
        
        #       print(S)
        if(conditional.test==T) {
          S <- S[which(S!=0)] #only count the configuration when there is at least one variant observed in the sibpair
        }
        c(observed=observed, mean=mean(S), var=sum((S-mean(S))^2)/length(S), n_carrier=n_carrier, n_unique_carrier=n_unique_carrier)
      }
    }
    )
    test.stat <- data.frame(do.call(rbind, test.stat))
    
    v <- test.stat$var
    se <- sqrt(sum(v))
    #e_avg <- mean(e) #overall expectation
    final.test.stat <- sum(test.stat$observed - test.stat$mean)/se
    p.value <- 2*pnorm(abs(final.test.stat), lower=F) #p-value
#     print(c(p.value=p.value, est=sum(test.stat$observed - test.stat$mean), var=se^2, final.test.stat=final.test.stat))
    c(p.value=p.value, est=sum(test.stat$observed - test.stat$mean), var=se^2)
#     est <- mean(test.stat$observed - test.stat$mean)
#     se <- sd(test.stat$mean)/sqrt(nrow(test.stat))
#     p.value <- 2*pnorm(abs(est/se), lower=F)
#     c(p.value=p.value, est=est, var=se^2)
  }
#family.test.trap.impute.founder
  
  #apply test with multiple imputation
  MI <- function(D) {
    #initialization
    diff <- NULL
    var <- NULL
    D <- D
    #run EM to estimate the allele frequency
    for(i in 1:D) {
      mi.result <- trap.impute()
      diff <- cbind(diff, mi.result[2])
      var <- cbind(var, mi.result[3])
    }
    
    TD <- mean(diff)
    VARD <- mean(var) + (1+1/D)*sum((diff-TD)^2)/(D-1)
    df <- (D-1)*(1+D*mean(var)/((D+1)*sum((diff-TD)^2)/(D-1)))^2    
#     c(2*pnorm(abs(TD/sqrt(VARD)), lower=F), 2*pt(abs(TD/sqrt(VARD)), df=df, lower=F))
    c(2*pnorm(abs(TD/sqrt(VARD)), lower=F))
  }
  result <- MI(D)
  result
}#family.test.trap.impute.founder


#pairwise trafic test aussuming the phase is known
trafic.pairwise.test <- function(data=family_generated_2g2c, f=risk.variant.id) {
  n_snp <- (ncol(data$data_family)-6)/2
  if(length(f)==1 & f[1] <1) { #support allele frequency or list of snps
    variant.id <-  which(snp$FREQ1 < f) # snp to look for
  } else(variant.id <-  f)
  n_variant <- length(variant.id) # no. of variant of interest
  ##initialization
  n_IBD_carrier <- 0
  n_nonIBD_carrier <- 0
  n_IBD <- 0
  n_nonIBD <- 0
  
  ##do by family family_idx <- 1
  for(a in 1:max(data$data_family[, "family"])) {
    ind_idx <- with(data$data_family, which(family==a & affect==2))
    ##extract the haplotype on variants of interest only on affected individual
    data_family <- data$data_family[ind_idx, c(1:6, 6+variant.id, 6+n_snp+variant.id)]
    tran_vec <- data$tran_vec[ind_idx, ]
    
    ##do every piar in the famiy
    n_affect <- length(ind_idx)
    pair_matrix <- combn(1:n_affect, 2) 
    for(b in 1:ncol(pair_matrix)) {#pair_idx <- 1
      ##assess the number of shared IBD chromosome region
      pair_data_family <- data_family[pair_matrix[, b], ]
      pair_tran_vec <- tran_vec[pair_matrix[, b], ]
      sibpair_haplotype <- as.data.frame(matrix(NA, nrow=4, ncol=6+n_variant, dimnames=list(NULL, colnames(pair_data_family)[1:(n_variant+6)])))
      sibpair_haplotype[seq(1,4, by=2),1:6] <- pair_data_family[,1:6]
      sibpair_haplotype[seq(2,4, by=2),1:6] <- pair_data_family[,1:6]
      sibpair_haplotype[seq(1,4, by=2),7:(n_variant+6)] <- pair_data_family[,7:(n_variant+6)]
      sibpair_haplotype[seq(2,4, by=2),7:(n_variant+6)] <- pair_data_family[,(n_variant+7):(6+2*n_variant)]
      
      sibpair_tran_vec <- as.data.frame(matrix(NA, nrow=4, ncol=1, dimnames=list(NULL, "Chrom")))
      sibpair_tran_vec[seq(1,4, by=2),1] <- pair_tran_vec[,2]
      sibpair_tran_vec[seq(2,4, by=2),1] <- pair_tran_vec[,3]
      
      ##allele frequency on the shared and non-shared chromosome
      S_hap <- rep(NA, nrow(sibpair_tran_vec)) #indicate the Shared (only count once) and non-shared chromosome
      sapply(LETTERS[1:4], function(x) {
          idx <- which(sibpair_tran_vec == x)
          if(length(idx)==2) S_hap[idx] <<- c("S", "D")
          if(length(idx)==1) S_hap[idx] <<- "NS"
        }
      )
      which(sibpair_tran_vec == sibpair_tran_vec[2,])
      n_IBD <- n_IBD + sum(S_hap=="S")
      n_nonIBD <- n_nonIBD + sum(S_hap=="NS")
      #count the true carrier haplotype and store them
      n_IBD_carrier <- n_IBD_carrier + sum(apply(2-sibpair_haplotype[which(S_hap=="S"), -(1:6)], 1, sum)>0)
      n_nonIBD_carrier <- n_nonIBD_carrier + sum(apply(2-sibpair_haplotype[which(S_hap=="NS"), -(1:6)], 1, sum)>0)
    }
  } 
  #trafic test assuming phase is known
  c(n_IBD, n_nonIBD, n_IBD_carrier, n_nonIBD_carrier)
  prop.test(x=c(n_IBD_carrier, n_nonIBD_carrier), n=c(n_IBD, n_nonIBD), correct=F)
}

trafic.test <- function(data=family_generated, f, true.only=F) {
  n_family <- max(data$data_family$family)
  
  if(length(f)==1 & f[1] <1) { #support allele frequency or list of snps
    variant.id <-  which(snp$FREQ1 < f) # snp to look for
  } else(variant.id <-  f)
  #convert sibpair data from Merlin format to each haplotype a row
  # data$data_family
  n_haplotype = 4*n_family
  sibpair_haplotype <- as.data.frame(matrix(-999, nrow=n_haplotype, ncol=6+nrow(snp), dimnames=list(NULL, colnames(data$data_family)[1:56])))
  sibpair_haplotype[seq(1,n_haplotype, by=2),1:6] <- data$data_family[,1:6]
  sibpair_haplotype[seq(2,n_haplotype, by=2),1:6] <- data$data_family[,1:6]
  sibpair_haplotype[seq(1,n_haplotype, by=2),7:56] <- data$data_family[,7:56]
  sibpair_haplotype[seq(2,n_haplotype, by=2),7:56] <- data$data_family[,57:106]
  
  sibpair_tran_vec <- as.data.frame(matrix(-999, nrow=n_haplotype, ncol=1, dimnames=list(NULL, "Chrom")))
  sibpair_tran_vec[seq(1,n_haplotype, by=2),1] <- data$tran_vec[,2]
  sibpair_tran_vec[seq(2,n_haplotype, by=2),1] <- data$tran_vec[,3]
  
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
  true.result <- fisher.test(table(S_hap[which(!S_hap=="D")], apply(2-as.matrix(sibpair_haplotype[, 6+variant.id]), 1, sum)[which(!S_hap=="D")]>0))
#   true.result <- prop.test(table(S_hap[which(!S_hap=="D")], apply(2-as.matrix(sibpair_haplotype[, 6+variant.id]), 1, sum)[which(!S_hap=="D")]>0), correct = F)
  
  if(true.only==T) {return(true.result)}
  
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

  c(trafic.true=true.result$p.value, trafic.haplotype=haplotype.result[3], trafic.genotype=genotype.result[3])
}


gene_case_control <- function(n_case_control_pair=1000, n_case=NA, n_control=NA) {
  if(is.na(n_case)) {
    n_case <- n_case_control_pair
    n_control <- n_case_control_pair
  }
  data_case_control <- matrix(NA, nrow=n_case+n_control, ncol=(2+n_snp*2))
  #generate cases first
  current_row <- 1 #which row is generating
  n_case.idx <- 1 #index of how many families have been generated
  while(n_case.idx <= n_case) {
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
  while(n_control.idx <= n_control) {
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
  n_case <- sum(data$affect==2)
  n_control <- nrow(data) - n_case
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
  
  carrier.case <- sum(carrier[,1:n_case])
  carrier.control <- sum(carrier[,(n_case+1):(nrow(data))])
  test.result <- fisher.test(matrix(c(carrier.case, carrier.control, 2*n_case - carrier.case, 2*n_control - carrier.control), nrow=2)) #turn off correct to avoid a conservative test 
  c(test.result$estimate[1],test.result$estimate[2], p.value=test.result$p.value)
}
# case_control.test()





#calculate the probability of obersing the family given the affected status and transmission pattern
#
#family.structure takes Merlin format: family.id individual.id father mother gender affected
#1 1 0 0 1 1
#1 2 0 0 2 2
#1 3 1 2 0 2
#1 4 1 2 0 2
#ex:
#family_strct.2g2c <- data.frame(family=c(1,1,1,1), person=c(1,2,3,4), father=c(0,0,1,1), 
#mother=c(0,0,2,2), sex=c(1,2,0,0), affect=c(1,2,2,2)) #1=male, 2=female, 1=unaffected, 2=affected
#trans.pattern keeps track how the haplotype is inhereted by offspring: family.id individual.id
#founder has its own number for founders haplotypes
#1 1 1 2
#1 2 3 4
#1 3 1 3
#1 4 2 3
CalculateFamilyProb <- function(family.structure, trans.pattern, mu, sigma2) {
  #no. of family member
  n.member <- nrow(family.structure)
  #read family.structure
  #read trans.pattern
  #calculate probability
}

#calculate all the possible configurations of transmission pattern given the affected status and family structure

#draw the family with certain transmission pattern with input of all possible probability configuration


