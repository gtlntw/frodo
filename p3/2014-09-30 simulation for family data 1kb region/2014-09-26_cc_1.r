### nohup R --vanilla < 2014-09-26_cc_1.r >! 2014-09-26_cc_1.txt &

source("function.r")
##the founder's haplotypes
#import the haplotypes generated by cosi
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
sum(snp$FREQ1 < 0.01) # number of snp with f < 0.05
sum(snp$FREQ1 == 0.0001) # number of snp with f < 0.05

##assign risk variants and the corresponding effect size (proportional to allele frequency)
n_haplo <- 10000
n_snp <- ncol(haplotype)-2
prevalence <- 0.05
b0_sqrt <- sqrt(prevalence)   #baseline
temp.idx <- sort(c(sample(which((snp$FREQ1 < 0.05)==T & snp$FREQ1 != 1/n_haplo ), 5), sample(which(snp$FREQ1 == 1/n_haplo), 5))) #randomly pick 5 variants with f < 0.05 and 5 singletons as casual 
b <- rep(1, length=(ncol(haplotype)-2))
b[temp.idx] <- abs(log10(snp$FREQ1[temp.idx])) #effect size is log10 of its allele frequency
#calculate the haplotype variants p(A|h)
haplotype.risk <- apply(2-haplotype[, -c(1:2)], 1, function(x) prod(b^x)*b0_sqrt)
#calculate the probability of drawing given the affected status of founder p(h|A)
# haplotype.draw_prob <- (haplotype.risk %o% haplotype.risk)/n_haplo^2/prevalence
mean(b) #mean relative risk
mean(apply(2-haplotype[, which((snp$FREQ1 < 0.01)==T)+2], 1, sum)) #carrier haplotype frequency


##simulation with 1000 replications
seed=1000
f=0.01
n_pair=1000
rep=1000
rep.idx <<- 1
sim_result <- replicate(rep, {
  print(rep.idx)
  rep.idx <<- rep.idx + 1
  case_control_generated <- gene_case_control(n=n_pair)
  case_control.test(data=case_control_generated, f=f)
})

write.csv(data.frame(f=f, prop.1=sim_result[1,], prop.2=sim_result[2,], p.value=sim_result[3,]),
          paste("res_",f,"_",n_pair,"_",seed,"_case_control",".csv",sep=""), row.names=FALSE)

# rep.idx <- 1
# sim_result <- replicate(n_rep, {
#   print(rep.idx)
#   family_generated <- gene_family(n=500) 
#   family.test(data=family_generated, f=0.01)
#   rep.idx <- rep.idx + 1
# })

