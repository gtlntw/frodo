source("create_dataset.R")
source("TRAFIC.r")

haplotype <- read.table("out_100k_10k_100kb.hap-1", header=F)
snp <-read.table("out_100k_10k_100kb.pos-1", header=T)
colnames(haplotype) <- c("HAP", "CHROM", paste(snp$CHROM, snp$CHROM_POS, sep=":"))
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
# null <- FALSE
n_haplo <- 10000
n_snp <- ncol(haplotype)-2
prevalence <- 0.05
b0_sqrt <- sqrt(prevalence)   #baseline


##gene drop simulation for three generations
family_strct.2g2c <- data.frame(family=c(1,1,1,1), person=c(1,2,3,4), father=c(0,0,1,1),
                                mother=c(0,0,2,2), sex=c(1,2,1,1), affect=c(1,2,2,2)) #1=male, 2=female, 1=unaffected, 2=affected

null<-TRUE

set.seed(0718)
# sim_result <- replicate(n_rep, {
##generate risk haplotypes
risk.variant.id <- sort(sample(which(snp$FREQ1 < 0.01), 10))
print(risk.variant.id)
risk.haplo.id <- which(apply(2-haplotype[, risk.variant.id+2], 1, sum)>0)
(rish.haplo.f <- mean(apply(2-haplotype[, risk.variant.id+2], 1, sum)>0)) #carrier haplotype frequency
# f= 0.0252
b <- rep(1, length=(ncol(haplotype)-2))
if(null==FALSE) b[risk.variant.id] <- abs(log10(snp$FREQ1[risk.variant.id])) #effect size is log10 of its allele frequency
#calculate the haplotype variants p(A|h)
haplotype.risk <<- apply(2-haplotype[, -c(1:2)], 1, function(x) prod(b^x)*b0_sqrt)
#calculate the probability of drawing given the affected status of founder p(h|A)
# haplotype.draw_prob <- (haplotype.risk %o% haplotype.risk)/n_haplo^2/prevalence
(r <<- mean(haplotype.risk[risk.haplo.id]/b0_sqrt)) #mean relative risk
(sigma2 <<- var(haplotype.risk[risk.haplo.id]/b0_sqrt)) #mean relative risk

##generate family data and index parents
family_generated <- gene_family(family_strct=family_strct.2g2c, n=1000)
temp.idx <- which(family_generated$data_family$person==1:2)
# family_generated$data_family <- family_generated$data_family[-temp.idx, ]
# family_generated$tran_vec <- family_generated$tran_vec[-temp.idx, ]

##convert to merlin format
data_family_temp <- family_generated$data_family
data_family_temp[, seq(7,n_snp*2+6, by=2)] <- family_generated$data_family[,7:(n_snp+6)]
data_family_temp[, seq(8,n_snp*2+6, by=2)] <- family_generated$data_family[,(n_snp+7):(n_snp*2+6)]
data_family_temp[, 7:(n_snp*2+6)] <- 3 - data_family_temp[, 7:(n_snp*2+6)] #minor allele:2 major allele:1
data_family_temp[temp.idx, -(1:6)] <- 0  ##replace father and mother genotype with missing 0
write.table(data_family_temp, "data_sibpair.ped", quote=F, sep=" ", col.names=F, row.names=F)
write.table(cbind(c("A", rep("M", n_snp)), c("disease", colnames(haplotype)[-(1:2)])), "data_sibpair.dat", quote=F, sep=" ", col.names=F, row.names=F)
write.table(cbind(snp$CHROM, colnames(haplotype)[-(1:2)], snp$CHROM_POS/10^6), "data_sibpair.map", quote=F, sep=" ", col.names=F, row.names=F)



#examine sibpair data
# family_generated$data_family
sibpair_haplotype <- as.data.frame(matrix(-999, nrow=4000, ncol=6+nrow(snp), dimnames=list(NULL, colnames(family_generated$data_family)[1:(n_snp+6)])))
sibpair_haplotype[seq(1,4000, by=2),1:6] <- family_generated$data_family[-temp.idx,1:6]
sibpair_haplotype[seq(2,4000, by=2),1:6] <- family_generated$data_family[-temp.idx,1:6]
sibpair_haplotype[seq(1,4000, by=2),7:(n_snp+6)] <- family_generated$data_family[-temp.idx,7:(n_snp+6)]
sibpair_haplotype[seq(2,4000, by=2),7:(n_snp+6)] <- family_generated$data_family[-temp.idx,(n_snp+7):(2*n_snp+6)]

sibpair_tran_vec <- as.data.frame(matrix(-999, nrow=4000, ncol=1, dimnames=list(NULL, "Chrom")))
sibpair_tran_vec[seq(1,4000, by=2),1] <- family_generated$tran_vec [-temp.idx,2]
sibpair_tran_vec[seq(2,4000, by=2),1] <- family_generated$tran_vec [-temp.idx,3]

#allele frequency on the shared and non-shared chromosome
S_hap <- rep(NA, nrow(sibpair_tran_vec)) #indicate the Shared (only count once) and non-shared chromosome
for(i in seq(1,nrow(sibpair_tran_vec), by=4)) { #S:shared  NS:non-shared  D: to dump shared
  S_hap[i] <- ifelse(sibpair_tran_vec[i,1]==sibpair_tran_vec[i+2,1], "S", "NS")
  S_hap[i+1] <- ifelse(sibpair_tran_vec[i+1,1]==sibpair_tran_vec[i+1+2,1], "S", "NS")
  S_hap[i+2] <- ifelse(sibpair_tran_vec[i+2,1]==sibpair_tran_vec[i,1], "D", "NS")
  S_hap[i+3] <- ifelse(sibpair_tran_vec[i+3,1]==sibpair_tran_vec[i+1,1], "D", "NS")
}
table(S_hap) #no. of shared and non-shared chromosome
tapply(apply(2-sibpair_haplotype[, 6+risk.variant.id], 1, sum), S_hap, sum) #allele on shared and non-shared
tapply(apply(2-sibpair_haplotype[, 6+risk.variant.id], 1, sum)>0, S_hap, sum) #number of risk haplotype based on collaping framework

#test statistics
t.test(apply(2-sibpair_haplotype[, 6+risk.variant.id], 1, sum)[which(!S_hap=="D")] ~ S_hap[which(!S_hap=="D")], var.equal=TRUE)

#create S vector
S_sibpair <- (sibpair_tran_vec[seq(1,length(S_hap), by=4),] == sibpair_tran_vec[seq(3,length(S_hap), by=4),]) +
  (sibpair_tran_vec[seq(2,length(S_hap), by=4),] == sibpair_tran_vec[seq(4,length(S_hap), by=4),])
S_sibpair
write.table(S_sibpair, "S_sibpair.ibd", , quote=F, sep=" ", col.names=F, row.names=F)
