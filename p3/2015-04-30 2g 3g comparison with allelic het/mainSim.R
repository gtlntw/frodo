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
print(r)
print(n_family) #number of family
print(p_dis) #prevelance

## ... your simualtion code
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
sum(snp$FREQ1 < 0.01) # number of snp with f < 0.01
sum(snp$FREQ1 == 0.0001) # number of singletons

##assign risk variants and the corresponding effect size (proportional to allele frequency)
# null <- FALSE
n_haplo <- 10000
n_snp <- ncol(haplotype)-2
prevalence <- p_dis
b0_sqrt <- sqrt(prevalence)   #baseline

#set up causual SNPs
#generate risk haplotypes
# risk.variant.id <- sort(sample(which(snp$FREQ1 < 0.20), 10))
risk.variant.id <- c(18,19,22,23,25,34,39,40,45,46)
risk.haplo.id <- which(apply(2-haplotype[, risk.variant.id+2], 1, sum)>0)
(rish.haplo.f <- mean(apply(2-haplotype[, risk.variant.id+2], 1, sum)>0)) #carrier haplotype frequency
haplotype.risk <- rep(1, length=nrow(haplotype))
#assign mean relative risk calculate the haplotype variants p(A|h)
haplotype.risk[risk.haplo.id] <- r
mean(haplotype.risk[risk.haplo.id]) #mean relative risk
haplotype.risk <<- haplotype.risk*b0_sqrt

##gene drop simulation for two generations
family_strct.2g.3a.1u <- data.frame(family=c(1,1,1,1), person=c(1,2,3,4), father=c(0,0,1,1), 
                                    mother=c(0,0,2,2), sex=c(1,2,1,1), affect=c(1,2,2,2)) #1=male, 2=female, 1=unaffected, 2=affected
family_strct.3g.3a.4u <- data.frame(family=c(1,1,1,1,1,1,1), person=c(1,2,3,4,5,6,7), father=c(0,0,1,1,0,4,4), 
                                    mother=c(0,0,2,2,0,5,5), sex=c(1,2,1,1,2,1,1), affect=c(1,2,1,2,1,2,1)) #1=male, 2=female, 1=unaffected, 2=affected
family_strct.3g.2a.5u <- data.frame(family=c(1,1,1,1,1,1,1), person=c(1,2,3,4,5,6,7), father=c(0,0,1,1,0,4,4), 
                                    mother=c(0,0,2,2,0,5,5), sex=c(1,2,1,1,2,1,1), affect=c(1,2,1,1,1,2,1)) #1=male, 2=female, 1=unaffected, 2=affected
family_strct.3g.3a.5u <- data.frame(family=c(1,1,1,1,1,1,1,1), person=c(1,2,3,4,5,6,7,8), father=c(0,0,1,1,0,4,4,4), 
                                    mother=c(0,0,2,2,0,5,5,5), sex=c(1,2,1,1,2,1,1,1), affect=c(1,2,1,2,1,2,1,1)) #1=male, 2=female, 1=unaffected, 2=affected

rep.idx <<- 1
sim_result <- replicate(n_rep, {
  print(rep.idx)
  rep.idx <<- rep.idx + 1
  
  #simulation for two-generation families
  family_generated_2g.3a.1u <- gene_family(family_strct=family_strct.2g.3a.1u, n=n_family, haplotype.risk=haplotype.risk) 
  result.trap.2g.3a.1u <- family.test(data=family_generated_2g.3a.1u, f=risk.variant.id)
  result.trafic.ext.2g.3a.1u <- family.test.trafic.ext(data=family_generated_2g.3a.1u, f=risk.variant.id)
  
  #simulation for three-generation families
  #three affected and four unaffected
  family_generated_3g.3a.4u <- gene_family(family_strct=family_strct.3g.3a.4u, n=n_family, haplotype.risk=haplotype.risk) 
  result.trap.3g.3a.4u <- family.test(data=family_generated_3g.3a.4u, f=risk.variant.id)
  result.trafic.ext.3g.3a.4u <- family.test.trafic.ext(data=family_generated_3g.3a.4u, f=risk.variant.id)
  #two affected and five unaffected
  family_generated_3g.2a.5u <- gene_family(family_strct=family_strct.3g.2a.5u, n=n_family, haplotype.risk=haplotype.risk) 
  result.trap.3g.2a.5u <- family.test(data=family_generated_3g.2a.5u, f=risk.variant.id)
  result.trafic.ext.3g.2a.5u <- family.test.trafic.ext(data=family_generated_3g.2a.5u, f=risk.variant.id)
  #three affected and five unaffected
  family_generated_3g.3a.5u <- gene_family(family_strct=family_strct.3g.3a.5u, n=n_family, haplotype.risk=haplotype.risk) 
  result.trap.3g.3a.5u <- family.test(data=family_generated_3g.3a.5u, f=risk.variant.id)
  result.trafic.ext.3g.3a.5u <- family.test.trafic.ext(data=family_generated_3g.3a.5u, f=risk.variant.id)
  #only report p.value
  c(result.trap.2g.3a.1u[7], result.trafic.ext.2g.3a.1u, 
    result.trap.3g.3a.4u[7], result.trafic.ext.3g.3a.4u,
    result.trap.3g.2a.5u[7], result.trafic.ext.3g.2a.5u,
    result.trap.3g.3a.5u[7], result.trafic.ext.3g.3a.5u)
})
## Write out your results to a csv file
result.df <- as.data.frame(t(sim_result))
colnames(result.df) <- c("result.trap.2g.3a.1u","result.trafic.ext.2g.3a.1u",
                         "result.trap.3g.3a.4u","result.trafic.ext.3g.3a.4u",
                         "result.trap.3g.2a.5u","result.trafic.ext.3g.2a.5u",
                         "result.trap.3g.3a.5u","result.trafic.ext.3g.3a.5u")
result.df <- cbind(seed,r,p_dis,n_family,result.df)
write.csv(result.df, paste("res_",r,"_",n_family,"_",seed,".csv",sep=""), row.names=FALSE)
## R.miSniam
## End(Not run)
## Not run: