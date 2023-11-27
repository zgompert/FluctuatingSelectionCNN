## simple simulation of data = initial, allele freqs, environment, ne, betas
## consider 100 possible causal variants

L<-100 ## 100 snps
N<-10 ## 10 pops.
G<-11 ## 10 gens

## sample p from beta, assumes no population structure
p<-rbeta(L,.5,.5)
## convert to maf
p[p>.5]<-1-p[p>.5]

## write results
p0<-matrix(rep(round(p,4),N),nrow=L,ncol=N,byrow=FALSE)

## write Ne, 100 samples from posterior, same beta for each
u<-600
l<-400
ne<-matrix(round(rbeta(100*N,10,10) * (u-l) + l,2),nrow=100,ncol=N)

## random environment from standard normal
env<-matrix(round(rnorm((G-1) * N,0,1),3),nrow=G-1,ncol=N)

## everything above is held constant, focus on different genetic architectures
## will set vz to 1, so vBV < 1 and is = h2
## h2 = .3 or .7
## known SNPs or pip = .1

## assume no uncertainty in SNPs in model
sd<-0.12
beta<-rnorm(L,0,sd)
vBV<-2 * sum(beta^2 * p * (1-p))
vBV
# 0.313816

write.table(cbind(rep(1,L),beta),file="sim_trait_1.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

# sd<-0.17
# beta<-rnorm(L,0,sd)
# vBV<-2 * sum(beta^2 * p * (1-p))
# vBV
# # 0.7065904
# 
# write.table(cbind(rep(1,L),beta),file="sim_trait_2.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
# 
# ## now model averaged version and pip = .1
# pip<-.1
# 
# sd<-0.35
# beta<-rnorm(L,0,sd)
# vBV<-2 * sum(pip*beta^2 * p * (1-p))
# vBV
# # 0.3055152
# 
# write.table(cbind(rep(1,L),beta),file="sim_trait_3.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
# 
# 
# pip<-.1
# 
# sd<-0.45
# beta<-rnorm(L,0,sd)
# vBV<-2 * sum(pip*beta^2 * p * (1-p))
# vBV
# # 0.705306
# 
# write.table(cbind(rep(1,L),beta),file="sim_trait_4.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
# 
# ## probably want to try one more with more SNPs and pip =.1 say 1000 snps, this would keep individual effects lower
# 
# 
