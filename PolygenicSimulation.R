#Simulate Necessary files for simulation

L<-100 ## 100 snps
Npop<-10 ## 10 pops.
Gen<-11 ## 10 gens
ss<-matrix(NA,nrow=40*2,ncol=2)
N<-1000
Npop<-10
Ngen<-11
Nsnp<-100

## sample p from beta, assumes no population structure
p<-rbeta(L,.5,.5)
## convert to maf
p[p>.5]<-1-p[p>.5]

## write results
p0<-matrix(rep(round(p,4),N),nrow=L,ncol=Npop,byrow=FALSE)

## write Ne, 100 samples from posterior, same beta for each
u<-600
l<-400
ne<-matrix(round(rbeta(100*N,10,10) * (u-l) + l,2),nrow=100,ncol=N)

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

trg<-read.table("sim_trait_1.txt",header=FALSE)

## sd of fitness function
#w<-c(10,40)
w<-c(2,5)
aGaus<-c(rep(sum(2*trg[,2]*p0[,1]),20),rep(sum(2*trg[,2]*p0[,1]),20)+.2)
bGaus<-rep(rep(c(0,.4),each=10),2)


for (sim in 2:250){
  ## random environment from standard normal
  env<-matrix(round(rnorm((Gen-1) * Npop,0,1),3),nrow=Gen-1,ncol=Npop)
  write.table(env, paste0("/uufs/chpc.utah.edu/common/home/gompert-group4/projects/fluctCNN/CNN_Polygenic_Data/Environment/DataSet_", sim, ".csv"),
              row.names=F, col.names = F)
  
  for(sel in 1:2){
    for(x in 1:40){
      
      
      ## snp selection 
      si<-array(dim=c(Ngen-1,Npop,Nsnp))
      
      ## summary stat containers
      dp<-matrix(NA,nrow=Ngen-1,ncol=Npop)
      
      
      ## sim loop
      p<-rep(NA,Nsnp)
      p1<-rep(NA,Nsnp)
      G<-matrix(NA,nrow=Nsnp,ncol=N)
      ph<-rep(NA,N)
      Wz<-rep(NA,N)
      S<-matrix(NA,nrow=Ngen-1,ncol=Npop)
      h2 <- 2*sum(trg[,2]^2 * p0[,1]* (1-p0[,1]))
      
      for(j in 1:Npop){
        p<-p0[,j]
        Ne<-floor(sample(ne[,j],size=1))
        for(k in 2:Ngen){
          for(l in 1:Nsnp){
            ## sample genotypes
            G[l,]<-rbinom(n=N,size=2,prob=p[l])
          }	
          ## compute phenotypes		
          ph<-t(G) %*% trg[,2]
          ph<-ph + rnorm(N,mean=0,sd=sqrt(1-h2))
          ## determine fitness function
          mn<-aGaus[x] + env[k-1,j] * bGaus[x]	
          Wz<-dnorm(x=ph,mean=mn,sd=w[sel])
          rWz<-Wz/sum(Wz)
          ## sample contributions to the next generation
          fitness<-sample(1:N,size=Ne,replace=TRUE,prob=rWz)
          p1<-apply(G[,fitness],1,mean)/2
          S[k-1,j]<-mean(ph[fitness])-mean(ph)
          ## calculate si
          si[k-1,j,]<-trg[,2] * as.numeric(S[k-1,j]/var(ph))
          ## calculate dp for bv
          dp[k-1,j]<-2 * sum((p1) * trg[,2])
          p<-p1
        }
      }
      write.table(dp, paste0("/uufs/chpc.utah.edu/common/home/gompert-group4/projects/fluctCNN/CNN_Polygenic_Data/Allele Frequency/DataSet_", sim, ".csv"),
                  row.name = F, col.names = F)
      
      #oo<-lm(as.numeric(S)~as.numeric(as.matrix(env)))
      #idx<-x + 40 * (sel-1)
      
      #params[idx,1:2]<-oo$coefficients
      #params[idx,3]<-mean(as.numeric(S))
      #params[idx,4]<-sd(as.numeric(S))
      #params[idx,5]<-mean(abs(as.vector(si)))
      #params[idx,6]<-sd(abs(as.vector(si)))
      
      #ss[idx,1]<-sum(as.vector(dp))
      #ss[idx,2]<-cov(as.vector(dp),as.vector(as.matrix(env)))
    }}
  print(sim)
}
