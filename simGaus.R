## ind. based simulation of evolution## stabalizing selection in fluct. environment
## approximaint a linear model for S

## read in genarch 1 initial conditions
#p0<-read.table("sim_p0.txt",skip=1,header=FALSE)
trg<-read.table("sim_trait_1.txt",skip=1,header=FALSE)
#ne<-read.table("sim_ne.txt",skip=1,header=FALSE)
#env<-read.table("sim_env.txt",skip=1,header=FALSE)

## row = sim no., col = alpha, beta, mean and SD S, mean and SD |si]
params<-matrix(NA,nrow=40*2,ncol=6)
## row = sim no., col=sum change, cov. change by env
ss<-matrix(NA,nrow=40*2,ncol=2)
N<-1000
Npop<-10
Ngen<-10
Nsnp<-100

## sd of fitness function
#w<-c(10,40)
w<-c(2,5)
aGaus<-c(rep(sum(2*trg[,2]*p0[,1]),20),rep(sum(2*trg[,2]*p0[,1]),20)+.2)
bGaus<-rep(rep(c(0,.4),each=10),2)

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
			dp[k-1,j]<-2 * sum((p1-p) * trg[,2])
			p<-p1
		}
	}

	oo<-lm(as.numeric(S)~as.numeric(as.matrix(env)))
	idx<-x + 40 * (sel-1)

	params[idx,1:2]<-oo$coefficients
	params[idx,3]<-mean(as.numeric(S))
	params[idx,4]<-sd(as.numeric(S))
	params[idx,5]<-mean(abs(as.vector(si)))
	params[idx,6]<-sd(abs(as.vector(si)))

	ss[idx,1]<-sum(as.vector(dp))
	ss[idx,2]<-cov(as.vector(dp),as.vector(as.matrix(env)))
}}

write.table(round(params,5),file="gausParams.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(round(ss,5),file="gausSs.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

