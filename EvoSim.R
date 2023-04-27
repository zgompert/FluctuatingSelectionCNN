## functions for biol 3020 project simulations

## gridSize = number of rows and columns, square to get number of populations
## popSize = number of dipliod individuals per population
## mig = migration rate (proportion) between neighboring populations
## nChrom = number of chromosomes
## nLoci = number of loci per chromosome (evenly spaced along chromosomes)
## Sigma = strength of stabilizing selection, larger values denote weaker selection
## spatial = model of spatial variation in seleciton: "none", "gradient", "random"
## temporal = model of temporal variation in selection: "none", "autocor", "random"
## gens = number of generations to simulate
## ssp = SD of selection in space
## sti = SD of selection in time
simevo<-function(gridSize=3, popSize=100, mig=0, nChrom=2, nLoci=10, Sigma=1,
		 spatial="none", temporal="none", gens=1000, ssp = 0.2, sti = 0.2){
	## initialize the populations, p = 0.5 for all pops and loci
	popGrid<-matrix(1:(gridSize^2),nrow=gridSize,ncol=gridSize,byrow=FALSE)
	nPop<-gridSize^2
	G<-vector("list",nPop)
	P<-matrix(NA,nrow=popSize,ncol=nPop)
	for(i in 1:gridSize){for(j in 1:gridSize){
		p<-popGrid[i,j]
		G[[p]]<-vector("list",popSize)
		for(n in 1:popSize){
			G[[p]][[n]]<-vector("list",nChrom)
			for(k in 1:nChrom){
				G[[p]][[n]][[k]]<-matrix(rbinom(n=nLoci*2,size=1,prob=0.5),nrow=nLoci,ncol=2)
				P[n,p]<-mean(unlist(G[[p]][[n]]))-0.5
			}
		}
	}}

	## initialize environment
	baseEnv<-matrix(0,nrow=gridSize,ncol=gridSize)
	if(spatial=="random"){
		baseEnv<-baseEnv+rnorm(nPop,mean=0,sd=ssp)
	}else if(spatial=="gradient"){
		baseEnv<-baseEnv+sort(rnorm(nPop,mean=0,sd=ssp))
	}	
		
	Fst<-rep(NA,gens)
	PropVar<-rep(NA,gens)
	mnHs<-rep(NA,gens)

	## main loop
	for(gg in 1:gens){
		if((gg %% 100)==0){
			cat("Generation",gg,"\n")
		}
		## set current environment
		if(temporal=="none"){
			curEnv<-baseEnv
		} else if(temporal=="random"){
			curEnv<-baseEnv+rnorm(nPop,mean=0,sd=sti)
		} else{
			if(gg==1){curEnv<-baseEnv+rnorm(nPop,mean=0,sd=sti)}
			if(gg>1){curEnv<-curEnv+rnorm(nPop,mean=0,sd=sti)}
		}
		
		## gene flow
		Gold<-G
		Pold<-P
		for(i in 1:gridSize){for(j in 1:gridSize){	
			## exchange to left x right
			if(j < gridSize){
				p1<-popGrid[i,j];p2<-popGrid[i,j+1]
				wmig<-which(rbinom(n=popSize,size=1,prob=mig)==1)
				G[[p1]][wmig]<-Gold[[p2]][wmig]			
				G[[p2]][wmig]<-Gold[[p1]][wmig]			
				P[wmig,p1]<-Pold[wmig,p2]
				P[wmig,p2]<-Pold[wmig,p1]
				Gold<-G
				Pold<-P
			}
			## exchange top bottom
			if(i < gridSize){
				p1<-popGrid[i,j];p2<-popGrid[i+1,j]
				wmig<-which(rbinom(n=popSize,size=1,prob=mig)==1)
				G[[p1]][wmig]<-Gold[[p2]][wmig]			
				G[[p2]][wmig]<-Gold[[p1]][wmig]			
				P[wmig,p1]<-Pold[wmig,p2]
				P[wmig,p2]<-Pold[wmig,p1]
				Gold<-G
				Pold<-P
			}
		}}
	
		## selection
		W<-P
		for(i in 1:gridSize){for(j in 1:gridSize){	
			p<-popGrid[i,j]
			W[,p]<-dnorm(x=P[,p],mean=curEnv[i,j],sd=Sigma)
			W[,p]<-W[,p]/sum(W[,p])
		}}
		## reproduction based on relative fitnesses
		for(p in 1:nPop){
			pa1<-sample(1:popSize,replace=TRUE,prob=W[,p])
			pa2<-sample(1:popSize,replace=TRUE,prob=W[,p])
			for(n in 1:popSize){
				rr<-sample(1:(nLoci-1),size=nChrom,replace=TRUE)
				for(k in 1:nChrom){
					if(runif(1)<0.5){
						G[[p]][[n]][[k]][,1]<-c(Gold[[p]][[pa1[n]]][[k]][1:rr[k],1],Gold[[p]][[pa1[n]]][[k]][(rr[k]+1):nLoci,2])
						G[[p]][[n]][[k]][,2]<-c(Gold[[p]][[pa2[n]]][[k]][1:rr[k],2],Gold[[p]][[pa2[n]]][[k]][(rr[k]+1):nLoci,1])
					} else{
						G[[p]][[n]][[k]][,1]<-c(Gold[[p]][[pa1[n]]][[k]][1:rr[k],2],Gold[[p]][[pa1[n]]][[k]][(rr[k]+1):nLoci,1])
						G[[p]][[n]][[k]][,2]<-c(Gold[[p]][[pa2[n]]][[k]][1:rr[k],1],Gold[[p]][[pa2[n]]][[k]][(rr[k]+1):nLoci,2])
					}
				}
				P[n,p]<-mean(unlist(G[[p]][[n]]))-0.5
			}
		}
		## compute summaries
		PP<-matrix(0,nrow=nChrom*nLoci,ncol=nPop)
		for(p in 1:nPop){
			for(n in 1:popSize){
				PP[,p]<-PP[,p] + unlist(lapply(G[[p]][[n]],apply,1,mean))
			}
			PP[,p]<-PP[,p]/popSize
		}
		Hs<-2 * PP * (1-PP)
		pbar<-apply(PP,1,mean)
		Ht<-2 * pbar * (1-pbar)
		Fst[gg]<-mean(Ht-Hs)/mean(Ht)
		PropVar[gg]<-mean(PP > 0 & PP < 1)
		mnHs[gg]<-mean(Hs)
	}
	
	## create output object
	out<-list(Hs=mnHs,PropVar=PropVar,Fst=Fst,PP=PP,Pheno=P)
	return(out)
}	
