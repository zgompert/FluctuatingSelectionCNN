#### Single locus allele frequency simulation

### Constant Selection ###
## Set parameters

# population number
nPol <- 10
# number of generations
nGen <- 10
# population size
popSize <- 200
# sample size between generations
sampleSize <- 30
# heterozygous effect
h <- 0.5
#selection coefficient
s <- 0.1

# Initialize the matrix
evoSim <- matrix(NA, nrow = nGen, ncol = 2*nPol)
# create a vector of odd numbers (for populations) and even numbers
# (environmental data)
popCol <- seq(from = 1, to = 19, by = 2)
envCol <- seq(from = 2, to = 20, by = 2)
# fill in row one with starting allele frequency and environmental data
evoSim[1, popCol] <- 0.5

selCof <- seq(from = -0.2, to = 0.2, by = 0.01)


## Simulation

for (pop in popCol) { #repeat the loop for each population
  for (gen in 2:nGen) { #repeat this loop for each generation
    freq <-  evoSim[gen - 1 , pop] #save the allele frequency for the previous generation
    expectedFreq <- (freq + (s*freq)*(1 - freq)*(freq + h*(1-2*freq))) #calculate expected frequency
    newFreq <- rbinom(1, (2*sampleSize), expectedFreq)/(2*sampleSize) #genetic drift
    evoSim[gen, pop] <- newFreq
  }
}

for (env in envCol) {
  for (gen in 1:nGen) {
    evoSim[gen, env] <- sample(selCof, 1) #generate environmental "white noise"
  }
}


### Fluctuating Selection ###
## Set parameters

# population number
nPol <- 10
# number of generations
nGen <- 10
# population size
popSize <- 200
# sample size between generations
sampleSize <- 30
# heterozygous effect
h <- 0.5

# Initialize the matrix
evoSim <- matrix(NA, nrow = nGen, ncol = 2*nPol)
# create a vector of odd numbers (for populations) and even numbers
# (environmental data)
popCol <- seq(from = 1, to = 19, by = 2)
envCol <- seq(from = 2, to = 20, by = 2)
# fill in row one with starting allele frequency and environmental data
evoSim[1, popCol] <- 0.5

# selection coefficient possibilities
selCof <- seq(from = -0.2, to = 0.2, by = 0.01)



##Simulation

for (env in envCol) {
  for (gen in 1:nGen) {
    evoSim[gen, env] <- sample(selCof, 1)
  }
}

for (pop in popCol){
  for (gen in 2:nGen) {
    s <- evoSim[gen - 1, pop + 1] #set selection equal to the environmental data from the previous generation
    freq <-  evoSim[gen - 1 , pop] #save the allele frequency for the previous generation
    expectedFreq <- (freq + (s*freq)*(1 - freq)*(freq + h*(1-2*freq))) #calculate expected frequency
    newFreq <- rbinom(1, (2*sampleSize), expectedFreq)/(2*sampleSize) #genetic drift
    evoSim[gen, pop] <- newFreq
  }
}

print(evoSim)
print("Hello")