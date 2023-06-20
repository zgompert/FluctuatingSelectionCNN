#### Single locus allele frequency simulation

### Constant Selection ###
## Set parameters
Constant_Selection <- function(populations = 10, generations = 11, popsize = 200,
                               samplesize = 30, h = 0.5, s = runif(1, -0.2, 0.2)) {
# population number
nPol <- populations
# number of generations
nGen <- generations
# population size
popSize <- popsize
# sample size between generations
sampleSize <- samplesize
# heterozygous effect
h <- h
#selection coefficient
s <- s

# Initialize the matrix
evoSim <- matrix(NA, nrow = nGen, ncol = 2*nPol)
# create a vector of odd numbers (for populations) and even numbers
# (environmental data)
popCol <- seq(from = 1, to = 2*nPol - 1, by = 2)
envCol <- seq(from = 2, to = 2*nPol, by = 2)
# fill in row one with starting allele frequency and environmental data
evoSim[1, popCol] <- 0.5



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
    evoSim[gen, env] <- runif(1, 0, 1)  #generate environmental "white noise"
  }
}
evoSimNew <- matrix(NA, nrow = nGen - 1, ncol = 2*nPol)

for (env in envCol) {
  for (gen in 1:nGen - 1) {
    evoSimNew[gen, env] <- evoSim[gen, env]
  }
}


for (pop in popCol){
  for (gen in 1:nGen-1){
    evoSimNew[gen, pop] <- (evoSim[gen + 1, pop]) - evoSim[gen, pop]
  }
}

return(evoSimNew)
}

evoSim <- Constant_Selection()

### Fluctuating Selection ###
## Set parameters
Fluctuating_Selection <- function(populations = 10, generations = 11, popsize = 200,
                                  samplesize = 30, h = 0.5){
# population number
nPol <- populations
# number of generations
nGen <- generations
# population size
popSize <- popsize
# sample size between generations
sampleSize <- samplesize
# heterozygous effect
h <- h

# Initialize the matrix
evoSim <- matrix(NA, nrow = nGen, ncol = 2*nPol)
# create a vector of odd numbers (for populations) and even numbers
# (environmental data)
popCol <- seq(from = 1, to = 2*nPol - 1, by = 2)
envCol <- seq(from = 2, to = 2*nPol, by = 2)
# fill in row one with starting allele frequency and environmental data
evoSim[1, popCol] <- 0.5

# selection coefficient possibilities



##Simulation

for (env in envCol) {
  for (gen in 1:nGen) {
    evoSim[gen, env] <- runif(1, 0, 1)
  }
}

a <- runif(1, -0.1, 0.1)
b <- runif(1, -0.3, 0.3)

for (pop in popCol){
  for (gen in 2:nGen) { #set selection equal to a function of environmental data from the previous generation
    s <- a + (b*((evoSim[gen - 1, pop + 1]) -0.5 ))
    freq <-  evoSim[gen - 1 , pop] #save the allele frequency for the previous generation
    expectedFreq <- (freq + (s*freq)*(1 - freq)*(freq + h*(1-2*freq))) #calculate expected frequency
    newFreq <- rbinom(1, (2*sampleSize), expectedFreq)/(2*sampleSize) #genetic drift
    evoSim[gen, pop] <- newFreq
  }
}

evoSimNew <- matrix(NA, nrow = nGen - 1, ncol = 2*nPol)

for (env in envCol) {
  for (gen in 1:nGen - 1) {
    evoSimNew[gen, env] <- evoSim[gen, env]
  }
}


for (pop in popCol){
  for (gen in 1:nGen-1){
    evoSimNew[gen, pop] <- (evoSim[gen + 1, pop]) - evoSim[gen, pop]
  }
}

return(evoSimNew)
}
#generate 5000 data sets for constant selection
for (i in 1:5000){
  data <- Constant_Selection()
  write.table(data, paste0("/uufs/chpc.utah.edu/common/home/gompert-group4/projects/fluctCNN/CNNTrainingData/NewType1DataSet_", i,".csv"),
            col.names = F, row.names = F)
}
#generate 5000 data sets for fluctuating selection
for (i in 1:5000){
  data <- Fluctuating_Selection()
  write.table(data, paste0("/uufs/chpc.utah.edu/common/home/gompert-group4/projects/fluctCNN/CNNTrainingData/NewType2DataSet_", i,".csv"),
            col.names = F, row.names = F)
}


evo <- Fluctuating_Selection()
evo
