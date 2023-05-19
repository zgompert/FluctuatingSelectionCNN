######  Simulation of single locus for multiple generations and multiple populations for different scenarios  #######


# no of populations
nPop <- 10
# no of generations
nGen <- 10
# population size
popSize <- 100
# heterozygous effect
h<-0.5

###################################################  Scenario 1:     Constant selection coefficient       #####################################


# selection coefficient of 0.02: Hypothetical
selCof<-0.02

# Making an epty matrix and Assigning allele frequency of 0.5 in the first generation of all 10 populations
Popsim<-matrix(NA, nrow=nPop, ncol=nGen)
Popsim[ ,1 ]<-0.5


# simulate allele frequencies over the generations and populations

# number of simulations/matrices required
nSil <-2
ConstSimOut<-list() ## ConstSimOut = simulation outputs for a constant selection coefficient for all populations and generations 

# Main loop
for (ns in 1:nSil){
  for (generation in 2:nGen) {
    for (population in 1:nPop) {
      
      
      # Calculate the frequency of the selected allele in the next generation
      ExpAlelFreq <-Popsim[population, generation-1] + selCof*Popsim[population, generation-1]*(1-Popsim[population, generation-1])*(Popsim[population, generation-1]+0.05*(1-2*(Popsim[population, generation-1])))
      allele_freq<- rbinom(prob=ExpAlelFreq, size =100, 2*popSize)/(2*popSize)
      
      # Store the allele frequency in the matrix
      Popsim[population, generation] <- mean(allele_freq)
    }
  }
  ConstSimOut[[ns]] <- Popsim
}
# Print the allele frequency matrix
print(ConstSimOut)


###################################################   Scenario 2:    Constant selection coefficient with selection coefficient NEXT to the allele frequency      #####################################


# no of populations
nPop <- 10
# no of generations
nGen <- 10
# population size
popSize <- 100
# heterozygous effect
h<-0.5
# selection coefficient of 0.02
selCof<-0.02

# Making an epty matrix and Assigning allele frequency of 0.5 in the first generation of all 10 populations
Popsim<-matrix(NA, nrow=nPop, ncol=2*nGen+1)
Popsim[ ,1 ]<-0.5

# number of simulations/matrices required
nSil <-2

ConstSimOutAdjac<-list() ## ConstSimOut = simulation outputs for a constant selection coefficient for all populations and generations (Allele frequency placed adjacent to the selection coefficient)

# Main loop
for (ns in 1:nSil){
  for (generation in 2:nGen) {
    for (on in seq(1, by=2, length=10)) {
      for (population in 1:nPop){
        # Calculate the frequency of the selected allele in the next generation
        ExpAlelFreq <-Popsim[population, on] + selCof*Popsim[population, on]*(1-Popsim[population, on])*(Popsim[population, on]+0.05*(1-2*(Popsim[population, on])))
        allele_freq<- rbinom(prob=ExpAlelFreq, size =100, 2*popSize)/(2*popSize)
        
        # Store the allele frequency in the matrix
        Popsim[population, on+2] <- mean(allele_freq)
        Popsim[population, on+1] <- 0.02
      }
    }
  }
  ConstSimOutAdjac[[ns]] <- Popsim[ ,1:20]
}
# Print the allele frequency matrix
print(ConstSimOutAdjac)


###################################################   Scenario 3:    VARIABLE selection coefficient with selection coefficient NEXT to the allele frequency      #####################################


# selection coefficient of 0.02
selCof<-seq(-.2, .2, by=0.001)
# Making an epty matrix and Assigning allele frequency of 0.5 in the first generation of all 10 populations
Popsim<-matrix(NA, nrow=nPop, ncol=2*nGen+1)
Popsim[ ,1 ]<-0.5

# number of simulations/matrices required
nSil <-2

VariableSimOut<-list() ## ConstSimOut = simulation outputs for a constant selection coefficient for all populations and generations 

# Main Loop
for (ns in 1:nSil){
  for (generation in 2:nGen) {
    for (on in seq(1, by=2, length=10)) {
      for (population in 1:nPop){
        # Calculate the frequency of the selected allele in the next generation
        ExpAlelFreq <-Popsim[population, on] + sample(selCof, 1)*Popsim[population, on]*(1-Popsim[population, on])*(Popsim[population, on]+0.05*(1-2*(Popsim[population, on])))
        allele_freq<- rbinom(prob=ExpAlelFreq, size =100, 2*popSize)/(2*popSize)
        
        # Store the allele frequency in the matrix
        Popsim[population, on+2] <- mean(allele_freq)
        Popsim[population, on+1] <- sample(selCof, 1)
      }
    }
  }
  VariableSimOut[[ns]] <- Popsim[ ,1:20]
}
# Print the list of allele frequency matrices
print(VariableSimOut)


