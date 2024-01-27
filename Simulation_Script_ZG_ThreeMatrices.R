#### Single locus allele frequency simulation

## This function simulates constant selection (a != 0), fluctuating selection (b != 0) and neutral evolution (a and b = 0)
## populations = number of populations
## generations = total number of generations, including initial
## p0 = initial allele frequency
## popsize = population size
## samplesize = sample size for observed allele freq, should be <= popsize
## h = het. effect
## a = constant selection parameter
## b = fluctuating selection parameter
## model is s = a + b*env

########################################################################################################################
# Function to generate a separate matrix for allele frequency difference
calc_AllelFreqDiffMatrix <- function(AllelFreqMatrx) {
  n_rows<-nrow(AllelFreqMatrx)
  n_cols<-ncol(AllelFreqMatrx)
  # new empty matrix
  result_sub_matrix<-matrix(NA, nrow = n_rows, ncol = n_cols)
  # looping through each column
  for (col in 1:n_cols) {
    result_sub_matrix[1, col] <-AllelFreqMatrx[1, col]
    # loop through all the rows and subtract values of the previous rows
    for (row in 2:n_rows) {
      result_sub_matrix[row, col] <-abs(AllelFreqMatrx[row, col] - AllelFreqMatrx[row-1, col])
    }
  }
  return(result_sub_matrix)
}

########################################################################################################################
Sim <- function(populations = 10, generations = 11, environments =10, p0 = runif(1, 0.2, 0.8), popsize = 200,
                samplesize = 30, h = 0.5, a = 0, b = 0, output = 1){
  # population number
  nPol <- populations
  # number of generations
  nGen <- generations
  # no. of environments
  nEnv <- environments
  cat("a =",a,";b = ",b,"\n")
  
  # Initialize the matrix
  ### for sepMatrix
  EnvMat <- matrix(NA, nrow = nGen, ncol = nEnv)
  FreqMat <- matrix(NA, nrow = nGen, ncol = nPol)
  
  
  # fill in row one with starting allele frequency and environmental data
  FreqMat[1, ] <- p0
  # 
  PopCol <- seq(from = 1, to = nPol, by = 1)
  EnvCol <- seq(from = 1, to = nEnv, by = 1)
  
  
  ## sample the environment
  for (env in EnvCol) {
    for (gen in 1:nGen) {
      EnvMat[gen, env] <- runif(1, 0, 1)
    }
  }
  
  
  for (pop in PopCol){
    for (gen in 2:nGen) { #set selection equal to a function of environmental data from the previous generation
      s <- a + (b*((EnvMat[gen - 1, pop]) -0.5 ))## center environment for simulations
      freq <-  FreqMat[gen - 1 , pop] #save the allele frequency for the previous generation
      expectedFreq <- freq + (s*freq*(1 - freq)*(freq + h*(1-2*freq))) #calculate expected frequency
      newFreq <- rbinom(1, (2*popsize), expectedFreq)/(2*popsize) #genetic drift based on population size
      FreqMat[gen, pop] <- newFreq
    }
    ## sample size perturbs allele frequences from true values
    FreqMat[,pop]<-rbinom(nGen,(2*samplesize),FreqMat[,pop])/(2*samplesize)
    ## arcsine square root transformation
    ## evoSim[,pop]<-2 * asin(sqrt(evoSim[,pop]))
  }
  
  
  if (output == 1){
    ## recode as change in p, but scaled to be between 0, and 1
    EnvMatNew <- matrix(NA, nrow = nGen - 1, ncol = nEnv)
    FreqMatNew <- matrix(NA, nrow = nGen - 1, ncol = nPol)
    
    for (env in EnvCol) {
      for (gen in 1:nGen - 1) {
        EnvMatNew[gen, env] <- EnvMat[gen, env]
      }
    }
    
    for (pop in PopCol){
      for (gen in 1:nGen - 1){ 
        ## pi is range of p, so 2 pi is range of delta p
        ##evoSimNew[gen, pop] <- (evoSim[gen + 1, pop] - evoSim[gen, pop])/(2*pi) + .5 
        ## SD standardize and scale to 0 to 1 
        FreqMatNew[gen, pop] <- (FreqMat[gen + 1, pop] - FreqMat[gen, pop])/sqrt(FreqMat[gen, pop] * (1-FreqMat[gen, pop]))
        FreqMatNew[gen, pop] <- (FreqMatNew[gen, pop] + 1)/2
      }
    }
  }
  if (output == 2){
    #record allele frequency and environment, but shift the frequency values up one 
    EnvMatNew <- matrix(NA, nrow = nGen - 1, ncol = nEnv)
    FreqMatNew <- matrix(NA, nrow = nGen - 1, ncol = nPol)
    
    for (env in EnvCol) {
      for (gen in 1:nGen - 1) {
        EnvMatNew[gen, env] <- EnvMat[gen, env]
      }
    }
    
    for (pop in PopCol){
      for (gen in 1:nGen - 1){
        FreqMatNew[gen, pop] <- FreqMat[gen + 1, pop]
      }
    }
  }
  return(list(EnvMatNew = EnvMatNew, FreqMatNew = FreqMatNew))
  AllelFreqDiffMatrix<-calc_AllelFreqDiffMatrix(FreqMatNew)
}
#######################################################################################################################################################


#generate 5000 data sets for strong, constant selection, samples a ~ +- .2
for (i in 1:5000){
  data <- Sim(10,11,10,0.5,popsize=200,samplesize=200,h=0.5,a=rbeta(1,10,40)*sample(c(-1,1),1),b=0)
  environmentMatrix <- data$EnvMatNew
  alleleFrequencyMatrix <- data$FreqMatNew
  
  # Write matrices to CSV files
  write.table(environmentMatrix, paste0("environment_matrix_Type1_", i, ".csv"), col.names = F, row.names = F,sep = ",")
  write.table(alleleFrequencyMatrix, paste0("allele_frequency_matrix_Type1_", i, ".csv"), col.names = F, row.names = F,sep = ",")
  write.table(AllelFreqDiffMatrix, paste0("DIFF_allele_frequency_matrix_Type1_", i, ".csv"), col.names = F, row.names = F,sep = ",")
}

#generate 5000 data sets for strong, fluctuating selection, a = 0 (no directional), b ~ +- .3
for (i in 1:5000){
  data <- Sim(10,11,10,0.5,popsize=200,samplesize=200,h=0.5,a=0,b=rbeta(1,15,35)*sample(c(-1,1),1))
  environmentMatrix <- data$EnvMatNew
  alleleFrequencyMatrix <- data$FreqMatNew
  
  # Write matrices to CSV files
  write.table(environmentMatrix, paste0("environment_matrix_Type2_", i, ".csv"), col.names = F, row.names = F,sep = ",")
  write.table(alleleFrequencyMatrix, paste0("allele_frequency_matrix_Type2_", i, ".csv"), col.names = F, row.names = F,sep = ",")
  write.table(AllelFreqDiffMatrix, paste0("DIFF_allele_frequency_matrix_Type2_", i, ".csv"), col.names = F, row.names = F,sep = ",")
}
#generate 5000 data sets for pure drift
for (i in 1:5000){
  data <- Sim(10,11,10,0.5,popsize=200,samplesize=200,h=0.5,a=0,b=0)
  environmentMatrix <- data$EnvMatNew
  alleleFrequencyMatrix <- data$FreqMatNew
  
  # Write matrices to CSV files
  write.table(environmentMatrix, paste0("environment_matrix_Type3_", i, ".csv"), col.names = F, row.names = F,sep = ",")
  write.table(alleleFrequencyMatrix, paste0("allele_frequency_matrix_Type3_", i, ".csv"), col.names = F, row.names = F,sep = ",")
  write.table(AllelFreqDiffMatrix, paste0("DIFF_allele_frequency_matrix_Type3_", i, ".csv"), col.names = F, row.names = F,sep = ",")
}

