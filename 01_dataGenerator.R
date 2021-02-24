#### Header ####
## Name: Acacia Ackles; Kate Skocelas; Julia [??]
## Purpose: This script generates the "data" for the project. 
## See our project explanation for more details at [filename]


#### Number of Datapoints ####
# Assume a balanced design
# 2 kinds of genes - A and B
# Arbitrary number of replicates performed
n.reps <- 50

#Total number of data points
n <- n.types*n.reps

# For a binomial distribution, we need to know the number of organisms in each replicate
# In our work we usually use the same number of organsims in each replicate
n.orgs <- rep(100,n)


#### Data Simulation ####

# Specify a categorical variable which indicates genotype
genA <- factor(rep(c(0,1), each=n/2)) # produces 0000...1111....
genB <- factor(rep(rep(c(0,1), each=n/4),2)) # produces 00...11...00...11...

# We also need a vector to indicate mass for each organism
# Assume this is a continuous variable between 1 and 10 grams
mass <- round(runif(n, min=1, max=10), digits=2)

# Chose the values for the parameters (logit transformed)
# Labeled to make it easier to think about
# Each element in beta.vec.names indicates the effect of that variable or interaction of variables
beta.vec.names <- c("genA", "genB", "genA:genB", "mass")
beta.vec <- c(-0.4, 0.1, 0.2, 0.2) # we can freely edit these
names(beta.vec) <- beta.vec.names


#### Model Matrix Creation ####
# Build the design matrix of the interactive combination of genotype and mass
Xmat = model.matrix(~type*mass)

#### Create Stochastic Data ####
#Generate the linear predictor (mutliple Xmat by the beta.vec)
lin.pred =  Xmat%*%beta.vec

#Transform the data with an inverse logit to get the expected proportion of cancerous samples
exp.p <- exp(lin.pred)/(1+exp(lin.pred))

# Add binomial noise
cancer.counts <- rbinom(n=n, size=n.orgs, prob=exp.p)

#### Combine Data ####
# Combine type data, mass data, and cancer counts
df <- data.frame(type, mass, cancer.counts)


#### Export ####
# Export the data to a csv file
filename <- "cancer_data.csv"
write.csv(df, file=filename, row.names = FALSE)
