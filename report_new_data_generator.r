---
title: "Project"
author: "Acacia, Kate, Julia"
date: "3/23/2021"
output: 
      html_document: default
      pdf_document: default
---

```{r setup, include=FALSE}
set.seed(1)
library(knitr)
opts_chunk$set(tidy.opts=list(width.cutoff=80),tidy=TRUE)            
```

## Scenario
Genes A and B have both been shown to be related to lifetime chance of cancer. We want to know if an organism with gene A has the same chance of cancer as an organism with gene B, and how that differs when an organism has both genes A and B, or neither gene. We need to include the covariate of the organism’s mass, because we know this also affects the lifetime chance of cancer. 

We sample N organisms of each of the four genotypes, and we record their body mass and whether or not cancer is present in that organism. 

### Story

Your lab is studying cancer-related genes using feline mammary carcinomas (FMCs) to both validate these tumors as models for human breast cancer (HBC) studies and to improve small animal veterinary practice. FMCs have been emerging as valuable models for human breast cancer, and the domestic cat is highly affected by spontaneous mammary tumors (Ferreira et al., 2019).
 
The cancer-related genes A and B are conserved between cat and human. Both genes have been shown to be related to lifetime incidence of HBC, but they have not been studied in cats. You want to know if a cat with gene A has the same lifetime chance of developing a FMC as a cat with gene B, and how that differs when an organism has both genes A and B, or neither gene (wild type).
 
A cat’s breed, sex, age, and whether or not it is intact all greatly influence their lifetime risk of FMCs (Ohio State University Veterinary Medical Center, 2021). To control for these factors, you use only intact female domestic shorthair cats that are 10 years and older (which show the highest incidence of FMCs) as study subjects. A cat’s mass affects their lifetime chance of cancer as well. 

Working with the MSU Small Animal Clinic, you recruit 1000 participants for the study: 250 with gene A, 250 with gene B, 250 with both genes A and B, and 250 with neither gene. For each cat, you record its genotype, mass at time of death and whether or not it had a FMC in its lifetime. Your null hypothesis is that the lifetime chance of developing a FMC is equal across all four genotypes. Your alternative hypothesis is that the lifetime chance of developing a FMC is different between genotypes.


Ferreira D, Martins B, Soares M, Correia J, Adega F, et al. (2019) Gene expression association study in feline mammary carcinomas. PLOS ONE 14(8): e0221776. https://doi.org/10.1371/journal.pone.0221776

Ohio State University Veterinary Medical Center (2021) Feline Mammary Tumors. Retrieved from https://vet.osu.edu/vmc/companion/our-services/oncology-and-hematology/common-tumor-types/feline-mammary-tumors


### Hypothesis

_Null Hypothesis_: The probability of cancer is equal across all genotypes.

_Alternative Hypothesis 1 (no interaction effect)_: A, B, and AB all show a higher probability of cancer than WT, where A+B = AB.

_Alternative Hypothesis 2 (interaction effect increases probability)_: A, B, and AB all show a higher probability of cancer than WT, where A + B > AB.

_Alternative Hypothesis 3 (interaction effect decreases probability)_: A, B, and AB all show a higher probability of cancer than WT, where A + B < AB.



### Variables

_Genotypes_: WT (neither A nor B), A, B, and AB. A and B each have an incidence rate for cancer.

_Mass_: a continuous variable between 7 and 15 lbs

_Genotype A x Genotype B_: Interaction effects

_Response variable_: Presence or absence of cancer in the organism (0/1) 

_Predictor variables_: Genotype A, Genotype B, interaction between Genotypes A & B, and organism mass 


## Modeling and Justification


**Binomial Distribution**

$y \sim Bin(p,N)$

- $p$ is probability of cancer
- $N$ is total number of individuals

_Justification_: We measure cancer as either present/ absent in each organism. The data simulation evolves N multicellular individuals. We have fixed N number of trials that we can repeat _ad infinitum_. This is a frequentists' approach. Hence, we chose the Binomial distribution.

**Deterministic Function (model) & Joint Probability (likelihood)**

$counts = \alpha + \beta_1 * genotype_A + \beta_2 * genotype_B + \delta * genotype_A * genotype_B + \gamma * mass + \epsilon$

- $counts$ = expected number of organisms with cancer present
- $\alpha$ = (intercept) -- expected incidence rate of cancer in wildtype reference
- $\beta_1$ = constant for $genotype_A$ (how much gene A influences incidence of cancer)
- $\beta_2$ = constant for $genotype_B$ (how much gene B influences incidence of cancer)
- $\delta$ = slope term indicating the interaction between genes A and B (how much genes A and B combined influence incidence of cancer)
- $\gamma$ = constant for mass (how much organism mass influences incidence of cancer)
- $\epsilon$ = randomness pulled from binomial distribution
  + $\epsilon \sim Bin(p,N)$

---

## Simulate the Data
```{r packages}
# Load packages
library(dplyr) # for data
library(tidyr) # also data
library(ggplot2) # for plotting
library(GGally) # for ggpairs function, pairwise plotting of multiple variables
library(lme4) # for GLMMs
library(AICcmodavg) # for AIC comparison
```
```{r simulate}
### Same random seed
set.seed(12345)

#### Number of Datapoints ####
n.orgs <- 20      # cats fitting our criteria per veterinary hospital
n.groups <- 100   # veterinary hospitals participating in our study
n <- n.orgs * n.groups # number of cats in our study

#### Data Simulation ####
# Specify a categorical variable which indicates genotype
#presence_genA <- rep(rep(c(0,1), each=n/4),2) # produces 00...11...00...11...
#presence_genB <- rep(c(0,1), each=n/2) # produces 0000...1111....
#types <- c("WT", "A", "B", "AB")
#genotype <- data.frame(presence_genA, presence_genB, named=factor(rep(types, each=n/4)))

#presence_genA <- rep(rep(c(0,1), each=n/4),2) # produces 00...11...00...11...
#presence_genB <- rep(c(0,1), each=n/2) # produces 0000...1111....
#types <- c("WT", "A", "B", "AB")
#genotype_names <- data.frame(presence_genA, presence_genB, named=factor(rep(types, each=n/4)))
#shuffle_genotype <- genotype_names[sample(1:nrow(genotype_names)), ]   # shuffle the rows in genotype to mix up what genotypes are found in which vet hospital
#genotype <- data.frame(hospital_id=rep(seq(1,n.groups), each=4), shuffle_genotype)

presence_genA <- rep(rep(c(0,1), each=n.groups/4),2) # produces 00...11...00...11...
presence_genB <- rep(c(0,1), each=n.groups/2) # produces 0000...1111....
types <- c("WT", "A", "B", "AB")
genotype_names <- data.frame(presence_genA, presence_genB, named=factor(rep(types, each=n.groups/4)))
shuffle_genotype <- genotype_names[sample(1:nrow(genotype_names)), ]   # shuffle the rows in genotype to mix up what genotypes are found in which vet hospital
genotype <- data.frame(hospital_id=rep(seq(1,n.groups), each=1), shuffle_genotype)


### Predictor variables ###
#mass <- round(runif(n, min=7, max=15), digits=2)    # Assume mass is a continuous variable between 7 and 15 pounds
#temp_genA <- round(runif(n, min=0, max=1), digits=4)   # gene A functional protein in blood nanograms/mL
#temp_genB <- round(runif(n, min=0, max=1), digits=4)   # gene B functional protein in blood nanograms/mL

mass <- round(runif(n.groups, min=7, max=15), digits=2)    # Assume mass is a continuous variable between 7 and 15 pounds
temp_genA <- round(runif(n.groups, min=0, max=10), digits=4)   # gene A functional protein in blood nanograms/mL
temp_genB <- round(runif(n.groups, min=20, max=40), digits=4)   # gene B functional protein in blood nanograms/mL

# multiply by our shuffle order genotypes to obtain randomized samples
genA <- shuffle_genotype$presence_genA * temp_genA  
genB <- shuffle_genotype$presence_genB * temp_genB

# Chose the values for the parameters (logit transformed)
# Labeled to make it easier to think about
# Each element in beta.vec.names indicates the effect of that variable or interaction of variables
beta.vec.names <- c("WT", "A", "B", "AB", "mass")
# this is nonlinear transformation on the logit scale, so 0.02 does NOT mean 2% increase in prob. of cancer!!!
# small change in these params will cause MAJOR changes because on logit scale
# A:B is the amount of EXTRA change on top of A+B, so anything over zero is an interaction effect
beta.vec <- c(0, 0.3, 0.2, 0.05, 0.05) 
names(beta.vec) <- beta.vec.names
#### Model Matrix Creation ####
# Build the design matrix of the interactive combination of genotype and mass
Xmat = model.matrix(~genA*genB + mass)
#### Create Stochastic Data ####
#Generate the linear predictor (mutliple Xmat by the beta.vec)
lin.pred =  Xmat[,]%*%beta.vec
#Transform the data with an inverse logit to get the expected proportion of cancerous samples
exp.p <- exp(lin.pred)/(1+exp(lin.pred))
# Add binomial noise
cancer <- rbinom(n=n.groups, size=n.orgs, prob=exp.p)  
no_cancer <- rep(n.orgs, n.groups) - cancer
percent_with_cancer <- cancer/rep(n.orgs, n.groups)
#### Combine Data ####
# Combine type data, mass data, and cancer counts
df <- data.frame(hospital_id=genotype$hospital_id,genA, genB, mass, cancer, no_cancer, percent_with_cancer, genotype=genotype$named)
df

ggpairs(df)
plot(percent_with_cancer ~ mass, df)
plot(percent_with_cancer ~ genA, filter(df, genotype=='A'))
plot(percent_with_cancer ~ genB, filter(df, genotype=='B'))
plot(percent_with_cancer ~ genA*genB, filter(df, genotype=='AB'))
#### Export ####
# Export the data to a csv file
#filename <- "cancer_data.csv"
#write.csv(df, file=filename, row.names = FALSE)
```

---

## Parameter estimation method

```{r models}
# non.cancer <- n-cancer
# binom.prop <- cbind(cancer, non.cancer)
# null.model <- glm(binom.prop ~ mass, data=df,binomial(link="logit"))
# no.interaction.model <- glm(binom.prop  ~ genA + genB + mass, data=df,binomial(link="logit"))
# interaction.model <- glm(binom.prop  ~ genA*genB + mass, data=df,binomial(link="logit"))
# summary(null.model)
# summary(no.interaction.model)
# summary(interaction.model)
null.model <- glm(cancer ~ 1, data=df,binomial(link="logit"))
mass.model <- glm(cancer ~ mass, data=df,binomial(link="logit"))
no.interaction.model <- glm(cancer  ~ genA + genB + mass, data=df,binomial(link="logit"))
interaction.model <- glm(cancer  ~ genA*genB + mass, data=df,binomial(link="logit"))
summary(null.model)
summary(mass.model)
summary(no.interaction.model)
summary(interaction.model)
```
---
## Statistical Test

```{r comparison}
models <- list(null.model, mass.model, no.interaction.model, interaction.model)
model.names <- c("null", "mass","additive", "interaction")

tab = AICctab(null.model, mass.model, no.interaction.model, interaction.model, base=TRUE, delta=TRUE, weights=TRUE, logLik=TRUE)
class(tab) = 'data.frame'
plot(logLik ~ df, data=tab)
plot(AICc ~ df, data=tab)

AICtable<-aictab(cand.set = models, modnames = model.names)

# Generate a sequence of new x values, evenly distributed across the range
x1new = rbinom(n=n, size=1, prob=exp.p)#seq(0, 1, length=100)

# Calculate predictions by hand
p_null = exp( null.model$coefficients[["(Intercept)"]] + 0 * x1new )      # null predictions
p_mass = exp( mass.model$coefficients[["(Intercept)"]] + 
            mass.model$coefficients[["mass"]] * x1new )      # mass predictions
p_add = exp( no.interaction.model$coefficients[["(Intercept)"]] + 
            no.interaction.model$coefficients[["genA1"]] * x1new +
            no.interaction.model$coefficients[["genB1"]] * x1new +
            no.interaction.model$coefficients[["mass"]] * x1new)  # additive predictions
p_inter = exp( interaction.model$coefficients[["(Intercept)"]] + 
            interaction.model$coefficients[["genA1"]] * x1new +
            interaction.model$coefficients[["genB1"]] * x1new +
            interaction.model$coefficients[["mass"]] * x1new + 
            interaction.model$coefficients[["genA1:genB1"]] * x1new)  # interaction predictions

# Calculate model averaged predictions
pAVG = p_null * AICtable$Cum.Wt[1] + 
  p_mass * AICtable$Cum.Wt[4] + 
  p_add * AICtable$Cum.Wt[3] + 
  p_inter * AICtable$Cum.Wt[2]

# Plot observations and the three sets of predictions
plot(cancer ~ mass, data=df)
lines(p_null ~ x1new, lwd=2, col="red")
lines(p_mass ~ x1new, lwd=2, col="blue")
lines(p_add ~ x1new, lwd=2, col="black")
lines(p_inter ~ x1new, lwd=2, col="yellow")
lines(pAVG ~ x1new, lwd=2, col='royalblue')


```
---
## Plots

```{r, echo=FALSE}
cancer.table <- df %>%
  group_by(genotype) %>%
  count(cancer) %>%
  spread(cancer, n) %>%
  rename(none="0", cancer="1")
cancer.table
# 
ggplot(data=df, aes(x=mass, y=percent_with_cancer)) + geom_point()
```

---
## Analysis of Effect & Uncertainty Sizes

```{r}
# W3_execise_key
df$mass_scaled = scale(df$mass)
df$genA_scaled = scale(df$genA)
df$genB_scaled = scale(df$genB)

mf = glmer(cbind(cancer, no_cancer) ~ genA + genB + mass + (1|hospital_id),
	family='binomial', data=df)
summary(mf)

mf_scaled = glmer(cbind(cancer, no_cancer) ~ genA_scaled + genB_scaled + mass_scaled + (1|hospital_id),
	family='binomial', data=df)
summary(mf)


```

---
## Interpretation of Model Results
