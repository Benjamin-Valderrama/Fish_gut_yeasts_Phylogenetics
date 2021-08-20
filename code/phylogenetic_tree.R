# Set up ------------------------------------------------------------------
# Libraries we will use
lib <- c("seqinr", "ape", "tidyverse", "phangorn", "parallel")
# Installing those we haven't installed yet
inst <- !lib %in% installed.packages()
sapply(lib[inst], install.packages, character.only = T)
# Loading the libraries
sapply(lib, require, character.only = T)


# Importing data ----------------------------------------------------------
MSA_path <- paste0(getwd(),"/data/processed/MSA_outgroup.fas")
MSA <- read.FASTA(file = MSA_path)
MSA_phydat <- phyDat(MSA, type = "DNA", levels = NULL)


# Basic tree (distance methods) -------------------------------------------
# We test all the models to find the best
model_test <- modelTest(MSA_phydat)
# We keep the best model (the one with the lowest AICc score)
best_model <- model_test$Model[which.min(model_test$AICc)]

# Distance matrix
dist_matrix_JC <- dist.ml(x = MSA_phydat)
treeUPGMA<- upgma(dist_matrix_JC)
treeUPGMA <- ladderize(treeUPGMA)

# Basic tree plot
plot(treeUPGMA)


# Sesion information ------------------------------------------------------
sessionInfo()
