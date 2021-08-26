# Set up ------------------------------------------------------------------
# Libraries we will use
lib <- c("seqinr", "ape", "tidyverse", "phangorn")
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
treeNJ<- NJ(dist_matrix_JC)

treeUPGMA <- ladderize(treeUPGMA)
treeNJ<- ladderize(treeNJ)

# Basic tree plot
plot(treeUPGMA)
plot(treeNJ)

# Bootstrap ---------------------------------------------------------------
set.seed(1)
fit <- pml(tree = treeUPGMA, data = MSA_phydat, model = best_model)
fitGTR <- optim.pml(fit)
bs <- bootstrap.pml(fitGTR, bs=1000, optNni=TRUE,
                   control = pml.control(trace = 0))

# Plot provisional
plotBS(tree = midpoint(fitGTR$tree), 
       BStrees = bs,
       p = 50,
       type= "phylogram",
       bs.col = "gray50")


# Metadata ----------------------------------------------------------------

# Changing the tip names to codes of the study and binomial names
my.tip.label <- fitGTR$tree$tip.label

# Home made function to change the names
rename.label <- function(x){
  # All tip labels are splited by its spaces
  splited_names <- strsplit(x = my.tip.label, split = " ")
  
  final.names <- c()
  # For each splited name, we ask whether they have just 2 elements or more
  for (i in 1:length(splited_names)) {
    # If the have just 2 elements, this sequence is part of the initial study
    if (lengths(splited_names[i]) == 2) {
      # We save the first element as part of the final list of tip names
      study.name <- splited_names[[i]][1]
      final.names <- c(final.names, study.name)
    } # If not, they are from the NCBI
    else{
      # We keep the second and third element as part of the final list of tip names
      ncbi.names <- paste(splited_names[[i]][2], splited_names[[i]][3])
      final.names <- c(final.names, ncbi.names)
    }
  }
  return(final.names)
} 

# Changing the names using the already made function
fitGTR$tree$tip.label <- rename.label(my.tip.label)

# If we look at the provisional plot, now tha names are changed
plotBS(tree = midpoint(fitGTR$tree), 
       BStrees = bs,
       p = 50,
       type= "phylogram",
       bs.col = "gray50")

# Check... plot.phylo


# Plot personalization ----------------------------------------------------



# Sesion information ------------------------------------------------------
sessionInfo()
