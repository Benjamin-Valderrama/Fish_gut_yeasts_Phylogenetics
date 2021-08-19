# Set up ------------------------------------------------------------------

# Libraries we will use
lib <- c("seqinr", "ape")
# Installing those we haven't installed yet
inst <- !lib %in% installed.packages()
sapply(lib[inst], install.packages, character.only = T)
# Loading the libraries
sapply(lib, require, character.only = T)


# Importing data ----------------------------------------------------------
# Import Multiple Sequence Alignment
alignment_path <- paste0(getwd(),"/data/processed/sequences_outgroup_alignment.fas")
alignment <- read.alignment(file = alignment_path, format = "fasta")
alignment_dnabin <- as.DNAbin.alignment(alignment)

# Making the tree ---------------------------------------------------------
# Calculating distance matrix
distance_T92 <- dist.dna(x = alignment_dnabin,
                         model = "T92") # According to the file in the 'exploratory' folder

tree <- ladderize(nj(distance_T92))


# Ploting the tree --------------------------------------------------------


