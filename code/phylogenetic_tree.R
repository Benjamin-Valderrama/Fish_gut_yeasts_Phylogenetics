# Set up ------------------------------------------------------------------

# Libraries we will use from CRAN
.cran <- c("ape", "phangorn", "dplyr", "stringr")
# Installing those we haven't installed yet
.to_inst_cran <- !.cran %in% installed.packages()
sapply(.cran[.to_inst_cran], install.packages, character.only = T)


# Libraries we will use from Bioconductor
.bioconductor <- c("BiocManager", "ggtree")
# Installing those we haven't installed yet
.to_inst_bioconductor <- !.bioconductor %in% installed.packages()
sapply(.bioconductor[.to_inst_bioconductor], install.packages, character.only = T)


# Loading the libraries
sapply(c(.cran, .bioconductor), library, character.only = T)



# Importing data ----------------------------------------------------------

MSA_path <- paste0(getwd(),"/data/processed/MSA_outgroup.fas")
MSA <- read.FASTA(file = MSA_path)
MSA_phydat <- phyDat(MSA, type = "DNA", levels = NULL)



# Basic tree (distance methods) -------------------------------------------

# Distance matrix
dist_matrix_JC <- dist.ml(x = MSA_phydat)

# Joining algorithms
treeUPGMA<- upgma(dist_matrix_JC)
treeNJ<- NJ(dist_matrix_JC)

# Ladderization
treeUPGMA <- ladderize(treeUPGMA)
treeNJ<- ladderize(treeNJ)

# Basic tree plot
plot(treeUPGMA)
plot(treeNJ)



# Adding metadata ---------------------------------------------------------

# Converting the tree to tibble to add the metadata
treebble <- as_tibble(treeUPGMA)

#' Things to do:
#' 1. Change tip names
#' 2. Adding the source of the sequence
#' 3. Adding a color for each source 


# 1. Change tip names
length_of_labels <- sapply(str_split(treebble$label, " "), length)

treebble <- treebble %>% 
    # Label will have the ID of the sequence
    mutate(label = ifelse(length_of_labels == 2, 
                        word(treebble$label, 1, 1), # Experimental ID
                        word(treebble$label, 2,3))) # Genus and Specie ID from NCBI


#' Things to do:
#' 1. Changin tip names [DONE]
#' 2. Adding the source of the sequence
#' 3. Adding a color for each source 


# 2. Adding the source of the sequence
treebble <- treebble %>% 
  mutate(source = case_when(as.numeric(label) <= 300 ~ "S. violacea",
                            as.numeric(label) > 300 ~ "G. chilensis",
                            str_detect(label, "[^0-9]") ~ "NCBI",
                            TRUE ~ "NA"),
         source = ifelse(source == "NA", NA, source))


#' Things to do:
#' 1. Changin tip names [DONE]
#' 2. Adding the source of the sequence [DONE]
#' 3. Adding a color for each source 


# 3. Adding a color for each source 
treebble <- treebble %>% 
  mutate(color = case_when(as.numeric(label) <= 300 ~ "blue",
                           as.numeric(label) > 300 ~ "orange",
                           str_detect(label, "[^0-9]") ~ "black",
                           TRUE ~ "NA"),
         color = ifelse(source == "NA", NA, color))


#' Things to do:
#' 1. Changin tip names [DONE]
#' 2. Adding the source of the sequence [DONE]
#' 3. Adding a color for each source [DONE]


# Back convertion from tibble to tree
final_treeUPGMA<- as.treedata(treebble)



# Plot personalization ----------------------------------------------------

source <- final_treeUPGMA@data$source

ggtree(final_treeUPGMA) +
  
  geom_treescale() + 
  geom_tiplab(aes(color = source), size = 5, fontface = "bold", show.legend = FALSE) + 
  geom_tippoint(aes(color = source), size = 2, alpha = 0.6) +
  
  scale_color_manual(values = c("Black", "Orange", "Blue"), 
                     labels = c("NCBI", "S. violacea", "G. chilensis")) +
  
  scale_x_continuous(limits = c(0,0.19)) +
  
  theme(legend.position = "bottom",
        legend.text = element_text(face = "bold", size = 18),
        legend.title = element_blank(),
        legend.key.size = unit(2, "cm")) +
  
  guides(colour = guide_legend(override.aes = list(size = 8)))



# Sesion information ------------------------------------------------------

sessionInfo()




# Bootstrap ---------------------------------------------------------------

# Adding bootstrap to the plot