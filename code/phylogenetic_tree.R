
# Ensure reproducibility --------------------------------------------------

#' I the renv.lock file you can find the list of the packages 
#' (and their versions) used in this project.
#' 
#' To install and load all these packages, you can install the renv
#' package, load it with library() and call :
#' 
#' renv::restore()


# Set up ------------------------------------------------------------------

# Libraries we will use from CRAN
.cran <- c("ape", "phangorn", "dplyr", "stringr", "ggplot2", "renv")
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
dist_matrix_JC <- dist.ml(x = MSA_phydat, model = "JC69")

# Joining algorithms
treeUPGMA<- upgma(dist_matrix_JC)
treeNJ<- NJ(dist_matrix_JC)

# Ladderization
treeUPGMA <- ladderize(treeUPGMA)
treeNJ<- ladderize(treeNJ)

# Basic tree plot
plot(treeUPGMA)
plot(treeNJ)


# Bootstrap ---------------------------------------------------------------

# Calculating bootstrap values
fit <- pml(tree = treeUPGMA, data = MSA_phydat, model = "JC69") # Try with other models

# Optimization of the model parameters
fitJC69 <- optim.pml(fit)

# Bootstrap with 1000 repetitions
bs <- bootstrap.pml(fitJC69, bs=1000, optNni=TRUE,
                    control = pml.control(trace = 0))

# Adding the bootstrap values to the tree object
bootstraped_tree<- plotBS(tree = midpoint(fitJC69$tree), 
                          BStrees = bs,
                          p = 0,
                          type = "n")

# Getting the bootstrap values in an independent vector
bootstrap.vals <- bootstraped_tree$node.label


# Adding metadata ---------------------------------------------------------

# Converting the tree to tibble to add the metadata
treebble <- as_tibble(treeUPGMA)

#' Things to do:
#' 1. Change tip names
#' 2. Adding the source of the sequence
#' 3. Adding a color for each source 
#' 4. Adding bootstrap values
#' 5. Adding colors to bootstrap values based on their significance


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
#' 4. Adding bootstrap values
#' 5. Adding colors to bootstrap values based on their significance


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
#' 4. Adding bootstrap values
#' 5. Adding colors to bootstrap values based on their significance


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
#' 4. Adding bootstrap values
#' 5. Adding colors to bootstrap values based on their significance


#' 4. Adding bootstrap values 
#' 5. Adding colors to bootstrap values based on their significance
treebble <- treebble %>% 
  # Adding the bootrstap values to the tree
  mutate(bootstrap.vals = c(rep(NA, nrow(treebble)-length(bootstrap.vals)), bootstrap.vals),
         # Adding colors to the bootstrap values
         bootstrap.cols = ifelse(bootstrap.vals >= 70, "grey25", "grey60"))


#' Things to do:
#' 1. Changin tip names [DONE]
#' 2. Adding the source of the sequence [DONE]
#' 3. Adding a color for each source [DONE]
#' 4. Adding bootstrap values [DONE]
#' 5. Adding colors to bootstrap values based on their significance [DONE]


# Back convertion from tibble to tree
final_treeUPGMA<- as.treedata(treebble)


# Plot personalization ----------------------------------------------------

# Check the defaults plot with the bootstrap values using phangorn's plotBS 
plotBS(tree = midpoint(fitJC69$tree), 
       BStrees = bs,
       p = 0,
       type= "phylogram",
       bs.col = "gray50")


# Making our own tree with ggtree's functions for deeper customization : 

# Ploting the resulting tree
ggtree(final_treeUPGMA) +
  
  #geom_treescale() + 
  geom_tiplab(aes(color = final_treeUPGMA@data$source), 
              size = 5, fontface = "bold", show.legend = FALSE) + 
  
  geom_tippoint(aes(color = source), size = 2, alpha = 0.6) +
  
  geom_text(label = final_treeUPGMA@data$bootstrap.vals,
            color = final_treeUPGMA@data$bootstrap.cols,
            hjust = 1.2, vjust = -0.5) +
  
  scale_color_manual(values = c("Blue", "Black", "Orange"), 
                     labels = c("G. chilensis", "NCBI", "S. violacea")) +
  
  scale_x_continuous(limits = c(0,0.19)) +
  
  theme(legend.position = "bottom",
        legend.text = element_text(face = "bold", size = 18),
        legend.title = element_blank(),
        legend.key.size = unit(2, "cm")) +
  
  guides(colour = guide_legend(override.aes = list(size = 8)))


# Creating a folder (if it doesn't exist yet) to store the plots
if(!dir.exists("figures")){
  dir.create("figures")
}

# Saving the plot in the figures folder
ggsave(filename = "figures/ggtree.tiff",
       last_plot(), 
       width = 12, height = 12)


# Sesion information ------------------------------------------------------

sessionInfo()

# R version 4.2.0 (2022-04-22 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 17763)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=Spanish_Chile.1252  LC_CTYPE=Spanish_Chile.1252    LC_MONETARY=Spanish_Chile.1252
# [4] LC_NUMERIC=C                   LC_TIME=Spanish_Chile.1252    
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggtree_3.4.0        BiocManager_1.30.18 ggplot2_3.3.6       stringr_1.4.0      
# [5] dplyr_1.0.9         phangorn_2.9.0      ape_5.6-2           tidytree_0.3.9     
# 
# loaded via a namespace (and not attached):
#   [1] treeio_1.20.0      tidyselect_1.1.2   purrr_0.3.4        lattice_0.20-45   
# [5] ggfun_0.0.6        colorspace_2.0-3   vctrs_0.4.1        generics_0.1.2    
# [9] utf8_1.2.2         gridGraphics_0.5-1 rlang_1.0.2        pillar_1.7.0      
# [13] glue_1.6.2         withr_2.5.0        DBI_1.1.3          lifecycle_1.0.1   
# [17] munsell_0.5.0      gtable_0.3.0       ragg_1.2.2         codetools_0.2-18  
# [21] labeling_0.4.2     parallel_4.2.0     fansi_1.0.3        Rcpp_1.0.8.3      
# [25] scales_1.2.0       jsonlite_1.8.0     farver_2.1.0       systemfonts_1.0.4 
# [29] textshaping_0.3.6  fastmatch_1.1-3    digest_0.6.29      aplot_0.1.6       
# [33] stringi_1.7.6      grid_4.2.0         quadprog_1.5-8     cli_3.3.0         
# [37] tools_4.2.0        yulab.utils_0.0.5  magrittr_2.0.3     lazyeval_0.2.2    
# [41] patchwork_1.1.1    tibble_3.1.7       crayon_1.5.1       tidyr_1.2.0       
# [45] pkgconfig_2.0.3    ellipsis_0.3.2     Matrix_1.4-1       ggplotify_0.1.0   
# [49] assertthat_0.2.1   rstudioapi_0.13    R6_2.5.1           igraph_1.3.1      
# [53] nlme_3.1-157       compiler_4.2.0    