######### STEP 1: environment preparation ######### 

### De-comment and run the following three lines if your system does not feature the following required libraries
# BiocManager::install("phyloseq")
# BiocManager::install("DESeq2")
# install.packages("ggplot2")
# install.packages("ggh4x")

### load the libraries
library("phyloseq")
library("DESeq2")
library("ggplot2")
library("ggh4x")

source("../CircoTax_DESeq2.R")  # loads the function "CircoTax_DESeq2", or...
#source("../CircoTax_ALDEx2.R") # loads the function "CircoTax_ALDEx2"
# NB: this tutorial assumes that your current working directory is the "Tutorial" folder downloaded from GitHub while the auto_DA functions are in its parent directory (as in CircoTax's GitHub repository)


### loads the R data containing the example phyloseq object "data" from https://doi.org/10.1186/s13293-023-00523-w
load("data.RData")

###look at possible contrasts:
phyloseq::sample_data(data)

### in this tutorial we will reproduce the analysis of the above mentioned paper, hence we will search for differences between males (M) and females (F), specified in the SEX column of sample_data, in the "Healthy" tissues, specified in the DISEASE column
### accordingly, let's subset the data mainteining only the healthy samples ("H") using phyloseq's "subset_sample" function
data_to_analyse<- subset_samples(data, DISEASE=="H")



######### STEP 2: USING THE AUTO_DA FUNCTION #########

CircoTax_DESeq2(data_to_analyse, contrast = c("SEX","M","F"))  # analysing the contrast M vs F

# ... that's all! ;) The whole differential analysis at each taxonomic level and results plotting is already finished! Check in your working directory.



######### STEP 3: LET's CUSTOMISE ANALYSIS AND PLOTS! #########

CircoTax_DESeq2(data_to_analyse, contrast = c("SEX","M","F"), size_taxon_circo = 2.5)  # Are the circoplot labels too big? Let's set a lower label size than default (3)


dir.create("Results_here") # To save every table and plot in a specified folder ...
CircoTax_DESeq2(data_to_analyse, contrast = c("SEX","M","F"), save_path = "Results_here/") # NB: without the "/" this character will be used as a prefix instead


CircoTax_DESeq2(data_to_analyse, contrast = c("SEX","M","F"), B = 0, lfc= 0, tax_out= c("Phylum","Class","Genus"), remove_redundants = F) # To disable the default filters


CircoTax_DESeq2(data_to_analyse, contrast = c("SEX","M","F"), auto_save = F) # make tables and plots available in your R environment rather than exporting them to files


CircoTax_DESeq2(data_to_analyse, contrast = c("SEX","M","F"), tax_out= c("Phylum","Class","Genus") ) # To excude from the analysis the specified ranks


CircoTax_DESeq2(data_to_analyse, contrast = c("SEX","M","F"), auto_log = F) # To disable text reporting every setting at the end of each analysis



### LET'S TRY A DIFFERENT RESULTS SORTING USING THE ARGUMENT "SORT_CIRCO"

dir.create("Trying_sorting_options/")
CircoTax_DESeq2(data_to_analyse, contrast = c("SEX","M","F"), sort_circo="rank", remove_redundants = F, save_path = "Trying_sorting_options/rank_") # default
CircoTax_DESeq2(data_to_analyse, contrast = c("SEX","M","F"), sort_circo="fc", remove_redundants = F, save_path = "Trying_sorting_options/fc_")
CircoTax_DESeq2(data_to_analyse, contrast = c("SEX","M","F"), sort_circo="absfc", remove_redundants = F, save_path = "Trying_sorting_options/absfc_")
CircoTax_DESeq2(data_to_analyse, contrast = c("SEX","M","F"), sort_circo="alpha", remove_redundants = F, save_path = "Trying_sorting_options/alpha_")
# check the files with "rank", "fc", "absfc" and "alpha" prefixes in your working directory: the results does not change (see box plots) but the CircoTax aesthetic does



### the same settings are available also in the sister function "Circoplot_ALDEx2"
### see the CircoTax DA logs or the manuals to explore every possible customization of plots and analysis !
