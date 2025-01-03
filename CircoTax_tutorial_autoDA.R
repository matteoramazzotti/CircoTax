##### PREPARING THE ENVIRONMENT ####

source("CircoTax_DESeq2.R")  # loads the function "CircoTax_DESeq2", or...
#source("Functions/CircoTax_ALDEx2.R") # loads the function "CircoTax_ALDEx2"
load("data.RData")   # loads the R data containing the example phyloseq object "data"

#look at possible contrasts:
sample_data(data)
#in this case we search for differences between males (M) and females (F) stored in the SEX column of sample_data



##### CLASSICAL USE ######

CircoTax_DESeq2(data, contrast = c("SEX","M","F"))  # analysing the contrast M vs F, then plot the CircoTax and the abundances related boxplot

# NB: the function will load in R the required packages (see documentation), but won't install them if they are not already on user's machine



##### EXAMPLE OF ADVANCED OPTIONS ##### 

CircoTax_DESeq2(data, contrast = c("SEX","M","F"), size_taxon_circo = 2.5)  # Are the circoplot labels too big? Let's set a lower label size than default (3)


dir.create("Results_here") # To save every table and plot in a specified folder ...
CircoTax_DESeq2(data, contrast = c("SEX","M","F"), save_path = "Results_here/") # NB: without the "/" this character will be used as a prefix instead


CircoTax_DESeq2(data, contrast = c("SEX","M","F"), B = 0, lfc= 0, tax_out= c("Phylum","Class","Genus"), remove_redundants = F) # To disable the default filters


CircoTax_DESeq2(data, contrast = c("SEX","M","F"), auto_save = F) # make tables and plots available in your R environment rather than exporting them to files


CircoTax_DESeq2(data, contrast = c("SEX","M","F"), tax_out= c("Phylum","Class","Genus") ) # To excude from the analysis the specified ranks


CircoTax_DESeq2(data, contrast = c("SEX","M","F"), auto_log = F) # To disable text reporting every setting at the end of each analysis


# the same settings are available also in the sister function "Circoplot_ALDEx2"
# see the CircoTax DA logs or the manuals to explore every possible customization of plots and analysis
