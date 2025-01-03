# CircoTax
CircoTax is an R function for that exploits the radial capabilities of ggplot2 for the representation of differential abundance analyses on microbiota data.

The file CircoTax_tutorial.pdf provides a guide on how to use the CircoTax package.
The CircoTax_tutorial_AutoDA.R contains an explained working script using CircoTax.

CircoTax can load data from an R data.frame (see example_data.txt) or can receive input directly from two R functions
- CircoTax_DESeq2.R
- CircoTax_ALDEx2.R
that implement a differential abundance analysis (with two differnt R packages, namely DESeq2 and ALDEx2) from the data.RData file and automatically create the CircoTax plot.
