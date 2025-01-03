# CircoTax
CircoTax is an R function for that exploits the radial capabilities of ggplot2 for the representation of differential abundance analyses on microbiota data.

The file CircoTax_tutorial.pdf provides a guide on how to use the CircoTax package.

Specific flavors of CircoTax has been developed to facilitate its usage:

- CircoTax.R is the basic function that can manages data files such exmaple_complete_data.txt

- CircoTax_custom.R is a simpler version that manages data files such exmaple_custom_data.txt

- CircoTax_DESeq2.R and CircoTax_ALDEx2.R load data from files such as example_data.RData
  
The latter two functions implement a differential abundance analysis (with two differnt R packages, namely DESeq2 and ALDEx2) and automatically create the CircoTax plot.

The file CircoTax_tutorial_AutoDA.R contains an explained working script on how to set up the diferential abundance analyses implemneted in CircoTax.
