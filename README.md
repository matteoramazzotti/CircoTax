# CircoTax
CircoTax is an R function implementing a specialized *ggplot2* graph to represent in a radial form a rank-aware collection of differentially abundant taxa. Specifically, the CircoTax radial bar plot shows 6 or 7 sectors that encode the taxonomy depth (from kingdom to genus) and, departing from the center, a number of radial bars that reach the appropriate sector and whose color and transparency are proportional to the log fold change intensity and direction.
A basic Quickstart guide is provided here, please refer to the Documentation for additional information and tutorials.

<p align="center">
	<img src="docs/img/CircoTax_plot_default.png" width="49%"  title="CircoTax Plot obtained from a complete taxonomy matrix.">
	<img src="docs/img/CircoRax_plot_custom.png" width="49%"  title="CircoTax Plot obtained from a custom R data.frame composed of three columns, namely taxon name, taxonomic rank name and value to display (e.g. FC).">
</p>

## Quickstart

Specific flavors of CircoTax has been developed to facilitate its usage:

- ***CircoTax()*** is the basic function to plot a complete *phyloseq* taxonomy matrix or even custom data frames.

- ***CircoTaxAutoDA()*** implement a differential abundance analysis using either the *DESeq2* or *ALDEx2* R packages. 
  

### How to install

To install with *devtools* :

```
	library("devtools")
	install_github("matteoramazzotti/CircoTax")
```
> [!NOTE]  
> Required dependencies should be installed automatically. In case of errors, please refer to the [Wiki](https://github.com/matteoramazzotti/CircoTax/wiki/Usage#installation).


### How to plot a CircoTax using a complete taxonomy matrix

```
	library("circotax")
	file_format_1 <- system.file("extdata", "circotax_example_1.txt", package = "circotax")
	data = read.delim(file_format_1, sep = "\t")
	CircoTax(data, tax_col = 2:7, fc_col = 1)
```

<p align="center">
	<img src="docs/img/CircoTax_plot_default.png" width="50%"  title="CircoTax Plot obtained from a complete taxonomy matrix.">
</p>


In this example the CircoTax plot is representing 62 differentially abundant ranks, among which 27 are at genus level resolution, 18 at family level resolution, 10 at order level resolution, 3 at the class level resolution and 4 at phylum level resolution.

CircoTax can take a complete taxonomy matrix as an input, here we provide a subset from our tutorial dataset [*circotax_example_1*](https://raw.githubusercontent.com/matteoramazzotti/CircoTax/refs/heads/main/inst/extdata/circotax_example_1.txt):


| **log2FoldChange** | **Domain** |         **Phylum**        |      **Class**      |     **Order**     |     **Family**     |    **Genus**    |
|:------------------:|:----------:|:-------------------------:|:-------------------:|:-----------------:|:------------------:|:---------------:|
|       -0.741       |  Bacteria  |       Proteobacteria      |          NA         |         NA        |         NA         |        NA       |
|        5.294       |  Bacteria  |        Chloroflexi        |          NA         |         NA        |         NA         |        NA       |
|       -1.841       |  Bacteria  |        Tenericutes        |          NA         |         NA        |         NA         |        NA       |
|       -1.984       |  Bacteria  | Cyanobacteria/Chloroplast |          NA         |         NA        |         NA         |        NA       |
|       -1.221       |  Bacteria  |       Proteobacteria      | Gammaproteobacteria |         NA        |         NA         |        NA       |
|        2.146       |  Bacteria  |       Bacteroidetes       |    Flavobacteriia   |         NA        |         NA         |        NA       |
|       -1.726       |  Bacteria  |         Firmicutes        |    Negativicutes    |         NA        |         NA         |        NA       |
|        2.370       |  Bacteria  |       Proteobacteria      | Alphaproteobacteria |  Caulobacterales  |  Caulobacteraceae  |  Brevundimonas  |
|       -3.611       |  Bacteria  |       Proteobacteria      | Gammaproteobacteria | Enterobacteriales | Enterobacteriaceae |     Pantoea     |
|       -3.204       |  Bacteria  |       Proteobacteria      | Gammaproteobacteria |   Pasteurellales  |   Pasteurellaceae  |   Haemophilus   |
|       -3.242       |  Bacteria  |       Proteobacteria      | Gammaproteobacteria |   Pasteurellales  |   Pasteurellaceae  | Aggregatibacter |


### How to plot a CircoTax using a custom data frame

The *CircoTax()* function allows to obtain a plot from an R data.frame composed of three columns, namely taxon name, taxonomic rank name and value to display (e.g. FC).

```
	library("circotax")
	file_format_2 <- system.file("extdata", "circotax_example_2.txt", package = "circotax")
	data = read.delim(file_format_2, sep = "\t")
	CircoTax(data, names = 1, tax_col = 2, fc_col = 3)
```

<p align="center">
	<img src="docs/img/CircoRax_plot_custom.png" width="49%"  title="CircoTax Plot obtained from a custom R data.frame composed of three columns, namely taxon name, taxonomic rank name and value to display (e.g. FC).">
</p>

An input example is included in the file [*circotax_example_2*](https://raw.githubusercontent.com/matteoramazzotti/CircoTax/refs/heads/main/inst/extdata/circotax_example_2.txt), which was used for the above CircoTax plot:

|     **Tax_name**    | **Rank_name** |  **FC**  |
|:-------------------:|:-------------:|:--------:|
| Gammaproteobacteria |     Class     |  1.06736 |
|    Tannerellaceae   |     Family    | 2.398398 |
|    Marinifilaceae   |     Family    | 3.136581 |
|    Sutterellaceae   |     Family    | 3.808097 |
|   Lactobacillaceae  |     Family    | 3.624271 |
|    Capnocytophaga   |     Genus     | -5.55057 |
|    Alloprevotella   |     Genus     | 3.172175 |
|    Muribaculaceae   |     Genus     |    10    |
|    Butyricimonas    |     Genus     |  4.38531 |
|      Clostridia     |     Genus     | 4.235196 |



> [!NOTE]  
> The function *CircoTax()* can also take optional arguments. More information regarding this can be found in the [Wiki](https://github.com/matteoramazzotti/CircoTax/wiki/Usage#circoTax).


### How to plot a CircoTax using Auto DA functions

***CircoTaxAutoDA()*** perform the differential analysis at each taxonomic level between two groups using either *DESeq2* or *ALDEx2*. It then filters the significant results and finally generate a tsv table, a box plot of percent abundances and a CircoTax plot of the results.

```
	library("circotax")
	library("phyloseq")

	data <- subset_samples(circotax::circotax_phyloseq_example, DISEASE == "H")
	df = CircoTaxAutoDA(data, mode = "DESeq2", contrast = c("SEX","M","F"))
	CircoTax(df) 
```

The above code is performed taking a *phyloseq* object named *circotax_phyloseq_example* which is bundled with our library and obtained by [*Elena Niccolai et al.*](https://doi.org/10.1186/s13293-023-00523-w) and results in the following CircoTax plot:

<p align="center">
	<img src="docs/img/CircoTax_plot_DESeq2.png" width="49%"  title='CircoTax Plot obtained using CircoTaxAutoDA(data, mode = "DESeq2", contrast = c("SEX","M","F"))'>
	<!-- <img src="docs/img/DESeq2_Boxplot.png" width="49%" title="CircoTax BoxPlot obtained using CircoTax_DESeq2()"> -->
</p>

> [!NOTE]  
> Functions *CircoTaxAutoDA()* accepts additional arguments to customise the analysis. More information regarding this can be found in the [Wiki](https://github.com/matteoramazzotti/CircoTax/wiki/Usage#circotaxautoda).