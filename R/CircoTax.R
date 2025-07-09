#' CircoTax Plot
#' 
#' CircoTax is an R function implementing a specialized ggplot2 graph to represent in a radial form 
#' a rank-aware collection of differentially abundant taxa.
#'
#' @param input_table Input DataFrame
#' @param title Title display on the plot. Default: CircoTax plot
#' @param fill_text Text to display above the color legend. Default: log2FC
#' @param ramp Character vector to represent a color gradient for the log2FC scale of values, from min to max. Default: Defaults to c("blue","white","orange")
#' @param tax_col Column index or range of column indexes holding IDs for the ranks. Default: 2:length(colnames(input_table))
#' @param fc_col Column index holding log2FC values. Default: 1 when tax_col argument is a range of indexes else 3 
#' @param names Column index holding IDs for the taxon, used when tax_col argument is a single index. Default: 1
#' @param size_taxon_circo Size of plot annotation. Default: 4
#' @param sort Sorting logic: "no" (no sorting), "rank" (by taxonomic rank), "fc" (by fold change), "absfc" (by absolute fold change), "alpha" (alphabetic order), "alpharank" (). Default: "rank".
#' @param rank_vector Character vector holding the ranking system ordered from Highest rank to lowest. Default: c("Domain","Phylum","Class","Order","Family","Genus","Species")
#' @param decrease_order Boolean to invert order in sorting. Default: FALSE
#' @param collapse_ranks Boolean to use only Ranks represented in the table. Default: TRUE
#' @param save_path Save path for plot image. Default: NULL
#' @param image_dim Pixel value for height and width of the plot.  Default: 4000
#' @param save_format "png" or "pdf". Default: "png"
#' @return CircoTax plot
#' @export
CircoTax = function (
  input_table,
  title = "CircoTax plot",
  fill_text = "log2FC",
  ramp = c("blue","white","orange"),
  tax_col = 2:length(colnames(input_table)),
  fc_col = ifelse(length(tax_col)>1,1,3),
  names = 1,
  size_taxon_circo=4,
  sort = c("rank","no","fc","absfc","alpha","alpharank"),
	rank_vector = c("Domain","Phylum","Class","Order","Family","Genus","Species"),
  decrease_order = FALSE,
	collapse_ranks = FALSE,
	save_path = NULL,
	image_dim = 4000,
	save_format = "png"
) {
    
  # library("ggplot2")
  # library("ggh4x")
  # library("stringi")
	# library("here")


	ow <- options("warn")
	options(warn = 0)
	if (!length(ramp) == 3) {
		warning("Default is 3 colors")
		ramp = c("blue","white","orange")
	}
	
	input_table = CircoTaxTableFormatter(
		input_table,
		names_col = names,
		tax_col = tax_col,
		fc_col = fc_col,
		rank_vector = rank_vector,
		collapse_ranks = collapse_ranks
	)
	tax_col = 2:length(colnames(input_table))
	fc_col = 1
	
	gplot_rank_labels <- unlist(
		lapply(
			colnames(input_table[,tax_col]),
			function(x) unlist(strsplit(x,""))[1]
		)
	)

	
  #build the taxa-related variables (tax index and label)
  #ranks are assumed to be in decreasing order from kinkdom(domain) to species
  #y represent the height (from the center of the circle) of the bars => domain=1 (min) species=7 (max)

	settings <- CircoTaxSorter(
		input_table = input_table,
		sort = sort[1],
		tax_col = tax_col,
		fc_col = fc_col,
		decrease_order = decrease_order
	)

  #builds the data.frame for ggplot2 

	df = data.frame(
		id = 1:dim(input_table)[1],
		name = settings$labels,
		rank = settings$y,
		FC = settings$fc
	)

  # plot settings
  circoTitle = ggplot2::ggtitle(title) # plot Title
  circoBars = ggplot2::geom_col(  # bars logic
    ggplot2::aes(fill=FC),
    position = "dodge"
  )
  circoGradient = ggplot2::scale_fill_gradient2(
    low = ramp[1],
    mid = ramp[2],
    high = ramp[3]
  )

  circoAnnotation = ggplot2::annotate(
    "text",
    label = gplot_rank_labels,
    x = rep(0,length(tax_col)),
    y = 1:length(tax_col),
    size = 4.5,
    vjust = 0.58
  )

  circoAnnotationTrim = ggh4x::geom_text_aimed(  # it allows to calculate the ideal angle to have labels perfectly perpendicular to the axis
    ggplot2::aes(
      x = id,
      y = length(tax_col) + 1,
      label = df[,"name"]
    ),
    color = "black",
    fontface = "bold",
    alpha = 0.6,
    size = size_taxon_circo,
  )

  circoLines = ggplot2::geom_hline( 
    yintercept= 1:length(tax_col),
    color="grey"
  ) 
  circoClip = ggplot2::coord_polar(start = 0, clip="off")
  circoLegendTitle = ggplot2::labs(fill=fill_text)

  #the plot starts here
	plot(
		ggplot2::ggplot(
			df,
			ggplot2::aes(x = id, y = rank)
		) +
		circoTitle +
		circoBars + 
		circoGradient +
		circoClip +
		circoLines +
		circoAnnotation +
		circoAnnotationTrim +
		circoLegendTitle + 
		ggplot2::theme(
			axis.text = ggplot2::element_blank(),
			axis.ticks = ggplot2::element_blank(),
			axis.title = ggplot2::element_blank(),
			panel.background = ggplot2::element_rect(fill = NA),
			plot.margin = grid::unit(rep(1,4), "cm"),      # It adjusts the margins of the plot to avoid the trimming of the labels!
			plot.title = ggplot2::element_text(hjust = 0.5,face="bold", size=18)
		)
	)
	if (!is.null(save_path)) {
		scale = 4000 /image_dim
		ggplot2::ggsave(filename = paste0(save_path,"CircoTax_plot.",save_format), device=save_format, width = image_dim, height = image_dim, units = "px", dpi = 300,scale = scale)	
	}
}