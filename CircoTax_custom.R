CircoTax_custom=function(input_table,
                          title="CircoTax plot",
                          fill_text="fold change",
                          ramp=c("orange","white","blue"),
                          names=1,
                          tax_col=2,
                          size_taxon_circo=3,
                          fc_col=3) {
  
  # version: 16/06/2025  (previous version: 17/12/2024)
  

  library("ggplot2")
  library("ggh4x")
  
  circular_names<- factor( input_table[,names] , unique(input_table[,names]) )
  ranks_names<- unique(input_table[,tax_col])
  ranks_ordered <- factor( input_table[,tax_col], ranks_names )
  ranks_numeric <- as.integer(ranks_ordered)
  
  number_point <- gsub(",",".",input_table[, fc_col])
  fc= as.numeric(number_point)
  
  #builds the data.frame for ggplot2 
  df=data.frame("id"=circular_names, "y_names"=ranks_ordered,"rank"= ranks_numeric, "FC"=fc)
  
  #the plot starts here
  ggplot(df, aes(x = id, y = rank)) +
    ggtitle(title) +
    geom_col(
      aes(fill=FC),
      position = "dodge"
    ) + 
    scale_fill_gradient2(
      low="blue",
      mid="white",
      high="orange"
    ) +
    annotate(
      "text",label=ranks_ordered, x=rep(0,length(ranks_numeric)), y=ranks_numeric, size=4.5, vjust=0.58
    ) +
    # geom_text( 
    geom_text_aimed(  # it allows to calculate the ideal angle to have labels perfectly perpendicular to the axis
      aes(
        x= id,
        y= length(ranks_names) + 1, # height of the name
        label=id),
      color="black",
      fontface="bold",
      alpha=0.6,
      size= size_taxon_circo ,
    ) +
    geom_hline( 
      yintercept=1:length(ranks_names),
      color="grey"
    ) +
    coord_polar(start = 0, clip="off") +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.background = element_rect(fill = NA),
      plot.margin = unit(rep(1,4), "cm"),      # It adjusts the margins of the plot to avoid the trimming of the labels!
      plot.title = element_text(hjust = 0.5,face="bold", size=18)
    ) +
    labs(fill=fill_text)
}
