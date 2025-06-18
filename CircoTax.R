CircoTax=function(input_table,
                  title="CircoTax plot",
                  fill_text="fold change",
                  ramp=c("orange","white","blue"),
                  tax_col=2:length(colnames(input_table)),
                  fc_col=1,
                  size_taxon_circo=3,
                  sort=c("no","rank","fc","absfc","alpha"),
                  sort_dir="d" 
                  ) {
  
  # version: 16/06/2025  (previous version: 17/12/2024)
  
  
  library("ggplot2")
  library("ggh4x")
  
  sort_dir=ifelse(sort_dir == "d",FALSE,TRUE)
  
  if( "Domain" %in% colnames(input_table) | "Kingdom" %in% colnames(input_table) ){
    file <- gsub("Kingdom","Domain",colnames(input_table))
    if(length(tax_col) == 5) {
      gplot_labels= unlist(strsplit("DPCOF",""))
    }
    if(length(tax_col) == 6) { #KPCOFG as in RDP
      gplot_labels= unlist(strsplit("DPCOFG",""))
    }
    if(length(tax_col) == 7) { #KPCOFGS as in silva
      gplot_labels= unlist(strsplit("DPCOFGS",""))
    }
  } else {
    if(length(tax_col) == 5) {
      gplot_labels= unlist(strsplit("PCOFG",""))
    }
    if(length(tax_col) == 6) { #PCOFGS
      gplot_labels= unlist(strsplit("PCOFGS",""))
    }
    # gplot_labels <<- gplot_labels
  }
    
  #build the taxa-related variables (tax index and label)
  #ranks are assumed to be in decreasing order from kinkdom(domain) to species
  #y represent the height (from the center of the circle) of the bars => domain=1 (min) species=7 (max)
  y=apply(input_table,1,function(x) length(tax_col)-sum(ifelse(is.na(x[tax_col]),1,0)))
  if(sort[1] == "no") {
    #y=apply(input_table,1,function(x) length(tax_col)-sum(ifelse(is.na(x[tax_col]),1,0)))
    synth=apply(input_table[,tax_col],1,function(x) paste0(x,collapse="-"))
    labels=unlist(lapply(strsplit(gsub("-NA.*","",perl=T,synth),"-"),function(x) rev(x)[1]))
    y=y
    fc=input_table[,fc_col]
  }
  if(sort[1] == "rank") {
    #y=apply(input_table,1,function(x) length(tax_col)-sum(ifelse(is.na(x[tax_col]),1,0)))
    o=order(y,decreasing=sort_dir)
    synth=apply(input_table[o,tax_col],1,function(x) paste0(x,collapse="-"))
    labels=unlist(lapply(strsplit(gsub("-NA.*","",perl=T,synth),"-"),function(x) rev(x)[1]))
    y=y[o]
    fc=input_table[o,fc_col]
  }
  if(sort[1] == "alphalin") {
    synth=apply(input_table[,tax_col],1,function(x) paste0(x,collapse="-"))
    o=order(synth,decreasing=sort_dir)
    synth=synth[o]
    labels=unlist(lapply(strsplit(gsub("-NA.*","",perl=T,synth),"-"),function(x) rev(x)[1]))
    y=y[o]
    fc=input_table[o,fc_col]
  }
  if(sort[1] == "alpha") {
    synth=apply(input_table[,tax_col],1,function(x) paste0(x,collapse="-"))
    labels=unlist(lapply(strsplit(gsub("-NA.*","",perl=T,synth),"-"),function(x) rev(x)[1]))
    o=order(labels)
    labels=labels[o]
    y=y[o]
    fc=input_table[o,fc_col]
  }
  if(sort[1] == "fc") {
    o=order(input_table[,fc_col],decreasing=sort_dir)
    synth=apply(input_table[,tax_col],1,function(x) paste0(x,collapse="-"))
    synth=synth[o]
    labels=unlist(lapply(strsplit(gsub("-NA.*","",perl=T,synth),"-"),function(x) rev(x)[1]))
    y=y[o]
    fc=input_table[o,fc_col]
  }
  if(sort[1] == "absfc") {
    o=order(abs(input_table[,fc_col]),decreasing=sort_dir)
    synth=apply(input_table[o,tax_col],1,function(x) paste0(x,collapse="-"))
    labels=unlist(lapply(strsplit(gsub("-NA.*","",perl=T,synth),"-"),function(x) rev(x)[1]))
    y=y[o]
    fc=input_table[o,fc_col]
  }
  
  #builds the data.frame for ggplot2 
  df=data.frame("id"=1:dim(input_table)[1],"name"=labels,"rank"=y,"FC"=fc)
  
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
      "text",label=gplot_labels, x=rep(0,length(tax_col)), y=1:length(tax_col),size=4.5, vjust=0.58
    ) +
    # geom_text( 
    geom_text_aimed(  # it allows to calculate the ideal angle to have labels perfectly perpendicular to the axis
      aes(
        x= id,
        y=7,
        label=name),
      color="black",
      fontface="bold",
      alpha=0.6,
      size= size_taxon_circo,
    ) +
    geom_hline( 
      yintercept=1:length(tax_col),
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
