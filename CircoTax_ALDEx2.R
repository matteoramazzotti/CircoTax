CircoTax_ALDEx2<- function(input_phyloseq, # REQUIRED
                             contrast, # REQUIRED
                             design="", # it is blank purposely (built later if not specified) 
                             eff_size = 0, # Effect size threshold for ALDEx2 results 
                             MCS = 128 , # Monte Carlo Samples, required for ALDEx2 call
                             p_val = 0.05,
                             p_adjustment= "BH", # Benj Holch
                             W = 10,
                             H = 7,
                             COLOR_A = "coral",
                             COLOR_B = "chartreuse",
                             format_image = "png",
                             auto_save= TRUE, # if FALSE then the objects are just returned to the environment
                             remove_redundants= TRUE,
                             remove_results="", # to manually discard bad results
                             tax_out="",
                             plot_circo= TRUE,
                             sort_circo = "rank", # it reorders the circoplot (see its functions)
                             plot_boxplot = TRUE,
                             sqrt_y_axis = TRUE, # it rescales the y axis in the boxplot
                             save_path = "",  # then no path will be pasted before the result/plots names by default
                             auto_log = TRUE
) {
  
  # version: 08/03/2024
  
  data=input_phyloseq
  contrast= contrast
  
  library("ggplot2")
  library("ALDEx2")
  library("phyloseq")
  library("ggh4x")
  
  
  
  ################## LOADING THE CIRCOTAX FUNCTION ##################
  
  CircoTax2=function(file,title="CircoTax plot",ramp=c("orange","white","blue"),tax_col=7:11,fc_col=2,sort=c("no","rank","fc","absfc","alpha"),sort_dir="d") {
    data=file
    sort_dir=ifelse(sort_dir == "d",FALSE,TRUE)
    if(length(tax_col) == 5) {
      gplot_labels= unlist(strsplit("PCOFG",""))
    }
    if(length(tax_col) == 6) { #KPCOFG as in RDP
      gplot_labels= unlist(strsplit("DPCOFG",""))
    }
    if(length(tax_col) == 7) { #KPCOFGS as in silva
      gplot_labels= unlist(strsplit("DPCOFGS",""))
    }
    #build the taxa-related variables (tax index and label)
    #ranks are assumed to be in decreasing order from kinkdom(domain) to species
    #y represent the height (from the center of the circle) of the bars => domain=1 (min) species=7 (max)
    y=apply(data,1,function(x) length(tax_col)-sum(ifelse(is.na(x[tax_col]),1,0)))
    if(sort[1] == "no") {
      #y=apply(data,1,function(x) length(tax_col)-sum(ifelse(is.na(x[tax_col]),1,0)))
      synth=apply(data[,tax_col],1,function(x) paste0(x,collapse="-"))
      labels=unlist(lapply(strsplit(gsub("-NA.*","",perl=T,synth),"-"),function(x) rev(x)[1]))
      y=y
      fc=data[,fc_col]
    }
    if(sort[1] == "rank") {
      #y=apply(data,1,function(x) length(tax_col)-sum(ifelse(is.na(x[tax_col]),1,0)))
      o=order(y,decreasing=sort_dir)
      synth=apply(data[o,tax_col],1,function(x) paste0(x,collapse="-"))
      labels=unlist(lapply(strsplit(gsub("-NA.*","",perl=T,synth),"-"),function(x) rev(x)[1]))
      y=y[o]
      fc=data[o,fc_col]
    }
    if(sort[1] == "alphalin") {
      synth=apply(data[,tax_col],1,function(x) paste0(x,collapse="-"))
      o=order(synth,decreasing=sort_dir)
      synth=synth[o]
      labels=unlist(lapply(strsplit(gsub("-NA.*","",perl=T,synth),"-"),function(x) rev(x)[1]))
      y=y[o]
      fc=data[o,fc_col]
    }
    if(sort[1] == "alpha") {
      synth=apply(data[,tax_col],1,function(x) paste0(x,collapse="-"))
      labels=unlist(lapply(strsplit(gsub("-NA.*","",perl=T,synth),"-"),function(x) rev(x)[1]))
      o=order(labels)
      labels=labels[o]
      y=y[o]
      fc=data[o,fc_col]
    }
    if(sort[1] == "fc") {
      o=order(data[,fc_col],decreasing=sort_dir)
      synth=apply(data[,tax_col],1,function(x) paste0(x,collapse="-"))
      synth=synth[o]
      labels=unlist(lapply(strsplit(gsub("-NA.*","",perl=T,synth),"-"),function(x) rev(x)[1]))
      y=y[o]
      fc=data[o,fc_col]
    }
    if(sort[1] == "absfc") {
      #y=apply(data,1,function(x) length(tax_col)-sum(ifelse(is.na(x[tax_col]),1,0)))
      o=order(abs(data[,fc_col]),decreasing=sort_dir)
      synth=apply(data[o,tax_col],1,function(x) paste0(x,collapse="-"))
      synth=synth[order(abs(data[,fc_col]),decreasing=sort_dir)]
      labels=unlist(lapply(strsplit(gsub("-NA.*","",perl=T,synth),"-"),function(x) rev(x)[1]))
      fc=data[o,fc_col]
    }
    
    #builds the data.frame for ggplot2 
    df=data.frame("id"=1:dim(data)[1],"name"=labels,"rank"=y,"FC"=fc)
    
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
        size=3
      ) +
      geom_hline( 
        yintercept=1:6,
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
      labs(fill="log2FC")
  }
  
  
  
  ################### CHECKING THE MANDATORY FLAGS ##################
  
  
  # checking if the object in input is a phyloseq object
  if ( is.null(sample_names(data)) ){
    stop( "\n\n Looks like the object in input has no 'sample names'. \nPlease check if this object is a normal phyloseq object (containing otu table, tax table and sample data) ")
  }
  
  
  # example of a proper contrast -->   contrast=c("Condition,Healthy,Tumor")
  if(length(contrast) != 3 ){ # then if badly written or it does not exist
    stop(paste(
      "\n\n Please define a proper 'contrast' argument \n(  'Column_name_Factor','reference_Level','Level_B',   for example:  c('Condition','Healthy','Tumor')   )\n\n",
      "You have the following column(s) available in your phyloseq sample data:",paste(colnames(sample_data(data)), collapse= "  "),"\n",
      "\nEach one of the levels in the specified contrast has to be a row of the target factor column in the phyloseq sample data\n\n")
    )
  }
  
  
  
  ################ CHECKING THE DEFAULT FLAGS STATUS ################
  
  ### the following list here allows a quick manual check of the function
  
  # design = ""
  # eff_size = 0
  # MCS = 5   # they should be 128 but with 5 is faster to test this function
  # W = 10
  # H = 7
  # p_val = 0.05
  # p_adjustment = "BH"
  # COLOR_A ="coral"
  # COLOR_B = "chartreuse"
  # format_image = "png"
  # auto_save=TRUE
  # remove_redundants=TRUE
  # remove_results=""
  # tax_out=""
  # plot_circo=TRUE
  # sort_circo = "rank"
  # plot_boxplot=TRUE
  # sqrt_y_axis=TRUE
  # save_path=""
  # auto_log=TRUE
  
  Defaults<-list(
    design = "",
    eff_size = 0,
    MCS = 128, 
    W = 10,
    H = 7, 
    p_val = 0.05,
    p_adjustment = "BH",
    COLOR_A ="coral",
    COLOR_B = "chartreuse",
    format_image = "png",
    auto_save=TRUE,
    remove_redundants=TRUE,
    remove_results="",
    # tax_out=""  # NB: tax_out is used as a 'stand alone' default and warning, then is disabled from this list
    plot_circo=TRUE,
    sort_circo = "rank", # it reorders the circoplot (see its functions)
    plot_boxplot=TRUE,
    sqrt_y_axis=TRUE , # it rescales the y axis in the boxplot
    save_path="",  # then no path will be pasted before the result/plots names
    auto_log=TRUE    # to print the help at the end of the function
  )
  
  
  Warnings_ls<-list(
    design = "\n * The argument design (statistical design of the analysis) is build only on the factor specified in contrast by default\n(to manually specify a design write it in a single string, for example  design='~ Gender+Condition' )",
    eff_size = "\n * The argument eff_size (effect size of the difference) has been set as 'eff_size=0' by default",
    W  = "\n * The argument W (width of the box plot) has been set as 'W=10' (inches) by default",
    H = "\n * The argument H (height of the box plot) has been set as 'H=7' (inches) by default",
    p_val = "\n * The argument p_val (significativity threshold) has been set as 'p_val=0.05' by default",
    p_adjustment = "\n * The argument p_adjustment (method to use to adjust the p-value) has been set as 'p_adjustment=BH' by default",
    MCS = "\n * The argument MCS (Monte Carlo Sample to be used in ALDEx2) has been set as 'MCS=128' by default",
    COLOR_A = "\n * The color of the second level has been set as COLOR_A='coral' by default",
    COLOR_B ="\n * The color of the reference level has been set as COLOR_B='chartreuse' by default",
    format_image = "\n * The argument format_image (image output) has been set as 'format_image='png'' by default\n(the possible choices are 'png' and 'pdf')",
    auto_save ="\n * The argument auto_save (automatic saving of the results) has been set as 'auto_save=TRUE' by default\n(set auto_save==FALSE to get the resulting objects in your environment)",
    remove_redundants ="\n * The argument remove_redundants (automatic removal of redundant results) has been set as 'remove_redundants=TRUE' by default\n(To disable it set 'remove_redundants=FALSE' )",
    remove_results ="\n * The argument remove_results (results that have not be displayed) has been set as blank by default\n(To manually specified which results should not be shown, provide a vector to 'remove_results=' argument)",
    # tax_out=""  # NB: tax_out is used as a 'stand alone' default and warning, then is disabled from this list
    plot_circo = "\n * The argument plot_circo (plotting through CircoTax plot) has been set as 'plot_circo=TRUE' by default",
    sort_circo = "\n * The argument sort_circo (sorting of the CircoPlot elements) has been set as 'sort_circo='rank'' by default\n(The possibile options are 'no' 'rank' 'fc' 'absfc' 'alpha') ",
    plot_boxplot= "\n * The argument plot_boxplot (plotting through boxplot) has been set as 'plot_boxplot=TRUE' by default",
    sqrt_y_axis= "\n * To allow the view of lower abundances, the y axis ticks in the boxplot follows a sqrt scale\n(To disable it add 'sqrt_y_axis=FALSE' to the function call) ",
    save_path="\n * The current working directory has been used as output path\n(to specify other paths add a string to 'save_path=' argument in the function call, e.g. 'save_path=folder_name/' )",
    auto_log="\n * To disable the display of this guide set 'auto_log=FALSE' "
  )
  
  Log_text<- "\n*** Analysis concluded! *** \n" # this will be always be printed at the end regardless 
  
  
  Warnings<- "The following unspecified options have been automatically set as described below:\n"
  for(i in 1:length(Defaults)){
    if( get( names(Defaults)[i] ) == Defaults[[i]] ){
      # "get" gets the object with that name from the environment (the one already set during the call, same name as the one in the list) and then this is compared with its corresponding default from the list
      Warnings<-c(Warnings, as.character( Warnings_ls[[i]] ))
    }
  }
  
  if(auto_log==TRUE){
    Log_text <- c(Log_text, Warnings)     # then the warnings rows are added to the log as help to customize the function
  }
  
  
  #################### CHECKING THE COMMON ERRORS #########################
  
  ####### spaces or symbols in the column names of the sample data causes errors
  
  if( ! grepl("~",design) & design!="" ){ # the design is blank by default --> if not empty then it is specified
    stop("\n\n The specified statistical design does not have required syntax \n( '~ Confounding_factor_name+Target_factor_name', e.g. '~Gender+Condition' ) \n\n")
    Sys.sleep(1)
  }
  
  
  ####### spaces or symbols in the column names of the sample data causes errors
  if( any(grepl(" ",contrast), grepl("-",contrast, fixed = T), grepl("+",contrast, fixed = T))  ){
    stop("\n\n The name specified factor and/or its levels contain spaces or symbols (e.g. - or +)! \nSince this may cause issues in R, please replace them with underscores or points.\n\n")
    Sys.sleep(1)
  }
  
  
  ####### Checking eventual invalid image format
  
  if(length(format_image)>1 | length(format_image)==0 | !(format_image %in% c("png","pdf")) & auto_save==TRUE){
    stop("\n\n Please choose one format between 'png' and 'pdf' \n\n")
    Sys.sleep(1)
  }
  
  
  ####### Checking eventual invalid p adjustment method
  
  if(length(p_adjustment)>1 | length(p_adjustment)==0 | !(p_adjustment %in% c("BH","HOLM")) ){
    stop("\n\n Please choose one multiple test p-value adjustment algorithm between 'BH' (Benjamini-Hochberg) and 'holm' (Holm) \n\n")
    Sys.sleep(1)
  }
  
  
  ######## if the logical are written as character then this silent error will not be visible in the log out file --> checking it there
  
  if(! all( is.logical(plot_circo),is.logical(auto_log),is.logical(sqrt_y_axis),is.logical(auto_save),is.logical(remove_redundants)) ){
    stop("\n\n One (or more) of the optional arguments which requires a logical input (TRUE or FALSE) has received a not logical input \n(e.g. not TRUE/FALSE or its input has been written between quotes) \n\n")
    Sys.sleep(1)
  }
  
  
  ######## if the numerical inputs are written as character then this silent error will causes issues (e.g. p_value<"Fuffolo" is FALSE but it is not an error)
  
  if(! all( is.numeric(eff_size),is.numeric(p_val),is.numeric(W),is.numeric(H), is.numeric(MCS) ) ){
    stop("\n\n One (or more) of the optional arguments which requires a numeric input (e.g. 'p_val=0.05' ) has received a not numeric input \n\n")
    Sys.sleep(1)
  }
  
  
  ######## every numeric input must have a length of 1 otherwise a silent error will causes issues
  
  if(! all( length(eff_size)==1,length(p_val)==1, length(W)==1, length(H)==1, length(MCS)==1) ){
    stop("\n\n One (or more) of the optional arguments which requires a SINGLE input (e.g. 'p_val=0.05' ) has received more than an input \n\n")
    Sys.sleep(1)
  }
  
  
  ######## checking the correctness of the factor name specified in 'contrast'
  
  if(! contrast[1] %in% colnames(sample_data(data)) ){
    stop( "\nThe main factor specified in the first string of 'contrast' have not been found among the columns of the phyloseq sample data!\n The following columns have been found: ", paste(colnames(sample_data(data)), collapse=" , " ) )
    Sys.sleep(1)
  }
  
  
  ######## checking the correctness levels specified in 'contrast'
  
  if(! contrast[2] %in% unique(sample_data(data)[[contrast[1]]]) | ! contrast[3] %in% unique(sample_data(data)[[contrast[1]]])){
    stop( "\nThe levels specified in 'contrast' have not been found in the related factor column!\n The following levels present in such column are: ", paste(unique(sample_data(data)[[contrast[1]]]), collapse=" , " ) )
  }
  
  
  ####### checking the atomic type of the target column in the sample_data
  
  if( is.character(sample_data(data)[[contrast[1]]])){
    sample_data(data)[[contrast[1]]]<-as.factor(sample_data(data)[[contrast[1]]])   # it would be converted by ALDEx2 anyways BUT with a lot of annoying warnings!
    cat("\nWarning: since it was a character vector, the specified sample data column has been coerced to factor automatically just for this analysis\n(It is suggested to manually coerce it before using this function) \n\n")
    Sys.sleep(2)
  }
  
  
  
  ####### checking the taxonomic level names of the phyloseq data
  
  colnames(tax_table(data))[colnames(tax_table(data))=="Kingdom"]<-"Domain"   # this is probable but easily resolvable if encountered, otherwise this row just won't change anything
  
  if( ! all( c("Domain","Phylum","Class","Order","Family","Genus") %in% colnames(tax_table(data)) ) ) {
    stop("\n\nThe taxonomic levels in your tax table have labeled as 'Domain,Phylum,Class,Order,Family,Genus' \n ")
  }
  
  
  ################## ADJUSTING AND CHECKING THE TAXONOMIC LEVELS AS NEEDED ######################
  
  tax_in <- c("Phylum","Class","Order","Family","Genus")
  tax_in <- tax_in [! tax_in %in% tax_out ]   # tax_out is empty by default
  cat("\nThe taxonomic levels that will be tested are:",tax_in)
  if( all(tax_out=="") & auto_log==TRUE ) {  # This warning alone is printed right here because in this way is much more clear to be understood
    cat("\n(To define which taxonomic levels should be excluded from the analysis, write them as input vector of 'tax_out=') \n ")
    Sys.sleep(2)
  }
  # tax_in <- "Family" # use this for a quick debug/check of the function
  
  if( length(tax_in) < 1 ){
    stop("\n\n You need to test at least one taxonomic level! \n\n  ")
    Sys.sleep(2)
  }
  
  
  ################ DEFINING AND CHECKING THE STATISTICAL DESIGN #########################
  
  # the object "design" is defined in the function call (e.g.  design="~ Age + Condition"  which is the common statistical design syntax according to R)
  
  # the object contrast is defining the levels of the main factor which have to be compared
  # contrast = c("Condition,Healthy,RCC") # this is an example of the object "contrast", with syntax factor, level A, level B
  
  
  if(design!=""){    # then it is manually specified by the user 
    
    there_it_is<-0    # to check if the specified main factor matches between the argument design and contrast
    split<-gsub("~","",design,fixed = T)
    split<-gsub(" ","",split,fixed = T)
    split<-unlist(strsplit(split, "+", fixed = T))
    for( i in 1:length(split) ){
      if( split[i] %in% contrast){
        there_it_is<-there_it_is + 1
      }
    }
    if(as.numeric(there_it_is) <1){
      stop("\n\n The main factor specified in contrast (e.g. 'Column_name,reference_level,level_B') is not included in the specified formula\n\n")
      Sys.sleep(1)
    }
    
  } else {    # if the design is not specified during the call (then it is still an empty string, as the default)...
    design <- paste("~",contrast[1])   # ...it takes it from the contrast (the first character is ALWAYS the targeted factor)
  }
  
  
  ### moreover, if the factor column in the input obj is not coerced as factor, or if it has more than 2 levels, then an error will be obtained where searching for the column names below!
  metadat <- as(sample_data(data),"data.frame")
  if( length( unique(metadat[[contrast[[1]]]]) ) > 2 ){ # if more than two levels in this column
    metadat[[contrast[1]]] <- factor(metadat[[contrast[1]]] , levels = c(contrast[[2]],  # reference/base
                                                                         contrast[[3]],  # higher level (target),
                                                                         metadat[[contrast[1]]] [! metadat[[contrast[1]]] %in% c(contrast[[2]], contrast[[3]])]
    )
    )
  }
  if ( class(metadat[[contrast[1]]])!="factor" ){ # if it has only 2 levels it won't be coerced by the function above (if it is not already a factor by itself)
    metadat[[contrast[1]]] <- factor(metadat[[contrast[1]]] , levels = c(contrast[[2]], contrast[[3]] ) )
  }
  
  
  
  ######################## STARTING DA WITH ALDEx2 ###########################
  
  
  ####### Preparing the objects for ALDEx2 loop
  
  suppressWarnings(rm(data_pruned, data.genus_pruned))
  # Trimming under sum of 10 (see ALDEx2 tutorial) and preparing new data (other ASV may be selected after glomming)
  suppressWarnings(data_pruned<- prune_taxa(taxa_sums(data) > 10, data)) # the warnings may derive from the tree (not important here!)
  
  Table_tot<-NULL
  Res_tot<-NULL
  
  for( t in tax_in ){
    cat("\nWorking on",t,"level...\n\n")
    suppressWarnings(rm(list=c("d", "d.prop", "Taxa.d", "res","DE", "target", "r", "r_level", "index","names_target", "logical_to_subset")))
    d <- tax_glom(data_pruned, taxrank = t, NArm = F)
    d.prop<- transform_sample_counts(d, function(x) x/sum(x)*100)
    
    if(t=="Genus"){ # updating missing names (NA and uncultured) but only for genus level
      taxa_temp<-as.data.frame(tax_table(d))
      for( x in 1: length(which(taxa_temp$Genus=="uncultured")) ) {
        taxa_temp$Genus[which(taxa_temp$Genus=="uncultured")[1]]<-paste("uncultured_ f",taxa_temp[which(taxa_temp$Genus=="uncultured")[1],"Family"])}
      for( x in 1: length(which(taxa_temp=="uncultured_ f uncultured")) ) {
        taxa_temp$Genus[ which(taxa_temp$Genus=="uncultured_ f uncultured")[1] ]<-paste("uncultured_ o",taxa_temp[which(taxa_temp$Genus=="uncultured_ f uncultured")[1],"Order"])}
      for( x in 1: length(which(is.na(taxa_temp$Genus))) ) {
        taxa_temp$Genus[ which(is.na(taxa_temp$Genus))[1] ]<-paste("NA_ f",taxa_temp[which(is.na(taxa_temp$Genus))[1],"Family"])}
      for( x in 1: length(which(taxa_temp=="NA_ f NA")) ) {
        taxa_temp$Genus[ which(taxa_temp$Genus=="NA_ f NA")[1] ]<-paste("NA_ o",taxa_temp[which(taxa_temp$Genus=="NA_ f NA")[1],"Order"])}
      for( x in 1: length(which(duplicated(taxa_temp$Genus[taxa_temp$Genus=="NA_ o NA"]))) ) {
        taxa_temp$Genus[ which(taxa_temp$Genus=="NA_ o NA")[1] ]<-paste("NA_ o NA",x+1) }
      tax_table(d)<-as.matrix(taxa_temp)
      tax_table(d.prop)<-as.matrix(taxa_temp)
      rm(taxa_temp) 
    }
    
    
    ### starting the ALDEx2 analysis
    
    # using the "manual" pipeline of Aldex2 to better specify the design ...
    mm <- model.matrix( formula(design), metadat ) # name of the sample data obj always for last!
    set.seed(1)
    aldx <- aldex.clr(otu_table(d),  mm,
                      mc.samples=MCS,   #"DMC"
                      # according to ALDEx2 vignette, the number of samples in the smallest group multiplied by the number of DMC be equal at least to 1000
                      denom="all", # every feature as denominator
                      verbose=T
    )
    
    aldx2 <- aldex.glm(aldx, verbose = T)
    #### THEN, the summary of the results
    aldx3 <- aldex.glm.effect(aldx, # calculates the effect size for every binary variable in the mm
                              CI = T # confidence interval of the effect size
    ) 
    aldx_final <- data.frame(aldx2,aldx3) 
    
    
    # Automatically searching the key columns in this dataframe ...
    column_with_Pval<- colnames(aldx_final) [ grepl(contrast[[1]],colnames(aldx_final)) &
                                                grepl(contrast[[3]],colnames(aldx_final)) &
                                                grepl("pval",colnames(aldx_final)) &
                                                !  grepl("adj",colnames(aldx_final)) ]
    column_with_eff_size<- colnames(aldx_final) [ grepl(contrast[[1]],colnames(aldx_final)) &
                                                    grepl(contrast[[3]],colnames(aldx_final)) &
                                                    grepl("effect$",colnames(aldx_final)) ]   # the dollar avoid to select the low and high CI column of the eff size
    
    
    resulted_p_val<-aldx_final[[column_with_Pval]]
    p_adj<-p.adjust(resulted_p_val, method= p_adjustment )
    
    resulted_eff_size<-aldx_final[[column_with_eff_size]]
    Enough_eff_size <- resulted_eff_size > eff_size    # logical vector
    res<- aldx_final
    res$p_adj_column<- p_adj
    res<- res[ res$p_adj_column<p_val & Enough_eff_size, ]
    results<-row.names(res)
    r<-as.data.frame(res)
    Taxa.d<-as.data.frame(tax_table(d))
    Taxa.d$ASV<-row.names(Taxa.d)
    r<- cbind.data.frame(Taxa.d[row.names(r), ] , r )
    r<-r[ r[,t] !="none", ] # removing the taxa labeled "none", they are artifact from glomming, composed of different taxa with no name for that level in NCBI
    if( length(r[,1])>0 ){ # if there are still results after removing the "none"
      r_level<-r
      r_level[, "Taxon"]<- rep(t)
      r_level <- r_level[r_level[,t] != "none" , ]
      Res_tot<-rbind.data.frame(Res_tot,r_level)
      
      ### beginning to build the basis for the boxplot (Table_DE)... 
      names_target <- as.character(r[,t])
      colnames(tax_table(d.prop))[colnames(tax_table(d.prop))==t]<-"Aimed_taxa"
      logical_to_subset <- as.data.frame(tax_table(d.prop))[["Aimed_taxa"]] %in% names_target
      # NB: if the function "subset_taxa" with %in% is used when launching this script as an function, then there will be a strange bug related to %in%  --> using prune_taxa with logical subset
      target <- suppressWarnings(prune_taxa(logical_to_subset, d.prop)) # cannot use t %in% target in this function, then it's scripted in this way
      # suppressWarnings avoids message related to the phylogenetic tree (useless here!) or to the column type (e.g. 'sample' in 'sample_sample') 
      Table_DE<-psmelt(target)
      colnames(Table_DE)[colnames(Table_DE)=="Aimed_taxa"]<-t # restored the original name
      
      ### appending to a unique box plot ...
      index<- which(colnames(Table_DE)=="Domain") : which(colnames(Table_DE)==t) # from : to
      index<- index[-length(index)] # removing the last index, regarding the taxa of interest
      Table_DE[,index]<-NULL
      Table_DE$Taxa<-t
      colnames(Table_DE)[colnames(Table_DE)==t]<-"Bacteria"
      Table_tot<-rbind.data.frame(Table_tot, Table_DE)
      # closing the nested "if r > 0" (after the "none" removal)
    } else {  # closing the initial "if r > 0"
      cat("Any results for the",t,"level\n")
      Sys.sleep(1)
    }
    
  }
  
  if( length(Res_tot[[1]]) > 1 & auto_save==TRUE ){ # if there are results
    write.csv(Res_tot, file=paste0(save_path,"ALDEx2_results.csv"), row.names = F)
  }
  
  
  
  ################## PLOTTING THE RESULTS (AN UNIQUE "IF" THAT STARTS FROM THERE) #####################
  
  if( length(Res_tot[[1]]) > 1){ # if there are results
    
    
    ###### to auto-detect redundant results
    
    auto_redund<-NULL
    if(remove_redundants==TRUE) {
      
      # 1) if redundant then same abundances
      auto_Max<-tapply(round(Table_tot$Abundance,0), Table_tot$Bacteria, max )
      auto_Min<-tapply(round(Table_tot$Abundance,0), Table_tot$Bacteria, min )
      auto_Median<-tapply(round(Table_tot$Abundance,0), Table_tot$Bacteria, median )
      # reordering the names (lower ranks first --> the redundants are those at higher levels)
      ordered_names <- c( unique(Table_tot[Table_tot$Taxa=="Genus", "Bacteria" ]),
                          unique(Table_tot[Table_tot$Taxa=="Family", "Bacteria" ]),
                          unique(Table_tot[Table_tot$Taxa=="Order", "Bacteria" ]),
                          unique(Table_tot[Table_tot$Taxa=="Class", "Bacteria" ]),
                          unique(Table_tot[Table_tot$Taxa=="Phylum", "Bacteria" ])
      )
      ordered_names<-ordered_names[!is.na(ordered_names)]
      auto_Max<-auto_Max[ ordered_names  ]
      auto_Min<-auto_Min[ ordered_names  ]
      auto_Median<-auto_Median[ ordered_names ] # the table it self follows the correct taxonomic order because it is built in this way
      # same abundances
      auto_redund<-names(auto_Max)[ duplicated(auto_Min) & duplicated(auto_Max) & duplicated(auto_Median)  ]
      
      # 2) if redundant then same ASV
      auto_ASV<-Res_tot$ASV[duplicated(Res_tot$ASV)]
      auto_ASV_Names<-unique(Table_tot$Bacteria[Table_tot$OTU %in% auto_ASV])
      auto_redund<-auto_redund[auto_redund %in% auto_ASV_Names]
      
      # 3) eventual redundant from the lower taxonomic level are not redundant!
      lower_level<-Res_tot$Taxon
      lower_level<-factor(Res_tot$Taxon)
      lower_level<-lower_level[match(lower_level, c("Phylum","Class","Order","Family","Genus")) ] # this row ensures that the order is ALWAYS correct!  
      lower_level<-as.character(lower_level[length(lower_level)]) # selects the last one
      if(! is.na(lower_level) ){ # it is NA if the target level is the Genus (no Species!) --> error in the row below
        auto_redund<-auto_redund[! auto_redund %in% Res_tot[,lower_level] ]  # NB: If lower level is still a factor then its string would be seen as number by R --> wrong column
      }
    } else {
      auto_redund<-c("")   # This object is required anyways!
    }
    
    
    ###### to add redundants manually
    
    if("remove_results" != "" ){ # then a vector of bacteria is specified during the function call
      if(! remove_results %in% as.vector(tax_table(data)) & length(remove_results)>1) {
        cat("\n Warning: the specified results to remove are not in the taxonomic table! \n\n")
        Sys.sleep(2)
      }
    }
    
    auto_redund<- c(auto_redund, remove_results)
    
    
    # Table_tot2 is the true "final" dataframe that is used to plot
    logical_selection<- Table_tot[,"Bacteria"] %in% auto_redund # to remove redundant results
    Table_tot2<-Table_tot[ ! logical_selection, ]
    logical_selection<- Table_tot2[ ,contrast[1]]  %in% contrast[-1] # if there are three (or four) groups in the factor then the third (which is excluded from the analysis) would be also plotted but gray and unlabeled!
    Table_tot2<-Table_tot2[ logical_selection, ]
    
    
    # cleaning the names
    Table_tot2$Bacteria<-gsub("_group","",Table_tot2$Bacteria) # e.g. from SILVA database
    Table_tot2$Bacteria<-gsub("[","",Table_tot2$Bacteria, fixed=T)
    Table_tot2$Bacteria<-gsub("]","",Table_tot2$Bacteria, fixed=T)
    Table_tot2$Bacteria<-gsub("_"," ",Table_tot2$Bacteria, fixed=T)
    
    # choosing colors is needed a named vector (creating it outside the plot is required to avoid errors)
    colors<- ifelse(as.character(unique(Table_tot2[,contrast[1]]))==contrast[2], COLOR_B, COLOR_A)
    names(colors)<- ifelse(unique(as.character(Table_tot2[,contrast[1]]))==contrast[2], contrast[2], contrast[3])
    
    # preparing the object with the plot
    DE_plot<- ggplot(Table_tot2, aes(x= Bacteria, y=Abundance,
                                     fill= get(contrast[1])
    )
    ) +   # NB: "get" is required to consider *the string from* variable contrast[1] as an variable name
      facet_grid(~factor(Taxa,levels = c("Phylum","Class","Order","Family","Genus")),
                 scales = "free_x", space="free") +
      geom_boxplot(width=0.8) + theme_bw(base_size = 15) +
      scale_fill_manual(values= colors ) +
      theme(strip.text.x=element_text(size=14,colour="black")) + 
      guides( fill=guide_legend(nrow=1) ) +
      theme(legend.margin=margin(-25, 0, 0, 0), legend.position="bottom", 
            legend.key.size=unit(0.8,"cm"), legend.text=element_text(size=18),
            axis.text.x = element_text(angle = 30, vjust=1, hjust=1, size=11.8), 
            axis.text.y = element_text(size=10.2), 
            plot.title= element_text(size=16),
            panel.grid.minor.y= element_blank(),
            plot.margin =  margin(t=5,r=5,b=5, l=5) ) +  
      scale_x_discrete(expand=c(-0.2, 1)) +
      labs(title= "Differently abundant Taxa", y="Percent Abundance", 
           fill=contrast[1],
           x="")
    
    
    if(sqrt_y_axis==TRUE){
      if(max(Table_tot2$Abundance)<2){
        DE_plot + scale_y_sqrt(breaks=c(0.1, 0.2, 0.5, 1, 1.5, 2))
      }
      if(max(Table_tot2$Abundance)> 2 & max(Table_tot2$Abundance)<4){
        DE_plot + scale_y_sqrt(breaks=c(0.1, 0.5, 1,2,2.5,3,3.5,4))
      }
      if(max(Table_tot2$Abundance)> 4 & max(Table_tot2$Abundance)<10){
        DE_plot + scale_y_sqrt(breaks=c(0.1, 0.5, 1,2,2.5,3,seq(4,max(Table_tot2$Abundance+2),2)))
      }
      if(max(Table_tot2$Abundance)> 10 & max(Table_tot2$Abundance)<28){
        DE_plot + scale_y_sqrt(breaks=c(0.1, 0.5, 1,seq(2,max(Table_tot2$Abundance+2),2)))
      }
      if( max(Table_tot2$Abundance)> 28 & max(Table_tot2$Abundance) < 51 ){
        DE_plot + scale_y_sqrt(breaks=c(0.2, 0.5, seq(2,24,2), seq(27,max(Table_tot2$Abundance+2),3)))
      }
      if( max(Table_tot2$Abundance)> 51 & max(Table_tot2$Abundance) < 101 ){
        DE_plot + scale_y_sqrt(breaks=c(0.2, 0.5, seq(2,24,2), seq(27,51,3), seq(55,max(Table_tot2$Abundance+2),5)))
      }
    } else {
      DE_plot   # without sqrt axis
    }
    
    
    if(format_image=="png" & plot_boxplot==TRUE & auto_save==TRUE){
      ggsave(filename = paste0(save_path,"ALDEx2_Boxplot.png"), device=format_image, width = W, height = H, units = "in", dpi=300)
    }
    if(format_image=="pdf" & plot_boxplot==TRUE & auto_save==TRUE){
      ggsave(filename = paste0(save_path,"ALDEx2_Boxplot.pdf"), device=format_image, width = W, height = H)
    }
    
    
    
    ################# CIRCOPLOT ALONE #####################
    
    ### preparing the objects
    Res_tot2<-Res_tot[! grepl("uncultured", Res_tot$Genus), ] # those uncultured does not have a name and here they are almost absent... better avoiding to plot those noises!
    ### removing the redundant results
    Res_tot2<-Res_tot2[!(Res_tot2$Family %in% auto_redund & Res_tot2$Taxon %in% "Family"), ] # NB: this structure is needed to avoid discarding a genus with the same name of its family (for example)
    Res_tot2<-Res_tot2[!(Res_tot2$Order %in% auto_redund & Res_tot2$Taxon %in% "Order"), ]
    Res_tot2<-Res_tot2[!(Res_tot2$Class %in% auto_redund & Res_tot2$Taxon %in% "Class"), ]
    Res_tot2<-Res_tot2[!(Res_tot2$Phylum %in% auto_redund & Res_tot2$Taxon %in% "Phylum"), ]
    Res_tot2$Genus<-gsub("-","_",Res_tot2$Genus, fixed=T) # to avoid errors
    # cleaning the names for the Circoplot
    for(t in tax_in){
      Res_tot2[,t]<-gsub("_group","",Res_tot2[,t])
      Res_tot2[,t]<-gsub("[","",Res_tot2[,t], fixed=T)
      Res_tot2[,t]<-gsub("]","",Res_tot2[,t], fixed=T)
      #Res_tot2[,t]<-gsub("_"," ",Res_tot2[,t], fixed=T)
    }
    
    if(length(Res_tot2$ASV)< 7 & plot_circo==TRUE){
      Log_text<-c(Log_text,"\n\nNB: Due to the low number of results, the CircoPlot may not match the expected aesthetic\n")
    }
    
    
    if(plot_circo==TRUE){   # preparing the object with the plot in it
      CircoTax2(Res_tot2,title="",ramp=c("orange","white","blue"),
                tax_col=2:7,
                fc_col=  which( colnames(Res_tot2)==column_with_eff_size ) ,
                sort=sort_circo) +
        labs(title=paste(contrast[3],"vs",contrast[2],"\n DA taxa\n"), fill = "Eff size")
    }
    
    
    # NB: IT IS NOT POSSIBLE TO USE THE PNG OR PDF FUNCTION INSIDE AN *UNIQUE* IF, THEM WOULD NOT WORK!
    
    if(format_image=="pdf" & auto_save==TRUE & plot_circo==TRUE){
      ggsave(filename = paste0(save_path,"CircoTax_plot.pdf"), device=format_image)
    }
    if(format_image=="png" & auto_save==TRUE & plot_circo==TRUE){
      ggsave(filename = paste0(save_path,"CircoTax_plot.png"), device=format_image, width = 6.2, height = 6.2, units = "in", dpi=300)
    }
  }
  
  
  ################# IF AUTO_SAVE OFF #####################
  
  if(length(Res_tot$ASV) > 1 & auto_save==FALSE){
    Results_ALDEx2<<-Res_tot2
    Boxplot_ALDEx2<<-Table_tot2
    Circoplot_ALDEx2 <<- CircoTax2(Res_tot2,title="",ramp=c("orange","white","blue"),
                                   tax_col=2:7,
                                   fc_col=  which( colnames(Res_tot2)==column_with_eff_size ) ,
                                   sort=sort_circo) +
      labs(title=paste(contrast[3],"vs",contrast[2],"\n DA taxa\n"), fill = "Eff size")
    Sys.sleep(1)
    cat("\nDue to auto_save==FALSE, the following objects have been saved in your environment:\nResults_ALDEx2 Boxplot_ALDEx2 Circoplot_ALDEx2\n ")
    Sys.sleep(1)
  }
  
  
  
  ###################### IF THERE ARE NO RESULTS ###################################
  
  if( length(Res_tot$ASV) < 1 ){ # (NB: not using "else", the condition is not quite the same)
    cat("\n\n Looks like there are not significant differences there! \n\n")
    Sys.sleep(2)
  }
  
  
  ############################# CLOSING #########################################
  
  cat(Log_text)
  
  #### if there are NAs at higher levels then every NA will be plotted as one --> being this an anomalous situation, a BIG warning is printed (it's no use to correct something which is wrong itself!) 
  if(auto_save==TRUE & length(which(is.na(Res_tot$Phylum)))> 0 ){   
    cat("\n\nSome of the results is unassigned (NA) at the phylum level! \nThis is *usually* an anomalous situation (e.g. host contaminat) and then may cause graphical issues.\nIt is strongly suggested to check the related OTU/ASV .")
    Sys.sleep(2)
  }
  
  cat("\n\n Take care, have a good work buddy! \n\n")
  
  Settings<-cbind.data.frame(statistical_design=design, # with cbind it's possible to specify the col names directly
                             Effect_size=eff_size,
                             p_value_threshold=p_val,
                             p_value_adj_method=p_adjustment,
                             MonteCarloSamples=MCS,
                             redundant_results_automatically_removed=remove_redundants,
                             manually_removed_results=remove_results,
                             taxa_analysed=paste(tax_in, collapse = " "), 
                             sqrt_y_axis=sqrt_y_axis,
                             output_format=format_image)
  
  if(auto_save==TRUE){
    write.table(x = t(Settings), 
                file=paste0(save_path,"ALDEx2_Circoplot_Settings.csv"),
                row.names = T, col.names = F, 
                sep=";", quote=F)
  }
}

