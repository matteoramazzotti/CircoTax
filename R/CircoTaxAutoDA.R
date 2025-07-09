#' CircoTax AutoDA
#' 
#' Starting from a phyloseq object, performs the differential analysis at each taxonomic level between two groups using either DESeq2 or ALDEx2. 
#' They then filter the significant results and generate a tsv table, ready to be displayed with CircoTax().
#'
#' @param data phyloseq input object
#' @param contrast Vector of three characters which are (in this order) the names of the factor of interest and two of its levels of which report the differences. e.g. c("SEX","M","F")
#' @param design String describing statistical design (within the limits of DESeq2 or ALDEx2) on which the analysis is conducted. e.g. '~ Gender+Condition'. Default: NULL
#' @param tax_out Vector whose characters are the taxonomic level to remove. Default: c()
#' @param mode Differential Analysis method to be used, either "DESeq2" or "ALDEx2". Default: "DESeq2"
#' @param MCS Monte-Carlo samples to generate during ALDEx2 computations. Default: 128.
#' @param p_adjustment Multiple test adjustment by adjusting (penalize) each p-value. The possible inputs are "BH" (Benjamini-Hochberg) and "holm" (Holm). Default: "BH".
#' @param eff_size ALDEx2 effect size threshold under which exclude results. Default: 0.
#' @param p_val p-value threshold (after the multiple test adjustment) to consider a result as significant. Default: 0.05
#' @param B DESeq2 base-mean threshold (mean after DESeq2 values transformation) under which exclude results. Default: 50.
#' @param lfc Log fold change threshold. Default: 1
#' @param lfc_shrink_method Shrinkage method of the log fold change, either "none", "apeglm" and "ashr". Default: "none".
#' @param remove_redundants Boolean to enable the automatic removal of redundant results (e.g. result repeated also in higher taxonomic levels being the only observation in that taxonomic clade). Default: TRUE.
#' @param remove_results Vector whose characters are the bacteria name to remove. Default: c()
#' @return Dataframe
#' @export
CircoTaxAutoDA <- function(
	data,
	contrast,
	design = NULL,
	tax_out = NULL,
	mode = c("DESeq2","ALDEx2"),
	MCS = 128,
	p_adjustment = "BH", # Benj Holch
	eff_size = 0, # Effect size threshold for ALDEx2 results
	p_val = 0.05,
	B = 50, # used as BaseMean Threshold in DESeq2 section
	lfc = 1, # log fold change (DESeq2)
	lfc_shrink_method= "none",
	remove_redundants= TRUE,
	remove_results="" # to manually discard bad results
) {

	# checking if the object in input is a phyloseq object
  if ( is.null(phyloseq::sample_names(data)) ){
    stop( "\n\n Looks like the object in input has no 'sample names'. \nPlease check if this object is a normal phyloseq object (containing otu table, tax table and sample data) ")
  }
	# example of a proper contrast -->   contrast=c("Condition,Healthy,Tumor")
  if(length(contrast) != 3 ){ # then if badly written or it does not exist
    stop(paste(
      "\n\n Please define a proper 'contrast' argument \n(  'Column_name_Factor','reference_Level','Level_B',   for example:  c('Condition','Healthy','Tumor')   )\n\n",
      "You have the following column(s) available in your phyloseq sample data:",paste(colnames(phyloseq::sample_data(data)), collapse= "  "),"\n",
      "\nEach one of the levels in the specified contrast has to be a row of the target factor column in the phyloseq sample data\n\n")
    )
  }

	if (is.character(phyloseq::sample_data(data)[[contrast[1]]])) {
    phyloseq::sample_data(data)[[contrast[1]]]<-as.factor(phyloseq::sample_data(data)[[contrast[1]]])
	}
	colnames(phyloseq::tax_table(data))[colnames(phyloseq::tax_table(data))=="Kingdom"]<-"Domain"
	tax_in <- c("Phylum","Class","Order","Family","Genus")
  tax_in <- tax_in [! tax_in %in% tax_out ]   # tax_out is empty by default
  
	if (!is.null(design)) {
		there_it_is<-0
    split<-gsub("~","",design,fixed = T)
    split<-gsub(" ","",split,fixed = T)
    split<-unlist(strsplit(split, "+", fixed = T))
    for( i in 1:length(split) ){
      if( split[i] %in% contrast){
        there_it_is<-there_it_is + 1
      }
    }
	} else {
		design <- paste("~",contrast[1])
	}

	if (tolower(mode[1]) == "aldex2") {
		metadat <- as.data.frame(phyloseq::sample_data(data))
		if( length( unique(metadat[[contrast[[1]]]]) ) > 2 ){ # if more than two levels in this column
			metadat[[contrast[1]]] <- factor(
				metadat[[contrast[1]]],
				levels = c(
					contrast[[2]],  # reference/base
					contrast[[3]],  # higher level (target),
					metadat[[contrast[1]]] [! metadat[[contrast[1]]] %in% c(contrast[[2]], contrast[[3]])]
				)
			)
		}
		if ( class(metadat[[contrast[1]]])!="factor" ){ # if it has only 2 levels it won't be coerced by the function above (if it is not already a factor by itself)
			metadat[[contrast[1]]] <- factor(metadat[[contrast[1]]] , levels = c(contrast[[2]], contrast[[3]] ) )
		}
	}
	suppressWarnings({
		rm(data_pruned, data.genus_pruned)
		data_pruned <- phyloseq::prune_taxa(phyloseq::taxa_sums(data) > 10, data)
	})
	Table_tot<-NULL
	Res_tot<-NULL
	resultsLength <- 0
	for (t in tax_in) {
		cat("\nWorking on",t,"level...\n\n")
    suppressWarnings(rm(list=c("d", "d.prop", "Taxa.d", "res","DE", "target", "r", "r_level", "index","names_target", "logical_to_subset")))
		d <- phyloseq::tax_glom(data_pruned, taxrank = t, NArm = F)
    d.prop <- phyloseq::transform_sample_counts(d, function(x) x/sum(x)*100)
		if(t == "Genus"){
			taxa_temp <- as.data.frame(phyloseq::tax_table(d))
			is_uncultured <- taxa_temp$Genus == "uncultured"
			taxa_temp$Genus[is_uncultured] <- paste0("uncultured_ f ", taxa_temp$Family[is_uncultured])
			fix_uf_u <- taxa_temp$Genus == "uncultured_ f uncultured"
			taxa_temp$Genus[fix_uf_u] <- paste0("uncultured_ o ", taxa_temp$Order[fix_uf_u])
			is_na_genus <- is.na(taxa_temp$Genus)
			taxa_temp$Genus[is_na_genus] <- paste0("NA_ f ", taxa_temp$Family[is_na_genus])
			fix_nf_na <- taxa_temp$Genus == "NA_ f NA"
			taxa_temp$Genus[fix_nf_na] <- paste0("NA_ o ", taxa_temp$Order[fix_nf_na])
			taxa_temp$Genus <- make.unique(taxa_temp$Genus)
			phyloseq::tax_table(d) <- as.matrix(taxa_temp)
			phyloseq::tax_table(d.prop) <- as.matrix(taxa_temp)
			rm(taxa_temp)
		}

		if (tolower(mode[1]) == "aldex2") {
			res = autoDAAldex2(
				data = d,
				metadata = metadat,
				design = design,
				contrast = contrast,
				MCS = MCS,
				p_adjustment = p_adjustment,
				p_val = 1,
				eff_size = eff_size
			)
			# results<-row.names(res)
			r<-as.data.frame(res)
			Taxa.d<-as.data.frame(phyloseq::tax_table(d))
			Taxa.d$ASV<-row.names(Taxa.d)
			r<- cbind.data.frame(Taxa.d[row.names(r), ] , r )
			r<-r[ r[,t] !="none", ] # removing artifacts from glomming, composed of different taxa with no name for that level in NCBI
			if (length(r[,1]) == 0 ) {
				cat("Any results for the",t,"level\n")
      	Sys.sleep(1)
				next
			}
			r_level<-r
      r_level[, "Taxon"]<- rep(t)
      r_level <- r_level[r_level[,t] != "none" , ]
      Res_tot<-rbind.data.frame(Res_tot,r_level)
			### beginning to build the basis for the boxplot (Table_DE)... 
      names_target <- as.character(r[,t])
      colnames(phyloseq::tax_table(d.prop))[colnames(phyloseq::tax_table(d.prop))==t]<-"Aimed_taxa"
      logical_to_subset <- as.data.frame(phyloseq::tax_table(d.prop))[["Aimed_taxa"]] %in% names_target
      # NB: if the function "subset_taxa" with %in% is used when launching this script as an function, then there will be a strange bug related to %in%  --> using prune_taxa with logical subset
      target <- suppressWarnings(phyloseq::prune_taxa(logical_to_subset, d.prop)) # cannot use t %in% target in this function, then it's scripted in this way
      # suppressWarnings avoids message related to the phylogenetic tree (useless here!) or to the column type (e.g. 'sample' in 'sample_sample') 
      Table_DE <- phyloseq::psmelt(target)
      colnames(Table_DE)[colnames(Table_DE)=="Aimed_taxa"]<-t # restored the original name
      ### appending to a unique box plot ...
      index<- which(colnames(Table_DE)=="Domain") : which(colnames(Table_DE)==t) # from : to
      index<- index[-length(index)] # removing the last index, regarding the taxa of interest
      Table_DE[,index]<-NULL
      Table_DE$Taxa<-t
      colnames(Table_DE)[colnames(Table_DE)==t]<-"Bacteria"
      Table_tot<-rbind.data.frame(Table_tot, Table_DE)
			resultsLength<-length(Res_tot[[1]])
		}
		if (tolower(mode[1]) == "deseq2") {
			res = autoDADESeq2(
				data = d,
				design = design,
				contrast = contrast,
				p_adjustment = p_adjustment,
				lfc_shrink_method = lfc_shrink_method,
				p_val = p_val,
				lfc = lfc,
				B = B
			)
			if (length(res[,1]) == 0 ) {
				cat("Any results for the",t,"level\n")
      	Sys.sleep(1)
				next
			}
			r<-as.data.frame(res)
      r$ASV<-row.names(r)
      Taxa.d<-as.data.frame(phyloseq::tax_table(d))
      # Taxa.d$ASV<-row.names(Taxa.d)
      Taxa.d <- Taxa.d[row.names(Taxa.d) %in% r$ASV, ]
      r<-cbind.data.frame( r[row.names(Taxa.d), ], Taxa.d )
      r$Domain<-NULL
      r$Species<-NULL 
      assign(paste(t,"results",sep="_"), r)
      r_level<-r
      r_level[, "Taxon"]<- rep(t)
      Res_tot<-rbind.data.frame(Res_tot,r_level)
      ### beginning to build the basis for the boxplot (Table_DE)...
      names_target <- as.character(r[,t])
      colnames(phyloseq::tax_table(d.prop))[colnames(phyloseq::tax_table(d.prop))==t]<-"Aimed_taxa"
      # NB: if the function "subset_taxa" with %in% is used when launching this script as an external function, then there will be a strange bug related to %in%  --> using prune_taxa with logical subset
      logical_to_subset <- as.data.frame(phyloseq::tax_table(d.prop))[["Aimed_taxa"]] %in% names_target
      target <- suppressWarnings(phyloseq::prune_taxa(logical_to_subset, d.prop)) # cannot use t %in% target in this function, then it's scripted in this way
      # NB: suppressWarnings avoids message related to the phylogenetic tree (useless here!) or to the column type (e.g. 'sample' in 'sample_sample') 
      Table_DE<- suppressWarnings(phyloseq::psmelt(target))
      colnames(Table_DE)[colnames(Table_DE)=="Aimed_taxa"]<-t # restored the original name
      Table_DE$ASV<-NULL
      # assign(paste("Table_DE_plot",t,sep="_"), Table_DE)
      # appending each box to unique box plot
      index<- which(colnames(Table_DE) %in% c("Domain","Kingdom")) : which(colnames(Table_DE)==t)
      index2<- index[-length(index)] # removing the last index, regarding the taxa of interest
      Table_DE[,index2]<-NULL
      Table_DE$Taxa<-t
      colnames(Table_DE)[colnames(Table_DE)==t]<-"Bacteria"
      Table_tot<-rbind.data.frame(Table_tot, Table_DE)
			resultsLength<-length(Res_tot$ASV)
		}
		
	}

	if (resultsLength == 0) {
		cat("\n\n Looks like there are not significant differences there! \n\n")
    Sys.sleep(2)
		return()
	}

	if (remove_redundants) {
		auto_redund <- autoDARedundancyHandler(
			table_tot = Table_tot,
			res_tot = Res_tot
		)
	} else {
		auto_redund <- c()
  }

	if(length(remove_results) > 1 & ! remove_results %in% as.vector(phyloseq::tax_table(data))) {
		cat("\n Warning: the specified results to remove are not in the taxonomic table! \n\n")
	}
 
	auto_redund <- c(auto_redund, remove_results)

	dfCirco = autoDACircoTaxDFMaker(
		res_tot = Res_tot,
		auto_redund = auto_redund,
		tax_in = tax_in
	)

	return(dfCirco)
}



autoDARedundancyHandler <- function(
	table_tot,
	res_tot
){
	# 1) if redundant then same abundances
	auto_Max<-tapply(round(table_tot$Abundance,0), table_tot$Bacteria, max )
	auto_Min<-tapply(round(table_tot$Abundance,0), table_tot$Bacteria, min )
	auto_Median<-tapply(round(table_tot$Abundance,0), table_tot$Bacteria, stats::median )
	# reordering the names (lower ranks first --> the redundants are those at higher levels)
	ordered_names <- c( 
		unique(table_tot[table_tot$Taxa=="Genus", "Bacteria" ]),
		unique(table_tot[table_tot$Taxa=="Family", "Bacteria" ]),
		unique(table_tot[table_tot$Taxa=="Order", "Bacteria" ]),
		unique(table_tot[table_tot$Taxa=="Class", "Bacteria" ]),
		unique(table_tot[table_tot$Taxa=="Phylum", "Bacteria" ])
	)
	ordered_names<-ordered_names[!is.na(ordered_names)]
	auto_Max<-auto_Max[ ordered_names  ]
	auto_Min<-auto_Min[ ordered_names  ]
	auto_Median<-auto_Median[ ordered_names ] # the table it self follows the correct taxonomic order because it is built in this way
	# same abundances
	auto_redund<-names(auto_Max)[ duplicated(auto_Min) & duplicated(auto_Max) & duplicated(auto_Median)  ]

	# 2) if redundant then same ASV
	auto_ASV<-res_tot$ASV[duplicated(res_tot$ASV)]
	auto_ASV_Names<-unique(table_tot$Bacteria[table_tot$OTU %in% auto_ASV])
	auto_redund<-auto_redund[auto_redund %in% auto_ASV_Names]
      
	# 3) eventual redundant from the lower taxonomic level are not redundant!
	lower_level<-res_tot$Taxon
	lower_level<-factor(res_tot$Taxon)
	lower_level<-lower_level[match(lower_level, c("Phylum","Class","Order","Family","Genus")) ] # this row ensures that the order is ALWAYS correct!  
	lower_level<-as.character(lower_level[length(lower_level)]) # selects the last one
	if(! is.na(lower_level) ){ # it is NA if the target level is the Genus (no Species!) --> error in the row below
		auto_redund<-auto_redund[! auto_redund %in% res_tot[,lower_level] ]  # NB: If lower level is still a factor then its string would be seen as number by R --> wrong column
	}

	return(auto_redund)
}



autoDACircoTaxDFMaker <- function(
	res_tot,
	auto_redund = c(),
	tax_in
) {
	### preparing the objects
	Res_tot2<-res_tot[! grepl("uncultured", res_tot$Genus), ] # those uncultured does not have a name and here they are almost absent... better avoiding to plot those noises!
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

	Res_tot3 <- data.frame(log2FC = Res_tot2[,2], padj = Res_tot2[,6]  ,Res_tot2[8:12])

	return(Res_tot3)
}

autoDADESeq2 <- function(
	data,
	design,
	contrast,
	lfc_shrink_method = NULL,
	p_val = 0.05,
	p_adjustment = "BH",
	lfc = 1,
	B = 50
){
	DEseq_data<-phyloseq::phyloseq_to_deseq2(data, design = stats::formula(design) )
	DE<-DESeq2::DESeq(DEseq_data)
	res = switch(
		lfc_shrink_method,
		apeglm = {
			coef_input <- DESeq2::resultsNames(DE)[ grepl(pattern = contrast[1], DESeq2::resultsNames(DE)) & grepl(pattern = contrast[2], DESeq2::resultsNames(DE))]
      DESeq2::lfcShrink(DE, coef= coef_input, type = lfc_shrink_method ) 
		},
		ashr = {
      DESeq2::lfcShrink(DE, contrast= contrast, type = lfc_shrink_method )
		},
		{
      DESeq2::results(DE, contrast= contrast) # In DESeq2 itself the contrast is specified as follows contrast=c("Condition","Healthy","RCC")
		}
	)
	res$padj<- stats::p.adjust( res$pvalue , method = p_adjustment )   # NB: is "BH" by default, but can be overwritten if otherwise specified 
	res = res[order(res$padj, na.last=NA), ]
	res<-res[(res$padj < p_val ) & abs(res$log2FoldChange)> lfc ,]
	res<-res[res$baseMean > B, ] # arbitrary threshold to avoid the most noisy result
  return(res)
}

autoDAAldex2 <- function(
	data,
	metadata,
	design,
	contrast,
	MCS,
	p_adjustment,
	p_val,
	eff_size = 0
) {
	mm <- stats::model.matrix( stats::formula(design), metadata )
	set.seed(1)
	aldx <- ALDEx2::aldex.clr(
		phyloseq::otu_table(data),
		mm,
		mc.samples=MCS,   #"DMC"
		# according to ALDEx2 vignette, the number of samples in the smallest group multiplied by the number of DMC be equal at least to 1000
		denom="all", # every feature as denominator
		verbose=T
	)
	aldx2 <- ALDEx2::aldex.glm(aldx, verbose = T)
	#### THEN, the summary of the results
	aldx3 <- ALDEx2::aldex.glm.effect(
		aldx, # calculates the effect size for every binary variable in the mm
		CI = TRUE # confidence interval of the effect size
	) 
	aldx_final <- data.frame(aldx2,aldx3)

	# Automatically searching the key columns in this dataframe ...
	search_contrast <- grepl(contrast[[2]],colnames(aldx_final))| grepl(contrast[[3]],colnames(aldx_final))
	search_factor <- grepl(contrast[[1]],colnames(aldx_final))
	search_pval_column <- grepl("pval",colnames(aldx_final)) & !grepl("adj",colnames(aldx_final))
	search_eff_column <- grepl("effect$",colnames(aldx_final))    # the symbol avoids the selection of the low and high CI columns with the eff size
	
	column_with_Pval<- colnames(aldx_final) [ search_factor & search_contrast & search_pval_column  ]
	column_with_eff_size<- colnames(aldx_final) [ search_factor & search_contrast & search_eff_column  ]
	
	resulted_p_val<-aldx_final[[column_with_Pval]]
	p_adj<- stats::p.adjust(resulted_p_val, method= p_adjustment )
	
	resulted_eff_size<-aldx_final[[column_with_eff_size]]
	Enough_eff_size <- resulted_eff_size > eff_size    # logical vector
	res<- aldx_final
	res$p_adj_column<- p_adj
	res<- res[ res$p_adj_column<p_val & Enough_eff_size, ]
}