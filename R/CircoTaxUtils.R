
CircoTaxTableFormatter <- function(
	input_table,
	names_col = 1,
	tax_col = 2,
	fc_col = 3,
	rank_vector = c("Domain","Phylum","Class","Order","Family","Genus","Species"),
	collapse_ranks = FALSE
) {
	if( any(c("Kingdom") %in% colnames(input_table)) ) {
		rank_vector = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
	}
	ranksN <- 1:length(rank_vector)

	if (length(tax_col) > 1) {
		if (collapse_ranks) {
			i<-0
			toRemove <- c()
			for (check in lapply(
				colnames(input_table[,tax_col]),
				function(x) all(is.na(input_table[[x]]))
			)) {
				i<-i+1
				if (check) {
					toRemove <- c(toRemove, i+1)
				}
			}
			if (length(toRemove) > 0) {
				input_table <- input_table[ -c(toRemove) ]
			}
			extra_cols = length(colnames(input_table)) - length(tax_col)
			uniqueRanksTagsInTable <- colnames(input_table[-c(1:extra_cols)])
		} else {
			uniqueRanksTagsInTable <- colnames(input_table[tax_col])
		}
	} else {
		if (collapse_ranks) {
			uniqueRanksTagsInTable <- unique(input_table[,tax_col])
		} else {
			uniqueRanksTagsInTable <- rank_vector
		}
	}
	uniqueRanksNumsInTable<-stringi::stri_replace_all_regex(uniqueRanksTagsInTable,rank_vector,ranksN,vectorize_all=FALSE)
	uniqueRanksNumsInTableO <- sort(uniqueRanksNumsInTable)
	uniqueRanksTagsInTableO<-stringi::stri_replace_all_regex(uniqueRanksNumsInTableO,ranksN,rank_vector,vectorize_all=FALSE)
	if (length(tax_col) > 1) {
		dataFrame <- data.frame(FC = input_table[,fc_col],input_table[uniqueRanksTagsInTableO])
		return(dataFrame)
	} else {
		dataFrame <- data.frame(FC = input_table[,fc_col])
		dataFrame[uniqueRanksTagsInTableO] = NA
	}

	names <- input_table[,names_col]
	for (i in 1:length(names)) {
		for (c in uniqueRanksTagsInTableO) {
			if (input_table[i,tax_col] == c) {
				dataFrame[[i,c]] = names[i]
				break
			} else {
				dataFrame[[i,c]] = "unk"
			}
		}
	} 

	return(dataFrame)
}


CircoTaxSorter <- function(
	input_table,
	sort = "rank",
	decrease_order = FALSE,
	tax_col=2:length(colnames(input_table)),
  fc_col= ifelse(length(tax_col)>1,1,3)
) {
	y = apply(input_table,1,function(x) length(tax_col)-sum(ifelse(is.na(x[tax_col]),1,0)))
	settings = switch(
		sort[1],
		rank = {
			o=order(y,decreasing=decrease_order)
			synth=apply(input_table[o,tax_col],1,function(x) paste0(x,collapse="-"))
			labels=unlist(lapply(strsplit(gsub("-NA.*","",perl=T,synth),"-"),function(x) rev(x)[1]))
			y=y[o]
			fc=input_table[o,fc_col]
			list(labels = labels, y = y, fc = fc)
		},
		alpharank = {
			synth=apply(input_table[,tax_col],1,function(x) paste0(x,collapse="-"))
			o=order(synth,decreasing=decrease_order)
			synth=synth[o]
			labels=unlist(lapply(strsplit(gsub("-NA.*","",perl=T,synth),"-"),function(x) rev(x)[1]))
			y=y[o]
			fc=input_table[o,fc_col]
			list(labels = labels, y = y, fc = fc)
		},
		alpha = {
			synth=apply(input_table[,tax_col],1,function(x) paste0(x,collapse="-"))
			labels=unlist(lapply(strsplit(gsub("-NA.*","",perl=T,synth),"-"),function(x) rev(x)[1]))
			o=order(labels)
			labels=labels[o]
			y=y[o]
			fc=input_table[o,fc_col]
			list(labels = labels, y = y, fc = fc)
		},
		fc = {
			o=order(input_table[,fc_col],decreasing=decrease_order)
			synth=apply(input_table[,tax_col],1,function(x) paste0(x,collapse="-"))
			synth=synth[o]
			labels=unlist(lapply(strsplit(gsub("-NA.*","",perl=T,synth),"-"),function(x) rev(x)[1]))
			y=y[o]
			fc=input_table[o,fc_col]
			list(labels = labels, y = y, fc = fc)
		},
		absfc = {
			o=order(abs(input_table[,fc_col]),decreasing=decrease_order)
			synth=apply(input_table[o,tax_col],1,function(x) paste0(x,collapse="-"))
			labels=unlist(lapply(strsplit(gsub("-NA.*","",perl=T,synth),"-"),function(x) rev(x)[1]))
			y=y[o]
			fc=input_table[o,fc_col]
			list(labels = labels, y = y, fc = fc)
		},
		{
			synth=apply(input_table[,tax_col],1,function(x) paste0(x,collapse="-"))
			labels=unlist(lapply(strsplit(gsub("-NA.*","",perl=T,synth),"-"),function(x) rev(x)[1]))
			y=y
			fc=input_table[,fc_col]
			list(labels = labels, y = y, fc = fc)
		}
	)
	return(settings)
}
