
CircoTaxTableFormatter <- function(
	input_table,
	names_col = 1,
	tax_col = 2,
	fc_col = 3,
	pval_col = NA,
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
		if (!is.na(pval_col)) {
			dataFrame <- data.frame(FC = input_table[,fc_col], PVAL = input_table[,pval_col], input_table[uniqueRanksTagsInTableO])
		} else {
			dataFrame <- data.frame(FC = input_table[,fc_col], PVAL = NA, input_table[uniqueRanksTagsInTableO])
		}
		return(dataFrame)
	} else {
		if (!is.na(pval_col)) {
			dataFrame <- data.frame(FC = input_table[,fc_col], PVAL = input_table[,pval_col])
		} else {
			dataFrame <- data.frame(FC = input_table[,fc_col], PVAL = NA)
		}
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
	tax_col = 2:length(colnames(input_table)),
	fc_col = ifelse(length(tax_col) > 1, 1, 3),
	pval_col = NA
) {
	tax <- input_table[, tax_col, drop = FALSE]

	# last non-NA name; "NA_ <code> <name>" placeholders are kept on purpose
	tax_label <- function(x) {
		x <- x[!is.na(x)]
		if (length(x)) x[length(x)] else NA_character_
	}

	y      <- apply(tax, 1, function(x) length(tax_col) - sum(is.na(x)))
	labels <- apply(tax, 1, tax_label)
	synth  <- apply(tax, 1, paste0, collapse = "-")

	o <- switch(
		sort[1],
		rank      = order(y, decreasing = decrease_order),
		alpharank = order(synth, decreasing = decrease_order),
		alpha     = order(labels),
		fc        = order(input_table[, fc_col], decreasing = decrease_order),
		absfc     = order(abs(input_table[, fc_col]), decreasing = decrease_order),
		pval      = order(input_table[, pval_col], decreasing = decrease_order),
		seq_len(nrow(input_table)) # unrecognized sort: keep the input order
	)

	list(
		labels = labels[o],
		y      = y[o],
		fc     = input_table[o, fc_col],
		pval   = input_table[o, pval_col]
	)
}
