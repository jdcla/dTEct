suppressPackageStartupMessages({
  library(edgeR)
  library(EnhancedVolcano)
  library(ggplot2)
  library(optparse)
  library(tidyverse)
  library(DESeq2)
  library(pheatmap)
  library(RColorBrewer)
  library(DEFormats)
  library(BiocParallel)
  library(svglite)
})
options(show.error.locations = TRUE)

# MDS -------------------------------------------------------------------------
eval_MDS <- function(dge.set, meta, cols_of_interest, prefix, suffix) {
  xy_dim <- max(dim(meta)[1]*0.10, 20)
  # MDS
  pcaData <- plotMDS(dge.set, ntop=Inf, returnData = TRUE, ntop=Inf)
  percentVar <- round(100 * pcaData$var.explained)
  plot_data <- data.frame(
    x = pcaData$x,
    y = pcaData$y,
    rep = meta$rep
  )
  for (col in cols_of_interest) {
    plot_data[[col]] <- meta[[col]]
  }
  for (i in seq_along(cols_of_interest)) {
    col <- cols_of_interest[[i]]
    plot <- ggplot(plot_data, aes(x = x, y = y, color=meta[[col]], shape=meta$rep)) +
      geom_point(size = 2, alpha = 0.8) +
      coord_fixed() + 
      geom_text(aes(label = meta$rep), size=4, color = "black") +
      xlab(paste0("PC1: ", percentVar[1], "% variance explained")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance explained")) +
      ggtitle(paste0("MDS on fitted ", suffix, " data (all counts/edgeR)")) +
      labs(color = col, shape="replicate")
    gfig(plot, str_c(prefix, "MDS_", i, "_", suffix), xy_dim, xy_dim)
    }
  }

# PCA -------------------------------------------------------------------------
eval_PCA <- function(dge.set, meta, cols_of_interest, prefix, suffix) {
  xy_dim <- max(dim(meta)[1]*0.10, 20)
  # PCA
  pcaData <- plotMDS(dge.set, ntop=Inf, returnData = TRUE, gene.selection="common")
  percentVar <- round(100 * pcaData$var.explained)
  plot_data <- data.frame(
    x = pcaData$x,
    y = pcaData$y,
    rep = meta$rep
  )
  for (col in cols_of_interest) {
    plot_data[[col]] <- meta[[col]]
  }
  for (i in seq_along(cols_of_interest)) {
    col <- cols_of_interest[[i]]
    plot <- ggplot(plot_data, aes(x = x, y = y, color=meta[[col]], shape=meta$rep)) +
      geom_point(size = 2, alpha = 0.8) +
      coord_fixed() + 
      geom_text(aes(label = meta$rep), size=4, color = "black") +
      xlab(paste0("PC1: ", percentVar[1], "% variance explained")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance explained")) +
      ggtitle(paste0(suffix, ": PCA on fitted data (all counts/edgeR)")) +
      labs(color = col, shape="replicate")
    gfig(plot, str_c(prefix, "PCA_", i, "_", suffix), xy_dim, xy_dim)
  }
}

# Distance Heatmap ------------------------------------------------------------

eval_heatmap <- function(norm_counts, meta, meta_cols, prefix, suffix) {
  xy_dim <- max(dim(meta)[1]*0.05, 5)
  sampleDists <- dist(t(norm_counts))
  sampleDistMatrix <- as.matrix(sampleDists)
  colnames(sampleDistMatrix) <- meta$counts_col 
  rownames(sampleDistMatrix) <- apply(meta[, meta_cols], 1, paste, collapse = "__")
  # colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  for (fmt in c("png", "svg")) {
    if (fmt == "png") {
      pngfig(str_c(prefix, "Heatmap_1_", suffix), xy_dim, xy_dim)
      pheatmap(sampleDistMatrix,
        clustering_distance_rows = sampleDists,
        clustering_distance_cols = sampleDists,
        # col = colors,
        main = paste0(suffix, ": Euclidean distance between counts of samples on fitted data.")
      )
      w()
    } else {
      svgfig(str_c(prefix, "Heatmap_1_", suffix), xy_dim, xy_dim)
      pheatmap(sampleDistMatrix,
        clustering_distance_rows = sampleDists,
        clustering_distance_cols = sampleDists,
        # col = colors,
        main = paste0(suffix, ": Euclidean distance between counts of samples on fitted data.")
      )
      w()
    }
  }
  # Sort the column and row names
  sorted_col_order <- order(colnames(sampleDistMatrix))
  sorted_row_order <- order(rownames(sampleDistMatrix))
  # Reorder the sampleDistMatrixrix based on the sorted column and row names
  sorted_sampleDistMatrix <- sampleDistMatrix[sorted_row_order, sorted_col_order]
  # Update column and row names to reflect the sorted order
  colnames(sorted_sampleDistMatrix) <- colnames(sampleDistMatrix)[sorted_col_order]
  rownames(sorted_sampleDistMatrix) <- rownames(sampleDistMatrix)[sorted_row_order]
  ## WRITE SAMPLE DISTANCES TO FILE
  write.table(
    cbind(sample = rownames(sorted_sampleDistMatrix), sorted_sampleDistMatrix),
    file=str_c(prefix, "sample.dists_", paste0(suffix, ".txt")),
    row.names=FALSE,
    col.names=TRUE,
    sep="\t",
    quote=FALSE
  )
}

# Gene Clustering -------------------------------------------------------------
eval_gene_clusters <- function(norm_counts, meta, meta_cols, prefix, suffix) {
  x_dim <- max(dim(meta)[1]*0.18, 7)
  k_genes <- max(dim(meta)[1]*0.7, 50)
  y_dim <- k_genes * 0.21
  for (i in 1:1) {
    topVarGenes <- head(order(rowVars(norm_counts), decreasing = TRUE), k_genes)
    mat <- norm_counts[topVarGenes, ]
    mat <- mat - rowMeans(mat)
    rownames(mat) <- feature2name[match(rownames(mat), feature2name[,"id"]), "name"]
    anno <- as.data.frame(meta[,c(meta_cols, "rep", "batch_date")])
    if (i == 1) {
      rownames(anno) <- colnames(mat)
    }
    for (fmt in c("png", "svg")) {
      if (fmt == "png") {
        pngfig(str_c(prefix, "GeneClusters_1_", suffix), x_dim, y_dim)
      } else {
        svgfig(str_c(prefix, "GeneClusters_1_", suffix), x_dim, y_dim)
      }
      pheatmap(mat, annotation_col = anno, main=suffix)
      w()
    }
  }
}

evaluate_unique_contrast <- function(meta, uniq_val, contrast_col, fit, outdir, log_file) {
  strat_string <- ""
  ribo_mask <- meta$seq_type == "Ribo"
  rna_mask <- meta$seq_type == "RNA"
  num_dots <- str_count(uniq_val, "\\.") + 1
  # Take into account combinatorial groups
  tmp_col <- paste0(contrast_col, ".", num_dots)
  # Check condition for _and_ in uniq_val
  if (!grepl("_and_", uniq_val) && any(grepl("_and_", meta[[tmp_col]]))) {
    # Split the entries by _and_
    split_entries <- strsplit(as.character(meta[[tmp_col]]), "_and_")
    # Determine which rows match
    uniq_mask <- sapply(split_entries, function(parts) {
      # Return TRUE if uniq_val matches either part
      uniq_val %in% parts
    })
  } else {
    uniq_mask <- meta[[tmp_col]] == uniq_val
  }
  ribo_grp_mask <- uniq_mask & ribo_mask
  rna_grp_mask <- uniq_mask & rna_mask
  n_rna <- sum(ribo_grp_mask)
  n_ribo <- sum(rna_grp_mask)
  if (sum(ribo_grp_mask) != 0 && sum(rna_grp_mask) != 0) {
    ribo_vals <- unique(meta[ribo_grp_mask, contrast_col])
    # Weigh factors according to contribution to supergroup
    weight_dict_ribo <- as.list(table(meta[ribo_grp_mask, contrast_col])/sum(ribo_grp_mask))
    rna_vals <- unique(meta[rna_grp_mask, contrast_col])
    weight_dict_rna <- as.list(table(meta[rna_grp_mask, contrast_col])/sum(rna_grp_mask))
    uniq_ribo <- construct_contrast_string(ribo_vals[[contrast_col]], "Ribo", meta, weight_dict_ribo)
    # TODO change disease_type entries to contain .NA up to the number of subclasses within the dataset!!!!!
    uniq_rna <- construct_contrast_string(rna_vals[[contrast_col]], "RNA", meta, weight_dict_rna)
    # TE
    uniq_contrast <- paste0("makeContrasts(", uniq_ribo, " - ", uniq_rna, ", levels=design)")
    if (!is.null(uniq_contrast)) {
      print(uniq_contrast)
    }
    contrast_id <- paste0(uniq_val, strat_string, "_TE")
    out_prefix <- paste0(outdir, "TE/", contrast_id)
    title <- paste0(contrast_id, "    (n = ", n_ribo, " / ", n_rna, ")")
    eval_contrast(fit, uniq_contrast, out_prefix, title, log_file)
  }
}

evaluate_combination_contrast <- function(meta, tup, contrast_col, fit, outdir, log_file) {
  ribo_mask <- meta$seq_type == "Ribo"
  rna_mask <- meta$seq_type == "RNA"
  grp_rna <- list()
  grp_ribo <- list()
  n_rna <- list()
  n_ribo <- list()
  order <- list()

  for (i in seq_along(tup)) {
    num_dots <- str_count(tup[[i]], "\\.") + 1
    # Take into account combinatorial groups
    tmp_col <- paste0(contrast_col, ".", num_dots)
    # Check condition for _and_ in tup[[i]]
    if (!grepl("_and_", tup[[i]]) && any(grepl("_and_", meta[[tmp_col]]))) {
      # Split the entries by _and_
      split_entries <- strsplit(as.character(meta[[tmp_col]]), "_and_")
      # Determine which rows match
      grp_mask <- sapply(split_entries, function(parts) {
        # Return TRUE if tup[[i]] matches either part
        tup[[i]] %in% parts
      })
    } else {
      grp_mask <- meta[[tmp_col]] == tup[[i]]
    }
    ribo_grp_mask <- grp_mask & ribo_mask
    n_ribo[[i]] <- sum(ribo_grp_mask)
    rna_grp_mask <- grp_mask & rna_mask
    n_rna[[i]] <- sum(rna_grp_mask)
    if (sum(ribo_grp_mask) != 0) {
      ribo_vals <- unique(meta[ribo_grp_mask, tmp_col])
      weight_dict_ribo <- as.list(table(meta[ribo_grp_mask, tmp_col])/sum(ribo_grp_mask))
      grp_ribo[[i]] <- construct_contrast_string(ribo_vals[[tmp_col]], "Ribo", meta, weight_dict_ribo)
      order[[i]] <- sum(meta[ribo_grp_mask, "control"] == "test")
    }
    if (sum(rna_grp_mask) != 0) {
      rna_vals <- unique(meta[rna_grp_mask, tmp_col])
      weight_dict_rna <- as.list(table(meta[rna_grp_mask, tmp_col])/sum(rna_grp_mask))
      grp_rna[[i]] <- construct_contrast_string(rna_vals[[tmp_col]], "RNA", meta, weight_dict_rna)
      order[[i]] <- sum(meta[rna_grp_mask, "control"] == "test")
    }
  }
  # Switch terms if first term has more "test" hits than second term
  if ((length(order) == 2) && (order[[1]] < order[[2]])) {
    if (length(grp_ribo) == 2) {
      grp_ribo <- grp_ribo[c(2, 1)]
    }
    if (length(grp_rna) == 2) {
      grp_rna <- grp_rna[c(2, 1)]
    }
    tup <- tup[c(2, 1)]
  }
  evaluate_contrasts(grp_rna, grp_ribo, tup, n_rna, n_ribo, fit, outdir, log_file)
}

evaluate_contrasts <- function(grp_rna, grp_ribo, tup, n_rna, n_ribo, fit, outdir, log_file) {
  strat_string <- ""
  if (length(unlist(grp_rna)) == 2) {
    # Evaluate RNA contrast
    grp_contrast <- paste0("makeContrasts(", paste(grp_rna, collapse = " - "), ", levels=design)")
    contrast_id <- paste0(tup[[1]], "__", tup[[2]], strat_string, "_RNA")
    out_prefix <- paste0(outdir, "RNA/", contrast_id)
    title <- paste0(contrast_id, "    (n = ", n_rna[[1]], " / ", n_rna[[2]], ")")
    if (!is.null(grp_contrast)) {
      print(grp_contrast)
    }
    eval_contrast(fit, grp_contrast, out_prefix, title, log_file)
  }
  if (length(unlist(grp_ribo)) == 2) {
    # Evaluate Ribo contrast
    grp_contrast <- paste0("makeContrasts(", paste(grp_ribo, collapse = " - "), ", levels=design)")
    contrast_id <- paste0(tup[[1]], "__", tup[[2]], strat_string, "_Ribo")
    out_prefix <- paste0(outdir, "Ribo/", contrast_id)
    title <- paste0(contrast_id, "    (n = ", n_ribo[[1]], " / ", n_ribo[[2]], ")")
    if (!is.null(grp_contrast)) {
      print(grp_contrast)
    }
    eval_contrast(fit, grp_contrast, out_prefix, title, log_file)
  }
  if ((length(unlist(grp_ribo)) == 2) && (length(unlist(grp_rna)) == 2)) {
    # Evaluate differential translation efficiency (dTE) contrast
    left_cmd <- paste0("(", grp_ribo[[1]], " - ", grp_rna[[1]], ")")
    right_cmd <- paste0("(", grp_ribo[[2]], " - ", grp_rna[[2]], ")")
    grp_contrast <- paste0("makeContrasts(", left_cmd, " - ", right_cmd, ", levels=design)")
    contrast_id <- paste0(tup[[1]], "__", tup[[2]], strat_string, "_dTE")
    out_prefix <- paste0(outdir, "dTE/", contrast_id)
    title <- paste0(contrast_id, "    (n = ", n_ribo[[1]], " / ", n_rna[[1]], " / ", n_ribo[[2]], " / ", n_rna[[2]], ")")
    if (!is.null(grp_contrast)) {
      print(grp_contrast)
    }
    eval_contrast(fit, grp_contrast, out_prefix, title, log_file)
  }
}

pngfig <- function(filename, width=7, height=7) {
  w <- width*90
  h <- height*90
  f_png <- str_c(filename, ".png")
  print(str_c("Outputting: ", f_png, " as PNG, ", w, "x", h))
  png(f_png, width=w, height=h, units="px")
}

svgfig <- function(filename, width=7, height=7) {
  w <- width*2
  h <- height*2
  f_svg <- str_c(filename, ".svg")
  print(str_c("Outputting: ", f_svg, " as SVG, ", w, "x", h))
  svg(f_svg, width=w, height=h)
}

w <- function() {
  dev.off()
}

gfig <- function(plot, filename, width=7, height=7) {
  w <- width*90
  h <- height*90
  # Always save as both PNG and SVG
  f_png <- paste0(filename, ".png")
  print(str_c("Outputting: ", f_png, " as PNG, ", w, "x", h))
  ggsave(f_png, plot=plot, width=w, height=h, units="px", dpi=250)
  f_svg <- paste0(filename, ".svg")
  print(str_c("Outputting: ", f_svg, " as SVG, ", width, "x", height))
  ggsave(f_svg, plot=plot, width=w*1.25, height=h*1.25, units="px")
}

# Function to get substring up to the n'th occurrence of a character
substr_to_nth_dot <- function(input_string, n) {
  # Split the string by the dot character
  parts <- strsplit(input_string, "\\.")[[1]]
  # Check if there are at least n occurrences
  if(length(parts) < n) {
    stop("There are fewer than n occurrences of '.'.")
  }
  # Concatenate the parts up to the n'th
  result <- paste(parts[1:n], collapse = ".")
  return(result)
}

construct_grp_mask <- function(meta, group_name, col) {
  # split id into subclasses by dot character
  group_name_split <- strsplit(group_name, "\\.")[[1]]
  # For every listed subclass, ensure exact match with sample
  masks <- list()
  for (i in seq_along(group_name_split)) {
    masks[[i]] <- convert_NA_to_false(meta[[paste0(col, ".", i)]] == group_name_split[[i]])
  }
  # Combine masks per subclass
  return (Reduce('&', masks))
}

construct_contrast_string <- function(contrast_vals , seq_type, meta, weight_dict) {
  # construct string for DE expression contrast grouping using IDs
  grp_strings <- list()
  for (i in seq_along(contrast_vals)) {
    split_val <- unlist(strsplit(as.character(contrast_vals[[i]]), "\\."))
    n_splits <- length(split_val)
    val_grps <- list()
    for (j in 1:n_splits) {
      val_grps[j] <- paste0(unlist(split_val)[1:j], collapse = ".")
    }
    key <- contrast_vals[[i]]
    grp_strings[i] <- paste0(weight_dict[[key]], "*group", paste0(val_grps, collapse = "__"), "__", seq_type)
  }
  grp_string <- paste0(grp_strings, collapse=" + ")
  return(paste0("(", grp_string, ")"))
}

eval_contrast <- function(fit, contrast, out_prefix, title, log_file) {
    cleaned_contrast <- sub("^makeContrasts\\(", "", sub(", levels=design\\)$", "", contrast))
    cat("Evaluating contrast ", title, " : ", cleaned_contrast, "\n", file = log_file, append = TRUE)
    # Evaluate the contrast
    lrt <- glmQLFTest(fit, contrast=eval(parse(text = contrast)))
    # Extract results
    res <- topTags(lrt, n=Inf)$table
    # Plot the top differentially translated genes
    res <- res |> tibble::rownames_to_column('gene_id')
    res <- (
        res
        |> arrange(PValue)
        |> mutate(gene_name = feature2name[match(res$gene_id, feature2name[,"id"]), "name"])
        |> select(gene_id, gene_name, everything())
    )
    write.csv(res, paste0(out_prefix, ".csv"), quote=FALSE, row.names=FALSE)
    ymax = max(-log10(res$PValue)[-log10(res$PValue) != Inf]) + 30
    plot <- EnhancedVolcano(res,
        lab = res$gene_name,
        x = "logFC",
        y = "FDR",
        pCutoff = 0.05,
        FCcutoff = 1,
        pointSize = 2.0,
        labSize = 3.5,
        ylim = c(0,ymax),
        title = title 
    )
    gfig(plot, paste0(out_prefix, "_Volcano"), 30,30)
}

parse_groups <- function(groups) {
  # Trim leading and trailing whitespace from the entire string
  groups <- trimws(groups)
  
  if (groups == "NA") {
    return("")
  } else if (grepl(",", groups)) {
    # Split by comma and trim each element of the resulting vector
    return(trimws(unlist(strsplit(groups, ","))))
  } else {
    # Trim and return the single group as a character vector
    return(c(trimws(groups)))
  }
}

convert_NA_to_false <- function(x) {
  replace(x, is.na(x), FALSE)
}

# Function to get valid combinations from a column. Valid combinations include 
# groups within the same supergroup (e.g. MBL.SHH vs MBL.MYC)
# .NA is not a valid subgroup
get_valid_combinations <- function(column_data) {
  data_df <- tibble(original = column_data) %>%
    mutate(
      supergroup = sapply(strsplit(original, "\\."), function(x) x[1]),
      levels_count = sapply(strsplit(original, "\\."), length)
    ) 
  same_levels <- function(a, b) {
    return(length(strsplit(a, "\\.")[[1]]) == length(strsplit(b, "\\.")[[1]]))
  }
  valid_combinations <- (
    data_df 
    |> group_by(supergroup)
    |> filter(n() > 1)
    |> distinct(original)
    |> reframe(combinations = list({
      # Only generate combinations if there are 2 or more distinct elements
      if (n() >= 2) {
        combs <- combn(original, 2, simplify = FALSE)
        Filter(function(x) same_levels(x[1], x[2]), combs)
      } else {
        list()  # Return an empty list if there are not enough elements
      }
    }))
    |> pull(combinations)
    |> unlist(recursive = FALSE)
  )
  # Additionally, check for direct combinations of unique entries
  unique_entries <- unique(data_df$original)
  if (length(unique_entries) > 1 && all(data_df$levels_count == 1)) {
    direct_combinations <- combn(unique_entries, 2, simplify = FALSE)
  } else {
    direct_combinations <- list()
  }
  all_combinations <- c(valid_combinations, direct_combinations)
  # Return unique combinations sorted
  return(unique(lapply(all_combinations, function(x) sort(x))))
}



option_list <- list(
    make_option(c("-m", "--metadata"  ), type="character", metavar="path"   , help="project sample sheet containing sample groups and ids."                        ),
    make_option(c("-i", "--rna_counts"    ), type="character", default=NULL                 , metavar="path"   , help="Count file matrix where rows are genes and columns are samples."                        ),
    make_option(c("-j", "--ribo_counts"   ), type="character", default=NULL                 , metavar="path"   , help="Count file matrix where rows are genes and columns are samples."                        ),
    make_option(c("-c", "--count_col"     ), type="integer"  , default=2                    , metavar="integer", help="First column containing sample count data."                                             ),
    make_option(c("-t", "--tx_table"      ), type="character", default=NULL                   , metavar="path",    help="Table linking transcript and gene IDs/names"                                          ),
    make_option(c("-f", "--tx_table_col"  ), type="character",  default="gene_id"        , metavar="string" , help="column name of tx_table to use as keys"                                                    ),
    make_option(c("-r", "--orf_tx_table"  ), type="character", default=NULL                   , metavar="path",    help="Table linking ORF and transcript IDs"                                          ),
    make_option(c("-d", "--id_col"        ), type="integer"  , default=1                    , metavar="integer", help="Column containing identifiers to be used."                                              ),
    make_option(c("-s", "--sep"           ), type="character", default=','                  , metavar="string" , help="Separator of input table."                                              ),
    make_option(c("-o", "--outdir"        ), type="character", default='out', metavar="path"   , help="Output directory."                                                                      ),
    make_option(c("-u", "--outer_join"    ), type="logical"  , default=FALSE                , metavar="boolean", help="Outer join RNA and Ribo reads, fill NAs with 0 expression."                                                     ),
    make_option(c("-l", "--cores"         ), type="integer"  , default=1                    , metavar="integer", help="Number of cores."                                                                       ),
    make_option(c("-a", "--contrast_cols" ), type="character", default="treatment_id"               , metavar="character", help="Column names from which contrasts are derived; separated by commas."          ),
    make_option(c("-e", "--no_batch_factor"), type="logical" , default=FALSE,               , metavar="boolean", help="Don't create factors for batch date. Can be necessary to achieve full rank for some settings")
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

### DEBUGING ### 
# opt$metadata <- "sample_sheet.csv"
# opt$rna_counts <- "out/rna/quants/salmon.STD_TS.NA/no_agg_NumReads_matrix.csv"
# opt$ribo_counts <- "out/ribo/quants/salmon.STD_TS.STD_TL/no_agg_NumReads_matrix.csv"
# opt$count_col <- 5
# opt$sep <- ","
# opt$outdir <- "test/"
# opt$tx_table <- "out/targets/STD_TS/STD_TL/tx_table.csv"
# opt$contrast_cols <- "disease_id"
# opt$cores <- 4
# opt$no_batch_factor <- TRUE

# Multi-core processing
register(MulticoreParam(opt$cores))
# Parsing multi-group contrasts
contrast_grps <- parse_groups(opt$contrast_cols)
contrast_col <- paste0(contrast_grps, collapse="__")
# Define the path to the log file
log_file <- paste0(opt$outdir, "run_info.txt")
# Clear existing content of the log file at the start of the script
# So the log file is fresh for each run
file.create(log_file) 

meta.table <- read_csv(
  opt$metadata,
  comment = "#",
  col_types = cols(
    replicate_num = col_character()
  )
) |> mutate(counts_col = ifelse(data_type %in% c("RNA_seq", "ONTRNA_seq"), paste0(smart_id, "_RNA"), paste0(smart_id, "_Ribo")))

if (!is.null(opt$ribo_counts)) {
  ribo_counts <- (
      read.csv(
          opt$ribo_counts,
          sep = opt$sep,
          check.names = FALSE
      ) %>%
      # Select the 'id_col' as row names
      tibble::column_to_rownames(var = names(read.csv(opt$ribo_counts, sep = opt$sep, nrows = 1))[opt$id_col]) %>%
      # Select only columns from `count_col` onwards, -1 to exclude the id_col
      dplyr::select((opt$count_col - 1):ncol(.)) %>%
      # Mutate across all columns to convert to integer
      mutate(across(
          .cols = everything(),
          .fns = ~ as.integer(.x)
      )) %>%
      # Rename columns to append '_Ribo'
      rename_with(.fn = function(s) { paste(s, "_Ribo", sep = "") })
  )
  # report the colnames of ribo_counts
  print(str_c("Ribo counts columns: ", paste(colnames(ribo_counts), collapse = ", ")))
  # Log the number of rows and columns in ribo_counts
  cat("Ribo counts matrix dimensions: ", nrow(ribo_counts), " rows, ", ncol(ribo_counts), " columns\n", file = log_file, append = TRUE)
  ribo_ids <- rownames(ribo_counts)
}
if (!is.null(opt$rna_counts)) {
  rna_counts <- read.csv(
      opt$rna_counts,
      sep=opt$sep,
      check.names=FALSE,
      ) %>%
      # Select the 'id_col' as row names
      tibble::column_to_rownames(var = names(read.csv(opt$rna_counts, sep = opt$sep, nrows = 1))[opt$id_col]) %>%
      # Select only columns from `count_col` onwards, -1 to exclude the id_col
      dplyr::select((opt$count_col - 1):ncol(.)) %>%
      # Mutate across all columns to convert to integer
      mutate(across(
          .cols = everything(),
          .fns = ~ as.integer(.x)
      )) %>%
      rename_with(.fn = function(s){paste(s, "_RNA", sep="")})
  # report the colnames of ribo_counts
  print(str_c("RNA counts columns: ", paste(colnames(rna_counts), collapse = ", ")))
  # Log the number of rows and columns in rna_counts
  cat("RNA counts matrix dimensions: ", nrow(rna_counts), " rows, ", ncol(rna_counts), " columns\n", file = log_file, append = TRUE)
  rna_ids <- rownames(rna_counts)
}

if (!is.null(opt$tx_table)) {
  # Get TX table to parse gene names
  tx.table <- read.csv(opt$tx_table)
  # Create dictionary that connects counts row id with <GENENAME:TRANSCRIPT_CDSTIS>
  feature2name <- (
      tx.table 
      |> select(opt$tx_table_col, "gene_name")
      |> distinct()
      |> dplyr::rename(id = opt$tx_table_col, name=gene_name)
      %>% {
          if (opt$tx_table_col == "transcript_id") 
              mutate(., name = paste(name, id, sep=":"))
          else .
      }
      |> mutate(
            name = if_else(is.na(name) | name == "" | is.nan(name), id, name),
          )
  )
} else {
  # Create dummy feature2name table from unique rownames of both ribo_counts and rna_counts
  if (!is.null(opt$rna_counts) && !is.null(opt$ribo_counts)) {
    feature2name <- data.frame(
        id = unique(c(rownames(ribo_counts), rownames(rna_counts))),
        name = unique(c(rownames(ribo_counts), rownames(rna_counts)))
    )
  } else if (!is.null(opt$rna_counts)) {
    feature2name <- data.frame(
        id = rownames(rna_counts),
        name = rownames(rna_counts)
    )
  } else {
    feature2name <- data.frame(
        id = rownames(ribo_counts),
        name = rownames(ribo_counts)
    )
  }
}

# Create CDS - transcript map
if (is.null(opt$orf_tx_table) || opt$tx_table_col == "gene_id") {
  ribo_rna_map <- data.frame(
    ORF_id = rownames(ribo_counts),
    transcript_id = rownames(ribo_counts)
  )
} else {
  ribo_rna_map <- read.csv(opt$orf_tx_table, sep = ",", header = TRUE)
  print(colnames(ribo_rna_map))
}
ribo_rna_map <- ribo_rna_map %>% tibble::column_to_rownames("ORF_id")

# Create dirs
dir.create(paste0(opt$outdir, 'dTE'), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(opt$outdir, 'Ribo'), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(opt$outdir, 'RNA'), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(opt$outdir, 'TE'), showWarnings = FALSE, recursive = TRUE)

# Intersect RNA and Ribo counts -------------------------------------------------
## canonical sequences

if (!is.null(opt$rna_counts) && !is.null(opt$ribo_counts)) {
  if (opt$outer_join) {
    print("FUNCTION IS NOT YET PROPERLY IMPLEMENTED")
    # Outer join: union of all rownames, fill missing with 0
    all_genes <- union(rownames(rna_counts), ribo_rna_map[rownames(ribo_counts),])
    ribo_counts_select <- ribo_counts[all_genes, , drop=FALSE]
    rna_counts_select <- rna_counts[all_genes, , drop=FALSE]
    # Fill NAs with 0
    ribo_counts_select[is.na(ribo_counts_select)] <- 0
    rna_counts_select[is.na(rna_counts_select)] <- 0
  } else {
    # Intersect (default)
    # Find transcripts and ORF counts with overlapping Names
    ribo_counts_select <- ribo_counts
    rna_counts_select <- rna_counts[ribo_rna_map[rownames(ribo_counts),],]
  }
  counts <- cbind(RNA=as.matrix(rna_counts_select), Ribo=as.matrix(ribo_counts_select))
  seq_types <- c("RNA", "Ribo")
  out_dirs <- c(dirname(opt$rna_counts), dirname(opt$ribo_counts))
} else if (!is.null(opt$rna_counts)) {
  counts <- rna_counts
  meta.table <- meta.table
  seq_types <- c("RNA")
  out_dirs <- c(dirname(opt$rna_counts))
} else {
  counts <- ribo_counts
  meta.table <- meta.table
  seq_types <- c("Ribo")
  out_dirs <- c(dirname(opt$ribo_counts))
}
######## DEBUG
# counts <- counts[1:2000,]
########

# Get rid of NA (can use more investigation how these appear)
mask <- apply(is.na(counts), 1, any)
if (sum(mask) > 0) {
  print(str_c("detected ", sum(mask), " NA rows"))
}
counts <- counts[mask == FALSE,]

# Import and process metadata --------------------------------------------------------------
meta.table <- (
  meta.table
  |> group_by(smart_id, data_type)
  |> mutate(
    seq_type = recode(data_type, "RNA_seq" = "RNA", "ONTRNA_seq" = "RNA", "Ribo_seq" = "Ribo"),
    across(everything(), ~ replace_na(as.character(.), "NA")),
  )
  |> slice_sample(n = 1)
  |> ungroup()
)
# Remove Subgroup info of samples with only a single instance of that contrast (stratified over R-seq data)
# Iterate samples from lowest level to highest
for (contrast in contrast_grps) {
  # Count the number of '.' in each string
  dot_counts <- sapply(meta.table[[contrast]], function(x) {
    count <- gregexpr("\\.", x)[[1]]
    
    # Check for no matches, where gregexpr returns -1
    ifelse(count[1] == -1, 0, length(count))
  })
  # Sort entries based on the number of '.' in descending order
  sorted_entries <- unique(meta.table[order(dot_counts, decreasing = TRUE),][[contrast]])
  for (entry in sorted_entries) {
    for (type in unique(meta.table$seq_type)) {
        mask <- (startsWith(meta.table[[contrast]], entry)) & (meta.table$seq_type %in% type)
        smart_ids <- unique(meta.table[mask,]$smart_id)
        if ((length(smart_ids) > 0) & (length(smart_ids) < 2)) {
          # If the entry contains a '.', rename the group to its supergroup part
          if (grepl("\\.", entry)) {
            # Logic to rename the group to its supergroup
            new_group <- sub("\\.[^.]*$", "", entry)  # Remove the last segment after the last dot
            meta.table[mask, contrast] <- new_group
            cat("Renaming group:", entry, "to supergroup:", new_group, "because of low sample count\n", file = log_file, append = TRUE)
          } else {
            # Otherwise, remove samples
            meta.table <- meta.table[!mask, ]
            # Identify the indices of the columns to remove
            indices_to_remove <- which(colnames(counts) %in% paste0(smart_ids, "_", type))
            counts <- counts[, -indices_to_remove]
            cat("Removing samples with group:", entry, ": ", smart_ids, "because of low sample count\n", file = log_file, append = TRUE)
          }
        }
    }
  }
}

# Dictionary to store combinations
comb_dict <- list()
uniq_dict <- list()
num_col_splits <- list()

meta_cols <- unique(c("location_id", "disease_id", "treatment_id", contrast_grps))
# Combine contrast group columns
if (length(contrast_grps) > 1) {
  meta.table[[contrast_col]] <- do.call(paste, c(meta.table[contrast_grps], sep=":"))
  meta_cols <- unique(c(meta_cols, contrast_grps, contrast_col))
} 
# Split columns into subclasses
for (column in meta_cols) {
  # For replace all : with . (only possible for treatment_id's)
  # TODO figure out how to use combinatorial effects more smartly
  if (!is.null(meta.table[[column]]) && length(meta.table[[column]]) > 0) {
    meta.table[[column]] <- gsub(":", "_and_", meta.table[[column]])
  } else {
    warning(sprintf("Column '%s' is missing or empty in meta.table. Skipping processing for this column.", column))
    next
  }
  # Get max splits
  max_splits <- max(str_count(meta.table[[column]], "\\.") + 1)
  num_col_splits[[column]] <- max_splits
  column_names <- paste0(column, ".", seq_len(max_splits))
  # Create new columns with increasing length of the original string
  for (i in seq_len(max_splits)) {
    meta.table[[column_names[i]]] <- sapply(strsplit(meta.table[[column]], "\\."), function(x) {
      paste(x[1:i], collapse = ".")
    })
  }
  # For the original column, ensure all elements have same number of subclasses by adding as many .NA
  meta.table[[column]] <- sapply(strsplit(meta.table[[column]], "\\."), function(x) {
    paste(c(x, rep("NA", max_splits - length(x))), collapse = ".")
  })
  # Get and store valid combinations
  combs <- list()
  for (i in seq_len(max_splits)) {
    vals <- unique(meta.table[[column_names[i]]])
    # Get rid of entries with ".NA"
    vals <- vals[!grepl("\\.NA", vals)]
    if (length(vals) > 1) {
      groups <- combn(vals, 2, simplify)
      # eval orthogonal combs (e.g. 1_and_a vs 2_and_b -> 1 vs 2, a vs. b)
      if (any(grep("_and_", meta.table[[column]]))) {
        # TODO hardcoded for a max of two combinatories factors, not for more
        orth_vals_1 <- unique(sapply(strsplit(meta.table[[column]], "_and_"), `[`, 1))
        orth_vals_2 <- unique(sapply(strsplit(meta.table[[column]], "_and_"), `[`, 2))
        orth_groups_1 <- combn(orth_vals_1, 2, simplify)
        orth_groups_2 <- combn(orth_vals_2, 2, simplify)
        groups <- cbind(groups, orth_groups_1, orth_groups_2)
      }
      for (j in 1:length(groups[1,])) {
        if (i > 1) {
          if (substr_to_nth_dot(groups[,j][1], i-1) == substr_to_nth_dot(groups[,j][2], i-1)) {
            combs <- append(combs, list(groups[,j]))
          }
        } else {
          combs <- append(combs, list(groups[,j]))
        }
      }
    }
  }
  comb_dict[[column]] <- combs
  uniq_dict[[column]] <- unique(unlist(combs))
  # Print warning if no combinations found
  # if (length(combs) == 0) {
  #   warning(sprintf("No valid combinations found for suggested column name '%s'.", column))
  # }
}
# Print the combinations and unique entries
cat("Combinations and unique entries for each column:\n", file = log_file, append = TRUE)
for (col in names(comb_dict)) {
  cat(str_c(col, ":\n"), file = log_file, append = TRUE)
  cat(paste0(" - ", paste(comb_dict[[col]], collapse = "\n - ")), "\n", file = log_file, append = TRUE)
  cat(str_c("Unique entries: ", paste(uniq_dict[[col]], collapse = ", ")), "\n", file = log_file, append = TRUE)
}
# TODO: validate presence of both data types
# Report on sample_ids filtered out because of missing data
# Identify filtered out sample_ids due to missing count data
filtered_sample_ids <- paste0(meta.table$smart_id, "_", meta.table$seq_type)
filtered_sample_ids <- filtered_sample_ids[!filtered_sample_ids %in% colnames(counts)]
if (length(filtered_sample_ids) > 0) {
  cat("Filtered out sample_ids due to missing count data:\n", file = log_file, append = TRUE)
  for (sid in filtered_sample_ids) {
    cat(" - ", sid, "\n", file = log_file, append = TRUE)
  }
}

# Generate column names according to the maximum number of splits
if (!"sample_type" %in% colnames(meta.table)) {
  meta.table$sample_type <- "NA"
}
if (!"batch_date" %in% colnames(meta.table)) {
  meta.table$batch_date <- "NA"
}
meta.samples <- (
  meta.table
  |> filter(
    counts_col %in% colnames(counts),
    )
  |> select(
    smart_id,
    sample_type,
    any_of(unlist(lapply(meta_cols, function(col) grep(paste0("^", col), colnames(meta.table), value = TRUE)))),
    rep = replicate_num,
    control = test_or_control,
    batch_date,
    counts_col,
    seq_type
  )
  |> mutate(
    across(where(is.character), as.factor),
    batch_date=as.factor(batch_date)
  )
  |> replace_na(list(rep = 1))
  |> unique()
)

# Create design formula and sample groups  ----------------------------
f = ("~0 + group")
# include non-unique contrast col names
contrast_cols <- list()
# Create all possible combinations of merged contrast groups
temp <- unlist(lapply(contrast_col, function(col) grep(paste0("^", col, "."), colnames(meta.samples), value = TRUE)))
# Check with meta.samples wether column exists
for (col in temp) {
  if (length(unique(meta.samples[[col]])) > 1) {
    contrast_cols <- append(contrast_cols, col[[1]])
  }
}
contrast_cols <- unlist(contrast_cols)
select_cols <- c("counts_col", contrast_cols, "seq_type")
# evaluate location_id, source_id, disease_id, treatment_id
common_cols <- c("source_id", "location_id", "disease_id", "treatment_id")
# exclude contrast_cols
common_cols <- common_cols[!common_cols %in% contrast_cols]
# include to selected columns if more than one unique value
for (col in common_cols) {
  if (col %in% colnames(meta.samples)) {
    if (length(unique(meta.samples[[col]])) > 1) {
      select_cols <- c(select_cols, col)
    }
  }
}


# ensure batch_date and seq_type are not perfectly correlated
corr_factor <- suppressWarnings({
  cf <- cor(as.numeric(meta.samples$seq_type), as.numeric(meta.samples$batch_date))
  if (is.na(cf)) 0 else cf
})
if (!opt$no_batch_factor && length(unique(meta.samples$batch_date)) > 1 && abs(corr_factor) < 0.99) {
  # Replace values that appear only once with "NaN"
  # batch_counts <- table(meta.samples$batch_date)
  # single_batches <- names(batch_counts[batch_counts == 1])
  # meta.samples$batch_date <- as.character(meta.samples$batch_date)
  # meta.samples$batch_date[meta.samples$batch_date %in% single_batches] <- "NaN"
  # meta.samples$batch_date <- as.factor(meta.samples$batch_date)
  select_cols <- c(select_cols, "batch_date")
  f = paste0(f, " + batch_date")
}
if (length(unique(meta.samples$sample_type)) > 1) {
  select_cols <- c(select_cols, "sample_type")
  f = paste0(f, " + sample_type")
}
# Select relevant columns 
meta.design <- (
  meta.samples 
  |> select(any_of(unlist(select_cols)))
  |> mutate(group = paste0(apply(meta.samples[, unlist(contrast_cols)], 1, paste, collapse = "__"), "__", seq_type))
  |> select(-all_of(unlist(contrast_cols)))
  |> column_to_rownames("counts_col")
)

# Rearrage counts according to meta.samples file
counts <- counts[,as.character(meta.samples$counts_col)]

# Combine contrast_grps and seq_type
design <- model.matrix(as.formula(f), data=data.frame(meta.design))

# write formula to text file
write(f, file = str_c(opt$outdir, "formula.txt"))
cat("\ngroup names:\n", colnames(design), "\n", file=log_file, append=TRUE)
cat("\nformula:\n", f, "\n", file=log_file, append=TRUE)

############### DEBUG
# counts_test <- counts[,]
# meta_test <- meta.design
# meta_test$batch_date <- meta_test$batch

# design <- model.marix(as.formula(f), data=data.frame(meta_test))
# # Compute the correlation matrix
# corr_matrix <- cor(design)
# # Find columns with high correlation
# corr_indices <- which(abs(corr_matrix) > 0.9, arr.ind = TRUE)
# high_corr_indices <- corr_indices[corr_indices[, "row"] < corr_indices[, "col"], ]
# print(high_corr_indices)

# fig(str_c('/home/clauwaer/data/data_processing/VanHeesch_Collab/debug2'), 20, 20)
# # Create the heatmap
# pheatmap(corr_matrix, cluster_rows = FALSE, cluster_cols = FALSE,
#          main = "Design Matrix Heatmap")
# w()
# fig(str_c('/home/clauwaer/data/data_processing/VanHeesch_Collab/debug'), 10, 10)
# # Create the heatmap
# pheatmap(design, cluster_rows = FALSE, cluster_cols = FALSE,
#          main = "Design Matrix Heatmap")
# w()
# # #####################

# create DGE object
dge <- DGEList(counts=counts, samples=data.frame(meta.design))
# Filter lowly expressed genes
keep <- filterByExpr(dge)
dge <- dge[keep, keep.lib.sizes = FALSE]
# Calculate normalization factors
dge <- calcNormFactors(dge, method="TMM")
# Estimate dispersion
dge <- estimateDisp(dge, design, min.row.sum=30)
# split normalized count data according to RNA and Ribo samples:
norm_counts_list <- list(RNA=dge[,meta.samples$seq_type == "RNA"], Ribo=dge[,meta.samples$seq_type == "Ribo"])
# Save normalized count data separate for RNA and Ribo samples within the input directory of both, respectively
for (i in seq_along(seq_types)) {
  seq_type <- seq_types[i]
  out_dir <- out_dirs[i]
  if (opt$tx_table_col == "transcript_id") {
    prefix <- "no_agg" 
  } else {
    prefix <- "gene_agg"
  }
  logcpm_counts <- cpm(norm_counts_list[[seq_type]], normalized.lib.sizes = TRUE, log=TRUE, prior.count=1)
  formatted_data <- rownames_to_column(data.frame(format(logcpm_counts, digits=3, nsmall = 3)), "Name")
  write.table(formatted_data, file=str_c(out_dir, "/", prefix, "_cpm_log_matrix.csv"), sep=",", row.names = FALSE, quote=FALSE)
}
# Fit the GLM model
fit <- glmQLFit(dge, design)

dge_idxs <- match(rownames(dge$samples), meta.samples$counts_col)
dge.meta <- meta.samples[dge_idxs,]

# Subset samples based on condition
if (length(seq_types) == 2) {
  seq_labels <- c("All", "Ribo", "RNA")
  seq_groups  <- list(c("RNA", "Ribo"), c("Ribo"), c("RNA"))
} else {
  seq_labels <- seq_types
  seq_groups <- seq_types
}

# Evaluate design to get metadata of interest
meta_of_interest <- colnames(meta.design)[!colnames(meta.design) %in% c("group", "seq_type")]
meta_of_interest <- c(contrast_grps, meta_of_interest)
for (i in seq_along(seq_labels)) {
  sample_mask <- dge.meta$seq_type %in% seq_groups[[i]]
  samples_to_keep <- rownames(dge$samples[sample_mask,])
  eval_MDS(dge[,samples_to_keep], dge.meta[sample_mask,], meta_of_interest, opt$outdir, seq_labels[i])
  eval_PCA(dge[,samples_to_keep], dge.meta[sample_mask,], meta_of_interest, opt$outdir, seq_labels[i])
  eval_heatmap(cpm(dge[,samples_to_keep], normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1), dge.meta[sample_mask,], meta_of_interest, opt$outdir, seq_labels[i])
  eval_gene_clusters(cpm(dge[,samples_to_keep], normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1), dge.meta[sample_mask,], meta_of_interest, opt$outdir, seq_labels[i])
}

bplapply(uniq_dict[[contrast_col]], function(val) {evaluate_unique_contrast(dge.meta, val, contrast_col, fit, opt$outdir, log_file)})
bplapply(comb_dict[[contrast_col]], function(val) {evaluate_combination_contrast(dge.meta, val, contrast_col, fit, opt$outdir, log_file)})
