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
  library(DEFormats)
  library(BiocParallel)
  if (requireNamespace("svglite", quietly = TRUE)) {
    library(svglite)
  }
})
options(show.error.locations = TRUE)

validate_inputs <- function(opt) {
  cat("Validating inputs...\n")
  errors <- c()
  warnings <- c()
  
  # 1. Check File Existence
  if (!is.null(opt$rna_counts) && !file.exists(opt$rna_counts)) errors <- c(errors, paste("RNA counts file not found:", opt$rna_counts))
  if (!is.null(opt$ribo_counts) && !file.exists(opt$ribo_counts)) errors <- c(errors, paste("Ribo counts file not found:", opt$ribo_counts))
  if (!is.null(opt$metadata) && !file.exists(opt$metadata)) errors <- c(errors, paste("Metadata file not found:", opt$metadata))
  if (!is.null(opt$tx_table_path) && !file.exists(opt$tx_table_path)) errors <- c(errors, paste("Tx Table file not found:", opt$tx_table_path))
  
  # 2. Check Logical Requirements
  if (is.null(opt$rna_counts) && is.null(opt$ribo_counts)) {
    errors <- c(errors, "You must provide at least one count matrix (--rna_counts or --ribo_counts).")
  }
  
  # 3. Validate Tx Table Columns
  if (!is.null(opt$tx_table_path) && file.exists(opt$tx_table_path)) {
    # Read just the header to check columns
    tx_header <- colnames(read.csv(opt$tx_table_path, nrows=1))
    
    # A. Check for gene_name (Recommended for plotting)
    if (!"gene_name" %in% tx_header) {
      warnings <- c(warnings, "Column 'gene_name' missing in tx_table. IDs will be used as names in plots.")
    }
    
    # B. TRANSCRIPT MODE REQUIREMENTS (Hardcoded columns)
    if (opt$feature_level == "transcript") {
        missing_cols <- setdiff(c("translon_id", "transcript_id"), tx_header)
        if (length(missing_cols) > 0) {
            errors <- c(errors, paste0("Transcript-level analysis (--feature_level transcript) requires specific columns in tx_table: ", 
                                       paste(missing_cols, collapse=", ")))
        }
    }
    
    # C. GENE MODE REQUIREMENTS (User defined col)
    if (opt$feature_level == "gene") {
        if (!opt$tx_table_col %in% tx_header) {
            errors <- c(errors, paste0("The specified ID column (--tx_table_col='", opt$tx_table_col, "') is not in the tx_table."))
        }
    }
  } else if (opt$feature_level == "transcript" && !is.null(opt$ribo_counts)) {
     # Transcript mode requested but no table provided
     errors <- c(errors, "Transcript-level dTE analysis requires --tx_table_path to map Translons.")
  }
  
  # Print Feedback
  if (length(warnings) > 0) {
    cat("WARNINGS:\n")
    for (w in warnings) cat(paste("  -", w, "\n"))
  }
  
  if (length(errors) > 0) {
    cat("CRITICAL ERRORS:\n")
    for (e in errors) cat(paste("  -", e, "\n"))
    stop("Input validation failed. Please check arguments.")
  }
  cat("Inputs validated successfully.\n")
}

get_hierarchical_cols <- function(levels) {
  # 1. Parse Supergroups
  levels <- as.character(levels)
  # Handle cases with no dots gracefully
  supergroups <- sapply(levels, function(x) {
      parts <- strsplit(x, "\\.")[[1]]
      if(length(parts) > 0) parts[1] else x
  })
  
  unique_super <- sort(unique(supergroups))
  n_super <- length(unique_super)
  
  # 2. Assign distinct Base Colors to Supergroups (Using Set1 or Dark2 for contrast)
  
  # --- CASE A: SINGLE SUPERGROUP (Full Spectrum for Subgroups) ---
  if (n_super == 1) {
    # If all groups belong to the same supergroup (e.g. brain.ctx, brain.hippo),
    # we treat them as distinct categorical variables using a full palette.
    subs <- sort(unique(levels))
    
    # Use a high-contrast palette (Set1 or Paired depending on N)
    pal_name <- if (length(subs) <= 8) "Set1" else "Paired"
    # Ensure palette is interpolated if n > palette_size
    final_cols <- colorRampPalette(brewer.pal(min(length(subs), 8), pal_name))(length(subs))
    names(final_cols) <- subs
    
    return(final_cols)
  }
  
  # --- CASE B: MULITPLE SUPERGROUPS (Hierarchical Gradients) ---
  base_palette <- colorRampPalette(brewer.pal(min(n_super, 8), "Set1"))(n_super)
  names(base_palette) <- unique_super
  
  final_cols <- c()
  
  # 3. Generate Base -> Light Gradients
  for (sup in unique_super) {
    subs <- sort(unique(levels[supergroups == sup]))
    n_sub <- length(subs)
    
    if (n_sub == 1) {
      final_cols[subs] <- base_palette[sup]
    } else {
      # Gradient: Base Color -> Lighter Version (mixed with 70% white)
      # We do NOT go all the way to white to keep it visible
      base_col <- base_palette[sup]
      light_col <- colorRampPalette(c(base_col, "white"))(10)[7] # 70% white mix
      
      # Generate ramp
      cols <- colorRampPalette(c(base_col, light_col))(n_sub)
      names(cols) <- subs
      final_cols <- c(final_cols, cols)
    }
  }
  return(final_cols)
}

# MDS -------------------------------------------------------------------------
eval_MDS <- function(dge.set, meta, cols_of_interest, prefix, suffix, plot_ids=FALSE) {
  xy_dim <- max(dim(meta)[1]*0.10, 20)
  # MDS
  pcaData <- plotMDS(dge.set, ntop=Inf, returnData = TRUE)
  percentVar <- round(100 * pcaData$var.explained)
  
  plot_data <- data.frame(
    x = pcaData$x,
    y = pcaData$y,
    rep = meta$rep,
    smart_id = meta$smart_id
  )
  
  for (col in cols_of_interest) {
    plot_data[[col]] <- meta[[col]]
  }
  
  # Select Label Column
  label_col <- if (plot_ids) "smart_id" else "rep"

  for (i in seq_along(cols_of_interest)) {
    col <- cols_of_interest[[i]]
    plot <- ggplot(plot_data, aes(x = x, y = y, color=meta[[col]], shape=meta$rep)) +
      geom_point(size = 3, alpha = 0.8) +
      coord_fixed() + 
      geom_text(aes(label = .data[[label_col]]), size=4, color = "black", vjust=1.5) +
      xlab(paste0("PC1: ", percentVar[1], "% variance explained")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance explained")) +
      ggtitle(paste0("MDS on fitted ", suffix, " data (all counts/edgeR)")) +
      scale_color_manual(values = get_hierarchical_cols(unique(as.character(meta[[col]])))) +
      labs(color = col, shape="replicate")
    gfig(plot, str_c(prefix, "MDS_", i, "_", suffix), xy_dim, xy_dim)
    }
}

# PCA -------------------------------------------------------------------------
eval_PCA <- function(dge.set, meta, cols_of_interest, prefix, suffix, plot_ids=FALSE) {
  xy_dim <- max(dim(meta)[1]*0.10, 20)
  # PCA
  pcaData <- plotMDS(dge.set, ntop=Inf, returnData = TRUE, gene.selection="common")
  percentVar <- round(100 * pcaData$var.explained)
  
  plot_data <- data.frame(
    x = pcaData$x,
    y = pcaData$y,
    rep = meta$rep,
    smart_id = meta$smart_id
  )
  
  for (col in cols_of_interest) {
    plot_data[[col]] <- meta[[col]]
  }
  
  label_col <- if (plot_ids) "smart_id" else "rep"

  for (i in seq_along(cols_of_interest)) {
    col <- cols_of_interest[[i]]
    plot <- ggplot(plot_data, aes(x = x, y = y, color=meta[[col]], shape=meta$rep)) +
      geom_point(size = 3, alpha = 0.8) +
      coord_fixed() + 
      geom_text(aes(label = .data[[label_col]]), size=4, color = "black", vjust=1.5) +
      xlab(paste0("PC1: ", percentVar[1], "% variance explained")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance explained")) +
      ggtitle(paste0(suffix, ": PCA on fitted data (all counts/edgeR)")) +
      scale_color_manual(values = get_hierarchical_cols(unique(as.character(meta[[col]])))) +
      labs(color = col, shape="replicate")
    gfig(plot, str_c(prefix, "PCA_", i, "_", suffix), xy_dim, xy_dim)
  }
}

# Distance Heatmap ------------------------------------------------------------
eval_heatmap <- function(norm_counts, meta, meta_cols, prefix, suffix) {
  xy_dim <- max(dim(meta)[1]*0.15, 5)
  # Calculate Distance
  sampleDists <- dist(t(norm_counts))
  sampleDistMatrix <- as.matrix(sampleDists)
  colnames(sampleDistMatrix) <- colnames(norm_counts)
  rownames(sampleDistMatrix) <- apply(meta[, meta_cols, drop=FALSE], 1, paste, collapse = "__")
  anno_df <- as.data.frame(meta[, meta_cols, drop=FALSE])
  rownames(anno_df) <- colnames(sampleDistMatrix)

  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  for (fmt in c("png", "svg")) {
    if (fmt == "png") {
      pngfig(str_c(prefix, "Heatmap_1_", suffix), xy_dim, xy_dim)
      pheatmap(sampleDistMatrix,
        clustering_distance_rows = sampleDists,
        clustering_distance_cols = sampleDists,
        main = paste0(suffix, ": Euclidean distance between counts of samples on fitted data.")
      )
      w()
    } else {
      svgfig(str_c(prefix, "Heatmap_1_", suffix), xy_dim, xy_dim)
      pheatmap(sampleDistMatrix,
        clustering_distance_rows = sampleDists,
        clustering_distance_cols = sampleDists,
        main = paste0(suffix, ": Euclidean distance between counts of samples on fitted data.")
      )
      w()
    }
  } 
  if (!is.null(colnames(sampleDistMatrix)) && !is.null(rownames(sampleDistMatrix))) {
      sorted_col_order <- order(colnames(sampleDistMatrix))
      sorted_row_order <- order(rownames(sampleDistMatrix))
      
      sorted_sampleDistMatrix <- sampleDistMatrix[sorted_row_order, sorted_col_order]
      
      write.table(
        cbind(sample = rownames(sorted_sampleDistMatrix), sorted_sampleDistMatrix),
        file=str_c(prefix, "sample.dists_", paste0(suffix, ".txt")),
        row.names=FALSE,
        col.names=TRUE,
        sep="\t",
        quote=FALSE
      )
  }
}


eval_gene_clusters <- function(norm_counts, meta, meta_cols, prefix, suffix) {
  x_dim <- max(dim(meta)[1]*0.18, 7)
  k_genes <- max(dim(meta)[1]*0.7, 50)
  y_dim <- k_genes * 0.21
  for (i in 1:1) {
    topVarGenes <- head(order(rowVars(norm_counts), decreasing = TRUE), k_genes)
    mat <- norm_counts[topVarGenes, ]
    mat <- mat - rowMeans(mat)
    
    # 1. Label Mapping with Fallback
    new_names <- feature2name[match(rownames(mat), feature2name[,"id"]), "name"]
    rownames(mat) <- ifelse(is.na(new_names), rownames(mat), new_names)
    # 2. Annotation Setup (Prevent Duplicates)
    cols_to_select <- unique(c(meta_cols, "rep", "batch_date", "sample_type"))
    # Ensure they actually exist in the metadata object
    cols_to_select <- cols_to_select[cols_to_select %in% colnames(meta)]
    
    anno <- as.data.frame(meta[, cols_to_select, drop=FALSE])
    
    if (i == 1) {
      rownames(anno) <- colnames(mat)
    }
    
    # 3. Custom Colors for Heatmap Annotations
    anno_colors <- list()
    for (col in cols_to_select) {
         # Only generate colors for categorical/character columns
         if(is.character(meta[[col]]) || is.factor(meta[[col]])) {
             vals <- unique(as.character(meta[[col]]))
             n_vals <- length(vals)
             # --- COLOR LOGIC SPLIT ---
             if (col %in% c("batch_date", "sample_type", "rep", "replicate_num")) {
                 # A. STANDARD PALETTE (Flat Variables)
                 std_pal <- colorRampPalette(brewer.pal(min(n_vals, 8), "Set2"))(n_vals)
                 names(std_pal) <- vals
                 anno_colors[[col]] <- std_pal
             } else {
                 # B. HIERARCHICAL GRADIENT (Grouping Variables)
                 anno_colors[[col]] <- get_hierarchical_cols(vals)
             }
         }
    }

    for (fmt in c("png", "svg")) {
      if (fmt == "png") {
        pngfig(str_c(prefix, "GeneClusters_1_", suffix), x_dim, y_dim)
      } else {
        svgfig(str_c(prefix, "GeneClusters_1_", suffix), x_dim, y_dim)
      }
      pheatmap(mat, 
        annotation_col = anno, 
        annotation_colors = anno_colors, 
        main=suffix
      )
      w()
    }
  }
}

evaluate_combination_contrast <- function(meta, tup, contrast_col, fit_paired, fit_rna, outdir, log_file) {
  rna_mask <- meta$seq_type == "RNA"
  ribo_mask <- meta$seq_type == "Ribo"
  
  # NEW: Detect Global Mode (Dual vs RNA-only)
  has_ribo_data <- any(ribo_mask)
  
  grp_rna <- list(); grp_ribo <- list(); n_rna <- list(); n_ribo <- list(); order <- list()

  for (i in seq_along(tup)) {
    num_dots <- str_count(tup[[i]], "\\.") + 1
    tmp_col <- paste0(contrast_col, ".", num_dots)
    if (!grepl("_and_", tup[[i]]) && any(grepl("_and_", meta[[tmp_col]]))) {
      split_entries <- strsplit(as.character(meta[[tmp_col]]), "_and_")
      grp_mask <- sapply(split_entries, function(parts) { tup[[i]] %in% parts })
    } else {
      grp_mask <- meta[[tmp_col]] == tup[[i]]
    }
    rna_grp_mask <- grp_mask & rna_mask
    n_rna[[i]] <- sum(rna_grp_mask)
    ribo_grp_mask <- grp_mask & ribo_mask
    n_ribo[[i]] <- sum(ribo_grp_mask)
    
    if (sum(rna_grp_mask) != 0) {
      rna_vals <- unique(meta[rna_grp_mask, contrast_col, drop=FALSE])
      weight_dict_rna <- as.list(table(meta[rna_grp_mask, contrast_col])/sum(rna_grp_mask))
      grp_rna[[i]] <- construct_contrast_string(rna_vals[,1], "RNA", meta, weight_dict_rna)
      order[[i]] <- sum(meta[rna_grp_mask, "control"] == "test")
    }
    if (sum(ribo_grp_mask) != 0) {
      ribo_vals <- unique(meta[ribo_grp_mask, contrast_col, drop=FALSE])
      weight_dict_ribo <- as.list(table(meta[ribo_grp_mask, contrast_col])/sum(ribo_grp_mask))
      grp_ribo[[i]] <- construct_contrast_string(ribo_vals[,1], "Ribo", meta, weight_dict_ribo)
      order[[i]] <- sum(meta[ribo_grp_mask, "control"] == "test")
    }
  }
  
  if ((length(order) == 2) && (order[[1]] < order[[2]])) {
    if (length(grp_rna) == 2) { grp_rna <- grp_rna[c(2, 1)]; n_rna <- n_rna[c(2, 1)] }
    if (length(grp_ribo) == 2) { grp_ribo <- grp_ribo[c(2, 1)]; n_ribo <- n_ribo[c(2, 1)] }
    tup <- tup[c(2, 1)]
  }
  
  # PASS has_ribo_data flag
  evaluate_contrasts(grp_rna, grp_ribo, tup, n_rna, n_ribo, fit_paired, fit_rna, has_ribo_data, outdir, log_file)
}

exec_one_vs_all <- function(meta, target_group, contrast_col, fit_paired, fit_rna, has_ribo_data, outdir, log_file) {
    # 1. Identify Target vs Other Samples
    # Note: We use startsWith to capture hierarchies (e.g. Target="MBL" includes "MBL.SHH")
    target_mask <- startsWith(as.character(meta[[contrast_col]]), target_group)
    other_mask  <- !target_mask
    
    # 2. Check Replicates for 'Other'
    # We need to ensure the 'Other' group exists and has sufficient data
    has_other_rna  <- sum(other_mask & meta$seq_type == "RNA") >= 2
    has_other_ribo <- sum(other_mask & meta$seq_type == "Ribo") >= 2
    
    # If we don't have enough "Other" samples, skip
    if (!has_other_rna) return(NULL)
    if (has_ribo_data && !has_other_ribo) return(NULL)

    # 3. Construct Strings
    # Helper to get string and N-count
    get_str <- function(mask, type) {
        sub_meta <- meta[mask & meta$seq_type == type, ]
        if (nrow(sub_meta) == 0) return(NULL)
        
        vals <- unique(sub_meta[[contrast_col]])
        # Weight by sample frequency (Micro-average)
        w_dict <- as.list(table(sub_meta[[contrast_col]]) / nrow(sub_meta))
        
        constr <- construct_contrast_string(vals, type, meta, w_dict)
        return(list(s=constr, n=nrow(sub_meta)))
    }

    # Target Strings
    t_rna <- get_str(target_mask, "RNA")
    t_ribo <- if(has_ribo_data) get_str(target_mask, "Ribo") else NULL
    
    # Other Strings
    o_rna <- get_str(other_mask, "RNA")
    o_ribo <- if(has_ribo_data) get_str(other_mask, "Ribo") else NULL
    
    # 4. Prepare Arguments for Standard Evaluator
    # We basically fake a "Tuple" input for evaluate_contrasts
    
    # Fake Tuple for naming: c(Target, "OTHER")
    tup_fake <- c(target_group, "OTHER")
    
    # RNA Inputs
    grp_rna <- list(t_rna$s, o_rna$s)
    n_rna <- list(t_rna$n, o_rna$n)
    
    # Ribo Inputs
    grp_ribo <- list()
    n_ribo <- list()
    if (has_ribo_data) {
        grp_ribo <- list(t_ribo$s, o_ribo$s)
        n_ribo <- list(t_ribo$n, o_ribo$n)
    }
    
    # Pass to existing function
    # It will handle the file creation, fitting, and plotting
    evaluate_contrasts(grp_rna, grp_ribo, tup_fake, n_rna, n_ribo, fit_paired, fit_rna, has_ribo_data, outdir, log_file)
}

evaluate_contrasts <- function(grp_rna, grp_ribo, tup, n_rna, n_ribo, fit_paired, fit_rna, has_ribo_data, outdir, log_file) {
  strat_string <- ""
  
  # --- 1. RNA Contrasts ---
  if (length(unlist(grp_rna)) == 2) {
    n_msg <- paste0("    (n = ", n_rna[[1]], " / ", n_rna[[2]], ")")

    if (has_ribo_data) {
        # --- CASE A: DUAL MODE (Ribo exists) ---
        # 1. Paired/Shared Model -> Suffix "_RNA"
        grp_contrast_paired <- paste0("makeContrasts(", paste(grp_rna, collapse = " - "), ", levels=design)")
        contrast_id_sub <- paste0(tup[[1]], "__", tup[[2]], strat_string, "_RNA")
        out_prefix_sub <- paste0(outdir, "RNA/", contrast_id_sub)
        title_sub <- paste0(contrast_id_sub, n_msg, " (Shared)")
        
        if (!is.null(grp_contrast_paired)) print(paste0("Evaluating RNA (Shared): ", grp_contrast_paired))
        eval_contrast(fit_paired, grp_contrast_paired, out_prefix_sub, title_sub, log_file)

        # 2. Independent Model -> Suffix "_RNA_full"
        grp_contrast_full <- paste0("makeContrasts(", paste(grp_rna, collapse = " - "), ", levels=design.rna)")
        contrast_id_full <- paste0(tup[[1]], "__", tup[[2]], strat_string, "_RNA_full")
        out_prefix_full <- paste0(outdir, "RNA/", contrast_id_full)
        title_full <- paste0(contrast_id_full, n_msg, " (All)")
        
        eval_contrast(fit_rna, grp_contrast_full, out_prefix_full, title_full, log_file)
        
    } else {
        # --- CASE B: RNA-ONLY MODE ---
        # Only output one file with suffix "_RNA", derived from the Independent (fit_rna) model.
        grp_contrast_full <- paste0("makeContrasts(", paste(grp_rna, collapse = " - "), ", levels=design.rna)")
        contrast_id_full <- paste0(tup[[1]], "__", tup[[2]], strat_string, "_RNA") # Simple Suffix
        out_prefix_full <- paste0(outdir, "RNA/", contrast_id_full)
        title_full <- paste0(contrast_id_full, n_msg, " (RNA Only)")
        
        eval_contrast(fit_rna, grp_contrast_full, out_prefix_full, title_full, log_file)
    }
  }

  # --- 2. Ribo Contrasts ---
  if (length(unlist(grp_ribo)) == 2) {
    n_msg <- paste0("    (n = ", n_ribo[[1]], " / ", n_ribo[[2]], ")")
    
    grp_contrast <- paste0("makeContrasts(", paste(grp_ribo, collapse = " - "), ", levels=design)")
    contrast_id <- paste0(tup[[1]], "__", tup[[2]], strat_string, "_Ribo")
    out_prefix <- paste0(outdir, "Ribo/", contrast_id)
    title <- paste0(contrast_id, n_msg, " (Shared/Paired)")
    
    if (!is.null(grp_contrast)) print(paste0("Evaluating Ribo: ", grp_contrast))
    eval_contrast(fit_paired, grp_contrast, out_prefix, title, log_file)
  }
  
  # --- 3. dTE Contrast ---
  if ((length(unlist(grp_ribo)) == 2) && (length(unlist(grp_rna)) == 2)) {
    left_cmd <- paste0("(", grp_ribo[[1]], " - ", grp_rna[[1]], ")")
    right_cmd <- paste0("(", grp_ribo[[2]], " - ", grp_rna[[2]], ")")
    grp_contrast <- paste0("makeContrasts(", left_cmd, " - ", right_cmd, ", levels=design)")
    contrast_id <- paste0(tup[[1]], "__", tup[[2]], strat_string, "_dTE")
    out_prefix <- paste0(outdir, "dTE/", contrast_id)
    title <- paste0(contrast_id, "    (n = ", n_ribo[[1]], " / ", n_rna[[1]], " / ", n_ribo[[2]], " / ", n_rna[[2]], ")")
    
    if (!is.null(grp_contrast)) print(grp_contrast)
    eval_contrast(fit_paired, grp_contrast, out_prefix, title, log_file)
  }
}


evaluate_unique_contrast <- function(meta, uniq_val, contrast_col, fit, outdir, log_file) {
  strat_string <- ""
  ribo_mask <- meta$seq_type == "Ribo"
  rna_mask <- meta$seq_type == "RNA"
  num_dots <- str_count(uniq_val, "\\.") + 1
  tmp_col <- paste0(contrast_col, ".", num_dots)
  
  if (!grepl("_and_", uniq_val) && any(grepl("_and_", meta[[tmp_col]]))) {
    split_entries <- strsplit(as.character(meta[[tmp_col]]), "_and_")
    uniq_mask <- sapply(split_entries, function(parts) { uniq_val %in% parts })
  } else {
    uniq_mask <- meta[[tmp_col]] == uniq_val
  }
  
  rna_grp_mask <- uniq_mask & rna_mask
  ribo_grp_mask <- uniq_mask & ribo_mask
  n_rna <- sum(rna_grp_mask)
  n_ribo <- sum(ribo_grp_mask)
  
  if (sum(ribo_grp_mask) != 0 && sum(rna_grp_mask) != 0) {
    ribo_vals <- unique(meta[ribo_grp_mask, contrast_col, drop=FALSE])
    weight_dict_ribo <- as.list(table(meta[ribo_grp_mask, contrast_col])/sum(ribo_grp_mask))
    rna_vals <- unique(meta[rna_grp_mask, contrast_col, drop=FALSE])
    weight_dict_rna <- as.list(table(meta[rna_grp_mask, contrast_col])/sum(rna_grp_mask))
    
    uniq_ribo <- construct_contrast_string(ribo_vals[,1], "Ribo", meta, weight_dict_ribo)
    uniq_rna <- construct_contrast_string(rna_vals[,1], "RNA", meta, weight_dict_rna)
    
    uniq_contrast <- paste0("makeContrasts(", uniq_ribo, " - ", uniq_rna, ", levels=design)")
    
    if (!is.null(uniq_contrast)) print(uniq_contrast)
    
    contrast_id <- paste0(uniq_val, strat_string, "_TE")
    out_prefix <- paste0(outdir, "TE/", contrast_id)
    
    # RESTORED: Title with sample counts
    title <- paste0(contrast_id, "    (n = ", n_ribo, " / ", n_rna, ")")
    
    eval_contrast(fit, uniq_contrast, out_prefix, title, log_file)
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
    # FIX IS HERE: Changed collapse="__" to collapse="_and_"
    grp_strings[i] <- paste0(weight_dict[[key]], "*group", paste0(val_grps, collapse = "_and_"), "__", seq_type)
  }
  grp_string <- paste0(grp_strings, collapse=" + ")
  return(paste0("(", grp_string, ")"))
}


eval_contrast <- function(fit, contrast, out_prefix, title, log_file) {
    cleaned_contrast <- sub("^makeContrasts\\(", "", sub(", levels=.*\\)$", "", contrast))
    cat("Evaluating contrast ", title, " : ", cleaned_contrast, "\n", file = log_file, append = TRUE)
    
    # 1. Evaluate
    lrt <- glmQLFTest(fit, contrast=eval(parse(text = contrast)))
    res <- topTags(lrt, n=Inf)$table
    res <- res |> tibble::rownames_to_column('row_id')
    
    has_tx <- exists("tx.table")
    
    # --- DEDUPLICATION LOGIC ---
    is_rna_subset <- grepl("_RNA", out_prefix) && grepl("_subset", out_prefix)
    # Check if we have necessary columns for deduplication
    can_dedup <- has_tx && all(c("translon_id", "transcript_id") %in% colnames(tx.table))
    
    if (opt$feature_level == "transcript" && is_rna_subset && can_dedup) {
        map_df <- tx.table[, c("translon_id", "transcript_id")]
        joined <- res %>% left_join(map_df, by=c("row_id"="translon_id"))
        
        if (sum(!is.na(joined$transcript_id)) > 0) {
             res <- joined %>%
                filter(!is.na(transcript_id)) %>%
                group_by(transcript_id) %>%
                slice_min(order_by = PValue, n = 1, with_ties = FALSE) %>%
                ungroup() %>%
                select(-transcript_id)
        }
    }

    # --- MAPPING LOGIC ---
    if (has_tx) {
        maps_list <- list()
        targets <- c(real_gene_id = "gene_id", real_gene_name = "gene_name")
        
        # 1. Map Translons
        if ("translon_id" %in% colnames(tx.table)) {
            maps_list$translon <- tx.table %>% select(key_id = translon_id, any_of(targets))
        }
        # 2. Map Transcripts
        if ("transcript_id" %in% colnames(tx.table)) {
            maps_list$transcript <- tx.table %>% select(key_id = transcript_id, any_of(targets))
        }
        # 3. Map Genes (User defined ID)
        if (opt$tx_table_col %in% colnames(tx.table)) {
             maps_list$gene <- tx.table %>% select(key_id = all_of(opt$tx_table_col), any_of(targets))
        }

        if (length(maps_list) > 0) {
            master_map <- bind_rows(maps_list) %>% 
                          filter(key_id != "" & !is.na(key_id)) %>%
                          distinct(key_id, .keep_all = TRUE)
            
            # Try exact match first
            res <- left_join(res, master_map, by=c("row_id" = "key_id"))
            
            # For unmatched rows, try stripping R's .N deduplication suffix
            # (e.g. TT_7_-_61798108.1 -> TT_7_-_61798108)
            unmatched <- is.na(res$real_gene_id) | res$real_gene_id == ""
            if (any(unmatched, na.rm=TRUE)) {
                stripped_ids <- sub("\\.[0-9]+$", "", res$row_id[unmatched])
                retry <- master_map[match(stripped_ids, master_map$key_id), , drop=FALSE]
                if (nrow(retry) > 0) {
                    res$real_gene_id[unmatched]   <- retry$real_gene_id
                    res$real_gene_name[unmatched] <- retry$real_gene_name
                }
            }
        }
        
        # Fill logic
        if ("real_gene_id" %in% colnames(res)) {
            res$gene_id <- ifelse(is.na(res$real_gene_id) | res$real_gene_id == "", res$row_id, res$real_gene_id)
        } else {
            res$gene_id <- res$row_id
        }
        
        if ("real_gene_name" %in% colnames(res)) {
            res$gene_name <- ifelse(is.na(res$real_gene_name) | res$real_gene_name == "", res$row_id, res$real_gene_name)
        } else {
            res$gene_name <- res$row_id
        }
        
    } else {
        res$gene_id   <- res$row_id
        res$gene_name <- res$row_id
    }
    
    # --- PLOT LABEL LOGIC ---
    if (opt$feature_level == "gene") {
        # Gene Mode: Clean Names only
        res$plot_label <- res$gene_name
    } else {
        # Transcript Mode: Composite Name to distinguish isoforms/ORFs
        # e.g. "TPM3 (ENST000...)"
        res$plot_label <- ifelse(res$gene_name == res$row_id, 
                                 res$gene_name, 
                                 paste0(res$gene_name, " (", res$row_id, ")"))
    }

    # 4. Final Formatting
    out_table <- res |> 
        select(gene_id, gene_name, row_id, logFC, logCPM, any_of("F"), PValue, FDR) |> 
        arrange(PValue)

    write.csv(out_table, paste0(out_prefix, ".csv"), row.names=FALSE)
    
    # 5. Plotting
    ymax_vals <- -log10(res$PValue)[-log10(res$PValue) != Inf]
    ymax <- if(length(ymax_vals) > 0) max(ymax_vals) + 1 else 10
    
    plot <- EnhancedVolcano(res,
        lab = res$plot_label, # Uses the logic defined above
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
    make_option(c("-t", "--tx_table_path"      ), type="character", default=NULL                   , metavar="path",    help="Table linking transcript and gene IDs/names"                                          ),
    make_option(c("-v", "--feature_level" ), type="character", default="gene"             , metavar="string" , help="Analysis level: 'gene' or 'transcript'. Determines expansion logic and naming."                                                  ),
    make_option(c("-f", "--tx_table_col"  ), type="character",  default="gene_id"        , metavar="string" , help="column name of tx_table to use as keys"                                                    ),
    make_option(c("-d", "--id_col"        ), type="integer"  , default=1                    , metavar="integer", help="Column containing identifiers to be used."                                              ),
    make_option(c("-s", "--sep"           ), type="character", default=','                  , metavar="string" , help="Separator of input table."                                              ),
    make_option(c("-o", "--outdir"        ), type="character", default='out', metavar="path"   , help="Output directory."                                                                      ),
    make_option(c("-u", "--outer_join"    ), type="logical"  , default=FALSE                , metavar="boolean", help="Outer join RNA and Ribo reads, fill NAs with 0 expression."                                                     ),
    make_option(c("-l", "--cores"         ), type="integer"  , default=1                    , metavar="integer", help="Number of cores."                                                                       ),
    make_option(c("-a", "--contrast_cols" ), type="character", default="treatment_id"               , metavar="character", help="Column names from which contrasts are derived; separated by commas."          ),
    make_option(c("-O", "--skip_one_vs_all"), action="store_true", default=FALSE, help="If set, skips One-vs-All contrasts (Group vs Rest). Default is to RUN them."),
    make_option(c("-p", "--plot_ids"), type="logical", default=FALSE, metavar="boolean", help="Plot Smart IDs instead of replicate numbers in PCA/MDS."),
    make_option(c("-e", "--no_batch_factor"), action="store_true", default=FALSE, help="Don't create factors for batch date. Can be necessary to achieve full rank for some settings"),
    make_option(c("-S", "--save_model"), type="character", default=NULL, metavar="path", help="Path to save the fitted model (RData file). Defaults to 'dTEct_model.RData' in outdir if not specified. Ignored if --load_model is used."),
    make_option(c("-L", "--load_model"), type="character", default=NULL, metavar="path", help="Path to load a pre-fitted model from (skips fitting)."),
    make_option(c("-k", "--skip_pairwise"), action="store_true", default=FALSE, help="If set, skips all pairwise contrasts and only runs One-vs-All (if enabled)."),
    make_option(c("-T", "--test_run"), action="store_true", default=FALSE, help="Run in test mode: subsets data to first valid contrast and fraction of genes for rapid debugging.")
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

# Validate Inputs (Skip some checks if loading model)
if (is.null(opt$load_model)) {
  validate_inputs(opt)
} else {
  if (!file.exists(opt$load_model)) {
    stop(paste("Model file not found:", opt$load_model))
  }
}


## DEBUG
# opt$metadata <- "sample_sheet.csv"
# opt$rna_counts <- "table.csv"
# # opt$ribo_counts <- "ribo/gene_agg_NumReads.csv"
# opt$count_col <- 5
# opt$sep <- ","
# # opt$tx_table <- "tx_table.csv"
# opt$cores <- 1
# opt$no_batch_factor <- TRUE
# opt$contrast_cols <- "disease_id"
# # opt$outdir <- "test/"

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
# Add trailing / to opt$outdir if missing
if (!grepl("/", opt$outdir)) {
  opt$outdir <- paste0(opt$outdir, "/")
}
# Check whether outdir exists and create it
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
# Define the path to the log file
log_file <- paste0(opt$outdir, "run_info.txt")

# --- LOAD MODEL LOGIC ---
if (!is.null(opt$load_model)) {
  cat(paste0("Loading model from ", opt$load_model, "...\n"))
  load(opt$load_model)
  # Update internal opt with new run parameters if needed (e.g. outdir change)
  # We generally want to keep the fitted objects but maybe change output location
  # For now, we assume the user provides consistent args or we rely on the loaded environment
  # EXCEPT for outdir and cores which might change
  
  # Ensure save_model is ignored when loading
  opt$save_model <- NULL
  
  # Re-register cores
  register(MulticoreParam(opt$cores))
  
  # Re-open log file since we just overwrote it or it might be closed
  log_file <- paste0(opt$outdir, "run_info.txt")
  cat("Model loaded. Skipping data processing and fitting.\n", file = log_file, append = TRUE)
  
} else {
  # --- NORMAL EXECUTION START ---
  
  # Set default save path if not provided
  if (is.null(opt$save_model)) {
      opt$save_model <- file.path(opt$outdir, "dTEct_model.RData")
      cat(paste0("No --save_model provided. Defaulting to: ", opt$save_model, "\n"))
  }
  
  # Clear existing content of the log file at the start of the script
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


if (!is.null(opt$tx_table_path)) {
  # Get TX table to parse gene names
  tx.table <- read.csv(opt$tx_table_path)

  if (opt$feature_level %in% c("transcript", "translon")) {
     # --- TRANSCRIPT/TRANSLON MODE ---
     cat("Generating feature2name map for TRANSCRIPT/TRANSLON level analysis...
", file = log_file, append = TRUE)
  


     
     # 1. Map for Translons (Used for dTE and Ribo)
     # We use the Translon ID as the key, and a composite name for the label
     map_translons <- (
       tx.table
       |> mutate(
            id = translon_id, 
            name = paste0(gene_name, " (", translon_id, ")")
          )
       |> select(id, name)
     )
     
     # 2. Map for Transcripts (Used for Full RNA model)
     # In this mode, the Full RNA model runs on Transcript IDs
     map_transcripts <- (
       tx.table
       |> mutate(
            id = transcript_id,
            # For RNA plots, we label with Gene Name + Transcript ID to be specific
            name = paste0(gene_name, " (", transcript_id, ")")
          )
       |> select(id, name)
     )
     
     # 3. Combine both maps
     feature2name <- rbind(map_translons, map_transcripts) |> distinct()
     
  } else {
    # --- GENE MODE ---
    cat("Generating feature2name map for GENE level analysis...\n", file = log_file, append = TRUE)
    
    # Standard 1:1 mapping using the user-defined column (usually gene_id)
    feature2name <- (
      tx.table 
      |> select(opt$tx_table_col, "gene_name")
      |> distinct()
      |> dplyr::rename(id = opt$tx_table_col, name=gene_name)
      |> mutate(
            name = if_else(is.na(name) | name == "" | is.nan(name), id, name),
          )
    )
  }
} else {
    # Fallback (Dummy Map)
    cat("No tx_table provided. Using IDs as Names.\n", file = log_file, append = TRUE)
    ids <- c()
    if(!is.null(opt$rna_counts)) ids <- c(ids, rownames(rna_counts))
    if(!is.null(opt$ribo_counts)) ids <- c(ids, rownames(ribo_counts))
    ids <- unique(ids)
    
    feature2name <- data.frame(id = ids, name = ids)
}


dir.create(paste0(opt$outdir, 'dTE'), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(opt$outdir, 'Ribo'), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(opt$outdir, 'RNA'), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(opt$outdir, 'TE'), showWarnings = FALSE, recursive = TRUE)

# Intersect RNA and Ribo counts -------------------------------------------------

if (!is.null(opt$rna_counts) && !is.null(opt$ribo_counts)) {

  # --- A. TRANSCRIPT MODE (Expansion Logic) ---
  if (opt$feature_level %in% c("transcript", "translon")) {

    
    if (is.null(opt$tx_table_path)) {
       stop("Error: 'transcript' mode requires --tx_table_path to map Translons to Transcripts.")
    }

    cat("Performing Translon-level expansion (Mapping RNA transcripts to multiple Translons)...\n", file = log_file, append = TRUE)
    
    # 1. Create Map: Translon ID -> Transcript ID
    # FIX: Use distinct() instead of slice() to avoid S4Vectors namespace conflict
    translon_rna_map <- tx.table %>% 
      select(translon_id, transcript_id) %>% 
      # Remove rows with missing/empty identifiers
      filter(!is.na(translon_id) & translon_id != "") %>%
      filter(!is.na(transcript_id) & transcript_id != "") %>%
      # Remove exact duplicates
      distinct() %>%
      # Ensure translon_id is unique (Key) by keeping the first occurrence
      distinct(translon_id, .keep_all = TRUE) %>%
      tibble::column_to_rownames("translon_id")

    # Save a global named vector: translon_id -> transcript_id
    # Used to remap IDs in RNA contrasts from the paired model
    translon_to_transcript_map <<- setNames(
        translon_rna_map$transcript_id,
        rownames(translon_rna_map)
    )

    # 2. Identify parent transcript for every Translon in Ribo matrix
    required_transcripts <- translon_rna_map[rownames(ribo_counts), "transcript_id"]
    
    # 3. Filter valid rows (Translon exists in Map AND Transcript exists in RNA Matrix)
    valid_mask <- !is.na(required_transcripts) & (required_transcripts %in% rownames(rna_counts))
    
    # 4. Subset Ribo counts
    ribo_counts_select <- ribo_counts[valid_mask, , drop=FALSE]
    
    # 5. Expand RNA Matrix
    # We grab the RNA row corresponding to the parent transcript of the Ribo row
    parent_ids_for_translons <- translon_rna_map[rownames(ribo_counts_select), "transcript_id"]
    rna_counts_expanded <- rna_counts[parent_ids_for_translons, , drop=FALSE]
    
    # Save with original transcript rownames for the independent RNA-only model
    # (transcript IDs must be preserved so RNA contrast output uses ST_.../ENST... IDs)
    rna_counts_original <- rna_counts_expanded

    # 6. Align Rownames (Critical for edgeR pairing)
    # The RNA matrix must now have Translon IDs as rownames to match the Ribo matrix
    rownames(rna_counts_expanded) <- rownames(ribo_counts_select)
    
    rna_counts_select <- rna_counts_expanded
    
    cat(sprintf("Expanded %d unique RNA transcripts to cover %d Translons.\n", length(unique(parent_ids_for_translons)), nrow(ribo_counts_select)), file = log_file, append = TRUE)
    
  } else {
    # --- B. GENE MODE (Simple Intersection) ---
    cat("Performing Gene-level intersection (1:1 Matching)...\n", file = log_file, append = TRUE)
    
    # We assume rownames are Gene IDs in both matrices
    common_ids <- intersect(rownames(rna_counts), rownames(ribo_counts))
    
    if (length(common_ids) == 0) {
        stop("Error: No common IDs found between RNA and Ribo counts. Check if files match 'feature_level' setting.")
    }
    
    ribo_counts_select <- ribo_counts[common_ids, , drop=FALSE]
    rna_counts_select <- rna_counts[common_ids, , drop=FALSE]
    
    cat(sprintf("Intersected %d genes common to both datasets.\n", length(common_ids)), file = log_file, append = TRUE)
  }
  
  # Combine into final counts matrix for the Paired Model
  counts <- cbind(RNA=as.matrix(rna_counts_select), Ribo=as.matrix(ribo_counts_select))
  seq_types <- c("RNA", "Ribo")
  
} else if (!is.null(opt$rna_counts)) {
  # RNA Only mode
  counts <- rna_counts
  seq_types <- c("RNA")
} else {
  # Ribo Only mode
  counts <- ribo_counts
  seq_types <- c("Ribo")
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

# 1. CLEAN & DEDUPLICATE METADATA
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

# -----------------------------------------------------------------------------
# 2. DICTIONARY GENERATION & COLUMN SPLITTING
# -----------------------------------------------------------------------------
cat("Parsing contrast groups and splitting hierarchies...\n", file=log_file, append=TRUE)

comb_dict <- list()
uniq_dict <- list()
num_col_splits <- list()

# Define columns to process: Input columns + Combined column (if multiple)
cols_to_process <- contrast_grps

# Combine contrast group columns if needed
if (length(contrast_grps) > 1) {
  # Legacy logic: Create combined column with ':' then replace with '_and_'
  meta.table[[contrast_col]] <- do.call(paste, c(meta.table[contrast_grps], sep=":"))
  cols_to_process <- unique(c(cols_to_process, contrast_col))
} 

# Split columns into subclasses and generate dictionary
for (column in cols_to_process) {
  
  if (!is.null(meta.table[[column]]) && length(meta.table[[column]]) > 0) {
    # Normalize separators
    meta.table[[column]] <- gsub(":", "_and_", meta.table[[column]])
    # Sanitize names (ensure valid R names)
    meta.table[[column]] <- make.names(meta.table[[column]], unique = FALSE)
  } else {
    warning(sprintf("Column '%s' is missing or empty. Skipping.", column))
    next
  }
  
  # Calculate depth of hierarchy
  max_splits <- max(str_count(meta.table[[column]], "\\.") + 1)
  num_col_splits[[column]] <- max_splits
  column_names <- paste0(column, ".", seq_len(max_splits))
  
  # Create split columns
  for (i in seq_len(max_splits)) {
    meta.table[[column_names[i]]] <- sapply(strsplit(meta.table[[column]], "\\."), function(x) {
      paste(x[1:i], collapse = ".")
    })
  }
  
  # Normalize main column (pad with .NA)
  meta.table[[column]] <- sapply(strsplit(meta.table[[column]], "\\."), function(x) {
    paste(c(x, rep("NA", max_splits - length(x))), collapse = ".")
  })
  
  # Log
  cat(paste0("  Column '", column, "': Max Depth = ", max_splits, "\n"), file=log_file, append=TRUE)

  # Build Combinations
  combs <- list()
  for (i in seq_len(max_splits)) {
    vals <- unique(meta.table[[column_names[i]]])
    vals <- vals[!grepl("\\.NA", vals)]
    
    if (length(vals) > 1) {
      groups <- combn(vals, 2, simplify)
      
      # Handle Orthogonal Combinations using _and_ logic
      if (any(grep("_and_", meta.table[[column]]))) {
        orth_vals_1 <- unique(sapply(strsplit(meta.table[[column]], "_and_"), `[`, 1))
        orth_vals_2 <- unique(sapply(strsplit(meta.table[[column]], "_and_"), `[`, 2))
        
        # Guard against single-value derived columns causing errors in combn
        orth_groups_1 <- if(length(orth_vals_1) > 1) combn(orth_vals_1, 2, simplify) else NULL
        orth_groups_2 <- if(length(orth_vals_2) > 1) combn(orth_vals_2, 2, simplify) else NULL
        
        if (!is.null(orth_groups_1)) groups <- cbind(groups, orth_groups_1)
        if (!is.null(orth_groups_2)) groups <- cbind(groups, orth_groups_2)
      }
      
      # Filter for Valid Comparisons (Shared Parent)
      for (j in 1:length(groups[1,])) {
        if (i > 1) {
          # Check if parents (up to i-1) are the same
          if (substr_to_nth_dot(groups[,j][1], i-1) == substr_to_nth_dot(groups[,j][2], i-1)) {
            combs <- append(combs, list(groups[,j]))
          }
        } else {
          combs <- append(combs, list(groups[,j]))
        }
      }
    }
  }
  comb_dict[[column]] <- unique(combs) # Ensure uniqueness
  uniq_dict[[column]] <- unique(unlist(combs))
}

# -----------------------------------------------------------------------------
# 3. INTELLIGENT GROUPING (Renaming sparse subgroups)
# -----------------------------------------------------------------------------
for (contrast in cols_to_process) {
  # Count depth of hierarchy
  dot_counts <- sapply(meta.table[[contrast]], function(x) {
    count <- gregexpr("\\.", x)[[1]]
    ifelse(count[1] == -1, 0, length(count))
  })
  
  # Process deepest hierarchies first
  sorted_entries <- unique(meta.table[order(dot_counts, decreasing = TRUE),][[contrast]])
  
  for (entry in sorted_entries) {
    # Check this group across the ENTIRE table (both RNA and Ribo)
    mask_all <- startsWith(meta.table[[contrast]], entry)
    
    # We only care if this group actually exists in the data
    if (sum(mask_all) > 0) {
      
      # FIX: Use dplyr::count explicitly to avoid namespace collision with plyr/matrixStats
      counts_by_type <- meta.table[mask_all, ] %>% dplyr::count(seq_type)
      
      # Identify if any present seq_type has < 2 replicates
      # Note: We check if 'n' is less than 2. 
      is_sparse <- any(counts_by_type$n < 2)
      
      if (is_sparse) {
        if (grepl("\\.", entry)) {
          # Collapse to parent group
          new_group <- sub("\\.[^.]*$", "", entry) 
          
          # Apply rename to ALL rows matching this entry (RNA AND Ribo)
          meta.table[mask_all, contrast] <- new_group
          
          cat("Optimization: Renaming sparse subgroup (consistent):", entry, "->", new_group, "\n", file = log_file, append = TRUE)
        } else {
          cat("Optimization: Group", entry, "has low replicates but is top-level. Keeping.\n", file = log_file, append = TRUE)
        }
      }
    }
  }
  
  # Regenerate hierarchy columns after renaming
  max_splits <- num_col_splits[[contrast]]
  column_names <- paste0(contrast, ".", seq_len(max_splits))
  
  for (i in seq_len(max_splits)) {
    meta.table[[column_names[i]]] <- sapply(strsplit(meta.table[[contrast]], "\\."), function(x) {
      paste(x[1:i], collapse = ".")
    })
  }
  # Re-normalize main column to ensure consistent NA padding
  meta.table[[contrast]] <- sapply(strsplit(meta.table[[contrast]], "\\."), function(x) {
    paste(c(x, rep("NA", max_splits - length(x))), collapse = ".")
  })
}

# -----------------------------------------------------------------------------
# 4. PREPARE FULL DATASET
# -----------------------------------------------------------------------------
# Select relevant columns 
meta.samples <- (
  meta.table
  |> filter(counts_col %in% colnames(counts))
  |> select(
    smart_id,
    any_of(c("sample_type", "batch_date")),
    # Capture all hierarchy columns correctly
    # Capture all hierarchy columns correctly
    any_of(unlist(lapply(cols_to_process, function(col) grep(paste0("^", col), colnames(meta.table), value = TRUE)))),
    # Capture all ID columns for plotting
    matches("_id$"),
    rep = replicate_num,
    control = test_or_control,
    counts_col,
    seq_type
  )
  |> mutate(across(where(is.character), as.factor))
  |> replace_na(list(rep = 1))
  |> distinct()
)

# Align Counts
counts <- counts[, as.character(meta.samples$counts_col), drop = FALSE]

# Create Grouping Column (FIXED)
if (length(contrast_grps) > 0) {
    
    grouping_cols <- unlist(lapply(contrast_grps, function(col) {
        # Find hierarchy columns (col.1, col.2, etc.)
        hier_cols <- grep(paste0("^", col, "\\.[0-9]+$"), colnames(meta.samples), value=TRUE)
        
        if (length(hier_cols) > 0) {
             # CASE A: Hierarchy Found
             # Sort numerically (.1, .2, .3)
             extract_num <- function(x) as.numeric(sub(paste0("^", col, "\\."), "", x))
             hier_cols <- hier_cols[order(extract_num(hier_cols))]
             
             # FIX: Return ONLY hierarchy columns. Do NOT append the main 'col'.
             return(hier_cols)
        } else {
             # CASE B: No Hierarchy (Flat group)
             # Use the original column name
             return(col)
        }
    }))
    
    # Filter to ensure they exist in this specific subset
    grouping_cols <- intersect(grouping_cols, colnames(meta.samples))
    
    # Construct Group ID
    meta.samples <- meta.samples %>% 
        mutate(group = paste0(apply(meta.samples[, grouping_cols, drop=FALSE], 1, paste, collapse = "_and_"), "__", seq_type))
        
    cat(paste0("Grouping columns used: ", paste(grouping_cols, collapse=", "), "\n"), file=log_file, append=TRUE)
    cat(paste0("Example Group ID: ", meta.samples$group[1], "\n"), file=log_file, append=TRUE)

} else {
    stop("No valid contrast columns found in metadata.")
}

# Define Design Formula
f = "~0 + group"

has_var <- function(vec) { length(unique(as.character(vec))) > 1 }

if (!opt$test_run && !opt$no_batch_factor && "batch_date" %in% colnames(meta.samples) && has_var(meta.samples$batch_date)) {
    # CHECK FOR COLLINEARITY
    # If batch_date is perfectly correlated with seq_type or group, it will cause rank deficiency.
    # We check if each batch contains only ONE level of the other factor.
    
    is_nested_in_type <- all(rowSums(table(meta.samples$batch_date, meta.samples$seq_type) > 0) == 1)
    is_nested_in_group <- all(rowSums(table(meta.samples$batch_date, meta.samples$group) > 0) == 1)
    
    if (is_nested_in_type || is_nested_in_group) {
        cat("WARNING: 'batch_date' is perfectly collinear with experimental factors. Dropping from design to avoid rank deficiency.\n", file=log_file, append=TRUE)
    } else {
        f = paste0(f, " + batch_date")
        cat("Adding 'batch_date' to design model.\n", file=log_file, append=TRUE)
    }
} else if (opt$test_run) {
    cat("TEST RUN MSG: Automatically skipping 'batch_date' to avoid rank deficiency on subset.\n", file=log_file, append=TRUE)
}

if (!opt$test_run && "sample_type" %in% colnames(meta.samples) && has_var(meta.samples$sample_type)) {
    f = paste0(f, " + sample_type")
    cat("Adding 'sample_type' to design model.\n", file=log_file, append=TRUE)
} else if (opt$test_run) {
    cat("TEST RUN MSG: Automatically skipping 'sample_type' to avoid rank deficiency on subset.\n", file=log_file, append=TRUE)
}


# Create Design Matrix
meta.design <- meta.samples %>% column_to_rownames("counts_col")
design <- model.matrix(as.formula(f), data=data.frame(meta.design))
# -----------------------------------------------------------------------------
# 5. NORMALIZE & QC
# -----------------------------------------------------------------------------
cat("Constructing DGE Object with ALL samples...\n", file=log_file, append=TRUE)
cat(paste0("Formula: ", f, "\n"), file=log_file, append=TRUE)

# --- TEST RUN SUBSETTING ---
if (opt$test_run) {
  cat("\n*** TEST RUN ENABLED ***\n", file = log_file, append = TRUE)
  cat("Subsetting data for rapid debugging...\n", file = log_file, append = TRUE)
  
  # 1. Subset Features (Top 2000)
  n_total <- nrow(counts)
  n_keep <- min(2000, n_total)
  counts <- counts[1:n_keep, , drop=FALSE]
  cat(sprintf("  - Kept top %d / %d features.\n", n_keep, n_total), file = log_file, append = TRUE)
  
  # 2. Subset Samples (First valid contrast only)
  # Ensure we keep the internal structure consistency (cols_to_process, etc.)
  # We parse the first requested contrast column
  target_col <- trimws(unlist(strsplit(opt$contrast_cols, ",")))[1]
  
  if (!target_col %in% colnames(meta.table)) {
      cat(sprintf("  - WARNING: Target contrast col '%s' not found. Using random 10 samples.\n", target_col), file = log_file, append = TRUE)
      keep_samples <- head(colnames(counts), 10)
  } else {
      # Find groups with >= 2 replicates (to allow for variance estimation)
      grp_counts <- table(meta.table[[target_col]])
      valid_grps <- names(grp_counts[grp_counts >= 2])
      
      if (length(valid_grps) >= 2) {
          # Keep samples from top 3 groups if available (to test One-vs-All), otherwise 2
          n_grps <- min(length(valid_grps), 3)
          grps_to_keep <- valid_grps[1:n_grps]
          
          # Collect sample IDs
          possible_ids <- rownames(meta.table)[meta.table[[target_col]] %in% grps_to_keep]
          
          # Intersect with available counts columns
          available_samples <- intersect(colnames(counts), possible_ids)
          
          # Further subset: Keep max 5 samples per group to keep run fast
          keep_samples <- c()
          for (g in grps_to_keep) {
              # Get samples for this group that are in the counts matrix
              g_samps <- rownames(meta.table)[meta.table[[target_col]] == g]
              g_samps <- intersect(g_samps, available_samples)
              
              # Take top 5
              if (length(g_samps) > 0) {
                  keep_samples <- c(keep_samples, head(g_samps, 5))
              }
          }
          
          if (length(keep_samples) < 4) {
             cat("  - Not enough samples found after intersection. Falling back to first 20 samples.\n", file = log_file, append = TRUE)
             keep_samples <- head(colnames(counts), 20)
          } else {
             cat(sprintf("  - Restricted to max 5 samples from groups: %s\n", paste(grps_to_keep, collapse=", ")), file = log_file, append = TRUE)
          }
          
      } else {
           keep_samples <- head(colnames(counts), 20)
           cat("  - Fewer than 2 valid groups found. Keeping first 20 samples.\n", file = log_file, append = TRUE)
      }
  }
  
  counts <- counts[, keep_samples, drop=FALSE]
  cat(sprintf("  - Kept %d samples.\n", ncol(counts)), file = log_file, append = TRUE)
  
  # CRITICAL FIX: Subset meta.design to match counts
  meta.design <- meta.design[colnames(counts), , drop=FALSE]
  
  cat("************************\n\n", file = log_file, append = TRUE)
}

# Pre-filtering (continue normal flow)
dge <- DGEList(counts=counts, samples=data.frame(meta.design))
keep <- filterByExpr(dge, design) # Use dge (with group info) or design if available, but design is now mismatched
# Actually, 'design' matrix (created line 1355) is also mismatched now!
# We must re-create the design matrix or subset it.
design <- design[colnames(counts), , drop=FALSE]

# Check if design has empty columns (factors dropped)
design <- design[, colSums(design != 0) > 0, drop=FALSE]

keep <- filterByExpr(dge, design)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method="TMM")

# --- SAVE NORMALIZED MATRICES ---
cat("Saving normalized LogCPM matrices...\n", file=log_file, append=TRUE)
seq_types_present <- intersect(c("RNA", "Ribo"), unique(dge$samples$seq_type))

for (st in seq_types_present) {
    # Subset DGE for this type
    sub_dge <- dge[, dge$samples$seq_type == st, keep.lib.sizes=FALSE]
    if (ncol(sub_dge) > 0) {
        logcpm_counts <- cpm(sub_dge, normalized.lib.sizes = TRUE, log=TRUE, prior.count=1)
        
        # For RNA in translon mode, remap rownames from translon IDs back to transcript IDs
        if (st == "RNA" && exists("translon_to_transcript_map")) {
            current_ids <- rownames(logcpm_counts)
            new_ids <- ifelse(current_ids %in% names(translon_to_transcript_map),
                              translon_to_transcript_map[current_ids],
                              current_ids)
            rownames(logcpm_counts) <- make.unique(new_ids)
        }
        
        # Use check.names=FALSE to prevent make.names() from mangling +/- in IDs
        formatted_data <- tibble::rownames_to_column(
            data.frame(format(logcpm_counts, digits=3, nsmall = 3), check.names=FALSE),
            "Name"
        )
        
        out_file <- paste0(opt$outdir, st, "_cpm_log_matrix.csv")
        write.table(formatted_data, file=out_file, sep=",", row.names = FALSE)
        cat(paste0("  Saved: ", out_file, "\n"), file=log_file, append=TRUE)
    }
}

# --- GENERATE QC PLOTS ---
cat("Generating QC Plots...\n", file=log_file, append=TRUE)
# Dynamic detection of all ID columns + contrast columns
qc_cols <- unique(c(contrast_grps, grep("_id$", colnames(dge$samples), value=TRUE)))
qc_cols <- qc_cols[qc_cols %in% colnames(dge$samples)]
# Filter to only keep columns with > 1 unique value
qc_cols <- qc_cols[sapply(qc_cols, function(col) length(unique(dge$samples[[col]])) > 1)]
# Explicitly exclude run_id and smart_id
qc_cols <- setdiff(qc_cols, c("run_id", "smart_id"))

present_types <- unique(dge$samples$seq_type)
qc_types <- list()
if ("RNA" %in% present_types) qc_types$RNA <- "RNA"
if ("Ribo" %in% present_types) qc_types$Ribo <- "Ribo"

# --- Per-modality RNA and Ribo plots (unchanged) ---
for (lbl in names(qc_types)) {
    mask <- dge$samples$seq_type %in% qc_types[[lbl]]
    if (sum(mask) >= 3) { 
        sub_dge <- dge[, mask]
        sub_meta <- dge$samples[mask, ]
        eval_MDS(sub_dge, sub_meta, qc_cols, opt$outdir, lbl, opt$plot_ids)
        eval_PCA(sub_dge, sub_meta, qc_cols, opt$outdir, lbl, opt$plot_ids)
        log_cpm <- cpm(sub_dge, log=TRUE, normalized.lib.sizes=TRUE)
        eval_heatmap(log_cpm, sub_meta, qc_cols, opt$outdir, lbl)
        eval_gene_clusters(log_cpm, sub_meta, qc_cols, opt$outdir, lbl)
    }
}

# --- Combined "All" plot: one point per biological sample ---
# Match RNA and Ribo by smart_id, concatenate their feature vectors so each
# biological sample is a single point with features = [RNA_features; Ribo_features].
if (length(present_types) > 1 && "RNA" %in% present_types && "Ribo" %in% present_types) {
    cat("Building concatenated-modality matrix for joint All QC plots...\n", file=log_file, append=TRUE)
    
    # Split metadata by seq_type
    meta_rna  <- dge$samples[dge$samples$seq_type == "RNA", , drop=FALSE]
    meta_ribo <- dge$samples[dge$samples$seq_type == "Ribo", , drop=FALSE]
    
    # Find biological samples present in BOTH modalities
    shared_ids <- intersect(as.character(meta_rna$smart_id),
                            as.character(meta_ribo$smart_id))
    
    if (length(shared_ids) >= 3) {
        cat(sprintf("  Found %d biological samples with both RNA and Ribo data.\n", length(shared_ids)),
            file=log_file, append=TRUE)
        
        # For each shared smart_id, get the corresponding counts_col (row in dge$samples)
        rna_col_lookup  <- setNames(rownames(meta_rna),  as.character(meta_rna$smart_id))
        ribo_col_lookup <- setNames(rownames(meta_ribo), as.character(meta_ribo$smart_id))
        
        rna_cols_matched  <- rna_col_lookup[shared_ids]
        ribo_cols_matched <- ribo_col_lookup[shared_ids]
        
        # Get the count sub-matrices
        rna_sub  <- as.matrix(dge$counts[, rna_cols_matched, drop=FALSE])
        ribo_sub <- as.matrix(dge$counts[, ribo_cols_matched, drop=FALSE])
        
        # Rename columns to the biological sample ID (smart_id)
        colnames(rna_sub)  <- shared_ids
        colnames(ribo_sub) <- shared_ids
        
        # Prefix row names to distinguish modalities
        rownames(rna_sub)  <- paste0("RNA:",  rownames(rna_sub))
        rownames(ribo_sub) <- paste0("Ribo:", rownames(ribo_sub))
        
        # Stack vertically: rows = [RNA features; Ribo features], cols = biological samples
        concat_mat <- rbind(rna_sub, ribo_sub)
        
        # Build metadata for the biological samples (use RNA metadata as base)
        meta_all <- meta_rna[rna_cols_matched, , drop=FALSE]
        rownames(meta_all) <- shared_ids
        # Override seq_type to "Both" since these are combined
        meta_all$seq_type <- factor("Both")
        
        # Build DGE and normalize
        dge_all <- DGEList(counts=concat_mat, samples=data.frame(meta_all))
        dge_all <- calcNormFactors(dge_all, method="TMM")
        
        eval_MDS(dge_all, meta_all, qc_cols, opt$outdir, "All", opt$plot_ids)
        eval_PCA(dge_all, meta_all, qc_cols, opt$outdir, "All", opt$plot_ids)
        log_cpm_all <- cpm(dge_all, log=TRUE, normalized.lib.sizes=TRUE)
        eval_heatmap(log_cpm_all, meta_all, qc_cols, opt$outdir, "All")
        # GeneClusters_All skipped — concatenated feature space not meaningful for top-variance gene heatmap
        cat("Joint All QC plots complete.\n", file=log_file, append=TRUE)
    } else {
        cat("  Fewer than 3 shared biological samples — skipping All plots.\n", file=log_file, append=TRUE)
    }
}


# -----------------------------------------------------------------------------
# 6. FILTER CONTRASTS & FIT MODELS
# -----------------------------------------------------------------------------
has_reps <- function(grp_name, meta, col_name, min_n=2) {
    num_dots <- str_count(grp_name, "\\.") + 1
    target_col <- paste0(col_name, ".", num_dots)
    
    if (target_col %in% colnames(meta)) {
        # Check if exact match exists
        is_exact <- grp_name %in% meta[[target_col]]
        # Check if column is a combined factor
        is_compound <- any(grepl("_and_", meta[[target_col]]))
        
        if (!is_exact && is_compound && !grepl("_and_", grp_name)) {
             # Orthogonal contrast case: Filter by component
             split_entries <- strsplit(as.character(meta[[target_col]]), "_and_")
             mask <- sapply(split_entries, function(parts) grp_name %in% parts)
             sub_meta <- meta[mask, ]
        } else {
             # Standard case: Exact match
             sub_meta <- meta[meta[[target_col]] == grp_name, ]
        }

        # Check counts PER SEQ TYPE. 
        # Both RNA and Ribo must have >= 2 for the model to work properly for dTE
        counts_by_type <- table(sub_meta$seq_type)
        if (length(counts_by_type) > 0 && all(counts_by_type >= min_n)) {
            return(TRUE)
        }
    }
    return(FALSE)
}

cat("Filtering contrast lists for sufficient replicates...\n", file=log_file, append=TRUE)

# Filter Combinations
valid_comb <- list()
for (tup in comb_dict[[contrast_col]]) {
    if (has_reps(tup[1], meta.design, contrast_col) && has_reps(tup[2], meta.design, contrast_col)) {
        valid_comb <- append(valid_comb, list(tup))
    } else {
        cat(paste0("Skipping contrast '", paste(tup, collapse=" vs "), "' (insufficient replicates)\n"), file=log_file, append=TRUE)
    }
}
comb_dict[[contrast_col]] <- valid_comb

# Filter Uniques
valid_uniq <- list()
for (val in uniq_dict[[contrast_col]]) {
    if (has_reps(val, meta.design, contrast_col)) {
        valid_uniq <- append(valid_uniq, list(val))
    }
}
uniq_dict[[contrast_col]] <- valid_uniq

# STOP Check
num_valid_contrasts <- length(comb_dict[[contrast_col]]) + length(uniq_dict[[contrast_col]])
if (num_valid_contrasts == 0) {
    cat("\n----------------------------------------------------------------\n", file=log_file, append=TRUE)
    cat("NO VALID CONTRASTS FOUND: All requested groups have insufficient replicates (n < 2 per SeqType).\n", file=log_file, append=TRUE)
    cat("QC plots have been generated. Exiting successfully.\n", file=log_file, append=TRUE)
    cat("----------------------------------------------------------------\n", file=log_file, append=TRUE)
    quit(save="no", status=0)
}

# Fit Models
cat("Valid contrasts found. Fitting GLM models...\n", file=log_file, append=TRUE)
dge <- estimateDisp(dge, design, min.row.sum=30)
fit_paired <- glmQLFit(dge, design)

# Fit Independent RNA Model
fit_rna <- NULL
meta.rna <- meta.design[meta.design$seq_type == "RNA", , drop=FALSE]
if (nrow(meta.rna) > 0) {
    meta.rna <- droplevels(meta.rna)
    f_rna <- "~0 + group"
    if (!opt$test_run && "batch_date" %in% colnames(meta.rna) && has_var(meta.rna$batch_date)) f_rna <- paste0(f_rna, " + batch_date")
    if (!opt$test_run && "sample_type" %in% colnames(meta.rna) && has_var(meta.rna$sample_type)) f_rna <- paste0(f_rna, " + sample_type")
    
    # Use rna_counts_original (with transcript IDs) if available (translon mode)
    # This ensures RNA contrasts report transcript IDs (ST_.../ENST...) not TM_...
    rna_mat_for_fit <- if (exists("rna_counts_original")) {
        cat("Using rna_counts_original (transcript IDs) for independent RNA model.\n", file=log_file, append=TRUE)
        rna_counts_original
    } else {
        counts[, rownames(meta.rna), drop=FALSE]
    }
    # Strip _RNA suffix from colnames to match meta.rna rownames (counts_col)
    rna_col_ids <- rownames(meta.rna)
    dge_rna <- DGEList(counts=rna_mat_for_fit[, rna_col_ids, drop=FALSE], samples=data.frame(meta.rna))
    design.rna <- model.matrix(as.formula(f_rna), data=data.frame(meta.rna))
    dge_rna <- estimateDisp(dge_rna, design.rna)
    fit_rna <- glmQLFit(dge_rna, design.rna)
}

# --- SAVE MODEL LOGIC ---
if (!is.null(opt$save_model)) {
    cat(paste0("Saving model to ", opt$save_model, "...\n"), file=log_file, append=TRUE)
    # List of objects to save to ensure the environment is reproducible
    # We save everything relevant for the contrast execution phase
    save(
        list = c("fit_paired", "fit_rna", "dge", "meta.samples", "design", "design.rna",
                 "contrast_col", "uniq_dict", "comb_dict", "qc_types", "opt", "log_file", 
                 "tx.table", "feature2name", "contrast_grps", "seq_types"), 
        file = opt$save_model
    )
    cat("Model saved successfully.\n", file=log_file, append=TRUE)
}

} # End of if (!is.null(opt$load_model)) else block

# Resume Execution
dge_idxs <- match(rownames(dge$samples), meta.samples$counts_col)
dge.meta <- meta.samples[dge_idxs,]

# -----------------------------------------------------------------------------
# 6. EXECUTE CONTRASTS
# -----------------------------------------------------------------------------
cat("Executing contrasts...\n", file=log_file, append=TRUE)

# Execute Unique Contrasts (TE mostly, uses Paired Fit)
bplapply(uniq_dict[[contrast_col]], function(val) {
    evaluate_unique_contrast(dge$samples, val, contrast_col, fit_paired, opt$outdir, log_file)
})

# Execute Combination Contrasts (DE & dTE, uses Paired + RNA Independent Fit)
# Execute Combination Contrasts (DE & dTE, uses Paired + RNA Independent Fit)
if (!opt$skip_pairwise) {
    bplapply(comb_dict[[contrast_col]], function(val) {
        evaluate_combination_contrast(dge$samples, val, contrast_col, fit_paired, fit_rna, opt$outdir, log_file)
    })
} else {
    cat("Skipping pairwise contrasts (--skip_pairwise is TRUE).\n", file=log_file, append=TRUE)
}

# -----------------------------------------------------------------------------
# 7. EXECUTE ONE-VS-ALL (Optional)
# -----------------------------------------------------------------------------
if (!opt$skip_one_vs_all) {
    cat("\nExecuting One-vs-All Contrasts...\n", file=log_file, append=TRUE)
    
    # We iterate over every valid unique group found
    # (These have already been filtered for replicates in Step 6)
    bplapply(uniq_dict[[contrast_col]], function(val) {
        # Check if we should run OvA for this group
        # (It must define a subset, i.e., not encompass the entire dataset)
        is_subset <- sum(startsWith(as.character(dge.meta[[contrast_col]]), val)) < nrow(dge.meta)
        
        if (is_subset) {
            exec_one_vs_all(
                dge.meta, 
                target_group = val, 
                contrast_col = contrast_col, 
                fit_paired = fit_paired, 
                fit_rna = fit_rna, 
                has_ribo_data = (length(qc_types$Ribo) > 0), # Check if Ribo exists
                outdir = opt$outdir, 
                log_file = log_file
            )
        }
    })
}

cat("Analysis complete.\n", file=log_file, append=TRUE)
