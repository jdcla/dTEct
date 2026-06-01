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
  if (requireNamespace("svglite", quietly = TRUE)) {
    library(svglite)
  }
})
options(show.error.locations = TRUE)

set_thread_guards <- function(per_process_threads = 1L) {
  per_process_threads <- max(1L, as.integer(per_process_threads))
  thread_vars <- c(
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "RCPP_PARALLEL_NUM_THREADS"
  )

  for (var in thread_vars) {
    do.call(Sys.setenv, as.list(stats::setNames(as.character(per_process_threads), var)))
  }

  invisible(thread_vars)
}

make_bpparam <- function(workers) {
  workers <- max(1L, as.integer(workers))
  BiocParallel::MulticoreParam(
    workers = workers,
    progressbar = FALSE,
    stop.on.error = TRUE,
    tasks = workers
  )
}

run_bplapply <- function(X, FUN, workers, ...) {
  if (length(X) == 0) return(list())
  bpworkers <- min(max(1L, as.integer(workers)), length(X))
  BiocParallel::bplapply(X, FUN, ..., BPPARAM = make_bpparam(bpworkers))
}


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
  if (isTRUE(opt$use_anota2) && (is.null(opt$rna_counts) || is.null(opt$ribo_counts))) {
    errors <- c(errors, "anota2seq mode requires both --rna_counts and --ribo_counts.")
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

calc_assay_norm_factors <- function(dge, assay_col = "seq_type", method = "TMM", log_file = NULL) {
    if (!assay_col %in% colnames(dge$samples)) {
        if (!is.null(log_file)) {
            cat("Assay column not found; falling back to global TMM normalization.\n", file=log_file, append=TRUE)
        }
        return(calcNormFactors(dge, method=method))
    }

    dge$samples$norm.factors <- 1
    assays <- sort(unique(as.character(dge$samples[[assay_col]])))
    if (!is.null(log_file)) {
        cat("Using assay-specific TMM normalization by ", assay_col, ".\n", sep="", file=log_file, append=TRUE)
    }

    for (assay in assays) {
        idx <- as.character(dge$samples[[assay_col]]) == assay
        if (sum(idx) < 2) {
            if (!is.null(log_file)) {
                cat("  Skipping TMM for assay '", assay, "' (fewer than 2 samples).\n", sep="", file=log_file, append=TRUE)
            }
            next
        }

        sub_dge <- dge[, idx, keep.lib.sizes=FALSE]
        sub_dge <- calcNormFactors(sub_dge, method=method)
        dge$samples$norm.factors[idx] <- sub_dge$samples$norm.factors
        if (!is.null(log_file)) {
            cat("  Normalized assay '", assay, "' with ", sum(idx), " samples.\n", sep="", file=log_file, append=TRUE)
        }
    }

    dge
}

build_combination_contrast_args <- function(meta, tup, contrast_col) {
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

  return(list(grp_rna=grp_rna, grp_ribo=grp_ribo, tup=tup, n_rna=n_rna, n_ribo=n_ribo, has_ribo_data=has_ribo_data))
}

evaluate_combination_contrast <- function(meta, tup, contrast_col, fit_paired, fit_rna, outdir, log_file) {
  args <- build_combination_contrast_args(meta, tup, contrast_col)
  # PASS has_ribo_data flag
  evaluate_contrasts(args$grp_rna, args$grp_ribo, args$tup, args$n_rna, args$n_ribo, fit_paired, fit_rna, args$has_ribo_data, outdir, log_file)
}

build_one_vs_all_contrast_args <- function(meta, target_group, contrast_col, has_ribo_data) {
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

    return(list(grp_rna=grp_rna, grp_ribo=grp_ribo, tup=tup_fake, n_rna=n_rna, n_ribo=n_ribo, has_ribo_data=has_ribo_data))
}

exec_one_vs_all <- function(meta, target_group, contrast_col, fit_paired, fit_rna, has_ribo_data, outdir, log_file) {
    args <- build_one_vs_all_contrast_args(meta, target_group, contrast_col, has_ribo_data)
    if (is.null(args)) return(NULL)

    # Pass to existing function
    # It will handle the file creation, fitting, and plotting
    evaluate_contrasts(args$grp_rna, args$grp_ribo, args$tup, args$n_rna, args$n_ribo, fit_paired, fit_rna, args$has_ribo_data, outdir, log_file)
}

evaluate_contrasts <- function(grp_rna, grp_ribo, tup, n_rna, n_ribo, fit_paired, fit_rna, has_ribo_data, outdir, log_file) {
  strat_string <- ""
  
  # --- 1. RNA Contrasts ---
  if (length(unlist(grp_rna)) == 2) {
    n_msg <- paste0("    (n = ", n_rna[[1]], " / ", n_rna[[2]], ")")

    if (has_ribo_data) {
        # --- CASE A: DUAL MODE (Ribo exists) ---
        contrast_id_sub <- paste0(tup[[1]], "__", tup[[2]], strat_string, "_RNA")
        out_prefix_sub <- paste0(outdir, "RNA/", contrast_id_sub)
        if (isTRUE(opt$use_anota2)) {
            title_sub <- paste0(contrast_id_sub, n_msg, " (anota2seq total mRNA)")
            eval_anota2seq_contrast(fit_anota2seq, grp_rna[[1]], grp_rna[[2]], "RNA", "total mRNA", out_prefix_sub, title_sub, log_file)
        } else {
            # 1. Paired/Shared Model -> Suffix "_RNA"
            grp_contrast_paired <- paste0("makeContrasts(", paste(grp_rna, collapse = " - "), ", levels=design)")
            title_sub <- paste0(contrast_id_sub, n_msg, " (Shared)")

            if (!is.null(grp_contrast_paired)) print(paste0("Evaluating RNA (Shared): ", grp_contrast_paired))
            eval_contrast(fit_paired, grp_contrast_paired, out_prefix_sub, title_sub, log_file, remap_to_transcript = TRUE)

            # 2. Independent Model -> Suffix "_RNA_full"
            grp_contrast_full <- paste0("makeContrasts(", paste(grp_rna, collapse = " - "), ", levels=design.rna)")
            contrast_id_full <- paste0(tup[[1]], "__", tup[[2]], strat_string, "_RNA_full")
            out_prefix_full <- paste0(outdir, "RNA/", contrast_id_full)
            title_full <- paste0(contrast_id_full, n_msg, " (All)")

            eval_contrast(fit_rna, grp_contrast_full, out_prefix_full, title_full, log_file)
        }
        
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

    contrast_id <- paste0(tup[[1]], "__", tup[[2]], strat_string, "_Ribo")
    out_prefix <- paste0(outdir, "Ribo/", contrast_id)
    if (isTRUE(opt$use_anota2)) {
        title <- paste0(contrast_id, n_msg, " (anota2seq translated mRNA)")
        eval_anota2seq_contrast(fit_anota2seq, grp_ribo[[1]], grp_ribo[[2]], "Ribo", "translated mRNA", out_prefix, title, log_file)
    } else {
        grp_contrast <- paste0("makeContrasts(", paste(grp_ribo, collapse = " - "), ", levels=design)")
        title <- paste0(contrast_id, n_msg, " (Shared/Paired)")

        if (!is.null(grp_contrast)) print(paste0("Evaluating Ribo: ", grp_contrast))
        eval_contrast(fit_paired, grp_contrast, out_prefix, title, log_file)
    }
  }
  
  # --- 3. dTE Contrast ---
  if ((length(unlist(grp_ribo)) == 2) && (length(unlist(grp_rna)) == 2)) {
    contrast_id <- paste0(tup[[1]], "__", tup[[2]], strat_string, "_dTE")
    out_prefix <- paste0(outdir, "dTE/", contrast_id)
    title <- paste0(contrast_id, "    (n = ", n_ribo[[1]], " / ", n_rna[[1]], " / ", n_ribo[[2]], " / ", n_rna[[2]], ")")

    if (isTRUE(opt$use_anota2)) {
        eval_anota2seq_contrast(fit_anota2seq, grp_ribo[[1]], grp_ribo[[2]], "Ribo", "translation", out_prefix, title, log_file)
    } else {
        left_cmd <- paste0("(", grp_ribo[[1]], " - ", grp_rna[[1]], ")")
        right_cmd <- paste0("(", grp_ribo[[2]], " - ", grp_rna[[2]], ")")
        grp_contrast <- paste0("makeContrasts(", left_cmd, " - ", right_cmd, ", levels=design)")

        if (!is.null(grp_contrast)) print(grp_contrast)
        eval_contrast(fit_paired, grp_contrast, out_prefix, title, log_file)
    }
  }
}

collect_anota2seq_jobs_from_args <- function(args, outdir) {
  jobs <- list()
  if (is.null(args)) return(jobs)

  grp_rna <- args$grp_rna
  grp_ribo <- args$grp_ribo
  tup <- args$tup
  n_rna <- args$n_rna
  n_ribo <- args$n_ribo
  has_ribo_data <- args$has_ribo_data
  strat_string <- ""

  if (has_ribo_data && length(unlist(grp_rna)) == 2) {
    n_msg <- paste0("    (n = ", n_rna[[1]], " / ", n_rna[[2]], ")")
    contrast_id <- paste0(tup[[1]], "__", tup[[2]], strat_string, "_RNA")
    jobs[[length(jobs) + 1]] <- list(
        analysis = "total mRNA",
        seq_type = "RNA",
        left_group = grp_rna[[1]],
        right_group = grp_rna[[2]],
        out_prefix = paste0(outdir, "RNA/", contrast_id),
        title = paste0(contrast_id, n_msg, " (anota2seq total mRNA)")
    )
  }

  if (length(unlist(grp_ribo)) == 2) {
    n_msg <- paste0("    (n = ", n_ribo[[1]], " / ", n_ribo[[2]], ")")
    contrast_id <- paste0(tup[[1]], "__", tup[[2]], strat_string, "_Ribo")
    jobs[[length(jobs) + 1]] <- list(
        analysis = "translated mRNA",
        seq_type = "Ribo",
        left_group = grp_ribo[[1]],
        right_group = grp_ribo[[2]],
        out_prefix = paste0(outdir, "Ribo/", contrast_id),
        title = paste0(contrast_id, n_msg, " (anota2seq translated mRNA)")
    )
  }

  if ((length(unlist(grp_ribo)) == 2) && (length(unlist(grp_rna)) == 2)) {
    contrast_id <- paste0(tup[[1]], "__", tup[[2]], strat_string, "_dTE")
    jobs[[length(jobs) + 1]] <- list(
        analysis = "translation",
        seq_type = "Ribo",
        left_group = grp_ribo[[1]],
        right_group = grp_ribo[[2]],
        out_prefix = paste0(outdir, "dTE/", contrast_id),
        title = paste0(contrast_id, "    (n = ", n_ribo[[1]], " / ", n_rna[[1]], " / ", n_ribo[[2]], " / ", n_rna[[2]], ")")
    )
  }

  return(jobs)
}

collect_anota2seq_jobs <- function(meta, dge_meta, contrast_col, combs, uniqs, outdir, skip_pairwise, skip_one_vs_all, log_file) {
  jobs <- list()

  if (!skip_pairwise) {
    for (tup in combs) {
      args <- build_combination_contrast_args(meta, tup, contrast_col)
      jobs <- c(jobs, collect_anota2seq_jobs_from_args(args, outdir))
    }
  }

  if (!skip_one_vs_all) {
    has_ribo_data <- any(dge_meta$seq_type == "Ribo")
    for (val in uniqs) {
      is_subset <- sum(startsWith(as.character(dge_meta[[contrast_col]]), val)) < nrow(dge_meta)
      if (is_subset) {
        args <- build_one_vs_all_contrast_args(dge_meta, val, contrast_col, has_ribo_data)
        jobs <- c(jobs, collect_anota2seq_jobs_from_args(args, outdir))
      }
    }
  }

  cat(sprintf("Collected %d anota2seq contrast outputs for batched evaluation.\n", length(jobs)), file=log_file, append=TRUE)
  return(jobs)
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
    # Group names must match the design/model column naming used elsewhere in
    # the script. Multi-level metadata groups are joined with _and_ before the
    # seq_type suffix, e.g. groupMBL_and_MBL.SHH__RNA.
    grp_strings[i] <- paste0(weight_dict[[key]], "*group", paste0(val_grps, collapse = "_and_"), "__", seq_type)
  }
  grp_string <- paste0(grp_strings, collapse=" + ")
  return(paste0("(", grp_string, ")"))
}

parse_anota2seq_group_weights <- function(group_string, seq_type) {
    # edgeR contrasts are represented as strings such as:
    #   (0.5*groupA__RNA + 0.5*groupB__RNA)
    # anota2seq needs a numeric matrix indexed by phenotype, so this parser
    # extracts the group names and their weights from those existing strings.
    pattern <- paste0("([0-9eE.+-]+)\\*group([^ +()]+)__", seq_type)
    matches <- gregexpr(pattern, group_string, perl=TRUE)
    tokens <- regmatches(group_string, matches)[[1]]

    if (length(tokens) == 0 || tokens[1] == "") {
        stop(paste0("Could not parse anota2seq contrast group: ", group_string))
    }

    weights <- c()
    for (token in tokens) {
        parts <- regmatches(token, regexec(pattern, token, perl=TRUE))[[1]]
        group_name <- parts[3]
        weight <- as.numeric(parts[2])
        weights[group_name] <- ifelse(group_name %in% names(weights), weights[group_name] + weight, weight)
    }
    return(weights)
}

normalize_anota2seq_weights <- function(weights, pheno_levels) {
    # anota2seq is built from paired RNA/Ribo samples only. During that setup we
    # may drop phenotype classes that do not have enough paired samples, but the
    # contrast strings were generated earlier from the full metadata table.
    #
    # This matters most for one-vs-all contrasts: the "OTHER" side can contain
    # many sparse/unpaired classes. Those classes are not present in the paired
    # anota2seq object, so they must be removed from the contrast. After removal,
    # the remaining side is renormalized to sum to 1 so the estimand remains:
    #   average(retained target classes) - average(retained comparison classes)
    weights <- weights[names(weights) %in% pheno_levels]
    if (length(weights) == 0) return(weights)

    weight_sum <- sum(weights)
    if (weight_sum == 0) return(weights)
    weights / weight_sum
}

force_anota2seq_zero_sum <- function(contrast) {
    # anota2seq requires every contrast column to sum to exactly zero, not just
    # numerically close to zero. Weighted one-vs-all contrasts often produce tiny
    # floating-point residuals such as 1e-17 after renormalization. Those are
    # mathematically zero, but anota2seq rejects them during input validation.
    #
    # Try adjusting the largest coefficients first, then fall back to all rows,
    # and keep only an adjustment that makes R's own colSums() exactly zero.
    # The correction is at machine precision and does not change the contrast in
    # any meaningful statistical sense.
    for (j in seq_len(ncol(contrast))) {
        col_sum <- colSums(contrast[,j,drop=FALSE])[[1]]
        if (col_sum == 0) next

        candidates <- unique(c(order(abs(contrast[,j]), decreasing=TRUE), seq_len(nrow(contrast))))
        for (idx in candidates) {
            adjusted <- contrast
            adjusted[idx,j] <- adjusted[idx,j] - colSums(adjusted[,j,drop=FALSE])[[1]]
            if (colSums(adjusted[,j,drop=FALSE])[[1]] == 0) {
                contrast <- adjusted
                break
            }
        }
    }
    contrast
}

make_anota2seq_base_contrast <- function(left_group, right_group, seq_type, pheno_levels) {
    left_weights <- parse_anota2seq_group_weights(left_group, seq_type)
    right_weights <- parse_anota2seq_group_weights(right_group, seq_type)

    # Project the metadata-derived contrast onto the phenotype classes that are
    # actually present in the paired anota2seq dataset. If either side disappears
    # completely, the contrast is not estimable and should be skipped.
    left_weights <- normalize_anota2seq_weights(left_weights, pheno_levels)
    right_weights <- normalize_anota2seq_weights(right_weights, pheno_levels)
    if (length(left_weights) == 0 || length(right_weights) == 0) {
        stop("anota2seq contrast has no paired samples for one side after filtering")
    }

    contrast <- matrix(0, nrow=length(pheno_levels), ncol=1)
    rownames(contrast) <- pheno_levels
    colnames(contrast) <- "contrast1"
    contrast[names(left_weights), 1] <- left_weights
    contrast[names(right_weights), 1] <- contrast[names(right_weights), 1] - right_weights
    return(force_anota2seq_zero_sum(contrast))
}

fill_anota2seq_contrast_matrix <- function(contrast, pheno_levels) {
    # anota2seq expects a full-rank contrast matrix with n_phenotypes - 1 columns.
    # Requested contrasts are placed first; filler columns are added solely to
    # satisfy the package's model-fitting requirements.
    needed_cols <- length(pheno_levels) - 1
    requested_cols <- ncol(contrast)
    if (needed_cols > 1) {
        current <- contrast
        group_pairs <- combn(pheno_levels, 2, simplify=FALSE)
        for (pair in group_pairs) {
            if (ncol(current) >= needed_cols) break
            filler_col <- matrix(0, nrow=length(pheno_levels), ncol=1)
            rownames(filler_col) <- pheno_levels
            filler_col[pair[1], 1] <- 1
            filler_col[pair[2], 1] <- -1
            trial <- cbind(current, filler_col)
            if (qr(trial)$rank > qr(current)$rank) {
                current <- trial
            }
        }
        if (ncol(current) < needed_cols || qr(current)$rank < needed_cols) {
            stop("could not construct full-rank anota2seq contrast matrix")
        }
        contrast <- current
        n_filler <- ncol(contrast) - requested_cols
        colnames(contrast) <- c(paste0("contrast", seq_len(requested_cols)), paste0("filler", seq_len(n_filler)))
    }
    return(force_anota2seq_zero_sum(contrast))
}

make_anota2seq_contrast_matrix <- function(left_group, right_group, seq_type, pheno_levels) {
    contrast <- make_anota2seq_base_contrast(left_group, right_group, seq_type, pheno_levels)
    return(fill_anota2seq_contrast_matrix(contrast, pheno_levels))
}

anota2seq_contrast_key <- function(analysis, contrast) {
    paste(analysis, paste(rownames(contrast), signif(as.numeric(contrast[,1]), 12), sep="=", collapse=";"), sep="|")
}

# anota2seq is deliberately handled more conservatively than edgeR below. In
# this container, forked R workers can fail while lazy-loading packages from the
# shared filesystem (for example cli/sysdata.rdb during plotting). Therefore the
# anota2seq model calls and final CSV/volcano writing run serially, while results
# are persisted to disk as soon as each batch completes.
format_elapsed <- function(start_time) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units="secs"))
    sprintf("%.1f min", elapsed / 60)
}

safe_file_label <- function(x) {
    gsub("[^A-Za-z0-9_.-]+", "_", x)
}

write_anota2seq_contrast_matrix <- function(contrast_matrix, file_path) {
    out <- data.frame(phenotype = rownames(contrast_matrix), contrast_matrix, check.names = FALSE)
    write.csv(out, file_path, row.names=FALSE)
}

write_anota2seq_batch_contrasts <- function(batch, analysis_name, batch_index, contrast_dir, log_file) {
    if (is.null(contrast_dir) || contrast_dir == "") return(invisible(NULL))

    dir.create(contrast_dir, recursive=TRUE, showWarnings=FALSE)
    analysis_label <- safe_file_label(analysis_name)
    batch_prefix <- file.path(contrast_dir, sprintf("%s_batch_%03d", analysis_label, batch_index))

    full_matrix_path <- paste0(batch_prefix, "_full_matrix.csv")
    write_anota2seq_contrast_matrix(batch$matrix, full_matrix_path)

    metadata <- data.frame(
        batch = batch_index,
        analysis = analysis_name,
        requested_contrast = seq_along(batch$jobs),
        title = vapply(batch$jobs, function(x) x$title, character(1)),
        out_prefix = vapply(batch$jobs, function(x) x$out_prefix, character(1)),
        key = vapply(batch$jobs, function(x) x$key, character(1)),
        stringsAsFactors = FALSE
    )
    write.csv(metadata, paste0(batch_prefix, "_metadata.csv"), row.names=FALSE)

    for (i in seq_along(batch$jobs)) {
        requested_path <- paste0(batch_prefix, sprintf("_requested_%03d.csv", i))
        write_anota2seq_contrast_matrix(batch$jobs[[i]]$contrast, requested_path)
    }

    cat(sprintf("Saved anota2seq contrast matrices for %s batch %d to %s\n", analysis_name, batch_index, contrast_dir), file=log_file, append=TRUE)
    invisible(NULL)
}

anota2seq_cache_key_file <- function(cache_dir, key) {
    file.path(cache_dir, paste0(digest::digest(key, algo="xxhash64"), ".rds"))
}

set_anota2seq_cache_result <- function(key, result) {
    cache <- anota2seq_contrast_cache
    cache$results[[key]] <- result
    anota2seq_contrast_cache <<- cache
    invisible(NULL)
}

load_anota2seq_disk_cache <- function(cache_dir, valid_jobs, log_file) {
    if (is.null(cache_dir) || cache_dir == "") return(0L)
    if (!dir.exists(cache_dir)) return(0L)

    loaded <- 0L
    for (job in valid_jobs) {
        cache_file <- anota2seq_cache_key_file(cache_dir, job$key)
        if (!file.exists(cache_file)) next
        set_anota2seq_cache_result(job$key, readRDS(cache_file))
        loaded <- loaded + 1L
    }

    if (loaded > 0) {
        cat(sprintf("Loaded %d anota2seq contrast results from disk cache %s.\n", loaded, cache_dir), file=log_file, append=TRUE)
    }
    loaded
}

save_anota2seq_disk_cache_result <- function(cache_dir, key, result) {
    if (is.null(cache_dir) || cache_dir == "") return(invisible(NULL))
    dir.create(cache_dir, recursive=TRUE, showWarnings=FALSE)
    saveRDS(result, anota2seq_cache_key_file(cache_dir, key))
    invisible(NULL)
}

evaluate_anota2seq_batch <- function(ads, batch, analysis_name, batch_index, total_batches, log_file) {
    start_time <- Sys.time()
    cat(sprintf(
        "Starting anota2seq %s batch %d/%d: %d requested outputs; contrast matrix %d x %d.\n",
        analysis_name, batch_index, total_batches, length(batch$jobs), nrow(batch$matrix), ncol(batch$matrix)
    ), file=log_file, append=TRUE)
    cat("  Note: only requested contrast columns are reported; filler columns are included only to satisfy anota2seq full-rank requirements.\n", file=log_file, append=TRUE)

    suppressMessages(suppressWarnings(capture.output({
        ads_res <- anota2seq::anota2seqAnalyze(
            Anota2seqDataSet = ads,
            contrasts = batch$matrix,
            analysis = analysis_name,
            useProgBar = FALSE,
            fileStem = tempfile("anota2seq_", tmpdir=tempdir())
        )
    }, file = NULL)))

    cat(sprintf(
        "Finished anota2seqAnalyze for %s batch %d/%d in %s; extracting %d outputs.\n",
        analysis_name, batch_index, total_batches, format_elapsed(start_time), length(batch$jobs)
    ), file=log_file, append=TRUE)

    results <- list()
    for (i in seq_along(batch$jobs)) {
        res <- as.data.frame(anota2seq::anota2seqGetOutput(
            object = ads_res,
            analysis = analysis_name,
            output = "full",
            selContrast = i,
            getRVM = TRUE
        ))
        results[[batch$jobs[[i]]$key]] <- res
        cat(sprintf(
            "  Extracted anota2seq %s batch %d/%d output %d/%d: %s\n",
            analysis_name, batch_index, total_batches, i, length(batch$jobs), batch$jobs[[i]]$title
        ), file=log_file, append=TRUE)
    }

    cat(sprintf(
        "Completed anota2seq %s batch %d/%d in %s.\n",
        analysis_name, batch_index, total_batches, format_elapsed(start_time)
    ), file=log_file, append=TRUE)
    results
}

finalize_contrast_result <- function(res, out_prefix, title, remap_to_transcript = FALSE) {
    # For RNA contrasts from the paired model, remap row_id from translon to transcript IDs
    # This ensures RNA outputs report TT_*/ENST... instead of ST_*/TM_...
    if (remap_to_transcript && exists("translon_to_transcript_map")) {
        res$row_id <- ifelse(res$row_id %in% names(translon_to_transcript_map),
                             translon_to_transcript_map[res$row_id],
                             res$row_id)
    }
    
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
    # Add transcript_id column: for translon-mode RNA outputs, maps ST_* -> TT_*
    # For RNA-only or gene-mode, transcript_id == row_id
    res$transcript_id <- if (exists("translon_to_transcript_map") && any(res$row_id %in% names(translon_to_transcript_map))) {
        ifelse(res$row_id %in% names(translon_to_transcript_map),
               translon_to_transcript_map[res$row_id],
               res$row_id)
    } else {
        res$row_id
    }

    out_table <- res |> 
        select(gene_id, gene_name, row_id, transcript_id, logFC, logCPM, any_of("F"), PValue, FDR) |> 
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

get_edgeR_contrast_cache <- function() {
    if (!exists("edgeR_contrast_cache", envir=.GlobalEnv, inherits=FALSE)) {
        edgeR_contrast_cache <<- new.env(parent=emptyenv())
    }
    edgeR_contrast_cache
}

make_edgeR_contrast_key <- function(fit, contrast_matrix) {
    contrast_matrix <- as.matrix(contrast_matrix)
    contrast_vals <- signif(as.numeric(contrast_matrix[,1]), 12)
    contrast_names <- rownames(contrast_matrix)
    if (is.null(contrast_names)) contrast_names <- names(contrast_vals)
    fit_cols <- colnames(fit$coefficients)
    if (is.null(fit_cols) && !is.null(fit$design)) fit_cols <- colnames(fit$design)

    paste(
        paste(dim(fit$coefficients), collapse="x"),
        paste(fit_cols, collapse=";"),
        paste(contrast_names, contrast_vals, sep="=", collapse=";"),
        sep="|"
    )
}

set_edgeR_contrast_cache_result <- function(cache, cache_key, res, max_entries = 64) {
    cache_keys <- attr(cache, "keys")
    if (is.null(cache_keys)) cache_keys <- character(0)

    if (!cache_key %in% cache_keys) {
        while (length(cache_keys) >= max_entries) {
            rm(list=cache_keys[[1]], envir=cache)
            cache_keys <- cache_keys[-1]
        }
        cache_keys <- c(cache_keys, cache_key)
    }

    assign(cache_key, res, envir=cache)
    attr(cache, "keys") <- cache_keys
}


eval_contrast <- function(fit, contrast, out_prefix, title, log_file, remap_to_transcript = FALSE) {
    cleaned_contrast <- sub("^makeContrasts\\(", "", sub(", levels=.*\\)$", "", contrast))
    cat("Evaluating contrast ", title, " : ", cleaned_contrast, "\n", file = log_file, append = TRUE)

    contrast_matrix <- eval(parse(text = contrast))
    cache <- get_edgeR_contrast_cache()
    cache_key <- make_edgeR_contrast_key(fit, contrast_matrix)

    if (exists(cache_key, envir=cache, inherits=FALSE)) {
        res <- get(cache_key, envir=cache, inherits=FALSE)
    } else {
        qlf <- glmQLFTest(fit, contrast=contrast_matrix)
        res <- topTags(qlf, n=Inf, sort.by="none")$table
        set_edgeR_contrast_cache_result(cache, cache_key, res)
    }

    res <- res |> tibble::rownames_to_column('row_id')
    finalize_contrast_result(res, out_prefix, title, remap_to_transcript)
}

prepare_anota2seq_contrast_cache <- function(ads, jobs, log_file) {
    anota2seq_contrast_cache <<- list(results=list())
    if (length(jobs) == 0) return(invisible(NULL))

    covariates <- anota2seq::anota2seqGetCovariates(ads)
    pheno_levels <- levels(factor(covariates$phenoVec))
    group_counts <- table(covariates$phenoVec)
    needed_cols <- length(pheno_levels) - 1
    valid_jobs <- list()
    seen_keys <- c()

    for (job in jobs) {
        contrast <- tryCatch(
            make_anota2seq_base_contrast(job$left_group, job$right_group, job$seq_type, pheno_levels),
            error = function(e) {
                cat("Skipping anota2seq contrast ", job$title, " (", conditionMessage(e), ")\n", file = log_file, append = TRUE)
                return(NULL)
            }
        )
        if (is.null(contrast)) next

        active_groups <- rownames(contrast)[contrast[,1] != 0]
        if (any(group_counts[active_groups] < 2)) {
            cat("Skipping anota2seq contrast ", job$title, " (insufficient paired samples)\n", file = log_file, append = TRUE)
            next
        }

        job$contrast <- contrast
        job$key <- anota2seq_contrast_key(job$analysis, contrast)
        if (job$key %in% seen_keys) next
        seen_keys <- c(seen_keys, job$key)
        valid_jobs[[length(valid_jobs) + 1]] <- job
    }

    if (length(valid_jobs) == 0) return(invisible(NULL))

    contrast_dir <- NULL
    cache_dir <- NULL
    if (exists("opt", envir=.GlobalEnv, inherits=FALSE) && !is.null(opt$outdir)) {
        contrast_dir <- file.path(opt$outdir, "anota2seq_contrasts")
        cache_dir <- file.path(opt$outdir, "anota2seq_cache")
    }
    cat(sprintf("Prepared %d valid anota2seq contrast outputs for evaluation.\n", length(valid_jobs)), file=log_file, append=TRUE)
    load_anota2seq_disk_cache(cache_dir, valid_jobs, log_file)

    for (analysis_name in unique(vapply(valid_jobs, function(x) x$analysis, character(1)))) {
        analysis_jobs <- valid_jobs[vapply(valid_jobs, function(x) x$analysis == analysis_name, logical(1))]
        batches <- list()
        current_jobs <- list()
        current_matrix <- NULL

        flush_batch <- function() {
            if (length(current_jobs) == 0) return(NULL)
            list(jobs=current_jobs, matrix=fill_anota2seq_contrast_matrix(current_matrix, pheno_levels))
        }

        for (job in analysis_jobs) {
            if (!is.null(anota2seq_contrast_cache$results[[job$key]])) next
            if (is.null(current_matrix)) {
                current_jobs <- list(job)
                current_matrix <- job$contrast
                next
            }

            trial <- cbind(current_matrix, job$contrast)
            if (ncol(current_matrix) < needed_cols && qr(trial)$rank > qr(current_matrix)$rank) {
                current_jobs[[length(current_jobs) + 1]] <- job
                current_matrix <- trial
            } else {
                batches[[length(batches) + 1]] <- flush_batch()
                current_jobs <- list(job)
                current_matrix <- job$contrast
            }
        }
        batches[[length(batches) + 1]] <- flush_batch()
        batches <- batches[!vapply(batches, is.null, logical(1))]
        if (length(batches) == 0) {
            cat(sprintf("All %d anota2seq %s requested outputs are already cached.\n", length(analysis_jobs), analysis_name), file=log_file, append=TRUE)
            next
        }

        cat(sprintf(
            "Running %d batched anota2seq %s jobs for %d requested outputs serially.\n",
            length(batches), analysis_name, length(analysis_jobs)
        ), file=log_file, append=TRUE)

        for (i in seq_along(batches)) {
            cat(sprintf(
                "Queued anota2seq %s batch %d/%d: %d requested outputs; contrast matrix %d x %d; outputs: %s\n",
                analysis_name, i, length(batches), length(batches[[i]]$jobs), nrow(batches[[i]]$matrix), ncol(batches[[i]]$matrix),
                paste(vapply(batches[[i]]$jobs, function(x) x$title, character(1)), collapse=" | ")
            ), file=log_file, append=TRUE)
            write_anota2seq_batch_contrasts(batches[[i]], analysis_name, i, contrast_dir, log_file)
        }

        completed_outputs <- 0L
        for (i in seq_along(batches)) {
            batch_result <- evaluate_anota2seq_batch(ads, batches[[i]], analysis_name, i, length(batches), log_file)
            for (key in names(batch_result)) {
                set_anota2seq_cache_result(key, batch_result[[key]])
                save_anota2seq_disk_cache_result(cache_dir, key, batch_result[[key]])
                completed_outputs <- completed_outputs + 1L
            }
            cat(sprintf(
                "Persisted anota2seq %s batch %d/%d to memory and disk cache (%d/%d outputs cached for this analysis).\n",
                analysis_name, i, length(batches), completed_outputs, length(analysis_jobs)
            ), file=log_file, append=TRUE)
        }
        cat(sprintf(
            "Cached %d/%d anota2seq %s requested outputs.\n",
            completed_outputs, length(analysis_jobs), analysis_name
        ), file=log_file, append=TRUE)
    }
    invisible(NULL)
}

eval_anota2seq_contrast <- function(ads, left_group, right_group, seq_type, analysis, out_prefix, title, log_file) {
    covariates <- anota2seq::anota2seqGetCovariates(ads)
    pheno_levels <- levels(factor(covariates$phenoVec))
    contrast <- tryCatch(
        make_anota2seq_base_contrast(left_group, right_group, seq_type, pheno_levels),
        error = function(e) {
            cat("Skipping anota2seq contrast ", title, " (", conditionMessage(e), ")\n", file = log_file, append = TRUE)
            return(NULL)
        }
    )
    if (is.null(contrast)) return(NULL)
    active_groups <- rownames(contrast)[contrast[,1] != 0]
    group_counts <- table(covariates$phenoVec)

    if (any(group_counts[active_groups] < 2)) {
        cat("Skipping anota2seq contrast ", title, " (insufficient paired samples)\n", file = log_file, append = TRUE)
        return(NULL)
    }

    key <- anota2seq_contrast_key(analysis, contrast)
    if (exists("anota2seq_contrast_cache") && !is.null(anota2seq_contrast_cache$results[[key]])) {
        res <- as.data.frame(anota2seq_contrast_cache$results[[key]])
    } else {
        cat("Evaluating anota2seq ", analysis, " contrast ", title, "\n", file = log_file, append = TRUE)
        suppressMessages(suppressWarnings(capture.output({
            ads_res <- anota2seq::anota2seqAnalyze(
                Anota2seqDataSet = ads,
                contrasts = fill_anota2seq_contrast_matrix(contrast, pheno_levels),
                analysis = analysis,
                useProgBar = FALSE,
                fileStem = tempfile("anota2seq_", tmpdir=tempdir())
            )
        }, file = NULL)))

        res <- as.data.frame(anota2seq::anota2seqGetOutput(
            object = ads_res,
            analysis = analysis,
            output = "full",
            selContrast = 1,
            getRVM = TRUE
        ))
    }
    res <- res |> tibble::rownames_to_column('row_id')
    res$logFC <- res$apvEff
    logcpm <- if (exists("anota2seq_logcpm") && analysis %in% names(anota2seq_logcpm)) anota2seq_logcpm[[analysis]] else NULL
    res$logCPM <- if (!is.null(logcpm)) logcpm[match(res$row_id, names(logcpm))] else NA_real_
    res$F <- res$apvRvmF
    res$PValue <- res$apvRvmP
    res$FDR <- res$apvRvmPAdj

    finalize_contrast_result(res, out_prefix, title)
}

build_anota2seq_dataset <- function(dge, log_file) {
    meta_rna <- dge$samples[dge$samples$seq_type == "RNA", , drop=FALSE]
    meta_ribo <- dge$samples[dge$samples$seq_type == "Ribo", , drop=FALSE]
    shared_ids <- intersect(as.character(meta_rna$smart_id), as.character(meta_ribo$smart_id))

    if (length(shared_ids) < 4) {
        stop("anota2seq mode requires at least four paired RNA/Ribo biological samples.")
    }

    rna_col_lookup <- setNames(rownames(meta_rna), as.character(meta_rna$smart_id))
    ribo_col_lookup <- setNames(rownames(meta_ribo), as.character(meta_ribo$smart_id))
    rna_cols <- rna_col_lookup[shared_ids]
    ribo_cols <- ribo_col_lookup[shared_ids]

    dataT <- as.matrix(dge$counts[, rna_cols, drop=FALSE])
    dataP <- as.matrix(dge$counts[, ribo_cols, drop=FALSE])
    colnames(dataT) <- shared_ids
    colnames(dataP) <- shared_ids

    phenoVec <- sub("__RNA$", "", as.character(meta_rna[rna_cols, "group"]))
    pheno_counts <- table(phenoVec)
    keep_pairs <- phenoVec %in% names(pheno_counts[pheno_counts >= 2])
    if (sum(!keep_pairs) > 0) {
        cat(sprintf("Dropping %d paired samples from anota2seq setup because their groups have n < 2.\n", sum(!keep_pairs)), file=log_file, append=TRUE)
        dataT <- dataT[, keep_pairs, drop=FALSE]
        dataP <- dataP[, keep_pairs, drop=FALSE]
        shared_ids <- shared_ids[keep_pairs]
        rna_cols <- rna_cols[keep_pairs]
        ribo_cols <- ribo_cols[keep_pairs]
        phenoVec <- phenoVec[keep_pairs]
    }
    if (length(unique(phenoVec)) < 2) {
        stop("anota2seq mode requires at least two paired sample classes with n >= 2.")
    }
    batchVec <- NULL
    if ("sample_type" %in% colnames(meta_rna) && length(unique(as.character(meta_rna[rna_cols, "sample_type"]))) > 1) {
        batchVec <- as.character(meta_rna[rna_cols, "sample_type"])
        cat("Using sample_type as anota2seq batchVec.\n", file=log_file, append=TRUE)
    } else if ("batch_date" %in% colnames(meta_rna) && "batch_date" %in% colnames(meta_ribo)) {
        rna_batch <- as.character(meta_rna[rna_cols, "batch_date"])
        ribo_batch <- as.character(meta_ribo[ribo_cols, "batch_date"])
        if (all(rna_batch == ribo_batch) && length(unique(rna_batch)) > 1) {
            batchVec <- rna_batch
            cat("Using batch_date as anota2seq batchVec.\n", file=log_file, append=TRUE)
        }
    }

    input_features <- nrow(dataP)
    cat(sprintf("Building anota2seq dataset from %d paired samples and %d input features.\n", length(shared_ids), input_features), file=log_file, append=TRUE)
    ads <- anota2seq::anota2seqDataSetFromMatrix(
        dataP = dataP,
        dataT = dataT,
        phenoVec = phenoVec,
        batchVec = batchVec,
        dataType = "RNAseq",
        normalize = TRUE,
        transformation = "TMM-log2",
        filterZeroGenes = TRUE,
        varCutOff = NULL
    )

    norm_data <- anota2seq::anota2seqGetNormalizedData(ads)
    retained_features <- nrow(norm_data$dataP)
    cat(sprintf(
        "anota2seq retained %d/%d features after zero-count filtering and normalization.\n",
        retained_features, input_features
    ), file=log_file, append=TRUE)
    anota2seq_logcpm <<- list(
        "translated mRNA" = rowMeans(norm_data$dataP, na.rm=TRUE),
        "total mRNA" = rowMeans(norm_data$dataT, na.rm=TRUE),
        "translation" = rowMeans(cbind(norm_data$dataP, norm_data$dataT), na.rm=TRUE)
    )
    return(ads)
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

parse_cli_bool <- function(x, option_name) {
  if (is.null(x)) return(NULL)
  value <- tolower(trimws(as.character(x)))
  if (value %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (value %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop(sprintf("Invalid value for %s: %s", option_name, x))
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
    make_option(c("-u", "--outer_join"), action="store_true", default=FALSE, help="Outer join RNA and Ribo reads, fill NAs with 0 expression."),
    make_option(c("-l", "--cores"         ), type="integer"  , default=1                    , metavar="integer", help="Number of cores."                                                                       ),
    make_option(c("-a", "--contrast_cols" ), type="character", default="treatment_id"               , metavar="character", help="Column names from which contrasts are derived; separated by commas."          ),
    make_option(c("-O", "--skip_one_vs_all"), action="store_true", default=FALSE, help="If set, skips One-vs-All contrasts (Group vs Rest). Default is to RUN them."),
    make_option(c("--one_vs_all"), type="character", default=NULL, metavar="TRUE/FALSE", help="Legacy compatibility flag. TRUE runs One-vs-All, FALSE skips it. Prefer --skip_one_vs_all to disable."),
    make_option(c("-p", "--plot_ids"), action="store_true", default=FALSE, help="Plot Smart IDs instead of replicate numbers in PCA/MDS."),
    make_option(c("-e", "--no_batch_factor"), action="store_true", default=FALSE, help="Don't create factors for batch date. Can be necessary to achieve full rank for some settings"),
    make_option(c("-S", "--save_model"), type="character", default=NULL, metavar="path", help="Path to save the fitted model (RData file). Defaults to 'dTEct_model.RData' in outdir if not specified. Ignored if --load_model is used."),
    make_option(c("-L", "--load_model"), type="character", default=NULL, metavar="path", help="Path to load a pre-fitted model from (skips fitting)."),
    make_option(c("-k", "--skip_pairwise"), action="store_true", default=FALSE, help="If set, skips all pairwise contrasts and only runs One-vs-All (if enabled)."),
    make_option(c("-T", "--test_run"), action="store_true", default=FALSE, help="Run in test mode: subsets data to first valid contrast and fraction of genes for rapid debugging."),
    make_option(c("-A", "--use_anota2"), action="store_true", default=FALSE, help="Use anota2seq for RNA, Ribo, and dTE contrasts. TE output is not generated in this mode.")
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

legacy_one_vs_all <- parse_cli_bool(opt$one_vs_all, "--one_vs_all")
if (!is.null(legacy_one_vs_all)) {
  opt$skip_one_vs_all <- !legacy_one_vs_all
}

# Prevent each worker from spawning its own BLAS/OpenMP thread team. The script
# uses process-level parallelism via BiocParallel; letting each process also use
# many math-library threads is what typically drives CPU usage far above the
# requested --cores value.
set_thread_guards(1L)

# Validate Inputs (Skip some checks if loading model)
if (is.null(opt$load_model)) {
  validate_inputs(opt)
} else {
  if (!file.exists(opt$load_model)) {
    stop(paste("Model file not found:", opt$load_model))
  }
}
if (isTRUE(opt$use_anota2) && !requireNamespace("anota2seq", quietly=TRUE)) {
  stop("Package 'anota2seq' is required for --use_anota2.")
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

# Configure a bounded BiocParallel backend once up front. We still pass an
# explicit BPPARAM at each parallel call so the effective worker count can be
# capped to the number of outstanding jobs.
register(make_bpparam(opt$cores))
thread_guard_vars <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "RCPP_PARALLEL_NUM_THREADS")
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
  runtime_opt <- opt
  cat(paste0("Loading model from ", runtime_opt$load_model, "...\n"))
  load(runtime_opt$load_model)
  loaded_use_anota2 <- opt$use_anota2
  if (is.null(loaded_use_anota2)) loaded_use_anota2 <- FALSE

  # A loaded model may contain the original CLI options used during fitting.
  # Preserve the runtime options from the current invocation for settings that
  # users reasonably expect to override on rerun, especially --cores and --outdir.
  opt$load_model <- runtime_opt$load_model
  opt$outdir <- runtime_opt$outdir
  opt$cores <- runtime_opt$cores
  opt$skip_pairwise <- runtime_opt$skip_pairwise
  opt$skip_one_vs_all <- runtime_opt$skip_one_vs_all
  opt$one_vs_all <- runtime_opt$one_vs_all
  opt$plot_ids <- runtime_opt$plot_ids
  opt$save_model <- NULL
  opt$use_anota2 <- loaded_use_anota2

  # Re-register cores using the runtime override, not the saved model value.
  register(make_bpparam(opt$cores))
  
  # Recreate output subdirectories because load_model skips the normal setup path.
  dir.create(paste0(opt$outdir, 'dTE'), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0(opt$outdir, 'Ribo'), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0(opt$outdir, 'RNA'), showWarnings = FALSE, recursive = TRUE)
  if (!isTRUE(opt$use_anota2)) {
    dir.create(paste0(opt$outdir, 'TE'), showWarnings = FALSE, recursive = TRUE)
  }

  # Re-open log file using the current output directory.
  log_file <- paste0(opt$outdir, "run_info.txt")
  cat("Model loaded. Skipping data processing and fitting.\n", file = log_file, append = TRUE)
  cat(sprintf("Using runtime overrides after model load: outdir=%s, cores=%d\n", opt$outdir, opt$cores), file = log_file, append = TRUE)
  cat(sprintf("Effective contrast flags after model load: skip_pairwise=%s, skip_one_vs_all=%s, use_anota2=%s\n", opt$skip_pairwise, opt$skip_one_vs_all, opt$use_anota2), file = log_file, append = TRUE)
  cat(sprintf("BiocParallel workers=%d; thread guards=%s\n", opt$cores, paste(sprintf("%s=%s", thread_guard_vars, Sys.getenv(thread_guard_vars, unset = NA_character_)), collapse=", ")), file = log_file, append = TRUE)
  
} else {
  # --- NORMAL EXECUTION START ---
  
  # Set default save path if not provided
  if (is.null(opt$save_model)) {
      opt$save_model <- file.path(opt$outdir, "dTEct_model.RData")
      cat(paste0("No --save_model provided. Defaulting to: ", opt$save_model, "\n"))
  }
  
  # Clear existing content of the log file at the start of the script
  file.create(log_file) 
  cat(sprintf("BiocParallel workers=%d; thread guards=%s\n", opt$cores, paste(sprintf("%s=%s", thread_guard_vars, Sys.getenv(thread_guard_vars, unset = NA_character_)), collapse=", ")), file = log_file, append = TRUE)
  
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
if (!isTRUE(opt$use_anota2)) {
  dir.create(paste0(opt$outdir, 'TE'), showWarnings = FALSE, recursive = TRUE)
}

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
    # (transcript IDs must be preserved so RNA contrast output uses TT_*/ENST... IDs)
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
dge <- calc_assay_norm_factors(dge, assay_col="seq_type", method="TMM", log_file=log_file)

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
if (isTRUE(opt$use_anota2)) {
    cat("Valid contrasts found. Building anota2seq model...\n", file=log_file, append=TRUE)
    fit_paired <- NULL
    fit_rna <- NULL
    design.rna <- NULL
    fit_anota2seq <- build_anota2seq_dataset(dge, log_file)
} else {
    cat("Valid contrasts found. Fitting edgeR quasi-likelihood GLM models...\n", file=log_file, append=TRUE)
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
        # This ensures RNA contrasts report transcript IDs (TT_*/ENST...) not ST_*/TM_...
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
        dge_rna <- calcNormFactors(dge_rna, method="TMM")
        dge_rna <- estimateDisp(dge_rna, design.rna)
        fit_rna <- glmQLFit(dge_rna, design.rna)
    }
}

# --- SAVE MODEL LOGIC ---
if (!is.null(opt$save_model)) {
    cat(paste0("Saving model to ", opt$save_model, "...\n"), file=log_file, append=TRUE)
    # List of objects to save to ensure the environment is reproducible
    # We save everything relevant for the contrast execution phase
    save_objs <- c("dge", "meta.samples", "design", "contrast_col", "uniq_dict", "comb_dict",
                   "qc_types", "opt", "log_file", "tx.table", "feature2name", "contrast_grps", "seq_types")
    if (isTRUE(opt$use_anota2)) {
        save_objs <- c(save_objs, "fit_anota2seq", "anota2seq_logcpm")
    } else {
        save_objs <- c(save_objs, "fit_paired", "fit_rna", "design.rna")
    }
    # Only include translon_to_transcript_map if it exists (translon mode only)
    if (exists("translon_to_transcript_map")) save_objs <- c(save_objs, "translon_to_transcript_map")
    save_objs <- save_objs[save_objs %in% ls()]
    save(list = save_objs, file = opt$save_model)
    cat("Model saved successfully.\n", file=log_file, append=TRUE)
}

} # End of if (!is.null(opt$load_model)) else block

# Resume Execution
dge_idxs <- match(rownames(dge$samples), meta.samples$counts_col)
dge.meta <- meta.samples[dge_idxs,]
contrast_workers <- if (isTRUE(opt$use_anota2)) 1L else opt$cores
if (isTRUE(opt$use_anota2)) {
    cat("Running final anota2seq contrast extraction/writing serially to avoid forked plotting/package lazy-load failures.\n", file=log_file, append=TRUE)
}

if (isTRUE(opt$use_anota2)) {
    anota2seq_jobs <- collect_anota2seq_jobs(
        dge$samples,
        dge.meta,
        contrast_col,
        comb_dict[[contrast_col]],
        uniq_dict[[contrast_col]],
        opt$outdir,
        opt$skip_pairwise,
        opt$skip_one_vs_all,
        log_file
    )
    prepare_anota2seq_contrast_cache(fit_anota2seq, anota2seq_jobs, log_file)
}

# -----------------------------------------------------------------------------
# 6. EXECUTE CONTRASTS
# -----------------------------------------------------------------------------
cat("Executing contrasts...\n", file=log_file, append=TRUE)

# Execute Unique Contrasts (TE mostly, uses Paired Fit)
if (!isTRUE(opt$use_anota2)) {
    run_bplapply(uniq_dict[[contrast_col]], function(val) {
        evaluate_unique_contrast(dge$samples, val, contrast_col, fit_paired, opt$outdir, log_file)
    }, workers = opt$cores)
} else {
    cat("Skipping TE contrasts in anota2seq mode.\n", file=log_file, append=TRUE)
}

# Execute Combination Contrasts (DE & dTE, uses Paired + RNA Independent Fit)
# Execute Combination Contrasts (DE & dTE, uses Paired + RNA Independent Fit)
if (!opt$skip_pairwise) {
    run_bplapply(comb_dict[[contrast_col]], function(val) {
        evaluate_combination_contrast(dge$samples, val, contrast_col, fit_paired, fit_rna, opt$outdir, log_file)
    }, workers = contrast_workers)
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
    run_bplapply(uniq_dict[[contrast_col]], function(val) {
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
    }, workers = contrast_workers)
}

cat("Analysis complete.\n", file=log_file, append=TRUE)
