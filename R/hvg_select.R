.run_hvg_select <- function(readers, writers, params, context) {
  annotations <- readers$feature_annotations()
  capability <- context$typed_slice_reads$matrix_chunks
  all_features <- as.character(capability$features)
  features <- all_features
  cells <- as.character(capability$cells)
  rows <- annotations$rows
  annotation_ids <- vapply(rows, `[[`, character(1), "feature_id")
  .rna_assert(identical(features, annotation_ids),
    "hvg_select: feature annotations are not matrix aligned")
  feature_type <- params$feature_type %||% NULL
  if (!is.null(feature_type)) {
    types <- vapply(rows, function(x) as.character(
      x$values$feature_type %||% NA_character_)[1], character(1))
    features <- features[!is.na(types) & types == as.character(feature_type)[1]]
  }
  if (!is.null(params$exclude_features))
    features <- setdiff(features, as.character(params$exclude_features))
  .rna_assert(length(features) > 0L,
    "hvg_select: no features available after filters")
  stats <- .rna_feature_stats(readers$matrix_chunks, all_features, cells,
    params$chunk_size %||% 256L)
  stats <- stats[stats$feature_id %in% features, , drop = FALSE]
  method <- as.character(params$method %||% "vst")[1]
  stats$dispersion <- stats$variance / pmax(stats$mean, .Machine$double.eps)
  stats$score <- if (identical(method, "variance")) stats$variance else if (
    identical(method, "dispersion")) stats$dispersion else
    log1p(pmax(stats$dispersion, 0))
  stats$method <- method
  matrix_id <- as.character(annotations$matrix_id)[1]
  stats$matrix_id <- matrix_id
  passes <- is.finite(stats$score)
  if (!is.null(params$min_mean)) passes <- passes & stats$mean >= params$min_mean
  if (!is.null(params$max_mean)) passes <- passes & stats$mean <= params$max_mean
  if (!is.null(params$min_dispersion)) passes <- passes &
    stats$score >= params$min_dispersion
  ord <- order(!passes, -stats$score, -stats$variance, stats$feature_id,
    na.last = TRUE)
  rank <- integer(nrow(stats)); rank[ord] <- seq_along(ord)
  selected <- ord[passes[ord]]
  selected <- selected[seq_len(min(length(selected),
    as.integer(params$n_features %||% 2000L)))]
  stats$selected <- FALSE; stats$selected[selected] <- TRUE; stats$rank <- rank
  stats <- stats[order(stats$rank), c("feature_id", "rank", "mean",
    "variance", "dispersion", "score", "detected_cells",
    "detected_fraction", "selected", "method", "matrix_id"), drop = FALSE]
  rownames(stats) <- NULL
  summary <- data.frame(matrix_id = matrix_id, method = method,
    n_cells = as.integer(length(cells)),
    n_input_features = as.integer(length(features)),
    n_passing_filters = as.integer(sum(passes, na.rm = TRUE)),
    n_selected_features = as.integer(sum(stats$selected, na.rm = TRUE)),
    min_mean = as.numeric(params$min_mean %||% NA_real_),
    max_mean = as.numeric(params$max_mean %||% NA_real_),
    min_dispersion = as.numeric(params$min_dispersion %||% NA_real_),
    batch_key = as.character(params$batch_key %||% NA_character_)[1],
    feature_type = as.character(params$feature_type %||% NA_character_)[1],
    chunk_size = as.integer(params$chunk_size %||% 256L),
    stringsAsFactors = FALSE)
  writers$write_terminal_primitive("selected_features",
    scworkbench::wb_prim_feature_table(stats,
      table_type = "feature_selection_table", meta = list(
        selection_kind = "highly_variable_features", matrix_id = matrix_id,
        method = method, n_selected_features =
          as.integer(sum(stats$selected, na.rm = TRUE)))))
  writers$write_terminal_primitive("fit_summary",
    scworkbench::wb_prim_obs_table(summary,
      table_type = "hvg_select_summary", entity_kind = "matrix",
      meta = list(matrix_id = matrix_id, method = method)))
  invisible(list(status = "ok"))
}
