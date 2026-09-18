.run_rna_qc_metrics <- function(readers, writers, params, context) {
  capability <- context$typed_slice_reads$matrix_chunks
  features <- as.character(capability$features)
  cells <- as.character(capability$cells)
  missing_policy <- as.character(params$missing_feature_policy %||% "error")[1]
  mito <- .rna_groups(features, params$mitochondrial_features,
    params$mitochondrial_pattern, isTRUE(params$pattern_ignore_case),
    missing_policy, "mitochondrial")
  ribo <- .rna_groups(features, params$ribosomal_features,
    params$ribosomal_pattern, isTRUE(params$pattern_ignore_case),
    missing_policy, "ribosomal")
  chunks <- scworkbench::wb_partition_axis(cells,
    as.integer(params$chunk_size %||% 2000L))
  matrix_id <- NULL
  values <- lapply(chunks, function(chunk) {
    primitive <- readers$matrix_chunks(features = features, cells = chunk)
    if (is.null(matrix_id)) matrix_id <<- as.character(primitive$matrix_id)[1]
    x <- primitive$values
    stored <- if (inherits(x, "sparseMatrix")) x@x else as.numeric(x)
    .rna_assert(all(is.finite(stored)) && all(stored >= 0),
      "rna_qc_metrics: matrix contains negative or non-finite values")
    tolerance <- as.numeric(params$integer_tolerance %||% 1e-8)
    if (isTRUE(params$require_integer_counts %||% TRUE))
      .rna_assert(all(abs(stored - round(stored)) <= tolerance),
        "rna_qc_metrics: matrix contains non-integer values; use a raw count matrix or explicitly set require_integer_counts=FALSE")
    total <- as.numeric(Matrix::colSums(x))
    out <- list(total_counts = total,
      detected_features = as.numeric(Matrix::colSums(x > 0)))
    if (length(mito$features)) {
      count <- as.numeric(Matrix::colSums(x[mito$features, , drop = FALSE]))
      out$mitochondrial_counts <- count
      out$mitochondrial_percent <- ifelse(total > 0, 100 * count / total, 0)
    }
    if (length(ribo$features)) {
      count <- as.numeric(Matrix::colSums(x[ribo$features, , drop = FALSE]))
      out$ribosomal_counts <- count
      out$ribosomal_percent <- ifelse(total > 0, 100 * count / total, 0)
    }
    out
  })
  keys <- names(values[[1]])
  metric_values <- stats::setNames(lapply(keys, function(key)
    unlist(lapply(values, `[[`, key), use.names = FALSE)), keys)
  metrics <- scworkbench::wb_prim_cell_vector_family(
    rep(keys, each = length(cells)), rep(cells, length(keys)),
    unlist(metric_values, use.names = FALSE), value_type = "continuous",
    meta = list(matrix_id = matrix_id, metric_keys = keys,
      percentage_unit = "percent", zero_count_percentage_convention = "zero",
      threshold_selection = "not_performed"))
  rows <- list()
  if (length(mito$features)) rows[[length(rows) + 1L]] <- data.frame(
    feature_id = mito$features, metric_group = "mitochondrial",
    selection_source = unname(mito$source), matrix_id = matrix_id,
    stringsAsFactors = FALSE)
  if (length(ribo$features)) rows[[length(rows) + 1L]] <- data.frame(
    feature_id = ribo$features, metric_group = "ribosomal",
    selection_source = unname(ribo$source), matrix_id = matrix_id,
    stringsAsFactors = FALSE)
  feature_table <- if (length(rows)) do.call(rbind, rows) else data.frame(
    feature_id = character(), metric_group = character(),
    selection_source = character(), matrix_id = character(),
    stringsAsFactors = FALSE)
  meta <- list(matrix_id = matrix_id,
    mitochondrial_pattern = mito$pattern,
    mitochondrial_missing_features = mito$missing,
    ribosomal_pattern = ribo$pattern,
    ribosomal_missing_features = ribo$missing,
    missing_feature_policy = missing_policy)
  zero <- sum(metric_values$total_counts == 0)
  summary <- do.call(rbind, lapply(keys, function(key) {
    value <- as.numeric(metric_values[[key]]); finite <- value[is.finite(value)]
    data.frame(matrix_id = matrix_id, metric = key,
      n_cells = as.integer(length(cells)),
      n_missing = as.integer(sum(!is.finite(value))),
      n_zero_denominator = if (grepl("_percent$", key)) as.integer(zero) else 0L,
      minimum = if (length(finite)) min(finite) else NA_real_,
      median = if (length(finite)) stats::median(finite) else NA_real_,
      mean = if (length(finite)) mean(finite) else NA_real_,
      maximum = if (length(finite)) max(finite) else NA_real_,
      stringsAsFactors = FALSE)
  }))
  writers$write_terminal_primitive("metrics", metrics)
  writers$write_terminal_primitive("qc_features",
    scworkbench::wb_prim_feature_table(feature_table,
      table_type = "rna_qc_feature_membership", meta = meta))
  writers$write_terminal_primitive("summary", scworkbench::wb_prim_obs_table(
    summary, table_type = "rna_qc_metrics_summary", entity_kind = "qc_metric",
    meta = c(meta, list(n_cells = length(cells), n_zero_count_cells = zero,
      zero_count_percentage_convention = "zero",
      threshold_selection = "not_performed"))))
  invisible(list(status = "ok"))
}
