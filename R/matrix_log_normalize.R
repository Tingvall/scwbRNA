.run_matrix_log_normalize <- function(readers, writers, params, context) {
  capability <- context$typed_slice_reads$matrix_chunks
  features <- as.character(capability$features)
  cells <- as.character(capability$cells)
  scale_factor <- as.numeric(params$scale_factor %||% 10000)
  .rna_assert(length(features) > 0L && length(cells) > 0L,
    "matrix_log_normalize: authenticated axes are empty")
  .rna_assert(is.finite(scale_factor) && scale_factor > 0,
    "matrix_log_normalize: scale_factor must be positive")
  handle <- writers$open_matrix_fragments("matrix")
  chunk_size <- as.integer(params$chunk_size %||% 2000L)
  library_sizes <- numeric(length(cells)); names(library_sizes) <- cells
  output_parts <- 0L
  for (source in handle$source_fragments) {
    start <- source$cell_start
    source_cells <- cells[seq.int(start,
      start + source$cell_count - 1L)]
    for (part_cells in scworkbench::wb_partition_axis(source_cells,
        chunk_size)) {
      primitive <- readers$matrix_chunks(features = features,
        cells = part_cells)
      values <- primitive$values
      library_size <- as.numeric(Matrix::colSums(values))
      names(library_size) <- part_cells
      denom <- library_size
      denom[!is.finite(denom) | denom <= 0] <- Inf
      normalized <- log1p(t(t(as.matrix(values)) *
        (scale_factor / denom)))
      dimnames(normalized) <- list(features, part_cells)
      handle$write_fragment(source$ordinal, normalized, features,
        part_cells)
      library_sizes[part_cells] <- library_size
      output_parts <- output_parts + 1L
    }
  }
  description <- paste0("Log-normalized matrix: log1p(",
    context$input_ids$matrix, " / per-cell library size * ",
    format(scale_factor, scientific = FALSE, trim = TRUE),
    "). Source matrix: ", context$input_ids$matrix, ".")
  handle$finish(list(source_matrix_id = context$input_ids$matrix,
    output_matrix_id = params$output_matrix_id,
    transform = "log1p_cp", scale_factor = scale_factor,
    source = "matrix_log_normalize", description = description))
  zero <- cells[!is.finite(library_sizes) | library_sizes <= 0]
  summary <- data.frame(
    output_matrix_id = as.character(params$output_matrix_id)[1],
    source_matrix_id = as.character(context$input_ids$matrix)[1],
    transform = "log1p_cp",
    formula = paste0("log1p(counts / library_size * ",
      format(scale_factor, scientific = FALSE, trim = TRUE), ")"),
    scale_factor = scale_factor, n_features = as.integer(length(features)),
    n_cells = as.integer(length(cells)),
    n_fragments = as.integer(output_parts),
    n_zero_library_cells = as.integer(length(zero)),
    min_library_size = suppressWarnings(min(library_sizes, na.rm = TRUE)),
    median_library_size = suppressWarnings(stats::median(library_sizes,
      na.rm = TRUE)),
    max_library_size = suppressWarnings(max(library_sizes, na.rm = TRUE)),
    storage_representation = "bpcells_uncompressed_double",
    stringsAsFactors = FALSE)
  writers$write_terminal_primitive("summary",
    scworkbench::wb_prim_obs_table(summary,
      table_type = "matrix_log_normalize_summary", entity_kind = "matrix",
      meta = list(source_matrix_id = context$input_ids$matrix,
        output_matrix_id = params$output_matrix_id,
        transform = "log1p_cp", scale_factor = scale_factor)))
  invisible(list(status = "ok"))
}
