.run_cell_suffix_label <- function(readers, writers, params, context) {
  axis <- readers$session_cells()
  cells <- as.character(unlist(axis$cell_id, use.names = FALSE))
  sep <- as.character(params$sep %||% "-")[1]
  .rna_assert(nzchar(sep), "sep must be a non-empty string")
  suffix <- sub(paste0("^.*", gsub("([.|()\\^{}+$*?]|\\[|\\]|\\\\)",
    "\\\\\\1", sep), ""), "", cells)
  label <- suffix
  mapping <- params$mapping %||% NULL
  if (!is.null(mapping)) {
    mapping <- unlist(mapping, use.names = TRUE)
    .rna_assert(is.character(mapping) && !is.null(names(mapping)),
      "mapping must be a named character vector")
    label <- unname(mapping[suffix])
    label[is.na(label) | !nzchar(label)] <-
      as.character(params$default_label %||% "UNKNOWN")[1]
  }
  vector <- scworkbench::wb_prim_cell_vector(stats::setNames(label, cells),
    value_type = "categorical")
  table <- scworkbench::wb_prim_obs_table(data.frame(cell_id = cells,
    suffix = suffix, label = label, stringsAsFactors = FALSE),
    table_type = "cell_suffix_label", entity_kind = "cell",
    meta = list(sep = sep,
      default_label = as.character(params$default_label %||% "UNKNOWN")[1]))
  writers$write_terminal_primitive("labels", vector)
  writers$write_terminal_primitive("label_table", table)
  invisible(list(status = "ok"))
}
