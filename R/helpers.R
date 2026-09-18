`%||%` <- function(x, y) if (is.null(x)) y else x

.rna_assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
  invisible(TRUE)
}

.rna_groups <- function(features, explicit, pattern, ignore_case,
    missing_policy, label) {
  features <- as.character(features)
  explicit <- unique(as.character(explicit %||% character()))
  explicit <- explicit[!is.na(explicit) & nzchar(explicit)]
  pattern <- as.character(pattern %||% "")[1]
  missing <- setdiff(explicit, features)
  if (length(missing) && identical(missing_policy, "error"))
    stop(label, " feature ids are absent from the matrix: ",
      paste(utils::head(missing, 20L), collapse = ", "), call. = FALSE)
  direct <- intersect(explicit, features)
  by_pattern <- character()
  if (nzchar(pattern)) {
    hit <- tryCatch(grepl(pattern, features, ignore.case = ignore_case,
      perl = TRUE), error = function(e) stop(label, " pattern is invalid: ",
        conditionMessage(e), call. = FALSE))
    by_pattern <- features[hit]
    .rna_assert(length(by_pattern) > 0L,
      paste0(label, " pattern matched no matrix features: ", pattern))
  }
  selected <- features[features %in% union(direct, by_pattern)]
  source <- vapply(selected, function(id) {
    if (id %in% direct && id %in% by_pattern) "explicit+pattern" else
      if (id %in% direct) "explicit" else "pattern"
  }, character(1))
  list(features = selected, source = source, missing = missing,
    pattern = pattern)
}

.rna_feature_stats <- function(reader, features, cells, chunk_size) {
  chunks <- scworkbench::wb_partition_axis(features, as.integer(chunk_size))
  pieces <- lapply(chunks, function(ids) {
    primitive <- reader(features = ids, cells = cells)
    x <- primitive$values
    n <- ncol(x); mu <- as.numeric(Matrix::rowMeans(x))
    ss <- as.numeric(Matrix::rowSums(x * x))
    variance <- if (n > 1L) pmax((ss - n * mu * mu) / (n - 1L), 0) else
      rep(NA_real_, length(mu))
    detected <- as.numeric(Matrix::rowSums(x > 0))
    data.frame(feature_id = ids, mean = mu, variance = variance,
      detected_cells = detected,
      detected_fraction = detected / max(n, 1L), stringsAsFactors = FALSE)
  })
  do.call(rbind, pieces)
}
