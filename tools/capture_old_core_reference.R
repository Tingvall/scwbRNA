args <- commandArgs(trailingOnly = TRUE)
value <- function(prefix) {
  hit <- args[startsWith(args, paste0("--", prefix, "="))]
  if (length(hit) != 1L) stop("exactly one --", prefix,
    "= argument is required")
  sub(paste0("--", prefix, "="), "", hit, fixed = TRUE)
}
core <- normalizePath(value("core-source"), mustWork = TRUE)
output <- value("output")
expected_commit <- "03accf6ac81af0a383cc549d2ae89351af55c7d0"
commit <- trimws(system2("git", c("-C", shQuote(core), "rev-parse", "HEAD"),
  stdout = TRUE))
dirty <- system2("git", c("-C", shQuote(core), "status", "--porcelain",
  "--untracked-files=all"), stdout = TRUE)
if (!identical(commit, expected_commit) || any(nzchar(dirty)))
  stop("old-core source must be clean and pinned to ", expected_commit)
devtools::load_all(core, quiet = TRUE)
root <- tempfile("scwbrna_old_core_")
bundle <- wb_init_bundle(root)
on.exit({ wb_close_bundle(bundle); unlink(root, recursive = TRUE) }, add = TRUE)
cells <- c("c2-A", "c1-B", "c3-A", "c4-Z")
features <- c("g2", "mt-1", "rpl1", "g1", "peak1")
values <- Matrix::Matrix(matrix(c(
  1, 0, 3, 2,
  0, 0, 1, 0,
  4, 2, 0, 1,
  2, 1, 5, 0,
  7, 0, 0, 3), nrow = length(features), byrow = TRUE,
  dimnames = list(features, cells)), sparse = TRUE)
wb_define_cells(bundle, cells)
wb_attach_matrix(bundle, "RNA", values, feature_type = "gene")
normalize_run <- wb_session_run_op(wb_session_new(), bundle,
  "matrix_log_normalize", inputs = list(matrix = "RNA", cell_set = "all",
    min_cells = 1L, output_matrix_id = "RNA_log", scale_factor = 1000,
    overwrite = FALSE))
normalize <- list(outputs = normalize_run$run$outputs,
  params = normalize_run$run$params)
identity <- function(x, identity_bundle = bundle) {
  out <- scworkbench:::.wb_prov_payload_identity(identity_bundle, x)
  list(scheme = out$content_identity_scheme, sha256 = out$content_sha256)
}
record <- function(result, identity_bundle = bundle) list(outputs = lapply(
    result$outputs, identity, identity_bundle = identity_bundle),
  params = result$params)
target <- output
dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
fixture <- jsonlite::read_json(target, simplifyVector = FALSE)
fixture$phase4c_source_commit <- expected_commit
fixture$methods$matrix_log_normalize <- record(normalize,
  wb_session_view(normalize_run$session))
writeChar(scworkbench:::.wb_operation_pack_canonical_json(fixture), target,
  eos = NULL, useBytes = TRUE)
