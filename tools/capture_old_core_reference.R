library(scworkbench)
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
suffix <- wb_op_cell_suffix_label(bundle, cells = cells, sep = "-",
  mapping = c(A = "sample_a", B = "sample_b"), default_label = "other")
qc <- wb_op_rna_qc_metrics(bundle, "RNA", cells,
  mitochondrial_pattern = "^mt-", ribosomal_pattern = "^rpl",
  missing_feature_policy = "drop", chunk_size = 2L)
hvg <- wb_op_hvg_select(bundle, "RNA", cells, n_features = 2L,
  method = "vst", min_mean = 0, exclude_features = "rpl1",
  feature_type = "gene", seed = 1L, chunk_size = 2L)
identity <- function(x) {
  out <- scworkbench:::.wb_prov_payload_identity(bundle, x)
  list(scheme = out$content_identity_scheme, sha256 = out$content_sha256)
}
record <- function(result) list(outputs = lapply(result$outputs, identity),
  params = result$params)
fixture <- list(schema = "scwbrna.old_core_reference.v1",
  source_commit = "ca3b4978a025b355aca6aad091c1f419de2cde12",
  cells = unname(as.list(cells)), features = unname(as.list(features)),
  values = unname(lapply(seq_len(nrow(values)), function(i)
    unname(as.list(as.numeric(values[i, ]))))),
  feature_type = list("gene", "gene", "gene", "gene", "gene"),
  methods = list(cell_suffix_label = record(suffix),
    rna_qc_metrics = record(qc), hvg_select = record(hvg)))
target <- file.path("tests", "fixtures", "old_core_reference.json")
writeChar(scworkbench:::.wb_operation_pack_canonical_json(fixture), target,
  eos = NULL, useBytes = TRUE)
