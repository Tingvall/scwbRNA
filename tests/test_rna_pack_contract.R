.rna_activate <- function() {
  try(scworkbench::wb_operation_pack_unload("scwbrna"), silent = TRUE)
  scworkbench::wb_operation_pack_activate(system.file("scworkbench",
    "operation-pack.json", package = "scwbRNA"))
}

.rna_fixture <- function() {
  reference <- jsonlite::read_json(file.path(.scwbrna_test_root, "fixtures",
    "old_core_reference.json"), simplifyVector = FALSE)
  root <- tempfile("scwbrna_contract_")
  bundle <- scworkbench::wb_init_bundle(root)
  cells <- as.character(unlist(reference$cells, use.names = FALSE))
  features <- as.character(unlist(reference$features, use.names = FALSE))
  values <- do.call(rbind, lapply(reference$values, function(x)
    as.numeric(unlist(x, use.names = FALSE))))
  dimnames(values) <- list(features, cells)
  scworkbench::wb_define_cells(bundle, cells)
  scworkbench::wb_attach_matrix(bundle, "RNA",
    Matrix::Matrix(values, sparse = TRUE), feature_type = "gene")
  list(root = root, bundle = bundle, reference = reference,
    cells = cells, features = features)
}

.rna_identity <- function(bundle, value) {
  x <- scworkbench:::.wb_prov_payload_identity(bundle, value)
  list(scheme = x$content_identity_scheme, sha256 = x$content_sha256)
}

.rna_committed_state <- function(bundle) {
  dag <- scworkbench::wb_prov_bundle_dag(bundle)
  list(scientific_state = scworkbench:::.wb_session_source_identity(
    bundle)$scientific_state_fingerprint,
    artifacts = nrow(scworkbench::wb_art_list(bundle,
      include_history = TRUE)),
    heads = nrow(scworkbench::wb_prov_list_heads(bundle)),
    runs = sum(as.character(dag$nodes$node_type) == "run"))
}

test_scwbrna_old_core_equivalence_contract <- function() {
  .rna_activate(); x <- .rna_fixture()
  on.exit({scworkbench::wb_close_bundle(x$bundle); unlink(x$root,
    recursive = TRUE, force = TRUE)}, add = TRUE)
  session <- scworkbench::wb_session_new()
  suffix <- scworkbench::wb_session_run_op(session, x$bundle,
    "scwbrna::cell_suffix_label", inputs = list(cell_set = "all",
      min_cells = 1L, sep = "-",
      mapping = list(A = "sample_a", B = "sample_b"),
      default_label = "other"))
  qc <- scworkbench::wb_session_run_op(suffix$session, x$bundle,
    "scwbrna::rna_qc_metrics", inputs = list(matrix = "RNA",
      cell_set = "all", min_cells = 1L, mitochondrial_pattern = "^mt-",
      ribosomal_pattern = "^rpl", missing_feature_policy = "drop",
      chunk_size = 2L))
  hvg <- scworkbench::wb_session_run_op(qc$session, x$bundle,
    "scwbrna::hvg_select", inputs = list(matrix = "RNA", cell_set = "all",
      min_cells = 1L, n_features = 2L, method = "vst", min_mean = 0,
      exclude_features = "rpl1", feature_type = "gene", seed = 1L,
      chunk_size = 2L))
  runs <- list(cell_suffix_label = suffix$run, rna_qc_metrics = qc$run,
    hvg_select = hvg$run)
  old_effective_params <- list(
    cell_suffix_label = list(sep = "-",
      mapping = list(A = "sample_a", B = "sample_b"),
      default_label = "other"),
    rna_qc_metrics = list(mitochondrial_features = NULL,
      mitochondrial_pattern = "^mt-", ribosomal_features = NULL,
      ribosomal_pattern = "^rpl", pattern_ignore_case = FALSE,
      missing_feature_policy = "drop", require_integer_counts = TRUE,
      integer_tolerance = 1e-8, chunk_size = 2L),
    hvg_select = list(n_features = 2L, method = "vst", min_mean = 0,
      max_mean = NULL, min_dispersion = NULL, batch_key = NULL,
      exclude_features = "rpl1", feature_type = "gene", seed = 1L,
      chunk_size = 2L)
  )
  for (method in names(runs)) {
    expected <- x$reference$methods[[method]]$outputs
    observed <- lapply(runs[[method]]$outputs, .rna_identity,
      bundle = x$bundle)
    stopifnot(setequal(names(observed), names(expected)),
      identical(observed[names(expected)], expected),
      identical(runs[[method]]$params, old_effective_params[[method]]))
  }
  TRUE
}

test_scwbrna_conditional_key_and_integer_contract <- function() {
  .rna_activate(); x <- .rna_fixture()
  on.exit({scworkbench::wb_close_bundle(x$bundle); unlink(x$root,
    recursive = TRUE, force = TRUE)}, add = TRUE)
  run <- function(...) scworkbench::wb_session_run_op(
    scworkbench::wb_session_new(), x$bundle, "scwbrna::rna_qc_metrics",
    inputs = c(list(matrix = "RNA", cell_set = "all", min_cells = 1L,
      chunk_size = 2L), list(...)))$run$outputs$metrics
  stopifnot(identical(unique(run()$key), c("total_counts",
    "detected_features")), identical(unique(run(
      mitochondrial_pattern = "^mt-")$key), c("total_counts",
        "detected_features", "mitochondrial_counts",
        "mitochondrial_percent")), identical(unique(run(
      ribosomal_pattern = "^rpl")$key), c("total_counts",
        "detected_features", "ribosomal_counts", "ribosomal_percent")),
    length(unique(run(mitochondrial_pattern = "^mt-",
      ribosomal_pattern = "^rpl")$key)) == 6L)
  counts <- run(mitochondrial_pattern = "^mt-")
  stopifnot(all(counts$value[counts$key == "mitochondrial_counts"] %% 1 == 0),
    all(is.finite(counts$value)), all(counts$value >= 0))
  TRUE
}

test_scwbrna_hvg_batch_key_rejection_contract <- function() {
  .rna_activate(); x <- .rna_fixture()
  on.exit({try(scworkbench::wb_close_bundle(x$bundle), silent = TRUE);
    unlink(x$root, recursive = TRUE, force = TRUE)}, add = TRUE)
  run <- function(batch_key) scworkbench::wb_session_run_op(
    scworkbench::wb_session_new(), x$bundle, "scwbrna::hvg_select",
    inputs = list(matrix = "RNA", cell_set = "all", min_cells = 1L,
      n_features = 2L, method = "vst", batch_key = batch_key,
      feature_type = "gene", chunk_size = 2L))

  # Null and omission remain the v3 representations of no batch key (omission
  # is also covered by the full old-core comparator). The generic input
  # boundary rejects empty strings before dispatch, so that old spelling is
  # not representable by the accepted v3 transport.
  null <- run(NULL)
  stopifnot(is.null(null$run$params$batch_key))
  scworkbench::wb_session_discard(null$session)
  before <- .rna_committed_state(x$bundle)
  error <- tryCatch({run("sample"); NULL}, error = identity)
  stopifnot(inherits(error, "error"),
    grepl("batch_key is not supported", conditionMessage(error),
      fixed = TRUE),
    identical(before, .rna_committed_state(x$bundle)))
  TRUE
}

test_scwbrna_rna_qc_failure_contract <- function() {
  .rna_activate(); x <- .rna_fixture()
  on.exit({try(scworkbench::wb_close_bundle(x$bundle), silent = TRUE);
    unlink(x$root, recursive = TRUE, force = TRUE)}, add = TRUE)
  cells <- x$cells[1:2]
  attach <- function(id, values) {
    dimnames(values) <- list(c("g1", "g2"), cells)
    scworkbench::wb_attach_matrix(x$bundle, id,
      Matrix::Matrix(values, sparse = TRUE), feature_type = "gene")
  }
  attach("RNA_fractional", matrix(c(1, 0.5, 0, 2), nrow = 2L))
  attach("RNA_negative", matrix(c(1, -1, 0, 2), nrow = 2L))
  attach("RNA_nonfinite", matrix(c(Inf, 0, 0, 2), nrow = 2L))
  run <- function(matrix_id, ...) scworkbench::wb_session_run_op(
    scworkbench::wb_session_new(), x$bundle, "scwbrna::rna_qc_metrics",
    inputs = c(list(matrix = matrix_id, cell_set = "all", min_cells = 1L,
      chunk_size = 2L), list(...)))
  expect_failure <- function(matrix_id, pattern, ...) {
    before <- .rna_committed_state(x$bundle)
    error <- tryCatch({run(matrix_id, ...); NULL}, error = identity)
    stopifnot(inherits(error, "error"),
      grepl(pattern, conditionMessage(error), fixed = TRUE),
      identical(before, .rna_committed_state(x$bundle)))
    invisible(error)
  }

  expect_failure("RNA_fractional", "matrix contains non-integer values")
  expect_failure("RNA_negative", "matrix contains negative or non-finite values")
  expect_failure("RNA_nonfinite", "finite",
    require_integer_counts = FALSE)
  expect_failure("RNA_fractional",
    "mitochondrial feature ids are absent from the matrix",
    require_integer_counts = FALSE,
    mitochondrial_features = "not_present")

  before_drop <- .rna_committed_state(x$bundle)
  dropped <- run("RNA_fractional", require_integer_counts = FALSE,
    mitochondrial_features = "not_present",
    missing_feature_policy = "drop")
  stopifnot(identical(unique(dropped$run$outputs$metrics$key),
      c("total_counts", "detected_features")),
    nrow(dropped$run$outputs$qc_features$table) == 0L,
    identical(dropped$run$outputs$qc_features$meta$missing_feature_policy,
      "drop"),
    identical(dropped$run$outputs$summary$meta$mitochondrial_missing_features,
      "not_present"))
  scworkbench::wb_session_discard(dropped$session)
  stopifnot(identical(before_drop, .rna_committed_state(x$bundle)))
  TRUE
}

test_scwbrna_promotion_reopen_contract <- function() {
  .rna_activate(); x <- .rna_fixture()
  on.exit({try(scworkbench::wb_close_bundle(x$bundle), silent = TRUE);
    unlink(x$root, recursive = TRUE, force = TRUE)}, add = TRUE)
  suffix <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    x$bundle, "scwbrna::cell_suffix_label", inputs = list(cell_set = "all",
      min_cells = 1L, sep = "-", default_label = "unknown"))
  run <- scworkbench::wb_session_run_op(suffix$session,
    x$bundle, "scwbrna::hvg_select", inputs = list(matrix = "RNA",
      cell_set = "all", min_cells = 1L, n_features = 2L,
      feature_type = "gene", chunk_size = 2L))
  promoted <- scworkbench::wb_session_promote_current(x$bundle, run$session,
    outputs = c("labels", "label_table", "selected_features", "fit_summary"),
    name_prefix = "rna",
    created_by = "scwbrna-test")
  ids <- as.character(promoted$promoted$artifact_id)
  transition_ids <- vapply(list(suffix$run, run$run), function(value)
    as.character(value$operation_pack_transition$transition_id)[1], character(1))
  for (transition_id in transition_ids) {
    reads <- scworkbench::wb_operation_pack_transition_registry_reads(
      x$bundle, transition_id)
    stopifnot(length(reads) == 1L,
      reads[[1L]]$kind %in% c("session_cell_axis", "matrix_feature_annotations"))
  }
  scworkbench::wb_finalize_recipient_bundle(x$bundle)
  scworkbench::wb_close_bundle(x$bundle)
  x$bundle <- scworkbench::wb_open_bundle(x$root, read_only = TRUE)
  stopifnot(length(ids) == 4L,
    all(vapply(ids, function(id) inherits(scworkbench::wb_art_load(
      x$bundle, id), "wb_primitive"), logical(1))),
    isTRUE(scworkbench::wb_validate_bundle(x$root)))
  TRUE
}

test_scwbrna_activation_rollback_contract <- function() {
  .rna_activate()
  before_ops <- scworkbench::wb_list_ops()
  before_packs <- scworkbench::wb_operation_packs()
  before_identity <- scworkbench::wb_operation_catalog_identity()
  source <- system.file("scworkbench", "operation-pack.json",
    package = "scwbRNA")
  root <- tempfile("scwbrna_bad_manifest_"); dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  manifest <- jsonlite::read_json(source, simplifyVector = FALSE)
  manifest$operations[[1L]]$dimension_policies[[1L]]$row_order <-
    "key_major_then_cell_v1"
  target <- file.path(root, "operation-pack.json")
  writeChar(jsonlite::toJSON(manifest, auto_unbox = TRUE, null = "null",
    digits = NA), target, eos = NULL, useBytes = TRUE)
  error <- try(scworkbench::wb_operation_pack_activate(target), silent = TRUE)
  stopifnot(inherits(error, "try-error"),
    identical(before_ops, scworkbench::wb_list_ops()),
    identical(before_packs, scworkbench::wb_operation_packs()),
    identical(before_identity, scworkbench::wb_operation_catalog_identity()))
  TRUE
}
