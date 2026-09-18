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
  normalize <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    x$bundle, "scwbrna::matrix_log_normalize", inputs = list(
      matrix = "RNA", cell_set = "all", min_cells = 1L,
      output_matrix_id = "RNA_log", scale_factor = 1000,
      overwrite = FALSE, chunk_size = 100L))
  expected <- x$reference$methods$matrix_log_normalize$outputs
  observed <- lapply(normalize$run$outputs, .rna_identity,
    bundle = scworkbench::wb_session_view(normalize$session))
  stopifnot(identical(observed[names(expected)], expected),
    identical(normalize$run$params$output_matrix_id, "RNA_log"),
    identical(normalize$run$params$scale_factor, 1000),
    identical(normalize$run$params$overwrite, FALSE),
    identical(normalize$run$params$chunk_size, 100L))
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

test_scwbrna_matrix_log_normalize_publication_contract <- function() {
  .rna_activate(); x <- .rna_fixture()
  on.exit({try(scworkbench::wb_close_bundle(x$bundle), silent = TRUE);
    unlink(x$root, recursive = TRUE, force = TRUE)}, add = TRUE)
  values <- as.matrix(scworkbench::wb_slice(x$bundle, "RNA",
    x$features, x$cells))
  first <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    x$bundle, "matrix_log_normalize", inputs = list(matrix = "RNA",
      cell_set = "all", min_cells = 1L, output_matrix_id = "RNA_log",
      scale_factor = 1000, chunk_size = 2L, overwrite = FALSE))
  staged <- scworkbench::wb_session_view(first$session)
  expected <- log1p(t(t(values) * (1000 / colSums(values))))
  observed <- as.matrix(scworkbench::wb_slice(staged, "RNA_log",
    x$features, x$cells))
  stopifnot(isTRUE(all.equal(observed, expected, tolerance = 1e-10,
      check.attributes = FALSE)),
    first$run$outputs$matrix$fragment_manifest$fragment_count >= 2L,
    first$run$outputs$summary$table$n_zero_library_cells[[1]] == 0L)
  promoted <- scworkbench::wb_session_promote_current(x$bundle,
    first$session, outputs = "summary")
  stopifnot(all(promoted$promoted$status == "ok"),
    "RNA_log" %in% scworkbench::wb_matrices(x$bundle)$matrix_id,
    isTRUE(scworkbench::wb_bundle_validation_report(x$root, "content")$ok))

  before <- as.matrix(scworkbench::wb_slice(x$bundle, "RNA_log",
    x$features, x$cells))
  options(scworkbench.operation_pack_v3_fragment_test_hook =
    function(event, fragment_ordinal, ...) {
      if (identical(event, "before_fragment_write") &&
          fragment_ordinal == 2L) stop("injected second fragment failure")
    })
  on.exit(options(scworkbench.operation_pack_v3_fragment_test_hook = NULL),
    add = TRUE)
  failed <- try(scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    x$bundle, "matrix_log_normalize", inputs = list(matrix = "RNA",
      cell_set = "all", min_cells = 1L, output_matrix_id = "RNA_log",
      scale_factor = 2000, chunk_size = 2L, overwrite = TRUE)), silent = TRUE)
  options(scworkbench.operation_pack_v3_fragment_test_hook = NULL)
  stopifnot(inherits(failed, "try-error"), identical(before,
    as.matrix(scworkbench::wb_slice(x$bundle, "RNA_log", x$features,
      x$cells))))

  replacement <- scworkbench::wb_session_run_op(
    scworkbench::wb_session_new(), x$bundle, "matrix_log_normalize",
    inputs = list(matrix = "RNA", cell_set = "all", min_cells = 1L,
      output_matrix_id = "RNA_log", scale_factor = 2000,
      chunk_size = 2L, overwrite = TRUE))
  replaced <- scworkbench::wb_session_promote_current(x$bundle,
    replacement$session, outputs = "summary")
  after <- as.matrix(scworkbench::wb_slice(x$bundle, "RNA_log",
    x$features, x$cells))
  stopifnot(all(replaced$promoted$status == "ok"), !identical(before, after),
    isTRUE(scworkbench::wb_bundle_validation_report(x$root, "content")$ok))
  TRUE
}

test_scwbrna_zero_library_fragment_partition_contract <- function() {
  .rna_activate()
  root <- tempfile("scwbrna_zero_library_fragment_")
  bundle <- scworkbench::wb_init_bundle(root, overwrite = TRUE)
  on.exit({try(scworkbench::wb_close_bundle(bundle), silent = TRUE)
    unlink(root, recursive = TRUE, force = TRUE)}, add = TRUE)
  features <- c("g1", "g2")
  cells <- c("positive_library", "zero_library")
  counts <- matrix(c(1, 2, 0, 0), nrow = 2L,
    dimnames = list(features, cells))
  scworkbench::wb_define_cells(bundle, cells)
  scworkbench::wb_attach_matrix(bundle, "RNA",
    Matrix::Matrix(counts, sparse = TRUE), feature_type = "gene")

  execute <- function(chunk_size) {
    output_id <- paste0("RNA_log_chunk", chunk_size)
    result <- scworkbench::wb_session_run_op(
      scworkbench::wb_session_new(), bundle,
      "scwbrna::matrix_log_normalize", inputs = list(matrix = "RNA",
        cell_set = "all", min_cells = 1L, output_matrix_id = output_id,
        scale_factor = 1000, overwrite = FALSE,
        chunk_size = as.integer(chunk_size)))
    summary <- result$run$outputs$summary$table
    manifest <- result$run$outputs$matrix$fragment_manifest
    stopifnot(summary$n_zero_library_cells[[1L]] == 1L,
      summary$n_fragments[[1L]] == length(cells) / chunk_size,
      manifest$fragment_count == length(cells) / chunk_size)
    promoted <- scworkbench::wb_session_promote_current(bundle,
      result$session, outputs = "summary",
      created_by = "zero-library-fragment-regression")
    stopifnot(all(promoted$promoted$status == "ok"))
    list(output_id = output_id, run_id = result$run$run_id,
      transition_id = result$run$operation_pack_transition$transition_id,
      summary = summary, manifest = manifest)
  }

  chunk_one <- execute(1L)
  chunk_two <- execute(2L)
  scworkbench::wb_close_bundle(bundle)
  bundle <- scworkbench::wb_open_bundle(root, read_only = TRUE)

  retained <- function(value) {
    matrix <- as.matrix(scworkbench::wb_slice(bundle, value$output_id,
      features, cells))
    fragments <- scworkbench::wb_matrix_fragments(bundle, value$output_id)
    binding_columns <- paste(c("direction", "binding_name", "entity_kind",
      "entity_id", "relationship_type", "content_hash",
      "content_identity_scheme", "source"), collapse = ",")
    source_bindings <- unique(DBI::dbGetQuery(bundle$con, paste0(
      "SELECT ", binding_columns,
      " FROM prov_run_entity_bindings WHERE run_id=? AND entity_kind='matrix' AND entity_id='RNA' ORDER BY direction,binding_name,relationship_type;"),
      params = list(value$run_id)))
    output_bindings <- DBI::dbGetQuery(bundle$con,
      "SELECT direction,binding_name,relationship_type,content_hash,content_identity_scheme FROM prov_run_entity_bindings WHERE run_id=? AND entity_kind='matrix' AND entity_id=? ORDER BY snapshot_id,direction,binding_name;",
      params = list(value$run_id, value$output_id))
    outputs <- DBI::dbGetQuery(bundle$con,
      "SELECT output_name,primitive_kind FROM prov_run_outputs WHERE run_id=? ORDER BY snapshot_id,output_name;",
      params = list(value$run_id))
    record <- scworkbench:::.wb_operation_pack_v2_record_validate(
      file.path(root, "operation_pack_v2", value$transition_id),
      load_payloads = TRUE)
    output_index <- match("matrix", vapply(record$request$outputs, `[[`,
      character(1), "name"))
    grant <- record$request$outputs[[output_index]]$descriptor$
      matrix_fragment_grant
    stopifnot(isTRUE(record$ok), nrow(source_bindings) == 1L,
      nrow(output_bindings) >= 1L,
      any(output_bindings$relationship_type %in% c("PRODUCES",
        "MATERIALIZES")), any(outputs$output_name == "matrix"),
      nrow(fragments) == value$manifest$fragment_count,
      identical(unlist(lapply(seq_len(nrow(fragments)), function(i)
        as.character(jsonlite::fromJSON(fragments$feature_ids_json[[i]]))),
        use.names = FALSE), rep(features, nrow(fragments))),
      identical(unlist(lapply(seq_len(nrow(fragments)), function(i)
        as.character(jsonlite::fromJSON(fragments$cell_ids_json[[i]]))),
        use.names = FALSE), cells))
    list(matrix = matrix,
      content = scworkbench:::.wb_operation_pack_v3_matrix_slice_content(
        matrix, features, cells),
      source_bindings = source_bindings,
      source_schedule_identity = grant$source_fragment_schedule$
        schedule_identity,
      source_fragment_identities = unique(vapply(
        grant$source_fragment_schedule$fragments,
        function(fragment) fragment$source_fragment_identity$sha256,
        character(1))))
  }

  one <- retained(chunk_one)
  two <- retained(chunk_two)
  expected <- log1p(t(t(counts) * (1000 / c(3, Inf))))
  dimnames(expected) <- list(features, cells)
  scientific_summary_fields <- c("source_matrix_id", "transform", "formula",
    "scale_factor", "n_features", "n_cells", "n_zero_library_cells",
    "min_library_size", "median_library_size", "max_library_size",
    "storage_representation")
  stopifnot(isTRUE(all.equal(one$matrix, expected,
      tolerance = 1e-12, check.attributes = TRUE)),
    identical(one$matrix, two$matrix), identical(one$content, two$content),
    identical(chunk_one$summary[scientific_summary_fields],
      chunk_two$summary[scientific_summary_fields]),
    chunk_one$summary$n_zero_library_cells[[1L]] == 1L,
    chunk_two$summary$n_zero_library_cells[[1L]] == 1L,
    identical(rownames(one$matrix), features),
    identical(colnames(one$matrix), cells),
    identical(one$source_bindings, two$source_bindings),
    identical(one$source_schedule_identity,
      two$source_schedule_identity),
    identical(one$source_fragment_identities,
      two$source_fragment_identities),
    isTRUE(scworkbench::wb_bundle_validation_report(root, "content")$ok))
  TRUE
}
