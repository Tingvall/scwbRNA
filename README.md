# scwbRNA (provisional)

This provisional operation pack contains `cell_suffix_label`,
`rna_qc_metrics`, `hvg_select`, and `matrix_log_normalize`. It is not the
complete planned RNA pack.
Installation is inert; explicit scworkbench activation is required.

The pack uses operation-pack ABI v3 only. `cell_suffix_label` receives an
authenticated, ordered session-cell axis without a fake matrix input;
`rna_qc_metrics` owns mitochondrial and ribosomal selection and emits one
atomic bounded metric family; `hvg_select` reads typed matrix chunks plus an
authenticated matrix-bound `feature_type` snapshot. Core rejects missing,
foreign, duplicate, or reordered cells, incomplete optional metric groups,
and mutated or cross-matrix annotation evidence. `matrix_log_normalize`
streams authenticated matrix slices into a closed, fragment-backed output;
core alone publishes or atomically replaces the matrix and its terminal
summary during promotion.

Activate and invoke it with:

```r
wb_operation_pack_activate(system.file(
  "scworkbench", "operation-pack.json", package = "scwbRNA"
))
wb_session_run_op(session, bundle, "scwbrna::hvg_select", inputs = list(
  matrix = "RNA", n_features = 2000L, feature_type = "gene"
))
wb_session_run_op(session, bundle, "scwbrna::matrix_log_normalize",
  inputs = list(matrix = "RNA", output_matrix_id = "RNA_log1p",
    scale_factor = 10000, overwrite = FALSE))
```

The four packed methods preserve the pinned former-core defaults, terminal
primitive identities, promotion closure, replay equivalence, and effective
parameter records. Batch-aware HVG selection remains explicitly unsupported;
the provisional pack makes no claim to contain PCA, UMAP, or biological
annotation methods.
