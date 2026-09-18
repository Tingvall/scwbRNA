#' Register the provisional scwbRNA operation pack
#' @return A validated scworkbench operation-pack registration.
#' @export
register_operation_pack <- function() {
  scworkbench::wb_operation_pack_register(
    system.file("scworkbench", "operation-pack.json", package = "scwbRNA"),
    list(`scwbrna::cell_suffix_label` = .run_cell_suffix_label,
      `scwbrna::rna_qc_metrics` = .run_rna_qc_metrics,
      `scwbrna::hvg_select` = .run_hvg_select,
      `scwbrna::matrix_log_normalize` = .run_matrix_log_normalize))
}
