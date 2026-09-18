library(scworkbench)
library(scwbRNA)
.scwbrna_test_root <- if (file.exists("test_rna_pack_contract.R")) "." else
  "tests"
source(file.path(.scwbrna_test_root, "test_rna_pack_contract.R"),
  local = FALSE)
tests <- sort(ls(pattern = "^test_", envir = .GlobalEnv))
failed <- character()
for (name in tests) {
  cat("RUN ", name, "\n", sep = "")
  result <- try(get(name, envir = .GlobalEnv)(), silent = TRUE)
  if (inherits(result, "try-error") || !isTRUE(result)) {
    cat("FAIL ", name, " -- ", as.character(result), "\n", sep = "")
    failed <- c(failed, name)
  } else cat("PASS ", name, "\n", sep = "")
}
cat("SCWBRNA TEST SUMMARY discovered=", length(tests), " executed=",
  length(tests), " passed=", length(tests) - length(failed),
  " skipped=0 failed=", length(failed), "\n", sep = "")
if (length(failed)) quit(status = 1L)
