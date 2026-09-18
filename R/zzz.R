.onAttach <- function(...) {
  conflicts <- cmstatrExt_conflicts()

  if (nrow(conflicts) > 0) {
    packageStartupMessage(cmstatrExt_conflict_message(conflicts))
  }
}
