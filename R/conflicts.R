
cmstatrExt_conflicts <- function() {  # nolint
  envs <- grep("^package:", search(), value = TRUE)
  envs <- rlang::set_names(envs)

  all_obs <- lapply(envs, function(e) ls(pos = e))
  all_obs <- utils::stack(all_obs)
  all_obs$ind <- as.character(all_obs$ind)

  cmstatr_ext_obs <- all_obs[all_obs$ind == "package:cmstatrExt", ]

  # ignore some generics
  ignore <- c("augment")
  cmstatr_ext_obs <- cmstatr_ext_obs[!cmstatr_ext_obs$values %in% ignore, ]

  other_obs <- all_obs[all_obs$ind != "package:cmstatrExt", ]

  is_conflict <- sapply(
    other_obs$values,
    function(x) sum(x == cmstatr_ext_obs$values) > 0
  )

  conflicts <- other_obs[is_conflict, ]

  conflicts
}

cmstatrExt_conflict_message <- function(conflicts) {  # nolint
  header <- cli::rule(
    left = cli::style_bold("Conflicts"),
    right = "cmstatrExt_conflicts"
  )

  funs <- conflicts$values
  pkgs <- sapply(conflicts$ind, function(x) gsub("^package:", "", x))

  bullets <- paste0(
    cli::col_red(cli::symbol$cross),
    " ",
    paste0(cli::col_blue("cmstatrExt::"), cli::col_green(funs), "()"),
    " masks ",
    paste0(cli::col_blue(pkgs), "::", funs, "()"),
    collapse = "\n"
  )

  paste0(
    header,
    "\n",
    bullets,
    "\n"
  )
}
