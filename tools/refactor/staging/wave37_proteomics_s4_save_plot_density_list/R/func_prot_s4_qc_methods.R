#' @export
savePlotDensityList <- function(input_list, prefix = "Density", suffix = c("png", "pdf"), output_dir) {
  list_of_filenames <- expand_grid(column = names(input_list), suffix = suffix) |>
    mutate(filename = paste0(prefix, "_", column, ".", suffix)) |>
    left_join(
      tibble(
        column = names(input_list),
        plots = input_list
      ),
      by = join_by(column)
    )

  purrr::walk2(
    list_of_filenames$plots,
    list_of_filenames$filename,
    \(.x, .y) {
      ggsave(plot = .x, filename = file.path(output_dir, .y))
    }
  )

  list_of_filenames
}

