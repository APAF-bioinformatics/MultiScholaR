# enrichedPathwayBarPlot
# ----------------------------------------------------------------------------
#' @export
enrichedPathwayBarPlot <- function( input_table, input_go_type = NA, remove_duplicted_entries = "merge", added_columns = "comparison") {

  duplicated_entries <- input_table %>%
    mutate(set_type = case_when( str_detect( gene_set, "positive") ~"positive",
                                 str_detect( gene_set, "negative") ~ "negative",
                                 TRUE ~ "neutral")) %>%
    group_by( across(c( any_of(added_columns), comparison, set_type, annotation_id) ) ) %>%
    dplyr::summarise( temp_qvalue = min(qvalue )) %>%
    ungroup() %>%
    dplyr::group_by( across(c( any_of(added_columns), comparison, annotation_id) ) ) %>%
    dplyr::summarise(counts = n(),
                     best_p_adj_value = min(temp_qvalue)) %>%
    ungroup() %>%
    dplyr::filter( counts > 1)

  if( remove_duplicted_entries == TRUE |
      remove_duplicted_entries == "delete" ) {
    input_table <- input_table  %>%
      anti_join( duplicated_entries, by =c("comparison" = "comparison",
                                           "annotation_id" = "annotation_id"))
  } else if( remove_duplicted_entries == "merge" ) {

    duplicates_tbl <- input_table %>%
      inner_join( duplicated_entries, by =c("comparison" = "comparison",
                                            "annotation_id" = "annotation_id")) %>%
      dplyr::filter( qvalue == best_p_adj_value ) %>%
      mutate( gene_set = "shared" )

    input_table <- input_table  %>%
      anti_join( duplicated_entries, by =c("comparison" = "comparison",
                                           "annotation_id" = "annotation_id")) %>%
      bind_rows( duplicates_tbl )

  }


  if (!is.na(input_go_type )) {
    if(! (input_go_type %in% (input_table %>% distinct(go_type) %>% dplyr::pull(go_type)) ) ) {
      stop( paste0( "input_go_type = ", input_go_type, ", is not in the input_table" ) )
    }

    filt_input_table <-   input_table %>%
      dplyr::filter( go_type == input_go_type )
  } else {
    filt_input_table <- input_table
  }


  bar_plot_data <- filt_input_table  %>%
    mutate( neg_log10_qvalue =  -log10(qvalue))  %>%
    arrange( neg_log10_qvalue) %>%
    mutate( gene_set = case_when( gene_set == "positive_list" ~ "Up-regulated"
                                  , gene_set == "negative_list" ~ "Down-regulated"
                                  , gene_set == "shared" ~ "Both")) %>%
    group_by(gene_set, term ) %>%
    summarise( neg_log10_qvalue = max(neg_log10_qvalue)) %>%
    ungroup()

  term_ordering <- bar_plot_data %>%
    dplyr::select( gene_set, neg_log10_qvalue, term) %>%
    distinct(gene_set, neg_log10_qvalue, term) %>%
    arrange( (gene_set), (neg_log10_qvalue), term ) %>%
    dplyr::pull(term )

  x_label <- "Pathway"
  if( !is.na(input_go_type) ) {
    x_label <- "GO Term"
  }

  bar_plot_data %>%
    mutate ( term = factor(term, levels= term_ordering)) %>%
    mutate ( gene_set = factor( gene_set, levels =c( "Up-regulated"
                                                     , "Down-regulated"
                                                     , "Both"))) %>%
    arrange( (gene_set), (neg_log10_qvalue) ) %>%
    ggplot( aes( term, neg_log10_qvalue, fill=gene_set  )) +
    geom_col() +
    scale_fill_manual(values = c( "Up-regulated" = "#F8766D", "Down-regulated" = "#00BFC4", "Both" = "grey")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))  +
    #facet_grid(  gene_set ~ ., scales = "free_y", space = "free_y") +
    coord_flip() +
    #theme(strip.text.y.right = element_text(angle = 0)) +
    labs(fill='Query Proteins') +
    ylab( expression("Significance, -log"[10]*"(q-value)") ) +
    xlab( x_label)

}


# ----------------------------------------------------------------------------
# enrichedGoTermBarPlot
# ----------------------------------------------------------------------------
#'@title Draw list of functional enrichment heatmaps
#'@description given input table, draw a bar plot representing the GO enrichment results.
#'The height of each bar represents the negative log (base 10) q-values of the query proteins.
#'@export
#' @param input_table,output_dir,analysis_type,file_suffix,width,height Runtime inputs used by this function; see the usage section for accepted values.
enrichedGoTermBarPlot <- function( input_table, output_dir,
                                   analysis_type = "GO", file_suffix, width=10, height = 7) {

  partial_go_term_bar_plot <- partial( enrichedPathwayBarPlot,
                                       input_table = input_table)

  list_of_go_type <- input_table %>%
    distinct( go_type) %>%
    arrange(go_type) %>%
    dplyr::pull(go_type)

  list_of_barplots <- purrr::map( list_of_go_type,
                                  ~partial_go_term_bar_plot(input_go_type = .))

  names( list_of_barplots) <- list_of_go_type

  suffix_list <- file_suffix

  plotGOBarPlotWithSuffix <- function(suffix, input_plot, plot_name, width, height) {
    ggsave(plot = input_plot,
           filename=file.path( output_dir,
                               paste0(analysis_type, "_",
                                      plot_name,
                                      "_bar_plot.",
                                      suffix )),
           width= width,
           height = height) }

  plotBarPlot <- function(plot_name, input_plot, width, height ) {
    purrr::walk( suffix_list,
                 ~plotGOBarPlotWithSuffix(., input_plot, plot_name, width, height) ) }

  purrr::walk2(names( list_of_barplots), list_of_barplots,
               ~plotBarPlot(.x, .y, width=width, height=height) )


}


# ----------------------------------------------------------------------------
# createWordCloudDataFrame
# ----------------------------------------------------------------------------
#'@title Create a word frequency distribution table for Word Cloud generation.
#'@description Create a word frequency distribution table for Word Cloud generation.
#'Based on article by Celine Van den Rul, How to Generate Word Clouds in R, Simple Steps on How and When to Use Them,
#' https://towardsdatascience.com/create-a-word-cloud-with-r-bde3e7422e8a (accessed 7th November 2022)
#'@export
#'@param text_list, a vector of text (e.g. a list of GO terms name)
createWordCloudDataFrame <- function( text_list) {

  docs <- tm::Corpus(tm::VectorSource(text_list))

  docs <- docs %>%
    tm::tm_map(tm::removeNumbers) %>%
    tm::tm_map(tm::removePunctuation) %>%
    tm::tm_map(tm::stripWhitespace)
  docs <- tm::tm_map(docs, tm::content_transformer(tolower))
  docs <- tm::tm_map(docs, tm::removeWords, tm::stopwords("english"))

  dtm <- tm::TermDocumentMatrix(docs)
  matrix <- as.matrix(dtm)
  words <- sort(rowSums(matrix),decreasing=TRUE)
  df <- data.frame(word = names(words),freq=words)
}


# ----------------------------------------------------------------------------
# cleanDuplicatesEnrichment
# ----------------------------------------------------------------------------
#'@export
cleanDuplicatesEnrichment <- function( input_table
                                       , pathway_column = term
                                       , fdr_column = qvalue
                                       , gene_set_column = gene_set) {

  duplicated_entries <- input_table |>
    group_by( across(c( {{gene_set_column}}, {{pathway_column}}) ) ) |>
    dplyr::summarise( temp_qvalue = min({{fdr_column}} )) |>
    ungroup() |>
    dplyr::group_by( across(c( {{pathway_column}}) ) ) |>
    dplyr::summarise(counts = n(),
                     best_p_adj_value = min(temp_qvalue)) |>
    ungroup() |>
    dplyr::filter( counts > 1)

  duplicates_tbl <- input_table |>
    inner_join( duplicated_entries, by =join_by( {{pathway_column}} )) |>
    dplyr::filter( qvalue == best_p_adj_value ) |>
    mutate( gene_set = "Both" )

  input_table_cln <- input_table  %>%
    anti_join( duplicated_entries, by = join_by( {{pathway_column}} ) ) |>
    bind_rows( duplicates_tbl )

  positive_label <- input_table_cln |>
    distinct({{gene_set_column}}) |>
    dplyr::filter( str_detect( {{gene_set_column}}, "[P|p]ositive")) |>
    dplyr::pull({{gene_set_column}})
  negative_label <- input_table_cln |>
    distinct({{gene_set_column}}) |>
    dplyr::filter( str_detect( {{gene_set_column}}, "[N|n]egative")) |>
    dplyr::pull({{gene_set_column}})
  both_label <- input_table_cln |>
    distinct({{gene_set_column}}) |>
    dplyr::filter( str_detect( {{gene_set_column}}, "[B|b]oth")) |>
    dplyr::pull({{gene_set_column}})

  proteomics_go_helper <- input_table_cln |>
    mutate( neg_log_10_fdr = -log10({{fdr_column}}))    |>
    mutate( {{gene_set_column}} := factor( {{gene_set_column}}, levels=c( positive_label
                                                                          , negative_label
                                                                          , both_label ))) |>
    arrange( {{gene_set_column}}, desc({{gene_set_column}}), desc(neg_log_10_fdr))


   proteomics_go_helper
}


# ----------------------------------------------------------------------------
# plotEnrichmentBarplot
# ----------------------------------------------------------------------------
#'@export
plotEnrichmentBarplot <- function( input_table
                                   , pathway_column = term
                                   , fdr_column = qvalue
                                   , gene_set_column = gene_set
                                   , xlab_string = expression(-log[10](FDR))
                                   , ylab_string = "Enriched Terms"
                                   , legend_title = "Gene Set"
                                   , legend_colours = c( "red", "blue", "black")) {

  proteomics_go_helper <- cleanDuplicatesEnrichment( input_table
                                                                , pathway_column = {{pathway_column}}
                                                                , fdr_column = {{fdr_column}}
                                                                , gene_set_column = {{gene_set_column}} )

  proteomics_go_helper2 <- proteomics_go_helper |>
    mutate( {{pathway_column}} := factor( {{pathway_column}}, levels = rev( proteomics_go_helper |> distinct({{pathway_column}}) |> dplyr::pull( {{pathway_column}}) )) )

  # proteomics_go_helper_2 <- proteomics_go_helper |>
  #   mutate( {{pathway_column}} := factor( {{pathway_column}}, levels = rev( unique( proteomics_go_helper |> dplyr::pull( {{pathway_column}}) ))  ) )

  output_barplot <- proteomics_go_helper2 |>
    ggplot( aes( neg_log_10_fdr,   {{pathway_column}},  fill={{gene_set_column}})) +
    geom_bar( stat="identity" ) +
    scale_fill_manual(legend_title, values = legend_colours) +
    xlab ( xlab_string) +
    ylab ( ylab_string) +
    theme_bw() +
    scale_x_continuous(limits = c(0,-log10( min ( proteomics_go_helper|> dplyr::pull({{fdr_column}}) ) )*1.05), expand = c(0, 0))
  #+
    #facet_grid( rows = vars({{gene_set_column}})  , scales="free_y", space = "free_y")

  output_barplot
}


# ----------------------------------------------------------------------------
# list2df
# ----------------------------------------------------------------------------
#' @export
list2df <- function(inputList) {
  # ldf <- lapply(1:length(inputList), function(i) {
  ldf <- lapply(seq_len(length(inputList)), function(i) {
    data.frame(categoryID=rep(names(inputList[i]),
                              length(inputList[[i]])),
               Gene=inputList[[i]])
  })

  do.call('rbind', ldf)
}


# ----------------------------------------------------------------------------
# list2graph
# ----------------------------------------------------------------------------
#' @importFrom igraph graph.data.frame V E
#' @export
list2graph <- function(inputList) {
  x <- list2df(inputList)
  g <- graph.data.frame(x, directed=FALSE)
  return(g)
}


# ----------------------------------------------------------------------------
# get_param_change_message
# ----------------------------------------------------------------------------
#' @export
get_param_change_message <- function(parameter, params_df) {
  paste0("Use '", params_df[parameter, "listname"],
         " = list(", params_df[parameter, "present"],
         " = your_value)' instead of '", params_df[parameter, "original"],
         "'.\n The ", params_df[parameter, "original"],
         " parameter will be removed in the next version.")
}


# ----------------------------------------------------------------------------
# node_add_alpha
# ----------------------------------------------------------------------------
#' @export
node_add_alpha <- function(p, hilight_category, hilight_gene, alpha_nohilight, alpha_hilight) {
  alpha_node <- rep(1, nrow(p$data))
  if (!is.null(hilight_category)) {
    alpha_node <- rep(alpha_nohilight, nrow(p$data))
    hilight_node <- c(hilight_category, hilight_gene)
    alpha_node[match(hilight_node, p$data$name)] <- alpha_hilight
  }
  p$data$alpha <- alpha_node
  return(p)
}


# ----------------------------------------------------------------------------
# get_enrichplot_color
# ----------------------------------------------------------------------------
#' @export
get_enrichplot_color <- function(n = 2) {
  colors <- getOption("enrichplot.colours")
  if (!is.null(colors)) return(colors)

  if (n != 2 && n != 3) stop("'n' should be 2 or 3")

  colors = c("#e06663", "#327eba")
  if (n == 2) return(colors)

  if (n == 3) return(c(colors[1], "white", colors[2]))
}


# ----------------------------------------------------------------------------
# set_enrichplot_color
# ----------------------------------------------------------------------------
#' @export
set_enrichplot_color <- function(colors = get_enrichplot_color(2),
                                 type = "color", name = NULL, .fun = NULL, ...) {

  type <- match.arg(type, c("color", "colour", "fill"))

  n <- length(colors)
  if (n < 2) {
    stop("'colors' should be of length >= 2")
  } else if (n == 2) {
    params <- list(low = colors[1], high = colors[2])
    fn_suffix <- "continuous"
  } else if (n == 3) {
    params <- list(low = colors[1], mid = colors[2], high = colors[3])
    fn_suffix <- "gradient2"
  } else {
    params <- list(colors = colors)
    fn_suffix <- "gradientn"
  }

  if (!is.null(.fun)) {
    if (n == 3) {
      # should determine parameter for user selected functions: 'gradient2' or 'gradientn'
      scale_parameters <- names(formals(.fun))
      fn_type <- if ("colours" %in% scale_parameters || "colors" %in% scale_parameters) {
        "gradientn"
      } else if ("mid" %in% scale_parameters) {
        "gradient2"
      } else {
        "continuous"
      }
      if (fn_type == "gradientn") {
        params <- list(colors = colors)
      } else {
        params <- list(low = colors[1], mid = colors[2], high = colors[3])
      }
    }
  } else {
    fn <- sprintf("scale_%s_%s", type, fn_suffix)
    .fun <- getFromNamespace(fn, "ggplot2")
  }

  params$guide <- guide_colorbar(reverse=TRUE, order=1)
  params$name <- name # no legend name setting by default as 'name = NULL'

  params <- modifyList(params, list(...))

  do.call(.fun, params)
}


# ----------------------------------------------------------------------------
# add_node_label
# ----------------------------------------------------------------------------
#' @importFrom ggraph geom_node_text geom_edge_arc geom_edge_link ggraph geom_node_point
#' @export
add_node_label <- function(p, data, label_size_node, cex_label_node, shadowtext) {
  # If use 'aes_(alpha =~I(alpha))' will put an error for AsIs object.
  # But I(alpha) is necessory, so use 'alpha = I(data$alpha)'.
  segment.size <- get_ggrepel_segsize()
  if (is.null(data)) {
    data <- p$data
  }
  if (shadowtext) {
    p <- p + geom_node_text(aes_(label=~name), data = data,
                            alpha = I(data$alpha),
                            size = label_size_node * cex_label_node, bg.color = "white",
                            repel=TRUE, segment.size = segment.size)
  } else {
    p <- p + geom_node_text(aes_(label=~name), data = data,
                            alpha = I(data$alpha),
                            size = label_size_node * cex_label_node, repel=TRUE,
                            segment.size = segment.size)
  }
  return(p)
}


# ----------------------------------------------------------------------------
# get_ggrepel_segsize
# ----------------------------------------------------------------------------
#' @export
get_ggrepel_segsize <- function(default = 0.2) {
  getOption("ggrepel.segment.size", default = default)
}


# ----------------------------------------------------------------------------
# cnetplotEdited
# ----------------------------------------------------------------------------
#' @export
cnetplotEdited <- function(
    geneSets,
    showCategory = 5,
    foldChange = NULL,
    layout = "kk",
    colorEdge = FALSE,
    circular = FALSE,
    node_label = "all",
    cex_category = 1,
    cex_gene = 1,
    cex_label_category = 1,
    cex_label_gene = 1,
    color_category = "#E5C494",
    color_gene = "#B3B3B3",
    shadowtext = "all",
    color.params=list(
      foldChange = NULL,
      edge = FALSE,
      category = "#E5C494",
      gene = "#B3B3B3"
    ),
    cex.params=list(
      category_node = 1,
      gene_node = 1,
      category_label = 1,
      gene_label = 1
    ),
    hilight.params=list(
      category = NULL,
      alpha_hilight = 1,
      alpha_no_hilight = 0.3
    ),
    ...) {

  label_size_category <- 5
  label_size_gene <- 5
  node_label <- match.arg(node_label, c("category", "gene", "all", "none"))

  params_df <- as.data.frame(rbind(
    c("foldChange", "color.params", "foldChange"),
    c("colorEdge", "color.params", "edge"),
    c("color_category", "color.params", "category"),
    c("color_gene", "color.params", "gene"),

    c("cex_category", "cex.params", "category_node"),
    c("cex_gene", "cex.params", "gene_node"),
    c("cex_label_category", "cex.params", "category_label"),
    c("cex_label_gene", "cex.params", "gene_label"))
  )
  colnames(params_df) <- c("original", "listname", "present")
  rownames(params_df) <- params_df$original

  default.color.params <- list(
    foldChange = NULL,
    edge = FALSE,
    category = "#E5C494",
    gene = "#B3B3B3"
  )
  default.cex.params <- list(
    category_node = 1,
    gene_node = 1,
    category_label = 1,
    gene_label = 1
  )
  default.hilight.params <- list(
    category = NULL,
    alpha_hilight = 1,
    alpha_no_hilight = 0.3
  )

  # use modifyList to change the values of parameter
  color.params <- modifyList(default.color.params, color.params)
  cex.params <- modifyList(default.cex.params, cex.params)
  hilight.params <- modifyList(default.hilight.params, hilight.params)
  params_list <- list( showCategory = showCategory,
                       foldChange = foldChange,
                       layout = layout,
                       colorEdge = colorEdge,
                       circular = circular,
                       node_label = node_label,
                       cex_category = cex_category,
                       cex_gene = cex_gene,
                       cex_label_category = cex_label_category,
                       cex_label_gene = cex_label_gene,
                       color_category = color_category,
                       color_gene = color_gene,
                       shadowtext = shadowtext,
                       color.params = color.params,
                       cex.params = cex.params,
                       hilight.params = hilight.params
  )

  # get all parameters value
  args <- as.list(match.call())
  removed_params <- intersect(params_df$original, names(args))
  if (length(removed_params) > 0) {
    for (i in removed_params) {
      params_list[[params_df[i, 2]]][[params_df[i, 3]]] <- get(i)
      warn <- get_param_change_message(i, params_df)
      warning(warn)
    }
  }

  color.params <- params_list[["color.params"]]
  cex.params <- params_list[["cex.params"]]
  hilight.params <- params_list[["hilight.params"]]

  foldChange <- color.params[["foldChange"]]
  colorEdge <- color.params[["edge"]]
  color_category <- color.params[["category"]]
  color_gene <- color.params[["gene"]]

  cex_category <- cex.params[["category_node"]]
  cex_gene <- cex.params[["gene_node"]]
  cex_label_category <- cex.params[["category_label"]]
  cex_label_gene <- cex.params[["gene_label"]]

  hilight_category <- hilight.params[["category"]]
  alpha_hilight <- hilight.params[["alpha_hilight"]]
  alpha_nohilight <- hilight.params[["alpha_no_hilight"]]

  if (circular) {
    layout <- "linear"
    geom_edge <- geom_edge_arc
  } else {
    geom_edge <- geom_edge_link
  }
  if (is.logical(shadowtext)) {
    shadowtext <- ifelse(shadowtext, "all", "none")
  }
  shadowtext_category <- shadowtext_gene <- FALSE
  if (shadowtext == "all") shadowtext_category <- shadowtext_gene <- TRUE
  if (shadowtext == "category") shadowtext_category <- TRUE
  if (shadowtext == "gene") shadowtext_gene <- TRUE

  g <- list2graph(geneSets)

  # if (!inherits(x,  "list")) {
  #     foldChange <- fc_readable(x, foldChange)
  # }

  size <- sapply(geneSets, length)
  V(g)$size <- min(size)/2
  n <- length(geneSets)
  V(g)$size[1:n] <- size
  node_scales <- c(rep(cex_category, n), rep(cex_gene, (length(V(g)) - n)))

  # add edge alpha
  hilight_category <- intersect(hilight_category, names(geneSets))

  if (!is.null(hilight_category) && length(hilight_category) > 0) {
    edges <- attr(E(g), "vnames")
    E(g)$alpha <- rep(alpha_nohilight, length(E(g)))
    hilight_edge <- grep(paste(hilight_category, collapse = "|"), edges)
    hilight_gene <- edges[hilight_edge]
    hilight_gene <- gsub(".*\\|", "", hilight_gene)
    E(g)$alpha[hilight_edge] <- min(0.8, alpha_hilight)
  } else {
    E(g)$alpha <- rep(0.8, length(E(g)))
  }

  show_legend <- c(FALSE, TRUE)
  names(show_legend) <- c("alpha", "color")
  if (colorEdge) {
    E(g)$category <- rep(names(geneSets), sapply(geneSets, length))
    edge_layer <- geom_edge(aes_(color = ~category, alpha = ~I(alpha)),
                            show.legend = show_legend)
  } else {
    edge_layer <- geom_edge(aes_(alpha = ~I(alpha)), colour='darkgrey',
                            show.legend = FALSE)
  }

  if (!is.null(foldChange)) {
    fc <- foldChange[V(g)$name[(n+1):length(V(g))]]
    V(g)$color <- NA
    # V(g)$color[1:n] <- color_category
    V(g)$color[(n+1):length(V(g))] <- fc
    show_legend <- c(TRUE, FALSE)
    names(show_legend) <- c("color", "size")
    p <- ggraph(g, layout=layout, circular = circular)
    p$data[-(1:n), "size"] <- 3 * cex_gene

    p <- node_add_alpha(p, hilight_category, hilight_gene, alpha_nohilight, alpha_hilight)

    alpha_category <- c(rep(1, n), rep(0, nrow(p$data)-n))
    alpha_gene <- c(rep(0, n), rep(1, nrow(p$data)-n))

    if (!is.null(hilight_category) && length(hilight_category) > 0) {
      alpha_category <- c(rep(alpha_nohilight, n), rep(0, nrow(p$data)-n))
      alpha_gene <- c(rep(0, n), rep(alpha_nohilight, nrow(p$data)-n))
      alpha_gene[match(hilight_gene, p$data$name)] <- alpha_hilight
      alpha_gene[match(hilight_category, p$data$name)] <- alpha_hilight
    }

    p <- p + edge_layer +
      geom_node_point(aes_(size=~size), color=I(color_category),
                      data = NULL, show.legend = show_legend,
                      alpha = I(alpha_category)) +
      ggnewscale::new_scale_color() +
      geom_node_point(aes_(color=~as.numeric(as.character(color)), size=~size),
                      data = NULL, alpha = I(alpha_gene)) +
      scale_size(range=c(3, 8) * cex_category) +
      # scale_colour_gradient2(name = "fold change") +
      set_enrichplot_color(colors = rev(get_enrichplot_color(3)), name = "fold change")


  } else {
    V(g)$color <- color_gene
    V(g)$color[1:n] <- color_category
    p <- ggraph(g, layout=layout, circular=circular)
    p$data[-(1:n), "size"] <- 3 * cex_gene
    p <- node_add_alpha(p, hilight_category, hilight_gene, alpha_nohilight, alpha_hilight)
    p <- p + edge_layer +
      geom_node_point(aes_(color=~I(color), size=~size, alpha=~I(alpha)))+
      scale_size(range=c(3, 8) * cex_category)

  }

  p <- p + theme_void()

  if (node_label == "category") {
    p$data[-c(1:n), "name"] <- NA
    p <- add_node_label(p = p, data = NULL, label_size_node = label_size_category,
                        cex_label_node = cex_label_category, shadowtext = shadowtext_category)
  } else if (node_label == "gene") {
    p$data[1:n, "name"] <- NA
    p <- add_node_label(p = p, data = NULL, label_size_node = label_size_gene,
                        cex_label_node = cex_label_gene, shadowtext = shadowtext_gene)
  } else if (node_label == "all") {
    p <- add_node_label(p = p, data = NULL,
                        label_size_node = c(rep(label_size_category, n), rep(label_size_gene, nrow(p$data)-n)),
                        cex_label_node = c(rep(cex_label_category, n), rep(cex_label_gene, nrow(p$data)-n)),
                        shadowtext = shadowtext_gene)
  }
  if (!is.null(foldChange)) {
    p <- p + guides(size  = guide_legend(order = 1),
                    color = guide_colorbar(order = 2))
  }
  return(p + guides(alpha = "none"))
}


# ----------------------------------------------------------------------------
# fc_readable
# ----------------------------------------------------------------------------
fc_readable <- function(x, foldChange = NULL) {
  if (is.null(foldChange))
    return(NULL)

  if (x@readable && x@keytype != "SYMBOL") {
    gid <- names(foldChange)
    if (is(x, 'gseaResult')) {
      ii <- gid %in% names(x@geneList)
    } else {
      ii <- gid %in% x@gene
    }
    gid[ii] <- x@gene2Symbol[gid[ii]]
    names(foldChange) <- gid
  }
  return(foldChange)
}


# ----------------------------------------------------------------------------
# update_n
# ----------------------------------------------------------------------------
#' @export
update_n <- function(x, showCategory) {
  if (!is.numeric(showCategory)) {
    if (inherits(x, 'list')) {
      showCategory <- showCategory[showCategory %in% names(x)]
    } else {
      showCategory <- intersect(showCategory, x$Description)
    }
    return(showCategory)
  }

  ## geneSets <- geneInCategory(x) ## use core gene for gsea result
  n <- showCategory
  if (inherits(x, 'list')) {
    nn <- length(x)
  } else {
    nn <- nrow(x)
  }
  if (nn < n) {
    n <- nn
  }

  return(n)
}


# ----------------------------------------------------------------------------
# extract_geneSets
# ----------------------------------------------------------------------------
#' @export
extract_geneSets <- function(x, n) {
  n <- update_n(x, n)

  if (inherits(x, 'list')) {
    geneSets <- x
  } else {
    geneSets <- clusterProfiler::geneInCategory(x) ## use core gene for gsea result
    y <- as.data.frame(x)
    geneSets <- geneSets[y$ID]
    names(geneSets) <- y$Description
  }
  if (is.numeric(n)) {
    return(geneSets[1:n])
  }
  return(geneSets[n]) ## if n is a vector of Description
}


# ----------------------------------------------------------------------------
