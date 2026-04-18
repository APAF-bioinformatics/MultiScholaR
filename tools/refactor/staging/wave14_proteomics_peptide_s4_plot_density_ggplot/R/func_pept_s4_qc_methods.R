setMethod(f="plotDensity"
          , signature="ggplot2::ggplot"
          , definition=function(theObject, grouping_variable, title = "", font_size = 8) {
            # First try to get data directly from the ggplot object's data element
            if (!is.null(theObject$data) && is.data.frame(theObject$data)) {
              pca_data <- as_tibble(theObject$data)
            } else {
              # Fall back to other extraction methods
              pca_data <- as_tibble(ggplot_build(theObject)$data[[1]])
              
              # If the data doesn't have PC1/PC2, try to extract from the plot's environment
              if (!("PC1" %in% colnames(pca_data) && "PC2" %in% colnames(pca_data))) {
                # Try to get the data from the plot's environment
                if (exists("data", envir = environment(theObject$mapping$x))) {
                  pca_data <- as_tibble(get("data", envir = environment(theObject$mapping$x)))
                } else {
                  stop("Could not extract PCA data from the ggplot object")
                }
              }
            }
            
            # Check if grouping variable exists in the data
            if (!grouping_variable %in% colnames(pca_data)) {
              stop(sprintf("grouping_variable '%s' not found in the data", grouping_variable))
            }

            # Create PC1 density plot
            pc1_density <- ggplot(pca_data, aes(x = PC1, fill = !!sym(grouping_variable), color = !!sym(grouping_variable))) +
              geom_density(alpha = 0.3) +
              theme_bw() +
              labs(title = title,
                   x = "PC1",
                   y = "Density") +
              theme(
                text = element_text(size = font_size),
                plot.margin = margin(b = 0, t = 5, l = 5, r = 5),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank()
              )

            # Create PC2 density plot
            pc2_density <- ggplot(pca_data, aes(x = PC2, fill = !!sym(grouping_variable), color = !!sym(grouping_variable))) +
              geom_density(alpha = 0.3) +
              theme_bw() +
              labs(x = "PC2",
                   y = "Density") +
              theme(
                text = element_text(size = font_size),
                plot.margin = margin(t = 0, b = 5, l = 5, r = 5),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank()
              )
            
            # Combine plots with minimal spacing
            combined_plot <- pc1_density / pc2_density +
              plot_layout(heights = c(1, 1)) +
              plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))
            
            return(combined_plot)
          }) 

