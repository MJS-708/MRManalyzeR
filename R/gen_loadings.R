gen_loadings = function(data, header, PC = 1, analyte_type = "Oxylipin"){

  PC = sprintf("PC%s", PC)

  # Build a named colour vector for the PUFA groups
  lvls <- unique(data[[header]])
  col_cols <- setNames(scales::hue_pal()(length(lvls)), lvls)

  # Map each feature to a hex colour
  axis_cols <- rev( col_cols[ data[[header]][order(data[[PC]], decreasing = T)  ] ] )

  p1 <- ggplot(
    data,
    aes(
      x = reorder(feature, .data[[PC]]),   # <- use column `feature`, reorder by PC
      y = .data[[PC]]
    )
  ) +
    geom_col(fill = "steelblue") +

    geom_point(
      aes(y = 0, colour = .data[[header]]),
      size = 0.00000002
    ) +

    scale_colour_manual(values = col_cols, name = header) +
    guides(colour = guide_legend(override.aes = list(size = 4))) +
    coord_flip() +
    labs(
      title = sprintf("Loadings for %s", PC),
      x = analyte_type,
      y = sprintf("Loading (%s)", PC)
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(colour = axis_cols),
      legend.position = "right"
    )

  return(p1)


}
