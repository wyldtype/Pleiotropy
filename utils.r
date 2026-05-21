plotExpressionProfile <- function(.counts, .info, .gene_idxs, .gene_names) {
  countdf <- map2(.gene_idxs, .gene_names, \(g, nm) {
    outdf <- tibble(sample = colnames(.counts)[-1],
                    gene_vec = as.numeric(.counts[.counts$gene_name == g, -1]))
    names(outdf) <- c("sample_name", nm)
    return(outdf)
  }) |> purrr::reduce(.f = left_join, by = "sample_name")
  plotdf <- countdf |> pivot_longer(cols = setdiff(names(countdf), "sample_name"),
                                    names_to = "gene", values_to = "expr") |> 
    left_join(y = .info, by = "sample_name")
  ggplot(plotdf, aes(x = factor(hour), y = expr)) +
    geom_point(aes(color = factor(gene)),
               size = 0.25) +
    geom_line(aes(color = factor(gene),
                  group = interaction(gene, replicate),
                  linetype = factor(replicate))) +
    #scale_color_brewer(palette = "Set1") +
    scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
    labs(color = "Cluster") +
    theme_classic() +
    ylab("Expression (tpm)") +
    xlab("Time (hours)") +
    facet_wrap(~factor(gene, levels = .gene_names, 
                       labels = .gene_names))
}

plotExpressionProfileAllEnvironments <- function(.gene_idxs, .gene_names,
                                                 .environments = c("immune", "dev", "oog"),
                                                 .color_vec = scales::hue_pal()(length(.gene_idxs))) {
  plotlist <- map(.environments, \(.e) {
    .counts <- counts[[.e]]
    .info <- infodf[[.e]]
    countdf <- map2(.gene_idxs, .gene_names, \(g, nm) {
      outdf <- tibble(sample = colnames(.counts)[-1],
                      gene_vec = as.numeric(.counts[.counts$gene_name == g, -1]))
      names(outdf) <- c("sample_name", nm)
      return(outdf)
    }) |> purrr::reduce(.f = left_join, by = "sample_name")
    plotdf <- countdf |> pivot_longer(cols = setdiff(names(countdf), "sample_name"),
                                      names_to = "gene", values_to = "expr") |> 
      left_join(y = .info, by = "sample_name")
    p <- ggplot(plotdf, aes(x = factor(hour), y = expr)) +
      geom_point(aes(color = factor(gene)),
                 size = 0.25) +
      geom_line(aes(color = factor(gene),
                    group = interaction(gene, replicate),
                    linetype = factor(replicate))) +
      scale_color_discrete(type = .color_vec, limits = .gene_names) +
      scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
      labs(color = "Cluster") +
      theme_classic() +
      ylab("Expression (tpm)") +
      xlab("Time (hours)") +
      facet_wrap(~factor(gene, levels = .gene_names, 
                         labels = .gene_names)) +
      ggtitle(.e)
    return(p)
  })
  names(plotlist) <- .environments
  return(plotlist)
}
