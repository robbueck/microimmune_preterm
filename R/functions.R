# Function to generate circular plots by sign
plot_circular_associations <- function(df_subset, title_text, text_size = 1, type_colors) {
  
  df_subset <- df_subset %>% 
    mutate(feature = as.character(feature),
           metaVariable = as.character(metaVariable))
  
  # Nodes dataframe
  nodes <- unique(
    rbind(
      df_subset %>% select(node = feature, type = feature_type, group = feature_group),
      df_subset %>% select(node = metaVariable, type = metavar_type, group = metavar_type)
    )
  )
  nodes <- nodes %>% 
    arrange(type)
  
  # Colors for each node type
  # unique_types <- unique(nodes$type)
  # type_colors <- viridisLite::turbo(length(unique_types))
  # names(type_colors) <- unique_types
  node_colors <- type_colors[nodes$type]
  names(node_colors) <- nodes$node
  
  node_groups <- nodes$group
  names(node_groups) <- nodes$node
  
  
  
  # Edges dataframe
  edges <- df_subset %>% 
    select(from = feature, to = metaVariable, Ds, status) %>%
    mutate(status = ifelse(status %in% c("AD", "OK_sd", "OK_nc", "OK_d"), yes = 1, no = 2))
  
  # Edge color gradient (positive: blue, negative: red)
  edge_col_fun <- if (all(df_subset$Ds > 0)) {
    colorRamp2(range(edges$Ds), c("#ffcccc", "#990000"))
  } else {
    colorRamp2(range(edges$Ds), c("#cce5ff", "#003399"))
  }
  
  # Initialize and plot circular diagram
  circos.clear()
  circos.par(gap.after = rep(1, nrow(nodes)), start.degree = 90, canvas.ylim = c(-1.5, 1.5),
             canvas.xlim = c(-1.5, 1.5))
  
  chordDiagram(edges[, c("from", "to")], 
               order = nodes$node,
               grid.col = node_colors,
               group = node_groups, 
               transparency = 0.4,
               col = edge_col_fun(edges$Ds),
               annotationTrack = "grid", 
               link.border = "black",
               link.lwd = 2,
               link.lty = edges$status,
               # directional = 1,
               # direction.type = "arrows",
               # link.arr.length = 0,
               # link.arr.col = edges[, c("from", "to")] %>% mutate(color = "black"), 
               preAllocateTracks = 1)
  
  circos.trackPlotRegion(track.index = 1, bg.border = NA, 
                         panel.fun = function(x, y) {
                           xlim <- get.cell.meta.data("xlim")
                           sector.name <- get.cell.meta.data("sector.index")
                           circos.text(mean(xlim), 0.1, 
                                       sector.name, 
                                       facing = "clockwise", 
                                       niceFacing = TRUE, 
                                       adj = c(0, 0.5),
                                       cex=text_size)
                         })
  
  title(title_text, line = 0)
  
  legend("bottomright", legend = names(type_colors), fill = type_colors,
         border = "black", bty = "n", cex = 1.2, title = "Data types")
  legend("bottom", legend = c("deconfounded", "confounded"), lty = c(1, 2),
         lwd = 2, bty = "n", cex = 1.2, title = "Confounding status")
  legend("bottomleft", legend = c(round(edges$Ds[which.min(abs(edges$Ds))], 2), 
                                  round(median(edges$Ds), 2), 
                                  round(edges$Ds[which.max(abs(edges$Ds))], 2)), 
         fill = if(all(df_subset$Ds > 0)) c("#ffcccc",edge_col_fun(median(edges$Ds)), "#990000") else c("#cce5ff", edge_col_fun(median(edges$Ds)), "#003399"),
         border = "black", bty = "n", cex = 1.2, title = "Effect size")
}




# plot a PCoA ordination
plot_pcoa <- function(ps_object, ordination,
                      color = NULL,
                      shape = NULL,
                      ellipses = F,
                      title = NULL,
                      alpha = 1,
                      size = 1,
                      axes = 1:2,
                      plot_dens = F){
  DF <- plot_ordination(ps_object, ordination, justDF = T, axes = axes)
  if (!is.null(color)) {
    if (!color %in% names(DF)) {
      warning("Color variable was not found in the available data you provided.", 
              "No color mapped.")
      color <- NULL
    }
  }
  if (!is.null(shape)) {
    if (!shape %in% names(DF)) {
      warning("Shape variable was not found in the available data you provided.", 
              "No shape mapped.")
      shape <- NULL
    }
  }
  x = colnames(DF)[1]
  y = colnames(DF)[2]
  if (ncol(DF) <= 2) {
    message("No available covariate data to map on the points for this plot `type`")
    ord_map = aes_string(x = x, y = y)
  } else {ord_map = aes_string(x = x, y = y, color = color, shape = shape, 
                               na.rm = TRUE)
  }
  p <- ggplot(DF, ord_map) + geom_point(na.rm = TRUE, alpha = alpha, size = size)
  if (!is.null(title)) {
    p = p + ggtitle(title)
  }
  if (length(ordination$values$Eigenvalues[axes]) > 0) {
    eigvec = ordination$values$Eigenvalues
    fracvar = eigvec[axes]/sum(eigvec)
    percvar = round(100 * fracvar, 1)
    strivar = as(c("PCO 1", "PCO 2"), "character")
    strivar = paste0(strivar, "   [", percvar, "%]")
    p = p + xlab(strivar[1]) + ylab(strivar[2])
  }
  if (ellipses) {
    p <- p + stat_ellipse(ord_map)
  }
  if (plot_dens) {
    xdens <- cowplot::axis_canvas(p, axis = "x")+
      geom_density(data = DF, aes(x = Axis.1, fill = !!sym(color)),
                   alpha = 0.7, size = 0.2)
    
    ydens <- cowplot::axis_canvas(p, axis = "y", coord_flip = TRUE)+
      geom_density(data = DF, aes(x = Axis.2, fill = !!sym(color)),
                   alpha = 0.7, size = 0.2)+
      coord_flip()
    
    p <- insert_xaxis_grob(p, xdens, grid::unit(.2, "null"), position = "top")
    p <- insert_yaxis_grob(p, ydens, grid::unit(.2, "null"), position = "right")
    p <- ggdraw(p)
  }
  return(p)
}



merge_taxa2_fixed <- function (x, taxa = NULL, pattern = NULL, name = "Merged")
{
  if (is.null(taxa) && is.null(pattern)) {
    return(x)
  }
  if (!is.null(pattern)) {
    if (!is.null(taxa)) {
      mytaxa <- taxa
    }
    else {
      mytaxa <- taxa(x)
    }
    if (length(grep(pattern, mytaxa)) == 0) {
      return(x)
    }
    mytaxa <- mytaxa[grep(pattern, mytaxa)]
  }
  else if (is.null(taxa)) {
    mytaxa <- taxa(x)
  }
  else {
    mytaxa <- taxa
  }
  x2 <- phyloseq::merge_taxa(x, mytaxa, 1)
  mytaxa <- gsub("\\)", "\\\\)", gsub("\\(", "\\\\(", mytaxa))
  taxa_names(x2) <- gsub(mytaxa[[1]], name, taxa_names(x2))
  tax_table(x2)[1, ] <- rep(name, ncol(tax_table(x2)))
  x2
}



BuildHeatmap_2 <- function (metaDeconfOutput, q_cutoff = 0.1, d_cutoff = 0.01,
                            cuneiform = FALSE, coloring = 0, showConfounded = TRUE,
                            intermedData = FALSE, featureNames = NULL, metaVariableNames = NULL,
                            d_range = "full", d_col = c("blue", "white", "red"), keepMeta = NULL,
                            keepFeature = NULL, trusted = c("OK_sd", "OK_nc", "OK_d",
                                                            "AD"), tileBordCol = "black", reOrder = "both", plotPartial = "Ds",
                            starSize = 2, starNudge_y = 0)
{
  if (length(d_col) != 3) {
    stop("wrong number of colors in d_col!\nSupply colors for c(min, middle, max)!")
  }
  if (!(d_range %in% c("fit", "full"))) {
    stop("d_range must be either \"fit\" or \"full\"!")
  }
  if (length(trusted) == 0) {
    stop("\"trusted\" must contain at least one trusted status label")
  }
  fromIntermed <- FALSE
  allLables <- c("OK_sd", "OK_nc", "OK_d", "AD", "NS")
  notTrusted <- allLables[!(allLables %in% trusted)]
  if (is(metaDeconfOutput, "list")) {
    effectSize <- reshape2::melt(data = metaDeconfOutput$Ds,
                                 varnames = c("feature", "metaVariable"), value.name = "Ds")
    effectSize$Ds[effectSize$Ds == Inf] <- 0
    effectSize$Ds[is.na(effectSize$Ds)] <- 0
    fdr <- reshape2::melt(data = metaDeconfOutput$Qs, varnames = c("feature",
                                                                   "metaVariable"), value.name = "Qs")
    fdr$Qs[is.na(fdr$Qs)] <- 1
    status <- reshape2::melt(data = metaDeconfOutput$status,
                             varnames = c("feature", "metaVariable"), value.name = "status")
  }
  else if ((sum(c("stars", "insignificant", "featureNames") %in%
                colnames(metaDeconfOutput))) == 3) {
    warning("treating input as 'intermedData = T' Buildheatmap output!!")
    fromIntermed <- TRUE
    effectSize <- metaDeconfOutput
  }
  else {
    if (!"groupingVar" %in% colnames(metaDeconfOutput)) {
      metaDeconfOutput$groupingVar <- "metadata"
    }
    effectSize <- metaDeconfOutput[, c("feature", "metaVariable",
                                       "Ds", "groupingVar")]
    if (plotPartial == "partial") {
      effectSize$Ds <- metaDeconfOutput$partial
      keepMeta <- unique(c(keepMeta, "maxRsq"))
    }
    else if (plotPartial == "partialRel") {
      effectSize$Ds <- metaDeconfOutput$partialRel
    }
    else if (plotPartial == "partialNorm") {
      effectSize$Ds <- metaDeconfOutput$partialNorm
    }
    fdr <- metaDeconfOutput[, c("feature", "metaVariable",
                                "Qs")]
    status <- metaDeconfOutput[, c("feature", "metaVariable",
                                   "status")]
    effectSize$Ds[effectSize$Ds == Inf] <- 0
    effectSize$Ds[is.na(effectSize$Ds)] <- 0
    fdr$Qs[is.na(fdr$Qs)] <- 1
  }
  if (is.null(effectSize$groupingVar)) {
    effectSize$groupingVar <- "metadata"
  }
  if (!fromIntermed) {
    insignificant <- unlist(lapply(strsplit(as.character(status$status),
                                            split = ", "), function(l) any(l %in% notTrusted)))
    trueDeconf <- unlist(lapply(strsplit(as.character(status$status),
                                         split = ", "), function(l) any(l %in% trusted)))
    effectSize$stars <- cut(fdr$Qs, breaks = c(-Inf, 0.001,
                                               0.01, 0.1, Inf), label = c("***", "**", "*", ""))
    effectSize$stars[insignificant] <- ""
    effectSize$status <- trueDeconf
    effectSize$stars <- as.character(effectSize$stars)
    effectSize$insignificant <- insignificant
    effectSize$trueDeconf <- !trueDeconf
    if (coloring == 1) {
      effectSize$Ds[effectSize$insignificant] <- 1e-06
    }
    if (coloring == 2) {
      effectSize$Ds[effectSize$trueDeconf] <- 1e-06
    }
    remove_metavariables <- vector()
    for (i in unique(effectSize$metaVariable)) {
      aMetaVariable <- fdr[fdr$metaVariable == i, ]
      aMetaVariableD <- effectSize[effectSize$metaVariable ==
                                     i, ]
      if (sum(na.exclude(abs(aMetaVariable$Qs)) > q_cutoff) ==
          length(na.exclude(aMetaVariable$Qs)) || sum(na.exclude(abs(aMetaVariableD$Ds)) <
                                                      d_cutoff) == length(na.exclude(aMetaVariableD$Ds)) ||
          all(aMetaVariableD$stars == "")) {
        remove_metavariables <- c(remove_metavariables,
                                  i)
      }
    }
    if (!is.null(keepMeta)) {
      remove_metavariables <- remove_metavariables[!(remove_metavariables %in%
                                                       keepMeta)]
    }
    effectSize <- effectSize[!(effectSize$metaVariable %in%
                                 remove_metavariables), ]
    remove <- vector()
    for (i in unique(effectSize$feature)) {
      aGenus <- fdr[fdr$feature == i, ]
      aGenusD <- effectSize[effectSize$feature == i, ]
      if (sum(na.exclude(abs(aGenus$Qs)) > q_cutoff) ==
          length(na.exclude(aGenus$Qs)) | sum(na.exclude(abs(aGenusD$Ds)) <
                                              d_cutoff) == length(na.exclude(aGenusD$Ds)) |
          all(aGenusD$stars == "")) {
        remove <- c(remove, i)
      }
    }
    remove <- as.vector(unique(remove))
    if (length(remove) > 0) {
      if (!is.null(keepFeature)) {
        remove <- remove[!(remove %in% keepFeature)]
      }
      if (length(remove) > 0) {
        effectSize <- effectSize[!(effectSize$feature %in%
                                     remove), ]
      }
    }
    effectSize <- droplevels(effectSize)
    if (length(unique(effectSize$metaVariable)) == 0) {
      stop("No associations pass current q_cutoff and/or d_cutoff filters!")
    }
    eff_cast <- reshape2::dcast(effectSize, effectSize[[1]] ~
                                  metaVariable, value.var = "Ds")
    rownames(eff_cast) <- eff_cast[[1]]
    eff_cast[[1]] <- NULL
    if ((reOrder %in% c("both", "feat")) & (nrow(eff_cast) >
                                            1)) {
      ord <- hclust(dist(eff_cast, method = "euclidean"),
                    method = "ward.D")$order
      effectSize$feature <- factor(as.factor(effectSize$feature),
                                   levels = levels(as.factor(effectSize$feature))[ord])
    }
    if ((reOrder %in% c("both", "meta")) & (ncol(eff_cast) >
                                            1)) {
      eff_cast <- scale(t(eff_cast))
      ord2 <- hclust(dist(eff_cast, method = "euclidean"),
                     method = "ward.D")$order
      effectSize$metaVariable <- droplevels(effectSize$metaVariable)
      effectSize$metaVariable <- factor(as.factor(effectSize$metaVariable),
                                        levels = levels(as.factor(effectSize$metaVariable))[ord2])
    }
    if (plotPartial == "partial") {
      effectSize$metaVariable <- factor(effectSize$metaVariable,
                                        levels = c(setdiff(levels(effectSize$metaVariable),
                                                           "maxRsq"), "maxRsq"))
    }
    effectSize$featureNames <- effectSize$feature
    effectSize$metaVariableNames <- effectSize$metaVariable
    if (!is.null(featureNames)) {
      if (!is(featureNames, "data.frame")) {
        warning("class(featureNames) was coerced to \"data.frame\"")
        featureNames <- as.data.frame(featureNames)
      }
      if (length(unique(featureNames[[2]])) != length(featureNames[[2]])) {
        featureNames[[2]] <- make.unique(featureNames[[2]])
        warning("non-unique human-readable feature names where made unique using base::make.unique")
      }
      map = stats::setNames(featureNames[[2]], featureNames[[1]])
      effectSize$featureNames <- map[as.vector(effectSize$feature)]
      effectSize$featureNames <- factor(as.factor(effectSize$featureNames),
                                        levels = map[levels(effectSize$feature)])
    }
    if (!is.null(metaVariableNames)) {
      if (!is(metaVariableNames, "data.frame")) {
        warning("class(metaVariableNames) was coerced to \"data.frame\"")
        metaVariableNames <- as.data.frame(metaVariableNames)
      }
      if (length(unique(metaVariableNames[[2]])) != length(metaVariableNames[[2]])) {
        metaVariableNames[[2]] <- make.unique(metaVariableNames[[2]])
        warning("non-unique human-readable metaVariable names where made unique using base::make.unique")
      }
      map = stats::setNames(metaVariableNames[[2]], metaVariableNames[[1]])
      effectSize$metaVariableNames <- map[as.vector(effectSize$metaVariable)]
      effectSize$metaVariableNames <- factor(as.factor(effectSize$metaVariableNames),
                                             levels = map[levels(effectSize$metaVariable)])
    }
  }
  if (intermedData == TRUE) {
    if (!any(c("*", "**", "***") %in% unique(effectSize$stars))) {
      warning("No unconfounded associations remain with the current cutoff values. ")
    }
    return(effectSize)
  }
  if (!any(c("*", "**", "***") %in% unique(effectSize$stars))) {
    stop("No unconfounded associations remain with the current cutoff values. Consider manually including categorical metaVariables into the Heatmap by listing them through the 'keepMeta' argument.")
  }
  lowerLim <- min(effectSize$Ds)
  upperLim <- max(effectSize$Ds)
  if (d_range == "full") {
    lowerLim <- -1
    upperLim <- 1
  }
  signifCol <- c("gray45", "black")
  signifMeaning <- c("confounded", "deconfounded")
  legendShapes <- c(1, 8)
  if (all(effectSize$status)) {
    signifCol <- c("black")
    signifMeaning <- c("deconfounded")
    legendShapes <- c(8)
  }
  if (cuneiform) {
    divShapes <- c()
    divShapesMeaning <- c()
    signs <- unique(sign(effectSize$Ds))
    if (-1 %in% signs) {
      divShapes <- c(divShapes, 25)
      divShapesMeaning <- c(divShapesMeaning, "negative association")
    }
    if (0 %in% signs) {
      divShapes <- c(divShapes, 23)
      divShapesMeaning <- c(divShapesMeaning, "no association/no data")
    }
    if (1 %in% signs) {
      divShapes <- c(divShapes, 24)
      divShapesMeaning <- c(divShapesMeaning, "positive association")
    }
    heatmapGGplot <- ggplot(effectSize, aes(x = .data$metaVariableNames,
                                            y = .data$featureNames)) + geom_point(aes(fill = .data$Ds,
                                                                                      shape = as.factor(sign(.data$Ds)), color = .data$status)) +
      scale_shape_manual(name = "Direction", values = divShapes,
                         labels = divShapesMeaning) + scale_fill_gradient2(low = d_col[1],
                                                                           mid = d_col[2], high = d_col[3], midpoint = 0, guide = guide_colorbar(display = "gradient"),
                                                                           limits = c(lowerLim, upperLim)) + scale_color_manual(name = "Confounding status",
                                                                                                                                values = signifCol, labels = signifMeaning) + guides(color = guide_legend(override.aes = list(shape = 24))) +
      facet_grid(cols = vars(.data$groupingVar), space = "free_x",
                 drop = T, scales = "free_x") + theme_classic() +
      theme(axis.text.x = element_text(size = 7, angle = 90,
                                       hjust = 1, vjust = 0.3), axis.text.y = element_text(size = 7,
                                                                                           angle = 0, hjust = 1, vjust = 0.35), plot.title.position = "plot",
            plot.title = element_text(hjust = 0), strip.background = element_blank()) +
      labs(title = "Summarizing cuneiform plot", x = "Metadata variables",
           y = "Omics features")
  }
  else {
    print("create heatmap")
    effectSize <- effectSize %>%
      arrange(status) %>%
      mutate(status = ifelse(status, yes = "Deconfounded\nassociation", no = "Confounded\nassociation"),
             status = ifelse(stars != "", yes = status, no = "Not significant"),
             status = factor(status, levels = c("Deconfounded\nassociation", "Confounded\nassociation", "Not significant")))
    print(unique(effectSize$status))
    heatmapGGplot <- ggplot(effectSize, aes(x = .data$metaVariableNames,
                                            y = .data$featureNames)) +
      geom_tile(aes(fill = .data$Ds, size = .data$status,
                    color = .data$status)) +
      scale_fill_gradient2(name = "Effect size",
                           low = d_col[1], mid = d_col[2], high = d_col[3],
                           midpoint = 0, guide = guide_colorbar(display = "gradient"),
                           limits = c(lowerLim, upperLim)) +
      geom_text(aes(label = .data$stars), size = starSize, key_glyph = "point",
                nudge_y = starNudge_y) +
      scale_size_manual(name = NULL,
                        values = c("Deconfounded\nassociation" = 2, "Confounded\nassociation" = 0.5,
                                   "Not significant" = 0.5),
                        breaks = "Deconfounded\nassociation") +
      scale_color_manual(name = NULL,
                         values = c("Deconfounded\nassociation" = "black", "Confounded\nassociation" = "grey",
                                    "Not significant" = "grey"),
                         breaks = "Deconfounded\nassociation") +
      guides(size = guide_legend(override.aes = list(fill = "white", color = c("black"))),
             color = "none") +
      facet_grid(~groupingVar,  scales = "free_x", drop = T, space = "free_x") +
      theme_classic() +
      theme(axis.text.x = element_text(size = 7, angle = 90,
                                       hjust = 1, vjust = 0.3),
            axis.text.y = element_text(size = 7,
                                       angle = 0, hjust = 1, vjust = 0.35), plot.title.position = "plot",
            plot.title = element_text(hjust = 0), plot.subtitle = element_text(size = 8),
            strip.background = element_blank()) +
      labs(title = "Summarizing heatmap",
           subtitle = "p.adjust-values: < 0.001 = ***, < 0.01 = **, < 0.1 = * ",
           x = "Metadata variables", y = "Omics features")
  }
  return(heatmapGGplot)
}
