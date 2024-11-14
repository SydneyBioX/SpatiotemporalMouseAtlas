shinyServer(function(input, output, session) {
  
  shinylogs::track_usage(storage_mode = shinylogs::store_json(path = "logs/"))
  
  values <- reactiveValues(meta = meta)
  
  gene_name_impGenerator = function() {
    if (input$gene_name_imp == "") {
      values$gene_name_imp_name <<- NA
      out = "Please type in a gene name"
      return(out)
    }
    
    ind = match(tolower(input$gene_name_imp), tolower(genes_imp))[1]
    if (is.na(ind)) {
      
      # genes_imp
      ind = grepl(tolower(input$gene_name_imp), tolower(genes_imp))
      
      # then plot for the first one
      values$gene_name_imp_name <<- genes_imp[ind][1]
      
      # display genes matching
      out = paste(c("Selecting first of multiple matching genes:",
                    sort(genes_imp[ind])[1:pmin(5,sum(ind))],
                    ifelse(sum(ind) > 5, "and more.", "")), collapse = " ")
    } else {
      values$gene_name_imp_name <<- genes_imp[ind][1]
      out = paste0("Selecting gene: ", values$gene_name_imp_name)
      return(out)
    }
    
    if (is.na(values$gene_name_imp_name)) {
      out = "no genes match input"
    }
    
    return(out)
  }
  
  
  output$gene_name_imp_parse_status <- renderText({
    out = gene_name_impGenerator()
    return(out)
  })
  
  
  
  # Normalising scales for multiple embryos
  norm = function(col){
    (col - min(col, na.rm = TRUE)) / (max(col, na.rm = TRUE) - min(col, na.rm = TRUE))
  }
  
  normaliseScale = function(data) {
    data = data |> 
      group_by(embryo) |> 
      mutate(dim1 = norm(dim1), 
             dim2 = 2*norm(dim2)) |> 
      ungroup()
    
    return(data)
  }
  
  # Syncs the sliders for point size for spatial plot and spatial plot (imputed)
  observeEvent(
    input$spatialPlotPointSize,
    updateSelectizeInput(session, "spatialPlotPointSizeImp", selected = input$spatialPlotPointSize),
  )
  
  observeEvent(
    input$spatialPlotPointSizeImp,
    updateSelectizeInput(session, "spatialPlotPointSize", selected = input$spatialPlotPointSizeImp),
    
  )
  
  

  # Generates spatial plot
  spatialPlotGenerator = function() {
    
    meta_sub = subset(
      meta,
      embryo %in% input$embryo_subset #& 
      # z %in% input$z_subset
    )
    
    
    if(input$embryo_same_size) {
      meta_sub = normaliseScale(meta_sub)
    }
    
    selectedNames = meta_sub$uniqueID

    
    # for the imputed plot, filter all NA cells
    # Set expression to 0 for NAs cell types when subsetting cell types
    # Plot by expression for log counts plots, dont grey out NAs
    
    
    if (input$tabs == "Spatial plot (Imputed)") {
      
      if(input$celltype_subset_all_SPImputed) {
        celltype_subset_imp = celltypes
      } else {
        celltype_subset_imp = input$celltype_subset_imp
      }
      
      add_imp()
      
      gene_name_impGenerator()
      
      validate(
        need(isolate(values$gene_name_imp_name) %in% genes_imp, "No valid gene selected.")
      )
      
      expr_imp = imp[isolate(values$gene_name_imp_name),]
      names(expr_imp) <- colnames(imp)
      expr <- expr_imp[rownames(meta)]
      
      pc_vals = expr[selectedNames]
      
      
      validate(
        need(length(pc_vals) > 0, "No cells are selected for subsetting, please tick at least one option for each of the embryo, z-slice and mapped cell type categories.")
      )
      
      # Remove NA cells, as they can't be imputed
      meta_sub = meta_sub |> 
        filter(!is.na(cellType))
      
      pc_vals = pc_vals[meta_sub$uniqueID]
      
      # Make the expression of NA cells 0 so they dont show up in front of selected cells.
      makeZero = meta_sub |> 
        filter(!cellType %in% celltype_subset_imp) |> 
        pull(uniqueID)
      
      pc_vals[makeZero] = 0
      
      pc_cols = colorRampPalette(colours)(100)[as.numeric(cut(pc_vals,100))]
      
      # If all expression is 0 then grey out all cells
      if(sum(pc_vals) == 0) {
        pc_cols = rep("grey95", length(pc_vals))
      }
      
      names(pc_cols) <- names(pc_vals)
      
      yl = " (Imputed expression)"
      
    } else {
      
      if(input$celltype_subset_all_SP) {
        celltype_subset = celltypes
      } else {
        celltype_subset = input$celltype_subset
      }
      
      if (input$colour_by == "Expression logcounts") {
        
        add_exprs_norm()
        
        pc_vals = exprs_norm[input$gene_name, selectedNames]
        
        validate(
          need(length(pc_vals) > 0, "No cells are selected for subsetting, please tick at least one option for each of the embryo, z-slice and mapped cell type categories.")
        )
        
        # Set the expression of subsetted cells to 0, so they don't show up in front, this wont affect violin plot as it uses separate data.frame
        pc_vals[rownames(meta_sub)[!meta_sub$cellType %in% celltype_subset]] = 0
        
        pc_cols = colorRampPalette(colours)(100)[as.numeric(cut(pc_vals,100))]
        
        # If all expression is 0 then grey out all cells
        if(sum(pc_vals) == 0) {
          pc_cols = rep("grey95", length(pc_vals))
        }
        
        names(pc_cols) <- names(pc_vals)
        
        yl = ""
        
      }
      if (input$colour_by == "Mapped cell type") {

        
        pc_vals = as.character(meta_sub$cellType)
        
        validate(
          need(length(pc_vals) > 0, "No cells are selected for subsetting, please tick at least one option for each of the embryo, z-slice and mapped cell type categories.")
        )
        
        
        pc_cols = celltype_colours[pc_vals]
        pc_cols[!pc_vals %in% celltype_subset] = "grey95"
        
        # Add ordering of cell types, so selected cell types show up first, make sure NA is the first factor so that its always plotted behind
        pc_vals = factor(pc_vals, levels = c(
          NA,
          unique(pc_vals[!pc_vals %in% c(celltype_subset, NA)]),
          unique(pc_vals[pc_vals %in% celltype_subset])
        ), exclude = NULL)
        
        names(pc_cols) <- rownames(meta_sub)
       
      }
      
    }

    # Order meta sub by pc_vals so cells of interest show up on top.
    meta_sub = meta_sub |>
      mutate(pc_vals = pc_vals) |>
      arrange(pc_vals)

    
    g_base = meta_sub |> 
      ggplot(
        aes(
          x = dim1,
          y = dim2,
          group = uniqueID,
          fill = uniqueID,
          tooltip = cellType,
          colour = uniqueID
        )
      ) + 
      scale_fill_manual(values = pc_cols, na.value = NA) +
      scale_colour_manual(values = pc_cols, na.value = NA) 
    
    g = g_base +
      geom_point_interactive(
        size = input$spatialPlotPointSize,
        shape = 21
      )
   
    g <- g +
      theme_classic() +
      coord_fixed() +
      xlab("") +
      # legend needs to be "none" because of colour is defined by uniqueID in aes not expression/celltype
      theme(legend.position = "none") + 
      theme(axis.line = element_blank()) + 
      theme(axis.text = element_blank()) + 
      theme(axis.title.y = element_blank()) +
      theme(axis.title.x = element_text(face = "italic")) +
      theme(axis.ticks = element_blank()) +
      theme(strip.text.y = element_blank()) +
      theme(strip.background.x = element_rect(colour = "white")) +
      NULL
    
    if (length(input$embryo_subset) > 1) {

      g <- g +
        facet_grid(~embryo, labeller = labeller(embryo = embryolabeller)) +
        NULL
      
    } else {
      g <- g +
        facet_grid(~embryo, labeller = labeller(embryo = embryolabeller)) +
        NULL
    }
    
    g_leg = NULL
    
    if (input$colour_by %in% c("Area-standardised gene expression", "Expression logcounts") | input$tabs == "Spatial plot (Imputed)") {
      
      heights = c(10,1)
      
      if (input$tabs == "Spatial plot (Imputed)") {
        gene_name_impGenerator()
        validate(
          need(isolate(values$gene_name_imp_name) %in% genes_imp, "No valid gene selected.")
        )
        g <- g + 
          xlab(paste0(isolate(values$gene_name_imp_name), yl)) + 
          NULL
        ttl = "Imputed gene expression"
      } else {
        g <- g + 
          xlab(paste0(input$gene_name, yl)) + 
          NULL
        ttl = input$colour_by
      }
      
      g_leg = g_leg_list[[ttl]]
    }
    
    if (input$colour_by == "Mapped cell type" & input$tabs != "Spatial plot (Imputed)") {
      
      heights = c(10,2)
      
      # get the legend
      g_leg <- g_leg_list[["Mapped cell type"]]
    }
    
    
    if (input$colour_by == "Mapped gut tube subtype" & input$tabs != "Spatial plot (Imputed)") {
      
      heights = c(10,1)
      
      # get the legend
      pc_df = data.frame(Expression = factor(unique(pc_vals), levels = names(AP_pseudo_colours)))
      
      g_leg_raw = ggplot(pc_df, aes(x = Expression, y = Expression)) +
        geom_point(aes(colour = Expression), size = 10) +
        theme_transparent() +
        scale_colour_manual(values = AP_pseudo_colours) +
        theme(legend.position = "bottom") +
        theme(legend.text = element_text(size = 15),
              legend.title = element_text(size = 20)) +
        guides(color = guide_legend(title.position = "top",
                                    title.hjust = 0.5,
                                    ticks = FALSE,
                                    title = "Mapped gut tube cell type")) +
        NULL
      
      g_leg = as_ggplot(get_legend(g_leg_raw))
    }
    
    
    list(g = g, g_leg = g_leg)
  }
  
  output$spatialPlot <- renderGirafe({
    
    
    g = spatialPlotGenerator()[["g"]]
    
    gi = girafe(code = print(g),
                options = list(
                  opts_hover(css = "fill:#FF3333;stroke:black;cursor:pointer;", reactive = TRUE),
                  opts_selection(type = "multiple", css = "fill:#FF3333;stroke:black;")
                ))
    gi
    
    return(gi)
  })
  
  
  output$spatialPlot_imp <- renderGirafe({
    
    g = spatialPlotGenerator()[["g"]]
    
    gi = girafe(code = print(g),
                options = list(
                  opts_hover(css = "fill:#FF3333;stroke:black;cursor:pointer;", reactive = TRUE),opts_selection(type = "multiple", css = "fill:#FF3333;stroke:black;")
                ))
    gi
    
    return(gi)
  })
  
  output$spatialPlot_leg <- renderPlot({
    return(spatialPlotGenerator()[["g_leg"]])
  })
  
  output$spatialPlot_leg_imp <- renderPlot({
    return(spatialPlotGenerator()[["g_leg"]])
  })
  
  
  
  
  violinPlotGenerator = function() {
    
    meta_sub = subset(
      meta,
      embryo %in% input$embryo_subset #& 
      # z %in% input$z_subset
    )
    selectedNames = rownames(meta_sub)
    
    if (input$tabs == "Spatial plot (Imputed)") {
      
      add_imp()
      
      # example gene
      # gene = "Pcdh19"
      # expr_imp = imp[input$gene_name_imp,]
      gene_name_impGenerator()
      validate(
        need(isolate(values$gene_name_imp_name) %in% genes_imp, "No valid gene selected.")
      )
      expr_imp = imp[isolate(values$gene_name_imp_name),]
      names(expr_imp) <- colnames(imp)
      expr <- expr_imp[rownames(meta)]
      
      pc_vals = expr[rownames(meta) %in% selectedNames]
      
      
      validate(
        need(length(pc_vals) > 0, "No cells are selected for subsetting, please tick at least one option for each of the embryo, z-slice and mapped cell type categories.")
      )
      
      yl = paste0(isolate(values$gene_name_imp_name), " (Imputed expression)")
      
    } else {
      
      add_exprs_norm()
      
      pc_vals = exprs_norm[input$gene_name, ][colnames(exprs_norm) %in% selectedNames]
      validate(
        need(length(pc_vals) > 0, "No cells are selected for subsetting, please tick at least one option for each of the embryo, z-slice and mapped cell type categories.")
      )
      
      yl = input$gene_name
      
      
    }
    
    ct_df = data.frame(pc_val = pc_vals, ct = meta[names(pc_vals), "extended_atlas_celltype"])
    ct_count = unclass(table(ct_df$ct))
    ct_df$ct_count = ct_count[as.character(ct_df$ct)]
    
    
    g_violin = ggplot(ct_df, aes(x = ct,
                                 y = pc_val,
                                 group = ct,
                                 fill = ct
    )) + 
      geom_text(aes(label = ct_count), angle = 60, size = 2, data = data.frame(ct_count = ct_count, ct = names(ct_count), pc_val = 0.5 + max(na.omit(ct_df$pc_val)))) +
      geom_violin(draw_quantiles = 0.5, colour = "black", scale = "width") + 
      theme_classic() + 
      scale_fill_manual(values = celltype_colours) + 
      xlab("") +
      ylab(yl) + 
      theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6)) + 
      theme(axis.title.y = element_text(face = "italic")) + 
      theme(legend.position = "none") +
      NULL
    
    g_violin
  }
  
  output$violinPlot <- renderGirafe({
    
    
    g = violinPlotGenerator()
    
    gi = girafe(code = print(g),
                options = list(
                  opts_hover(css = "fill:#FF3333;stroke:black;cursor:pointer;", reactive = TRUE),opts_selection(type = "multiple", css = "fill:#FF3333;stroke:black;")
                ))
    gi
    
    return(gi)
  })
  
  output$violinPlot_imp <- renderGirafe({
    
    
    g = violinPlotGenerator()
    
    gi = girafe(code = print(g),
                options = list(
                  opts_hover(css = "fill:#FF3333;stroke:black;cursor:pointer;", reactive = TRUE),opts_selection(type = "multiple", css = "fill:#FF3333;stroke:black;")
                ))
    gi
    
    return(gi)
  })
  
  
  
  # for Virtual Dissection section
  VirtualDissectionPlotGenerator = function() {
    
    meta2 <- isolate(values$meta)
    
    meta_sub = subset(meta2,
                      meta2$embryo %in% input$embryo_subset #& 
                      #meta2$z %in% input$z_subset
    )
  
    
    if(input$embryo_same_size) {
      meta_sub = normaliseScale(meta_sub)
    }
    
    
    validate(
      need(nrow(meta_sub) > 0,
           "No cells are selected for subsetting, please tick at least one option for each of the embryo and z-slice categories.")
    )
    
    
    if (input$virtual_colour_by == "Expression logcounts") {
      
      add_exprs_norm()
      
      pc_vals = exprs_norm[input$virtual_gene_name, meta_sub$uniqueID]
      
      validate(
        need(length(pc_vals) > 0, "No cells are selected for subsetting, please tick at least one option for each of the embryo, z-slice and mapped cell type categories.")
      )
      
      pc_cols = colorRampPalette(colours)(100)[as.numeric(cut(pc_vals,100))]
      names(pc_cols) <- names(pc_vals)
      
      yl = ""
      
    }
    if (input$virtual_colour_by == "Imputed expression") {
      
      
      add_imp()
      
      expr_imp = imp[input$virtual_gene_name_imputed,]
      names(expr_imp) <- colnames(imp)
      expr <- expr_imp[rownames(meta)]
      
      pc_vals = expr[meta_sub$uniqueID]
      
      
      validate(
        need(length(pc_vals) > 0, "No cells are selected for subsetting, please tick at least one option for each of the embryo, z-slice and mapped cell type categories.")
      )
      
      # Remove NA cells, as they can't be imputed
      meta_sub = meta_sub |>
        filter(!is.na(cellType))
      
      pc_vals = pc_vals[meta_sub$uniqueID]
      
      # Make the expression of NA cells 0 so they dont show up in front of selected cells.
      makeZero = meta_sub |>
        filter(is.na(cellType)) |>
        pull(uniqueID)

      pc_vals[makeZero] = 0
      
      pc_cols = colorRampPalette(colours)(100)[as.numeric(cut(pc_vals,100))]
      
      # If all expression is 0 then grey out all cells
      if(sum(pc_vals) == 0) {
        pc_cols = rep("grey95", length(pc_vals))
      }
      
      names(pc_cols) <- names(pc_vals)
      
      yl = " (Imputed expression)"
      
    }
    if (input$virtual_colour_by == "Mapped cell type") {
      
      pc_vals = as.character(meta_sub$cellType)
      
      validate(
        need(length(pc_vals) > 0, "No cells are selected for subsetting, please tick at least one option for each of the embryo, z-slice and mapped cell type categories.")
      )
      
      
      pc_cols = celltype_colours[pc_vals]
      
      # Add ordering of cell types, so selected cell types show up first, make sure NA is the first factor so that its always plotted behind
      pc_vals = factor(pc_vals, levels = c(
        NA,
        unique(pc_vals[!is.na(pc_vals)])
      ), exclude = NULL)
      
      names(pc_cols) <- rownames(meta_sub)
      
    } 
    
    if (input$virtual_colour_by == "AP axis") {
      pc_vals = meta_sub$AP
      pc_cols = colorRampPalette(APDV_colours)(100)[as.numeric(cut(pc_vals,100))]
      names(pc_cols) <- rownames(meta_sub)
    } 
    
    if (input$virtual_colour_by == "DV axis") {
      
      pc_vals = meta_sub$DV
      pc_cols = colorRampPalette(APDV_colours)(100)[as.numeric(cut(pc_vals,100))]
      names(pc_cols) <- rownames(meta_sub)
    }
    
    # Order meta sub by pc_vals so cells of interest show up on top.
    meta_sub = meta_sub |>
      mutate(pc_vals = pc_vals) |> 
      arrange(selected, !is.na(pc_vals), pc_vals)
    
    # Virtual dissection using Physical slice
    if (input$virtualdissect_choice == "Physical") {
      g = ggplot(meta_sub, 
                 aes(x = dim1,
                     y = dim2,
                     # tooltip = uniqueID,
                     tooltip = cellType,
                     data_id = uniqueID)) + 
        facet_grid(~embryo, labeller = labeller(embryo = embryolabeller))
      
    }
    
    # Virtual dissection using UMAP
    if (input$virtualdissect_choice == "UMAP") {
      
      g = ggplot(meta_sub, 
                 aes(x = UMAP1,
                     y = UMAP2,
                     # tooltip = uniqueID,
                     tooltip = cellType,
                     data_id = uniqueID))
      
      g <- g + geom_point(data = meta, colour = "grey95", size = 0.8, alpha = 1, show.legend = FALSE)
    }
    
    
    g = g  +
      geom_point_interactive( aes(colour = selected, fill = uniqueID),
                              shape = 21,
                              size = input$dissectionPointSize,
                              stroke = input$dissectionStrokeSize) +
      scale_colour_manual(values = c("Group A" = "red",
                                     "Unselected" = alpha("grey", 0.0),
                                     "Group B" = "blue"),
                          breaks = c("Group A", "Group B")) +
      scale_fill_manual(values = pc_cols) +
      theme_classic() + 
      coord_fixed() + 
      theme(legend.position = "bottom") + 
      theme(axis.line = element_blank()) + 
      theme(axis.text = element_blank()) + 
      theme(axis.title = element_blank()) +
      theme(axis.ticks = element_blank()) +
      theme(strip.background = element_rect(colour = "white")) +
      guides(colour = guide_legend(title = "",
                                   override.aes = list(size = 5)),
             fill = "none") +
      NULL
    
    gi = girafe(code = print(g),
                options = list(
                  opts_hover(css = "fill:#FF3333;stroke:black;cursor:pointer;", reactive = TRUE),
                  opts_selection(type = "multiple", css = "fill:#FF3333;stroke:black;"),
                  opts_toolbar(position = "top", fixed = TRUE)
                ))
    
    gi
  }
  
  updateVirtualDissection = reactive({
    paste(input$embryo_same_size,
          input$embryo_subset, 
          input$virtualdissect_choice,
          input$add,
          input$add2,
          input$remove,
          input$removeAll,
          input$dissectionPointSize,
          input$dissectionStrokeSize,
          input$virtualDissectionCelltypes,
          input$addCT,
          input$addCT2,
          input$removeCT,
          input$removeAllButCT,
          input$preselected,
          input$add_pre,
          input$add_pre2,
          input$remove_pre,
          input$virtual_colour_by,
          input$virtual_gene_name,
          input$virtual_gene_name_imputed)
  })
  
  VirtualDissectionPlot = eventReactive(updateVirtualDissection(),
                                        VirtualDissectionPlotGenerator()
  )
  
  
  output$virtualDissection <- renderGirafe({
    
    VirtualDissectionPlot = VirtualDissectionPlot()
    VirtualDissectionPlot
  })
  
  
  output$virtualDissection_leg <- renderPlot({
    g_leg_list[[input$virtual_colour_by]]
  })
  
  
  observeEvent(input$virtualDissection_reset, {
    session$sendCustomMessage(type = 'virtualDissection_set', message = character(0))
  })
  
  virtualDissection_selected_state <- reactive({
    input$virtualDissection_selected
  })
  
  addGenerator = function() {
    meta2 <- isolate(values$meta)
    meta2[meta2$uniqueID %in% virtualDissection_selected_state(),"selected"] <- "Group A"
    values$meta <<- meta2
    showNotification("Virtually dissected cells added to Group A!")
    return(meta)
  }
  
  
  observeEvent(input$add, {
    addGenerator()
  })
  
  add2Generator = function() {
    meta2 <- isolate(values$meta)
    meta2[meta2$uniqueID %in% virtualDissection_selected_state(),"selected"] <- "Group B"
    values$meta <<- meta2
    showNotification("Virtually dissected cells added to Group B!")
    return(meta)
  }
  
  observeEvent(input$add2, {
    add2Generator()
  })
  
  removeGenerator = function() {
    meta2 <- isolate(values$meta)
    meta2[meta2$uniqueID %in% virtualDissection_selected_state(),"selected"] <- "Unselected"
    values$meta <<- meta2
    showNotification("Virtually dissection cells removed from either group!")
    return(meta)
  }
  
  observeEvent(input$remove, {
    removeGenerator()
  })
  
  removeAllGenerator = function() {
    meta2 <- isolate(values$meta)
    meta2[,"selected"] <- "Unselected"
    values$meta <<- meta2
    
    showNotification("All cells removed from selection!")
    return(meta)
  }
  
  observeEvent(input$removeAll, {
    removeAllGenerator()
  })
  
  ###### Adding cell types
  
  addCTGenerator = function() {
    
    meta2 <- isolate(values$meta)
    selectedUniqueIDs = subset(meta2,
                               cellType
                               %in% input$virtualDissectionCelltypes & 
                                 embryo %in% input$embryo_subset# & 
                               # z %in% input$z_subset
    )[,"uniqueID"]
    meta2[meta2$uniqueID %in% selectedUniqueIDs, "selected"] <- "Group A"
    values$meta <<- meta2
    showNotification("Selected cell type cells added to Group A!")
    return(meta)
  }
  
  observeEvent(input$addCT, {
    addCTGenerator()
  })
  
  addCT2Generator = function() {
    meta2 <- isolate(values$meta)
    selectedUniqueIDs = subset(meta2,
                               cellType
                               %in% input$virtualDissectionCelltypes & 
                                 embryo %in% input$embryo_subset #& 
                               # z %in% input$z_subset
    )[,"uniqueID"]
    meta2[meta2$uniqueID %in% selectedUniqueIDs, "selected"] <- "Group B"
    values$meta <<- meta2
    showNotification("Selected cell type cells added to Group B!")
    return(meta)
  }
  
  observeEvent(input$addCT2, {
    addCT2Generator()
  })
  
  removeCTGenerator = function() {
    meta2 <- isolate(values$meta)
    selectedUniqueIDs = subset(meta2,
                               cellType %in% input$virtualDissectionCelltypes & 
                                 embryo %in% input$embryo_subset
    )[,"uniqueID"]
    meta2[meta2$uniqueID %in% selectedUniqueIDs, "selected"] <- "Unselected"
    values$meta <<- meta2
    showNotification("Selected cell type cells removed from either group!")
    return(meta)
  }
  
  observeEvent(input$removeCT, {
    removeCTGenerator()
  })
  
  
  removeAllButCTGenerator = function() {
    meta2 <- isolate(values$meta)
    selectedUniqueIDs = subset(meta2,
                                 embryo %in% input$embryo_subset &
                                 !(cellType %in% input$virtualDissectionCelltypes)
    )[,"uniqueID"]
    meta2[meta2$uniqueID %in% selectedUniqueIDs, "selected"] <- "Unselected"
    values$meta <<- meta2
    showNotification("All but selected cell type cells removed from either group!")
    return(meta)
  }
  
  observeEvent(input$removeAllButCT, {
    removeAllButCTGenerator()
  })
  
  ######### CT end
  
  addPreselectedGenerator = function() {
    selectedUniqueIDs = do.call(c,
                                sapply(input$preselected$datapath, 
                                       readRDS, simplify = FALSE)
    )
    meta2 <- isolate(values$meta)
    meta2[meta2$uniqueID %in% selectedUniqueIDs, "selected"] <- "Group A"
    values$meta <<- meta2
    showNotification("Preselected cells added to Group A!")
    return(meta)
  }
  
  observeEvent(input$add_pre, {
    addPreselectedGenerator()
  })
  
  addPreselected2Generator = function() {
    selectedUniqueIDs = do.call(c,
                                sapply(input$preselected$datapath, 
                                       readRDS, simplify = FALSE)
    )
    meta2 <- isolate(values$meta)
    meta2[meta2$uniqueID %in% selectedUniqueIDs, "selected"] <- "Group B"
    values$meta <<- meta2
    showNotification("Preselected cells added to Group B!")
    return(meta)
  }
  
  observeEvent(input$add_pre2, {
    addPreselected2Generator()
  })
  
  removePreselectedGenerator = function() {
    selectedUniqueIDs = do.call(c,
                                sapply(input$preselected$datapath, 
                                       readRDS, simplify = FALSE)
    )
    meta2 <- isolate(values$meta)
    meta2[meta2$uniqueID %in% selectedUniqueIDs, "selected"] <- "Unselected"
    values$meta <<- meta2
    showNotification("Preselected cells removed from either group!")
    return(meta)
  }
  
  observeEvent(input$remove_pre, {
    removePreselectedGenerator()
  })
  
  output$virtualDissectionCellDownload <- downloadHandler(
    filename = function() {
      paste("cells_groupA.Rds")
    },
    content = function(file) {
      meta2 <- isolate(values$meta)
      meta_sub = subset(meta2, selected == "Group A" &
                          embryo %in% input$embryo_subset #& 
                        # z %in% input$z_subset
      )
      saveRDS(sort(as.character(unique(meta_sub[,"uniqueID"]))), file)
    }
  )
  
  output$virtualDissectionCellDownload2 <- downloadHandler(
    filename = function() {
      paste("cells_groupB.Rds")
    },
    content = function(file) {
      meta2 <- isolate(values$meta)
      meta_sub = subset(meta2, selected == "Group B" &
                          embryo %in% input$embryo_subset #& 
                        # z %in% input$z_subset
      )
      saveRDS(sort(as.character(unique(meta_sub[,"uniqueID"]))), file)
    }
  )
  
  VirtualDissectionBarPlotGenerator = function() {
    
    meta2 <- isolate(values$meta)
    
    meta_subList = sapply(c("Group A", "Group B"), function(sel) {
      subset(meta2, selected == sel &
               embryo %in% input$embryo_subset #&
             # z %in% input$z_subset
      )
    }, simplify = FALSE)
    
    meta_subList <- meta_subList[unlist(lapply(meta_subList, nrow)) > 0]
    
    gList = lapply(meta_subList, function(meta_sub) {
      
      if (nrow(meta_sub) == 0) return(NULL)
      
      meta_sub_sum = data.frame(
        count = tapply(meta_sub$cellType, meta_sub$cellType, length))
      meta_sub_sum$cellType <- rownames(meta_sub_sum)
      meta_sub_sum = na.omit(meta_sub_sum)
      
      g = ggplot(meta_sub, 
                 aes(x=reorder(cellType,cellType,
                               function(x)-length(x)))) +    
        geom_bar(aes(fill = cellType)) + 
        geom_text(aes(y = count,
                      label = count), 
                  data = meta_sub_sum, 
                  vjust = "bottom",
                  size = 5) +
        scale_fill_manual(values = celltype_colours) + 
        theme_classic() + 
        theme(legend.position = "none") +
        theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size = 15)) +
        theme(axis.title.y = element_text(size = 15)) +
        theme(plot.title = element_text(size = 15)) +
        xlab("") +
        ylab("Number of cells selected") +
        ggtitle(paste0("Mapped cell types for ", 
                       sum(meta_sub_sum$count), 
                       " selected cells")) + 
        scale_x_discrete(expand = c(0.1, 0.15)) +
        scale_y_continuous(expand = c(0.1, 0.2)) +
        NULL
      return(g)
    })
    
    if (all(unlist(lapply(gList, is.null)))) return(NULL)
    
    gList_named <- sapply(names(gList), function(gname){
      gList[[gname]] + ylab(paste0("Number of ", gname, " cells"))
    }, simplify = FALSE)
    
    wrap_plots(gList_named, nrow = 1)
  }
  
  VirtualDissectionBarPlot = eventReactive(updateVirtualDissection(),
                                           VirtualDissectionBarPlotGenerator()
  )
  
  output$VirtualDissectionBarPlot <- renderPlot({
    VirtualDissectionBarPlot = VirtualDissectionBarPlot()
    
    validate(
      need(!is.null(VirtualDissectionBarPlot),
           "Need cells assigned to either Group A or Group B to be able to display barplot")
    )
    
    VirtualDissectionBarPlot
  })
  
  
  VirtualDissectionMAPlotGenerator = function() {
    
    
    meta2 <- isolate(values$meta)
    
    virtualDissectionCells = as.character(
      subset(meta2, 
             selected == "Group A" &
               embryo %in% input$embryo_subset #& 
             # z %in% input$z_subset
      )[,"uniqueID"])
    
    validate(
      need(length(virtualDissectionCells) > 0,
           "Need cells assigned to both Group A and Group B to be able to display MA-plot"
      )
    )
    
    
    
    virtualDissectionCells_not = as.character(
      subset(meta2, 
             selected == "Group B" &
               embryo %in% input$embryo_subset #& 
             # z %in% input$z_subset
      )[,"uniqueID"])
    
    validate(
      need(length(virtualDissectionCells_not) > 0,
           "Need cells assigned to both Group A and Group B to be able to display MA-plot")
    )
    
    
    # Log counts MA vals
    add_exprs()
    
    
    A_val = rowMeans(exprs_counts[, c(virtualDissectionCells, virtualDissectionCells_not)])
    M_val = rowMeans(exprs_counts[, virtualDissectionCells]) - rowMeans(exprs_counts[, virtualDissectionCells_not])
    MA_df_seq = data.frame(
      A_val = A_val,
      M_val = M_val,
      gene = rownames(exprs_counts),
      M_rank = rank(-M_val),
      M_rank_bottom = rank(M_val),
      type = "seqFISH"
    )
    
    
    # Imputed MA vals
    add_imp()
    
    # Need to filter NA cells for imputed data
    na_cells = meta2 |> 
      filter(is.na(cellType)) |> 
      pull(uniqueID)
    
    virtualDissectionCells_Imp = virtualDissectionCells[!(virtualDissectionCells%in% na_cells)]
    virtualDissectionCells_not_Imp = virtualDissectionCells_not[!(virtualDissectionCells_not%in% na_cells)]
    
    A_val_Imp = rowMeans(imp[, c(virtualDissectionCells_Imp, virtualDissectionCells_not_Imp)])
    M_val_Imp = rowMeans(imp[, virtualDissectionCells_Imp]) - rowMeans(imp[, virtualDissectionCells_not_Imp])
    
    
    MA_df_imp = data.frame(
      A_val = A_val_Imp,
      M_val = M_val_Imp,
      gene = rownames(imp),
      M_rank = rank(-M_val_Imp),
      M_rank_bottom = rank(M_val_Imp),
      type = "Imputed"
    )
    
    # Combine logcounts and imputed MA plots
    MA_df = rbind(MA_df_seq, MA_df_imp)
    MA_df$type <- factor(MA_df$type, levels = c("seqFISH", "Imputed"))
    
    
    
    # Plot
    g = ggplot(MA_df, aes(x = A_val, y = M_val, label = gene, tooltip = gene, data_id = gene)) + 
      geom_point_interactive(colour = "grey") + 
      geom_hline(yintercept = 0) + 
      geom_text_interactive(aes(label = gene, colour = "red"),
                            data = subset(MA_df, M_rank <= 5),
                            size = 5,
                            fontface = "italic") +
      geom_text_interactive(aes(label = gene, colour = "blue"),
                            data = subset(MA_df, M_rank_bottom<= 5),
                            size = 5,
                            fontface = "italic") +
      theme_classic() +
      xlab("Mean expression of Group A and Group B cells") + 
      ylab("Difference in means of Group A versus Group B cells") +
      ggtitle("MA-Plot of expression logcounts") +
      
      scale_color_identity(name = "",
                           breaks = c("red", "blue"),
                           labels = c("Higher in Group A", "Higher in Group B"),
                           guide = "legend") +
      theme(legend.position = "none") +
      theme(legend.text = element_text(size = 10)) +
      
      guides(
        colour = guide_legend(
          title = "",
          override.aes = aes(label = "X",
                             size = 8)
        )) +
      
      theme(axis.title.y = element_text(size = 10)) +
      theme(axis.title.x = element_text(size = 10)) +
      facet_wrap(~type, nrow = 1) +
      NULL
    g
    
  }
  
  VirtualDissectionMAPlot = eventReactive(updateVirtualDissection(),
                                          VirtualDissectionMAPlotGenerator()
  )
  
  # comment testing ggtips
  output$VirtualDissectionMAPlot <- renderGirafe({
    VirtualDissectionMAPlot = VirtualDissectionMAPlot()
    girafe(code = print(VirtualDissectionMAPlot),
           options = list(
             opts_hover(css = "fill:#000000;stroke:black;cursor:pointer;", reactive = TRUE),
             opts_selection(
               type = "none", css = "fill:#000000;stroke:black;")
           ))
  })
  
  mRNARegionPlot_dblclickGenerator = function(coord = "mRNARegionPlot_dblclick") {
    
    dblclick_df = nearPoints(meta,
                             coordinfo = input[[coord]],
                             xvar = "x_global_affine",
                             yvar = "y_global_affine",
                             maxpoints = 1
    )
    uniq = dblclick_df[,"uniqueID"]
    
    observeEvent(input$mRNARegionPlot_dblclick, {
      newval_x = meta[uniq,"x_global_affine"]# - (input$embryo1_centre_x + radius)
      newval_y = meta[uniq,"y_global_affine"]# - (-input$embryo1_centre_y + radius)
      updateSliderInput(session,
                        paste0(meta[as.character(uniq), "embryo"],
                               "_centre_x"),
                        value = newval_x,
                        min = embryo_coords_range_x[1],
                        max = embryo_coords_range_x[2]
      )
      updateSliderInput(session,
                        paste0(meta[as.character(uniq), "embryo"],
                               "_centre_y"),
                        value = newval_y,
                        min = embryo_coords_range_y[1],
                        max = embryo_coords_range_y[2]
      )
    })
    
    
    return(uniq)
    
  }
  
  
  # this is here so that double click in digital in situ works
  output$info <- renderText({
    xy_str <- function(e) {
      if(is.null(e)) return("NULL\n")
      paste0("x=", round(e$x, 1), " y=", round(e$y, 1), "\n")
    }
    
    paste0(
      "dblclick: ", xy_str(input$mRNARegionPlot_dblclick)
    )
    
    uniq = mRNARegionPlot_dblclickGenerator()
    paste0("")
  })
  
  
  
  
  
  
  ## AP DV plot
  
  
  APDVPlotGenerator = function() {
    
    # Cells to subset
    if(input$celltype_subset_all_apdv) {
      celltype_subset = celltypes
    } else {
      celltype_subset = input$apdv_plot_celltype
    }
    
    # Filter cell type and embryo
    meta_sub = meta %>% 
      filter(embryo %in% input$embryo_subset) %>% 
      filter(cellType %in% celltype_subset)
    
    
    # Add gene expression 
    add_imp()
    validate(need(input$apdv_plot_genename %in% rownames(imp), "No valid gene selected."))
    meta_sub = cbind(meta_sub[c("cellType", input$ap_or_dv)], t(imp)[meta_sub$uniqueID, input$apdv_plot_genename])
    colnames(meta_sub)[3:length(colnames(meta_sub))] = input$apdv_plot_genename
    
    meta_sub = meta_sub %>% 
      pivot_longer(!c("cellType", input$ap_or_dv), names_to = "genes", values_to = "expression")
    
    
    
    plots = lapply(input$apdv_plot_genename, function(geneName){
      
      data = meta_sub %>% 
        filter(genes == geneName)
      
      p <- ggplot(data, aes(x = .data[[input$ap_or_dv]], y = expression)) +
        geom_point(aes(text = paste0(input$ap_or_dv, " axis: ", round(.data[[input$ap_or_dv]], 2), "<br>",
                          genes," expression: ", round(expression, 2), "<br>",
                          "Cell Type: ", cellType, "<br>"))) +
        geom_smooth() +
        theme_classic() +  
        labs(x = paste(input$ap_or_dv, "Axis"),
             y = paste("Gene:", geneName))
      
      # Convert ggplot object to plotly
      return(ggplotly(p, tooltip = "text"))
      
    })
    
    final_plot = subplot(plots,
                         nrows = length(input$apdv_plot_genename),
                         shareX = TRUE,
                         titleY = TRUE)
    
    return(final_plot)
    
  }
  
  output$APDV_plot <- renderPlotly(
    
      APDVPlotGenerator()
  
  )
  

  
  
})





