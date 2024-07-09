library(ggplot2)
library(patchwork)
library(ggpubr)
library(DelayedArray)
library(HDF5Array)
library(shinycssloaders)
library(ggiraph)
library(Matrix)
library(SingleCellExperiment)

source("celltype_colours.R")

colours = c("grey85","cornflowerblue", "black")

# function to check if colours are legit colours
areColours <- function(x) {
    sapply(x, function(X) {
        tryCatch(is.matrix(col2rgb(X)), 
                 error = function(e) FALSE)
    })
}
subsetColours = function(x) {
    # given a vector of colours remove the ones that aren't colours
    # give defaults if none or only one are given
    if (x[1] == "jonny") return(c("gray75", "cornflowerblue", "black"))
    if (x[1] == "default") return(c("yellow", "red"))
    xx = x[areColours(x)]
    if (length(xx) >= 2) return(xx)
    if (length(xx) == 1) return(c(xx, "white"))
    if (length(xx) == 0) return(c("white", "white"))
}



# Reading data
seqFISH_spe = readRDS("data/seqFISH_spe_final.Rds")


meta = seqFISH_spe %>% 
  colData() %>% 
  data.frame()

# Creating cell type column
# meta$cellType = meta$extended_atlas_celltype
# celltypes = unique(meta$extended_atlas_celltype)
meta$cellType = meta$refined_annotation
celltypes = sort(unique(meta$refined_annotation))

meta$uniqueID = rownames(meta)
meta$selected = factor("Unselected", levels = c("Group A", "Group B", "Unselected"))


genes = sort(rownames(seqFISH_spe))
embryos = paste0("embryo", 1:7)
# zvals = c(2,5)

embryo_coords_range_x = range(meta$dim1, na.rm = TRUE)
embryo_coords_range_y = range(meta$dim2, na.rm = TRUE)


# commented may 2024
file_combined = "data/combined_compressed.h5"
cnames = readRDS("data/combined_cnames.Rds")
rnames = readRDS("data/combined_rnames.Rds")
genes_imp = rnames


## TEST UMAPDATA, MUST REMOVE
# set.seed(10)
# meta$UMAP1 = rnorm(length(meta$dim1))
# meta$UMAP2 = rnorm(length(meta$dim1))
meta$UMAP1 = meta$UMAP_1
meta$UMAP2 = meta$UMAP_2


##################

subsetCellsLogical = function(textToParse, input) {
    
    require(reshape)
    # meta, exprs, rnames, and imp are taken from the global environment
    # example inputs, remember quotes need to be escaped
    # textToParse = "Tbx18 < 0.5 & Shh > 0.5"
    # textToParse = "rank(x) < 10"
    # textToParse = "cluster == 4"
    # textToParse = "uniqueID %in% c(\"embryo3_Pos3_cell348_z2\")"
    
    # check if any are viable gene names
    textSplit = strsplit(textToParse, "[ \t\r\n]|=|<|>|!|\\(|\\)|&")[[1]]
    genesParse = intersect(textSplit, union(rnames, rownames(exprs)))
    print(paste(c("subsetting according to these genes:", genesParse), collapse = " "))
    
    # append these genes to the colData, prioritise measured over imputed
    colsToAddSeqFISH = as.matrix(exprs[intersect(genesParse, rownames(exprs)),,drop = FALSE])
    colsToAddImputed = as.matrix(imp[intersect(setdiff(genesParse, rownames(exprs)), rnames),,drop = FALSE])
    dfParse = cbind(meta, t(colsToAddSeqFISH), t(colsToAddImputed))
    
    dfParse_sub = subset(dfParse, embryo %in% input$embryo_subset)
    
    selectedCells <- tryCatch(rownames(subset(dfParse_sub, eval(parse(text = textToParse)))), error = function(e) e)
    
    if (any(class(selectedCells) == "error")) {
        showNotification("Error in logical statement, e.g. gene not found, please try again...")
        return(NULL)
    }
    return(selectedCells)
}

# NOTE: Switched the labels to match the paper.
embryolabeller = function(string) {
    newstring = string
    newstring[newstring == "embryo1"] <- "Embryo 7"
    newstring[newstring == "embryo2"] <- "Embryo 6"
    newstring[newstring == "embryo3"] <- "Embryo 5"
    newstring[newstring == "embryo4"] <- "Embryo 4"
    newstring[newstring == "embryo5"] <- "Embryo 3"
    newstring[newstring == "embryo6"] <- "Embryo 2"
    newstring[newstring == "embryo7"] <- "Embryo 1"
    return(newstring)
}

# commented may 2024
add_boundary_polygons = function() {
    if (!"boundary_polygons" %in% ls(envir = .GlobalEnv)) {
        showNotification("Loading cell segmentation...")
        boundary_polygons <<- meta #same info as meta
        showNotification("Loading cell segmentation... done!")
    }
}

add_exprs = function() {
    if (!"exprs_counts" %in% ls(envir = .GlobalEnv)) {
        showNotification("Loading expression...")
        exprs_counts <<- counts(seqFISH_spe)
        showNotification("Loading expression... done!")
    }
}

# commented out may 2024, dont have imputed data yet.
add_imp = function() {
    if (!"imp" %in% ls(envir = .GlobalEnv)) {
        showNotification("loading imputed data...")
        imp <<- HDF5Array(filepath = file_combined, name = "logcounts")
        rownames(imp) <<- rnames
        colnames(imp) <<- cnames
        showNotification("loading imputed data... done!")
    }
}

add_exprs_norm = function() {
    if (!"exprs_norm" %in% ls(envir = .GlobalEnv)) {
        showNotification("Loading expression logcounts...")
        exprs_norm <<- logcounts(seqFISH_spe)
        showNotification("Loading expression logcounts... done!")
    }
}

# commented 2024
# add_mRNA_df = function() {
#     if (!"mRNA_df" %in% ls(envir = .GlobalEnv)) {
#         showNotification("Loading mRNA data...")
#         mRNA_df <<- readRDS("data/mRNA.Rds")
#         showNotification("Loading mRNA data... done!")
#     }
# }

# commented 2024
g_leg_list = readRDS("data/g_leg_list.Rds")
