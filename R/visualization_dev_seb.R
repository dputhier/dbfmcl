#################################################################
##    Define the plot_heatmap_non_interactive
#################################################################

#' @title
#' plot_heatmap_non_interactive
#' @description
#' Plot the observed and simulated distance with the Kth nearest neighbors.
#' @param object A ClusterSet object.
#' @param center A logical to indicate whether to center row.. 
#' @param ceil A value for ceiling (NULL for no ceiling). Ceiling is performed centering.
#' @param floor A value for flooring (NULL for no flooring). Flooring is performed after centering.
#' @param cell_clusters A vector of cell clusters with cell barcodes as names.
#' @param show_dendro A logical to indicate whether to show column dendrogram.
#' @param gene_clusters A cluster id to plot. Default is NULL for plotting all cluster.
#' @param use_top_genes A logical to indicate whether to use highly similar genes in the slot top_genes of ClusterSet.
#' @param name A title for the heatmap.
#' @param xlab A title for the x axis.
#' @param ylab A title for the y axis.
#' @param colorbar_name A title for the colorbar.
#' @param show_legend A logical to indicate whether to show colorbar.
#' @param colors A vector of colors.
#' @param colors_cell_clusters A vector of colors for column annotations.
#' @param row_labels A logical to indicate whether to show row labels.
#' @param col_labels A logical to indicate whether to show col labels.
#' @param label_size A value for label font size.
#' @param line_size_vertical An integer for the size of horizontal white line which separate gene clusters.
#' @param line_size_horizontal An integer for the size of vertical white line  which separate cell clusters.
#'
#' @return pheatmap-class object.
#' @export
#'
#' @examples
#' m <- create_3_rnd_clust()
#' 
#' res <- find_gene_clusters(data=m,
#'                              distance_method="pearson",
#'                              inflation = 2,
#'                              k=75,
#'                              row_sum=-Inf,
#'                              highest=0.3,
#'                              min_nb_supporting_cell = 0,
#'                              fdr = 1e-8)
#' plot_heatmap_non_interactive(object = res)
#' plot_heatmap_non_interactive(object = res, cluster = "1")
#' 

#' @rdname plot_heatmap_non_interactive

plot_heatmap_non_interactive <- function(
  object,
  center = TRUE,
  ceil = 1,
  floor = -1,
  cell_clusters = NULL,
  show_dendro = TRUE,
  gene_clusters = NULL,
  use_top_genes = FALSE,
  name = NULL,
  xlab = NULL,
  ylab = NULL,
  colorbar_name = "Expression level",
  show_legend = TRUE,
  colors = c("#A9D6E5", "#2166AC", "#000000", "#B2182B", "#FFCA3A"),
  colors_cell_clusters = c("#9F1717", "#AE5B11", "#C48D00", "#517416", 
                           "#115C8A", "#584178", "#9D1C70", "#E96767", 
                           "#EC9342", "#FFCA3A", "#8AC926", "#4DADE8",
                           "#9579B9", "#E25CB4", "#DB2020", "#DA7316", 
                           "#F0AE00", "#6D9D1E", "#1882C0", "#71529A",
                           "#D02494", "#EF9292", "#F2B57D", "#FFDA77", 
                           "#B6E36A", "#7BC4EE", "#AD98C9", "#EA8AC9"),
  col_labels = FALSE,
  label_size = 9,
  line_size_vertical = 15,
  line_size_horizontal = 15) {
  
  m <- object@data
  
  # # Config
  if (is.null(gene_clusters)) {
    gene_clusters <- object@gene_clusters_metadata$cluster_id
  }
  
  # Cell order
  if (!is.null(cell_clusters)){
    print_msg("Ordering cells.", msg_type="INFO")
    m <- m[,names(sort(cell_clusters))]
  } else {
    dist_cells <- cor(m, method = "pearson")
    dist_cells <- as.dist((1-dist_cells)/2)
    hclust_cells <- hclust(dist_cells, method = "complete")
  }
  
  # Centering
  if(center){
    print_msg("Centering matrix.", msg_type="INFO")
    m <- t(scale(t(m), center = TRUE, scale = FALSE))
  }
  
  # Ceiling and flooring
  if(!is.null(ceil)){
    print_msg("Ceiling matrix.", msg_type="INFO")
    m[m > ceil] <- ceil
  }
  
  if(!is.null(floor)){
    print_msg("Flooring matrix.", msg_type="INFO")
    m[m < floor] <- floor
  }  
  
  
  # Reduce matrix to gene clusters in gene_clusters parameter
  if(!is.null(gene_clusters)){
    gene_cl_int <- unlist(object@gene_clusters[gene_clusters], use.names = FALSE)
    m <- m[gene_cl_int,]
  }
  
  # Reduce m rows to only keep genes from top_genes
  if(use_top_genes) {
    if (length(object@top_genes) == 0) {
      print_msg("The slot top_genes of the input ClusterSet object is empty. Be sure to run top_genes() before.", msg_type = "STOP")
    }
    
    if(is.null(gene_clusters)){
      genes_top <- unlist(object@top_genes, use.names = FALSE)
      m <- m[genes_top,]
    } else {
      genes_top <- unlist(object@top_genes[gene_clusters], use.names = FALSE)
      m <- m[genes_top,]
    }
  }
  
  
  # Add blank row to separate gene clusters in heatmap
  if(!length(gene_clusters) == 0){
    gap_rows <- head(as.numeric(cumsum(object@gene_clusters_metadata$size[gene_clusters])), -1)
  }else{
    gap_rows <- NULL
  }
  
  
  # Compute gaps to add on heatmap columns
  if(!is.null(cell_clusters)){
    #TODO stop when length(cell_cluster) > ncol(m) 
    if(!length(cell_clusters) == ncol(m)){
      print_msg(paste("The number of cells provided in argument",
                      "'cell_clusters' is not equal to the number of",
                      "cells in the ClusterSet object"), msg_type = "STOP")
    }
    # Create the annotation dataframe for pheatmap
    cell_clusters_ordered <- cell_clusters[hclust_cells$order]
    gap_order <- table(cell_clusters_ordered)[unique(cell_clusters_ordered)]
    gap_cols <- head(as.numeric(cumsum(gap_order)), -1)
    annotation_cols <- data.frame(cell_clusters = cell_clusters_ordered)
    annotation_cols$cell_clusters <- as.factor(annotation_cols$cell_clusters)
    annotation_cols[hclust_cells$order,]
    col_to_use <- colors_cell_clusters[1:max(cell_clusters)]
    names(col_to_use) <- unique(cell_clusters)
    ann_colors = list(cell_clusters = col_to_use)
  }else{
    gap_cols <- NULL
  }
  
  ####### Heatmap #######
  # Main heatmap
  print_msg("Plotting heatmap.", msg_type="INFO")
  #   if(!is.null(cell_clusters)){
  #     cell_clusters_anno <- cell_clusters[match(cell_names_blank, names(cell_clusters))]
  #     htmp <- htmp %>% add_col_annotation(as.factor(cell_clusters_anno), colors = list(colors_cell_clusters))
  #   }
  #   
  #   if(show_dendro & is.null(cell_clusters)) {
  #     htmp <- htmp %>% iheatmapr::add_col_dendro(hclust_cells, reorder = FALSE)
  #   }
  #   
  #   print_msg("Adding Titles.", msg_type="DEBUG")
  #   
  #   if(!is.null(ylab)){
  #     htmp <- htmp %>% add_row_title(ylab, side="right", font = list(size = 12))}
  #   if(!is.null(xlab)){
  #     htmp <- htmp %>% add_col_title(xlab, side="top", font = list(size = 12))}
  #   if(!is.null(name)){
  #     htmp <- htmp %>% add_col_title(name, side="top", font = list(size = 24))}
  # Reorder rows to get the same order as the interactive heatmap
  
  if(show_dendro){
    htmp <- pheatmap(mat = m, 
                     annotation_legend = show_legend, legend = show_legend,
                     color = colorRampPalette(colors)(50),
                     cluster_rows = F, cluster_cols = hclust_cells, 
                     show_rownames = F, show_colnames = col_labels, 
                     border_color = NA, scale = "none", na_col = "white", 
                     gaps_col = gap_cols, gaps_row = gap_rows, 
                     annotation_col = annotation_cols,
                     annotation_colors = ann_colors)
  }else{
    htmp <- pheatmap(mat = m[, hclust_cells$order], 
                     annotation_legend = show_legend, legend = show_legend,
                     color = colorRampPalette(colors)(50),
                     cluster_rows = F, cluster_cols = F, 
                     show_rownames = F, show_colnames = col_labels, 
                     border_color = NA, scale = "none", na_col = "white", 
                     gaps_col = gap_cols, gaps_row = gap_rows, 
                     annotation_col = annotation_cols,
                     annotation_colors = ann_colors)
  }
  
  return(htmp)
}

##### FOR TEST ONLY
# cell_clusters <- c(rep(1, 10), rep(2, 10))
# cell_clusters <- c(rep(1, 4), rep(2, 6), rep(3, 5), rep(4, 5))
# names(cell_clusters) <- colnames(res@data)
m <- create_4_rnd_clust()

object <- find_gene_clusters(data=m,
                             distance_method="pearson",
                             inflation = 2,
                             k=75,
                             row_sum=-Inf,
                             highest=0.3,
                             min_nb_supporting_cell = 0,
                             fdr = 1e-8)
center = TRUE
ceil = 3
floor = -3
show_dendro = TRUE
gene_clusters = NULL
use_top_genes = FALSE
name = NULL
xlab = NULL
ylab = NULL
colorbar_name = "Expression level"
show_legend = TRUE
colors = c("#A9D6E5", "#2166AC", "#000000", "#B2182B", "#FFCA3A")
colors_cell_clusters = c("#9F1717", "#AE5B11", "#C48D00", "#517416", 
                         "#115C8A", "#584178", "#9D1C70", "#E96767", 
                         "#EC9342", "#FFCA3A", "#8AC926", "#4DADE8",
                         "#9579B9", "#E25CB4", "#DB2020", "#DA7316", 
                         "#F0AE00", "#6D9D1E", "#1882C0", "#71529A",
                         "#D02494", "#EF9292", "#F2B57D", "#FFDA77", 
                         "#B6E36A", "#7BC4EE", "#AD98C9", "#EA8AC9")
col_labels = FALSE
label_size = 9
line_size_vertical = 15
line_size_horizontal = 15