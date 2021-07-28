################################################################
##        Main script of the R PACKAGE : DBFMCL
##
## Authors : BERGON A, J. BAVAIS
##  with the collaboration of LOPEZ F., TEXTORIS J. and PUTHIER D.
##
##
## R CMD SHLIB dbf.c -o dbf
#################################################################
##    UTILS (a set of useful function)
#################################################################
# 0 : No message
# 1 : info
# 2 : DEBUG

VERBOSITY_DBFMCL = 3

print_msg <- function(msg, msg_type="INFO"){
  if(msg_type == "INFO")
    if(VERBOSITY_DBFMCL > 0)
      cat(paste("|-- ", msg, "\n"))
  if(msg_type == "DEBUG")
    if(VERBOSITY_DBFMCL > 1)
      cat(paste("|-- ", msg, "\n"))
  if(msg_type == "WARNING")
      cat(paste("|-- ", msg, "\n"))
}
#################################################################
##    DEFINITION OF A SPECIFIC CLASS OBJECT : ClusterSet
#################################################################
 
library(Biobase)
library(ggplot2)
library(reshape2)
library(dplyr)
library(MCL)
library(igraph)
library(gprofiler2)

#' @title
#' ClusterSet
#' @description
#' This class is a representation of a partitioning algorithm and is intented to store gene clusters.
#' @slot name character. The original input file name (if applicable).
#' @slot data matrix. The matrix containing the filtered/partitionned data.
#' @slot cluster vector. Mapping of row/genes to clusters.
#' @slot size vector. The size of each cluster.
#' @slot center matrix. The centers of each clusters.
#' @slot parameters list. The parameter used.
#' @slot algorithm vector. The algorithm used to produce the clusters.
#' @slot cell_types vector. The cell types.
#' @slot cell_colors vector. The cell types to color mapping.
#' @slot cell_order vector. How cell should be ordered.
#' @slot cluster_annotations list. Functional annotation of clusters.
#' @return A ClusterSet object.
#' @export
#'
#' @examples
#' 
#' \dontrun{
#'   m <- matrix(rnorm(80000), nc=20)
#'    m[1:100,1:10] <- m[1:100,1:10] + 4
#'    m[101:200,11:20] <- m[101:200,11:20] + 3
#'    m[201:300,5:15] <- m[201:300,5:15] + -2
#'    res <- DBFMCL(data=m,
#'                  distance_method="pearson",
#'                  av_dot_prod_min = 0,
#'                  inflation = 1.2,
#'                  k=25,
#'                  fdr = 10)
#' plot_clust(res, ceil = 10, floor = -10)
#' plot_clust(res, type="tile", ceil = 10, floor = -10)
#' write_clust(res, filename_out = "ALL.sign.txt")
#'   is(res)
#' }
#'               
setClass("ClusterSet",
  representation = list(
    name = "character",
    data = "matrix",
    cluster = "vector",
    size = "vector",
    center = "matrix",
    parameters = "list",
    algorithm = "vector",
    cell_types = "vector",
    cell_colors = "vector",
    cell_order = "vector",
    cluster_annotations = "list"
  ),
  prototype = list(
    name = character(),
    data = matrix(nr = 0, nc = 0),
    cluster = numeric(),
    size = numeric(),
    center = matrix(nc = 0, nr = 0),
    parameters = list(),
    algorithm = character(),
    cell_types = vector(),
    cell_colors = vector(),
    cell_order = vector(),
    cluster_annotations = list()
  )
)

#################################################################
##    REDEFINE SHOW METHOD FOR CLASS OBJECT : ClusterSet
#################################################################

setMethod(
  "show", signature("ClusterSet"),

  function(object) {
    
    cat("\t\tAn object of class ClusterSet\n")
    cat("\t\tName:", slot(object, "name"), "\n")
    cat("\t\tMemory used: ", object.size(object), "\n")
    cat("\t\tNumber of cells: ", ncol(slot(object, "data")), "\n")
    cat(
      "\t\tNumber of informative genes: ",
      nrow(slot(object, "data")), "\n"
    )
    cat("\t\tNumber of clusters: ", length(slot(object, "size")), "\n")
    cat("\t\tThis object contains the following informations:\n")

    for (i in slotNames(object)) {
      cat("\t\t\t - ", i, "\n")
    }
    if (length(slot(object, "parameters")) > 0) {
      for (i in 1:length(slot(object, "parameters"))) {
        cat(
          "\t\t\t\t * ", names(slot(object, "parameters"))[[i]],
          " = ", slot(object, "parameters")[[i]], "\n"
        )
      }
    }
  }
)


#################################################################
##    REDEFINE dim/ncol/nrow METHOD FOR CLASS OBJECT : ClusterSet
#################################################################

#' @title
#' ncol
#' @description
#' The number of columns of a ClusterSet object.
#' @param x A ClusterSet object.
#' @return The number of columns.
#' @export
#'
#' @examples
#'
#' \dontrun{
#'   m <- matrix(rnorm(80000), nc=20)
#'    m[1:100,1:10] <- m[1:100,1:10] + 4
#'    m[101:200,11:20] <- m[101:200,11:20] + 3
#'    m[201:300,5:15] <- m[201:300,5:15] + -2
#'    res <- DBFMCL(data=m,
#'                 distance_method="pearson",
#'                  av_dot_prod_min = 0,
#'                  inflation = 1.2,
#'                  k=25,
#'                  fdr = 10)
#'   ncol(res)
#' }
#'
setMethod(
  "ncol", signature("ClusterSet"),

  function(x) {

      return(ncol(x@data))

      }
)

#' @title
#' nrow
#' @description
#' The number of rows of a ClusterSet object.
#' @param x A ClusterSet object.
#' @return The number of rows.
#' @export
#'
#' @examples
#'
#' \dontrun{
#'   m <- matrix(rnorm(80000), nc=20)
#'    m[1:100,1:10] <- m[1:100,1:10] + 4
#'    m[101:200,11:20] <- m[101:200,11:20] + 3
#'    m[201:300,5:15] <- m[201:300,5:15] + -2
#'    res <- DBFMCL(data=m,
#'                  distance_method="pearson",
#'                  av_dot_prod_min = 0,
#'                  inflation = 1.2,
#'                  k=25,
#'                  fdr = 10)
#'   nrow(res)
#' }
#'
setMethod(
  "nrow", signature("ClusterSet"),

  function(x) {

      return(nrow(x@data))

      }
)


#' @title
#' dim
#' @description
#' The number of rows/columns of a ClusterSet object.
#' @param x A ClusterSet object.
#' @return The number of rows/columns.
#' @export
#'
#' @examples
#'
#' \dontrun{
#'   m <- matrix(rnorm(80000), nc=20)
#'    m[1:100,1:10] <- m[1:100,1:10] + 4
#'    m[101:200,11:20] <- m[101:200,11:20] + 3
#'    m[201:300,5:15] <- m[201:300,5:15] + -2
#'   res <- DBFMCL(data=m,
#'                  distance_method="pearson",
#'                  av_dot_prod_min = 0,
#'                  inflation = 1.2,
#'                  k=25,
#'                 fdr = 10)
#'   dim(res)
#' }
#'
setMethod(
  "dim", signature("ClusterSet"),

  function(x) {

      return(dim(x@data))

      }
)


#################################################################
##    Define the load_seurat function for class ClusterSet
#################################################################

#' @title
#' load_seurat
#' @description
#' Load a seurat object into a DBFMCL. At the moment the objective is mainly
#' to store cell identity (i.e cell types/groups to barcode mapping) and
#' cell type to color mapping.
#' @param object A ClusterSet object.
#' @param seurat_obj A seurat object to extract cell to group mapping.
#' @param dimplot_obj To display cell clusters on the profile diagram with colors extracted from Dimplot output (provide also a Seurat object).
#
#' @return A ClusterSet object.
#' @export
#'
#' @examples
#' # see online examples
setGeneric("load_seurat",
    function(object,
            seurat_obj=NULL,
            dimplot_obj=NULL) {
      standardGeneric("load_seurat")
})

#' @rdname load_seurat
setMethod("load_seurat",
    signature(object = "ClusterSet"),
    function(object,
            seurat_obj=NULL,
            dimplot_obj=NULL) {
      if (!inherits(seurat_obj, "Seurat")) {
        stop("Please provide a Seurat and patchwork object.")
      }
      if (!inherits(dimplot_obj, "patchwork")) {
        stop("Please provide a Seurat and patchwork object.")
      }
      if(ncol(object) != ncol(seurat_obj)){
        stop("The number of cells is not the same in DBFMCL and Seurat Object.")
      }

      g <- ggplot2::ggplot_build(dimplot_obj)
      tmp_mat <- dplyr::distinct(as.data.frame(cbind(g$data[[1]]$colour,
                                              as.character(g$plot$data$ident))))
      cell_col_tmp <- tmp_mat[,1]
      cell_grp_tmp <- tmp_mat[,2]
      object@cell_colors <- stats::setNames(as.character(cell_col_tmp), cell_grp_tmp)

      tmp_mat <- dplyr::distinct(as.data.frame(cbind(rownames(g$plot$data), as.character(g$plot$data$ident))))
      object@cell_types <- stats::setNames(as.character(tmp_mat[,2]), tmp_mat[,1])

      object@cell_order <- rownames(g$plot$data)[order(g$plot$data$ident)]

      return(object)
  }
)

#################################################################
##    Define the plot_profile for class ClusterSet
#################################################################

#' @title
#' plot_clust
#' @description
#' Plot the results (heatmap or profiles) contained in a ClusterSet object.
#' @param object A ClusterSet object.
#' @param type The type of diagram ("line or "tile").
#' @param to_log2 Whether data should be transform in logarithm base 2 (+ 1 as a pseudocount).
#' @param average_only Only display the average gene expression profile (if type="line").
#' @param average_line_color The color of the average profile.
#' @param standardizing Whether rows should be divided by standard deviation.
#' @param ceil A value for ceiling (NULL for no ceiling). Ceiling is performed after log transformation, centering and standardization.
#' @param floor A floor value (NULL for no ceiling). flooring is performed after log transformation, centering and standardization.
#' @param centering Whether rows should be centered. 
#'
#' @return A ggplot diagram.
#' @export
#'
#' @examples
#' # see online examples
setGeneric("plot_clust",
           
           function(object,
                    type = c("line", "tile"),
                    to_log2 = FALSE,
                    average_only=FALSE,
                    average_line_color="skyblue4",
                    standardizing = FALSE,
                    ceil=1,
                    floor=-1,
                    centering = TRUE) {
             
              standardGeneric("plot_clust")
})

#' @rdname plot_clust
setMethod(
  "plot_clust",
  signature(object = "ClusterSet"),
  function(object,
           type = c("line", "tile"),
           to_log2 = FALSE,
           average_only=FALSE,
           average_line_color="skyblue4",
           standardizing = FALSE,
           ceil=1,
           floor=-1,
           centering = TRUE) {

    # The type of diagram
    type <- type[1]

    print_msg("getting matrix", msg_type="DEBUG")

    nb <- length(object@size)
    m <- object@data
    if(length(object@cell_order) != 0){
      m <- m[, object@cell_order]
    }

    if (to_log2) {
      m <- log2(m + 1)
    }

    ## median-centering of row
    if (centering) {
      print_msg("Median-centering rows.", msg_type="DEBUG")
      mean_row <- apply(m, 1, mean)
      m <- sweep(m, MARGIN = 1, STATS = mean_row, FUN = "-")
    }

    ## Standardizing row
    if (standardizing) {
      print_msg("Standardizing rows.", msg_type="DEBUG")
      sd_row <- apply(m, 1, sd)
      m <- sweep(m, MARGIN = 1, STATS = sd_row, FUN = "/")
    }

    ## Ceiling / flooring
    if(!is.null(ceil)){
      print_msg("Ceiling matrix.", msg_type="DEBUG")
      m[m > ceil] <- ceil
    }

    if(!is.null(floor)){
      print_msg("Flooring matrix.", msg_type="DEBUG")
      m[m < floor] <- floor
    }

    ## melting
    print_msg("Melting matrix.", msg_type="DEBUG")
    m_melt <- as.data.frame(m)
    m_melt$cluster <- object@cluster
    m_melt$gene <- row.names(object@data)

    m_melt <- melt(m_melt,
      id.vars = c("cluster", "gene"),
      variable.name = "samples"
    )

    print_msg("Storing cell types.", msg_type="DEBUG")
    ## Storing cell types:

    if(!is.null(object@cell_types)){
      m_melt$cell_types <- as.character(object@cell_types[as.character(m_melt$samples)])
    }else{
      print_msg("Warning: cell type is undefined.", msg_type="WARNING")
        m_melt$cell_types <- "unknown_cell_type"
    }

    ## plotting
    # Note that samples, value, gene, cluster
    # may appear as undefined variable to "R check" command.
    # A workaround is to define them as NULL first...
    samples <- value <- gene <- cluster <- cluster_mean <- NULL
    if (type == "line") {
      print_msg("Preparing diagram.", msg_type="DEBUG")
      p <- ggplot(data = m_melt, aes(
        x = samples,
        y = value
      ))

      ## displaying cell types:
      if(length(object@cell_types) > 0 & length(object@cell_colors) > 0){
        print_msg("Adding cell populations to the diagram.", msg_type="DEBUG")
        cell_types <- NULL # Avoid "no visible binding for global variable" inn R check.
        p <- p + ggplot2::geom_vline(aes(xintercept= samples, color=cell_types))
        p <- p + ggplot2::scale_color_manual(values=object@cell_colors,  guide = ggplot2::guide_legend(override.aes = list(size = 5)))
      }

      if(! average_only){
        print_msg("Adding gene profile.", msg_type="DEBUG")
        p <- p + ggplot2::geom_line(color = "azure3", aes(group = gene), size=0.1)
      }

      print_msg("Adding average profile.", msg_type="DEBUG")

      p <- p + ggplot2::geom_line(
        data = m_melt %>%
          group_by(cluster, samples) %>%
          summarise(cluster_mean = mean(value)),
        aes(
          x = samples,
          y = cluster_mean,
          group = cluster
        ),
        color = average_line_color,
    size=0.2
      )

      print_msg("Faceting.", msg_type="DEBUG")
        p <- p + facet_grid(cluster ~ .)

      print_msg("Theming.", msg_type="DEBUG")
      p <- p + theme_bw()
      p <- p + theme(
        strip.text.y = element_text(angle = 0),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )


    } else if (type == "tile") {
      print_msg("Preparing diagram (tile).", msg_type="DEBUG")
      p <- ggplot(
        data = m_melt,
        aes(
          x = samples,
          y = gene,
          fill = value
        )
      )

      print_msg("Preparing color palette.", msg_type="DEBUG")
      col <- unlist(strsplit("#67001f,#b2182b,#d6604d,#f4a582,#fddbc7,#f7f7f7,#d1e5f0,#92c5de,#4393c3,#2166ac,#053061", ","))
      color.ramp <- colorRampPalette(col)(10)
      p <- p + geom_tile()
      p <- p + theme_bw()
      p <- p + scale_fill_gradientn(
        colours = color.ramp,
        name = "Signal"
      )

      print_msg("Theming.", msg_type="DEBUG")
      p <- p + theme(
        strip.text.y = element_text(angle = 0),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank()
      )
    }

    print_msg("Adding facets.", msg_type="DEBUG")
    p <- p + facet_grid(cluster ~ ., scales = "free_y")

    return(p)
  }
)

#################################################################
##    Define the write function for class ClusterSet
#################################################################


#' @title
#' write_clust
#' @description
#' Write a ClusterSet into a flat file.
#' @param object ClusterSet. 
#' @param filename_out The outfile name.
#' @param out_path The path to the file.
#' @param nb_na_row Number of separating rows (containing NAs).
#' @return Write a file.
#' @export
#'
#' @examples
#' # see online help.
setGeneric("write_clust",
           
           function(object, 
                    filename_out = NULL,
                    out_path = ".",
                    nb_na_row=3) {
                standardGeneric("write_clust")
})


#' @rdname write_clust
setMethod(
  "write_clust",
  signature(object = "ClusterSet"),
  function(object,
           filename_out = NULL,
           out_path = ".",
           nb_na_row=5) {
           
    if (out_path == ".") out_path <- getwd()

    if (is.null(filename_out)) {
      filename_out <- "exprs.dataMods.txt"
    }

    data <- object@data
    nb <- 0
    dataT <- c("clusters", colnames(data))

    ## processing data
    for (i in 1:length(object@size)) {

      print_msg(paste("Cluster ", i, " --> ", object@size[i], " probes"), msg_type="INFO")

      subData <- data[object@cluster == i, ]
      subData <- cbind(rownames(subData), subData)
      if (nb_na_row > 0){
        intLine <- matrix(rep(NA, (ncol(data) + 1)*nb_na_row), nrow = nb_na_row)
        dataT <- rbind(dataT, subData, intLine)
      }
      nb <- nb + 1

    }


    ## exporting results
    print_msg("Exporting results", msg_type="DEBUG")
    write.table(dataT, file.path(out_path, filename_out),
      col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE
    )

      print_msg(paste("--> Creating file : ",
                      file.path(out_path, filename_out)),
                msg_type="DEBUG")
  }
)

#################################################################
##    DBF-MCL
#################################################################

#' @title
#' The "Density Based Filtering and Markov CLustering" algorithm (DBF-MCL).
#' @description
#' DBFMCL is a tree-steps adaptative algorithm that \emph{(i)} find elements
#' located in dense areas (DBF), \emph{(ii)} uses selected items to construct a
#' graph, \emph{(iii)} performs graph partitioning using the Markov CLustering
#' Algorithm (MCL).
#'
#' This function requires installation of the mcl program
#' (\url{http://www.micans.org/mcl}). See "Warnings" section for more
#' informations.
#'
#' When analyzing a noisy dataset, one is interested in isolating dense regions
#' as they are populated with genes/elements that display weak distances to
#' their nearest neighbors (i.e. strong profile similarities). To isolate these
#' regions DBF-MCL computes, for each gene/element, the distance with its kth
#' nearest neighbor (DKNN).In order to define a critical DKNN value that will
#' depend on the dataset and below which a gene/element will be considered as
#' falling in a dense area, DBF-MCL computes simulated DKNN values by using an
#' empirical randomization procedure. Given a dataset containing n genes and p
#' samples, a simulated DKNN value is obtained by sampling n distance values
#' from the gene-gene distance matrix D and by extracting the kth-smallest
#' value. This procedure is repeated n times to obtain a set of simulated DKNN
#' values S. Computed distributions of simulated DKNN are used to compute a FDR
#' value for each observed DKNN value. The critical value of DKNN is the one
#' for which a user-defined FDR value (typically 10\%) is observed. Genes with
#' DKNN value below this threshold are selected and used to construct a graph.
#' In this graph, edges are constructed between two genes (nodes) if one of
#' them belongs to the k-nearest neighbors of the other. Edges are weighted
#' based on the respective coefficient of correlation (\emph{i.e.}, similarity)
#' and the graph obtained is partitioned using the Markov CLustering Algorithm
#' (MCL).
#'
#' @param data a \code{matrix}, \code{data.frame} or \code{Seurat} object.
#' @param filename a character string representing the file name.
#' @param name a prefix for the names of the intermediary files created by DBF
#' and MCL.
#' @param path a character string representing the data directory where
#' intermediary files are to be stored. Default to current working directory.
#' @param mcl_cmd_line Boolean. Whether to use the fast MCL version through command line.
#' @param mcl_cmd_line_threads If mcl_cmd_line is TRUE, how many threads should be used (integer).
#' @param distance_method a method to compute the distance to the k-th nearest
#' neighbor. One of "pearson" (Pearson's correlation coefficient-based
#' distance), "spearman" (Spearman's rho-based distance), "euclidean".
#' @param av_dot_prod_min Any cluster with average dot product below this value is discarded. This allow to delete
#' clusters in which correlation is influenced/supported by very few samples (typically 1).
#' @param min_cluster_size Minimum number of element inside a cluster. MCL tend to create lots of clusters with
#' very few (e.g 2) objects.
#' @param silent if set to TRUE, the progression of distance matrix calculation
#' is not displayed.
#' @param k the neighborhood size.
#' @param random the number of simulated distributions S to compute. By default
#' \code{random = 3}.
#' @param memory_used size of the memory used to store part of the distance
#' matrix. The subsequent sub-matrix is used to computed simulated distances to
#' the k-th nearest neighbor (see detail section).
#' @param fdr an integer value corresponding to the false discovery rate
#' (range: 0 to 100).
#' @param inflation the main control of MCL. Inflation affects cluster
#' granularity. It is usually chosen somewhere in the range \code{[1.2-5.0]}.
#' \code{inflation = 5.0} will tend to result in fine-grained clusterings, and
#' whereas \code{inflation = 1.2} will tend to result in very coarse grained
#' clusterings. By default, \code{inflation = 2.0}. Default setting gives very
#' good results for microarray data when k is set between 70 and 250.
#' @param set.seed specify seeds for random number generator.
#' @return a ClusterSets class object.
#' @section Warnings: With the current implementation, this function only works
#' only on UNIX-like plateforms.
#'
#' MCL should be installed. One can used the following command lines in a
#' terminal:
#'
#' \code{# Download the latest version of mcl (the script has been tested
#' successfully with the 06-058 version).}
#' \code{wget http://micans.org/mcl/src/mcl-latest.tar.gz}
#' \code{# Uncompress and install mcl}
#' \code{tar xvfz mcl-latest.tar.gz}
#' \code{cd mcl-xx-xxx}
#' \code{./configure}
#' \code{make}
#' \code{sudo make install}
#' \code{# You should get mcl in your path}
#' \code{mcl -h}
#' @author Bergon A., Bavais J., Textoris J., Granjeaud S., Lopez F and Puthier
#' D.
#' @references Van Dongen S. (2000) A cluster algorithm for graphs. National
#' Research Institute for Mathematics and Computer Science in the 1386-3681.
#' @keywords manip
#' @examples
#'
#' \dontrun{
#'   m <- matrix(rnorm(80000), nc=20)
#'   m[1:100,1:10] <- m[1:100,1:10] + 4
#'    m[101:200,11:20] <- m[101:200,11:20] + 3
#'    m[201:300,5:15] <- m[201:300,5:15] + -2
#'    res <- DBFMCL(data=m,
#'                  distance_method="pearson",
#'                 av_dot_prod_min = 0,
#'                 inflation = 1.2,
#'                  k=25,
#'                 fdr = 10)
#' plot_clust(res, ceil = 10, floor = -10)
#' plot_clust(res, type="tile", ceil = 10, floor = -10)
#' write_clust(res, filename_out = "ALL.sign.txt")
#' }
#'
#' @export DBFMCL
DBFMCL <- function(data = NULL, 
                   filename = NULL, 
                   path = ".",
                   mcl_cmd_line=FALSE,
                   mcl_cmd_line_threads=1,
                   name = NULL,
                   distance_method = c("pearson", "spearman", "euclidean"),
                   av_dot_prod_min=2,
                   min_cluster_size=10,
                   silent = FALSE,
                   k = 50,
                   random = 3, 
                   memory_used = 1024,
                   fdr = 10,
                   inflation = 8,
                   set.seed = 123) {

  ## testing the system
  if (.Platform$OS.type == "windows") {
    stop("\t--> A unix-like OS is required to launch mcl and cluster programs.")
  }

  ## getting parameters
  data_source <- get_data_4_DBFMCL(data = data, filename = filename, path = path)
  data_matrix <- data_source$data

  # A simple function to create a random string
  create_rand_str <- function() {
      v = c(sample(LETTERS, 3, replace = TRUE),
            sample(0:9, 4, replace = TRUE),
            sample(letters, 3, replace = TRUE))
      return(paste0(sample(v),collapse = ""))
  }

  if (is.null(name)) name <- data_source$name
  if (is.null(name)) name <- create_rand_str()

  distance_method <- match.arg(distance_method)
  txt <- paste("\n\tInflation: ", inflation, sep = "")

  ## writting all parameters

  cat(
    "The following parameters will be used :",
    "\n\tWorking directory: ", getwd(),
    "\n\tName: ", name,
    "\n\tDistance method: ", distance_method,
    "\n\tMinimum average dot product for clusters: ", av_dot_prod_min,
    "\n\tMinimum cluster size: ", min_cluster_size,
    "\n\tNumber of neighbors: ", k,
    "\n\tNumber of randomizations: ", random,
    "\n\tFDR: ", fdr, "%", txt,
    "\n\tVisualize standard outputs from both mcl and cluster",
    "commands: ", silent,
    "\n\tMemory used : ", memory_used, "\n\n"
  )


  ## DBF algorithm, returns a ClusterSet object
  obj <- DBF(data_matrix,
             name,
             distance_method = distance_method,
             silent = silent,
             k = k,
             random = random,
             memory_used = memory_used,
             fdr = fdr,
             set.seed = set.seed
  )

  dbf_out_file <- paste0(name, ".dbf_out.txt")
  mcl_out_file <- paste0(name, ".mcl_out.txt")

  print_msg("DBF completed. Starting MCL step.", msg_type="DEBUG")

  if (length(readLines(dbf_out_file)) > 0) {

      ## Launching mcl

      if(mcl_cmd_line){
        print_msg("Running MCL through the command line for best performance.")
        mcl_system_cmd(name, inflation = inflation, silent = silent, threads=mcl_cmd_line_threads)
      }else{
        print_msg("You are using the R implementation of MCL.", msg_type="WARNING")
        print_msg("Use the command line version for best performance (mcl_cmd_line)", msg_type="WARNING")
        print_msg(paste("Reading : ", dbf_out_file), msg_type="WARNING")
        graph_tab <- read.csv(dbf_out_file, sep=" ", header=F)
        colnames(graph_tab) <- c("source", "dest", "weight")
        graph_igraph <-  igraph::graph_from_data_frame(graph_tab, directed=F)
        graph_adj <- igraph::as_adj(graph_igraph, attr='weight')
        mcl_res <- MCL::mcl(graph_adj,  expansion = 2, inflation = 8, allow1 = TRUE, addLoops = FALSE)
        cluster_to_genes <- split(rownames(graph_adj), mcl_res$Cluster)
        print_msg("Writing gene clusters", msg_type="INFO")
        print_msg(mcl_out_file, msg_type="DEBUG")

        for(i in 1:length(cluster_to_genes)){
          write.table(paste(cluster_to_genes[[i]], collapse="\t"),
                      file=mcl_out_file,
                      append=TRUE,
                      col.names=FALSE,
                      row.names=FALSE,
                      quote=FALSE)
        }
      }

      print_msg(paste0("Reading MCL output: ", mcl_out_file), msg_type="DEBUG")

      ## getting mcl results into the ClusterSet object
      mcl_cluster <- readLines(mcl_out_file)
      gene_list <- NULL
      clusters <- NULL
      size <- NULL
      nb <- 0
      nb_cluster_deleted <- 0

      for (i in 1:length(mcl_cluster)) {
        h <- unlist(strsplit(mcl_cluster[i], "\t"))
        cur_clust <- data_matrix[h,]
        cur_clust[cur_clust > 0 ] <- 1
        cur_dot_prod <- cur_clust %*% t(cur_clust)

        if(mean(cur_dot_prod) > av_dot_prod_min & length(h) > min_cluster_size){

          nb <- nb + 1
          gene_list <- c(gene_list, h)
          clusters <- c(clusters, rep(nb, length(h)))
          if (is.null(size)) {
            size <- length(h)
          }
          else {
            size <- c(size, length(h))
          }
        }else{
          nb_cluster_deleted <- nb_cluster_deleted + 1
        }
      }
      print_msg(paste(nb, " clusters conserved after MCL partitioning."),
                msg_type="INFO")
      print_msg(paste(nb_cluster_deleted,
                      " clusters filtered out from MCL partitioning (size and mean dot product)."),
                msg_type="INFO")


      ## build ClusterSet object
      if (nb > 0) {
        obj@name <- name
        obj@data <- as.matrix(data_matrix[gene_list, ])
        obj@cluster <- clusters
        obj@size <- size

        centers <- matrix(ncol = ncol(data_matrix), nrow = nb)
        ## calcul of the mean profils
        for (i in 1:nb) {
          centers[i, ] <- apply(obj@data[obj@cluster == i, ],
            2, mean,
            na.rm = TRUE
          )
        }
        obj@center <- centers

        ## add DBFMCL parameters used to build this object
        obj@parameters <- list(
          distance_method = distance_method,
          k = k,
          random = random,
          fdr = fdr,
          set.seed = set.seed,
          inflation = inflation
        )
      }
  }
  else {
    stop("\t--> There is no conserved gene.\n\n")
  }
  return(obj)

}

###############################################################
##    COMPUTE DBF algorithm
###############################################################

#' @title
#' DBF
#' @description
#' This function is an internal function used by \code{\link{DBFMCL}} to detect
#' informative elements (\emph{i.e.}, those that belong to dense regions). User
#' should not use this function. Instead they can use the \code{\link{DBFMCL}}
#' function with \code{clustering} argument set to \code{FALSE}.
#'
#' See \code{\link{DBFMCL}}
#'
#' @param data a matrix or data.frame
#' @param name a prefix for the file name
#' @param distance_method a method to compute the distance to the k-th nearest
#' neighbor. One of "pearson" (Pearson's correlation coefficient-based
#' distance), "spearman" (Spearman's rho-based distance) or "euclidean".
#' @param silent if set to TRUE (default), the progression of distance matrix
#' calculation is not displayed.
#' @param k the neighborhood size.
#' @param random the number of simulated distributions S to compute. By default
#' \code{random = FALSE}.
#' @param fdr a value for the false discovery rate.
#' @param memory_used size of the memory used to store part of the distance
#' matrix. The subsequent sub-matrix is used to computed simulated distances to
#' the k-th nearest neighbor (see detail section).
#' @param set.seed specify seeds for random number generator.
#' @section Warnings: Works only on UNIX-alikes platforms.
#' @author Bergon A., Bavais J., Textoris J., Granjeaud S., Lopez F and Puthier
#' D.
#' @seealso \code{\link{DBFMCL}}
#' @references Lopez F.,Textoris J., Bergon A., Didier G., Remy E., Granjeaud
#' S., Imbert J. , Nguyen C. and Puthier D. TranscriptomeBrowser: a powerful
#' and flexible toolbox to explore productively the transcriptional landscape
#' of the Gene Expression Omnibus database. PLoSONE, 2008;3(12):e4001.
#' @keywords manip
#' @export DBF
DBF <- function(data, name = NULL,
                distance_method = c("spearman", "pearson", "euclidean"),
                silent = FALSE,
                k = 100,
                random = 3,
                fdr = 10,
                memory_used = 1024,
                set.seed = 123) {

  ## testing the system
  if (.Platform$OS.type != "windows") {
    if (!is.null(data)) {
      ## getting data and parameters
      if (is.null(name)) name <- "exprs"
      data <- get_data_4_DBFMCL(data = data)$data
      row <- rownames(data)
      col <- colnames(data)
      distance_method <- match.arg(distance_method)

      ## transforming data into double
      data <- apply(data, 2, as.double)

      if (silent) {
        print_msg(paste0("Computing distances to the kth-nearest neighbors ",
                         "and associated FDR values... \n"),
                  msg_type = "INFO")

      }
      outfile <- paste(name, ".dbf_out.txt", sep = "")

      ## launching DBF
      a <- .C("DBF",
        data,
        as.integer(nrow(data)),
        as.integer(ncol(data)),
        row,
        col,
        distance_method,
        as.integer(k),
        as.integer(random),
        as.integer(!silent),
        as.integer(memory_used),
        as.integer(fdr),
        as.integer(!silent),
        m2 = vector(length = nrow(data), mode = "character"),
        outfile,
        as.integer(set.seed)
      )

      ## creation of the ClusterSet object
      obj <- new("ClusterSet")
      obj@algorithm <- "DBFMCL"

      informative <- a$m2[a$m2 != ""]
      if (length(informative) > 0) {
        obj@data <- as.matrix(data[a[[4]] %in% informative, ])
        obj@cluster <- rep(1, nrow(obj@data))
        obj@size <- nrow(obj@data)
        obj@center <- matrix(
          apply(obj@data[obj@cluster == 1, ],
            2,
            mean,
            na.rm = TRUE
          ),
          nrow = 1
        )
      }
      return(obj)
    } else {
      stop("\t--> Please provide a matrix...\n\n")
    }
  }
  else {
    stop("\t--> A unix-like OS is required to launch mcl and cluster programs.")
  }
}

##############################################################
##    MCL
##############################################################


#' @title
#' Invokes the command line version of Markov CLustering algorithm (MCL).
#' @description
#'  This function invokes the mcl system command. MCL is a clustering algorithm
#' for graphs that was developped by Stijn van Dongen (see references for
#' further informations).
#' @param name a character string corresponding to the file name.
#' @param inflation the main control of MCL. Inflation affects cluster
#' granularity. It is usually chosen somewhere in the range \code{[1.2-5.0]}.
#' \code{inflation = 5.0} will tend to result in fine-grained clusterings, and
#' whereas \code{inflation = 1.2} will tend to result in very coarse grained
#' clusterings. By default, \code{inflation = 2.0}. Default setting gives very
#' good results for microarray data when k is set around 100.
#' @param silent if set to TRUE, the progression of the MCL partitionning is
#' not displayed.
#' @param threads The number of threads to use.
#' @return Returns a file with the ".mcl\_out.txt" extension.
#' @section warning: Works only on UNIX-like plateforms. MCL should be
#' installed. The following command lines can be used for installation.
#'
#' \code{# Download the latest version of mcl (RTools4TB has been tested
#' successfully with the 06-058 version).}
#' \code{wget http://micans.org/mcl/src/mcl-latest.tar.gz}
#' \code{# Uncompress and install mcl}
#' \code{tar xvfz mcl-latest.tar.gz}
#' \code{cd mcl-xx-xxx}
#' \code{./configure}
#' \code{make}
#' \code{sudo make install}
#' \code{# You should get mcl in your path}
#' \code{mcl -h}
#' @author Bergon A., Lopez F., Textoris J., Granjeaud S. and Puthier D.
#' @references Stijn van Dongen. A cluster algorithm for graphs.  Technical
#' Report INS-R0010, National Research Institute for Mathematics and Computer
#' Science in the Netherlands, Amsterdam, May 2000.
#' \url{http://www.cwi.nl/ftp/CWIreports/INS/INS-R0010.ps.Z}
#' @keywords manip
#' @export mcl_system_cmd
mcl_system_cmd <- function(name, inflation = 2.0, silent = FALSE, threads=1) {
  ## testing the system
  if (.Platform$OS.type != "windows") {

    ## Testing mcl installation
    if (system("mcl --version | grep 'Stijn van Dongen'", intern = TRUE) > 0) {
      if (!silent) {
        cat("Running mcl (graph partitioning)... ")
        verb <- ""
      }
      else {
        verb <- "-V all "
      }
      if (inflation != 2) {
        i <- paste("-I ", as.character(round(inflation, 1)), sep = "")
      } else {
        i <- "-I 2.0"
      }
      threads <- paste("-te", threads, sep = " ")
      ## launching mcl program
      cmd <- paste0("mcl ",
                   name, ".dbf_out.txt ",
                   i,
                   " --abc -o ",
                   name, ".mcl_out.txt ",
                   verb,
                   threads)

      system(cmd)
      if (!silent) {
        print_msg("Done", msg_type="INFO")
        print_msg(paste0("creating file : ",
                        file.path(getwd(), paste(name, ".mcl_out.txt", sep = ""))),
                  msg_type="INFO")
      }
    } else {
      stop(
        "\t--> Please install mcl on your computer...\n",
        "\t--> You can download it from : 'http://www.micans.org/mcl/'\n\n"
      )
    }
  }
  else {
    stop("--> A unix-like OS is required to launch mcl and cluster programs.")
  }
}


#################################################################
##    Getting data
#################################################################

#' @title
#' Fetch an expression matrix from a file, dataframe or Seurat object.
##' @description
#' This function fetchs an expression matrix from a file, dataframe or Seurat object.
#' @param data A \code{matrix}, \code{data.frame} or \code{Seurat} object.
#' @param filename A character string representing the file name.
#' @param path A character string representing the data directory where
#' intermediary files are to be stored. Default to current working directory.
#'
#' @return A list containing a matrix and the filename (if filename argument is used).
#' @export get_data_4_DBFMCL
#'
#' @examples
#' 
#' \dontrun{
#' ## with an artificial dataset
#'
#' m <- matrix(rnorm(80000), nc = 20)
#' res <- get_data_4_DBFMCL(data=m)
#' }
#' 
get_data_4_DBFMCL <- function(data = NULL, filename = NULL, path = ".") {

  ## getting matrix (probesID vs SamplesID)
  if (!is.null(data)) {
    if (inherits(data, "Seurat")) {
      data <- as.matrix(data@assays$RNA@data)
    }
    else if (is.data.frame(data)) {
      data <- as.matrix(data)
    }
    if (!is.matrix(data)) {
      stop(
        "\t--> Please provide a Seurat Object, a data.frame",
        " or a matrix.\n"
      )
    }
    name <- NULL
  }
  else {
    if (!is.null(filename)) {
      data <- as.matrix(read.table(file.path(path, filename),
        sep = "\t", header = TRUE, row.names = 1, quote = ""
      ))
      name <- unlist(strsplit(filename, "\\."))[1]
    }
    else {
      stop(
        "\t--> Please provide an ExpressionSet, a data.frame, ",
        "a matrix or a tabular file\n"
      )
    }
  }
  ## adding dimnames if not provided
  if (is.null(rownames(data))) {
    print_msg("Row names not provided. Adding.", msg_type = "DEBUG")
    rownames(data) <- paste("gene", 1:nrow(data), sep = "")
  }
  if (is.null(colnames(data))) {
    print_msg("Colum names not provided. Adding.", msg_type = "DEBUG")
    colnames(data) <- paste("sample", 1:ncol(data), sep = "")
  }

  return(list(data = data, name = name))
}


#########################################################
##      END PACKAGE DBFMCL
#########################################################
