################################################################################
#' Select genes in co-expression patterns
#'
#' This function selects genes based on their correlation with other genes. 
#' It takes as input a normalized gene expression matrix and 
#' a set of parameters that affect the selection of genes.
#'
#' @param data a matrix, data.frame or Seurat object.
#' @param dist_threads an integer specifying the number of threads for the kendall correlation.
#' @param distance_method a character string indicating the method for computing distances (one of "pearson", "cosine", "euclidean" or "kendall").
#' @param highest During the process, genes will be ordered by their distance to their k nearest neighbors (dknn). 
#' This parameter controls the fraction of genes with high dknn (ie. noise) whose neighborhood (i.e associated distances) 
#' will be used to compute simulated values. 0 would mean to use all the genes.
#' @param k an integer specifying the k nearest neighbours to compute.
#' @param row_sum an integer specifying the minimum number of cells expressing a gene. 
#' Genes expresssed in less than the chosen values are filtered out from the analysis.To keep all genes, use -Inf. 
#' @param fdr a numeric value indicating the false discovery rate threshold (range: 0 to 100).
#' @param which_slot a character string indicating which slot to use from the input scRNA-seq object (one of "data", "sct" or "counts"). 
#' SCT is the recommended method from Seurat package when working with spatial transcriptomics data.
#' @param no_dknn_filter a logical indicating whether to skip the k-nearest-neighbors (KNN) filter. If FALSE, all genes are kept for the next steps.
#' @param seed an integer specifying the random seed to use.
#'
#' @return a ClusterSet class object
#' 
#' @author Julie Bavais, Sebastien Nin, Aurelie Bergon, Fabrice Lopez, Julien Textoris, Samuel Granjeaud, Lionel Spinelli and Denis Puthier
#' 
#' @references
#' - Van Dongen S. (2000) A cluster algorithm for graphs. National
#' Research Institute for Mathematics and Computer Science in the 1386-3681.
#' - Lopez F.,Textoris J., Bergon A., Didier G., Remy E., Granjeaud
#' S., Imbert J. , Nguyen C. and Puthier D. TranscriptomeBrowser: a powerful
#' and flexible toolbox to explore productively the transcriptional landscape
#' of the Gene Expression Omnibus database. PLoSONE, 2008;3(12):e4001.
#'
#' @examples
#' # Set verbosity to 1 to display info messages only.
#' set_verbosity(1)
#' 
#' # Create a matrix with 4 signatures
#' m <- create_4_rnd_clust()
#' 
#' # Select informative genes
#' res <- select_genes(m,
#'                     distance = "kendall",
#'                     k = 75,
#'                     highest = 0.3,
#'                     fdr = 1e-8,
#'                     row_sum = -Inf)
#' 
#' # Display selected genes
#' res@gene_clusters$`1`
#'
#' @export select_genes

select_genes <- function(data = NULL,
                         dist_threads = 1,
                         distance_method = c("pearson",
                                             "cosine",
                                             "euclidean",
                                             "kendall"),
                         highest = 0.00005,
                         k = 80,
                         no_anti_cor=F,
                         row_sum = 1,
                         fdr = 0.005,
                         which_slot = c("data", "sct", "counts"),
                         no_dknn_filter = FALSE,
                         seed = 123) {
  
  # ======================
  ## Set parameters
  # ======================
  
  ## set a seed for reproductibility
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Get normalized gene expression matrix
  which_slot <- match.arg(which_slot)
  data <- get_data_for_scigenex(data = data, which_slot=which_slot)
  
  # Check if distance provided match with one used by DBF()
  if (length(distance_method) > 1) {
    distance_method <- distance_method[1]
  }
  distance_method <- match.arg(distance_method, c("pearson",
                                                  "cosine",
                                                  "euclidean",
                                                  "kendall"))
  
  # Check if highest is between 0 and 1
  if (highest < 0 | highest > 1) {
    print_msg("highest argument should be >= 0 and <= 1.",
              msg_type = "STOP")              
  }
  
  
  
  # ======================
  #### Compute distance matrix ####
  # ======================
  
  #################### Correlation and distance matrices
  # Remove genes with 0 values for all cells
  
  if(!inherits(data, "dgCMatrix")){
    genes_to_keep <- rowSums(data) > row_sum
  }else{
    genes_to_keep <- Matrix::rowSums(data) > row_sum
  }
  
  select_for_correlation <- data[genes_to_keep, ]
  
  # Compute gene-gene correlation/distance matrix
  
  print_msg(
    paste0(
      "Computing distances using selected method: ",
      distance_method
    ),
    msg_type = "INFO"
  )
  print_msg(paste0("Number of selected rows/genes (row_sum): ", nrow(select_for_correlation)), 
            msg_type = "DEBUG")
  
  # Compute correlations and corresponding 
  # distance matrix. Note that for pearson
  # and cosine, 0 < distance < 2.
  if (distance_method == "pearson") {
    dist_matrix <- qlcMatrix::corSparse(t(select_for_correlation))
    if(no_anti_cor)
      dist_matrix[dist_matrix < 0] <- 0
    dist_matrix <- 1 - dist_matrix
  }else  if (distance_method == "kendall") {
    dist_matrix <- as.matrix(amap::Dist(select_for_correlation, method="kendall", nbproc=dist_threads))
    if(no_anti_cor)
      dist_matrix[dist_matrix < 0] <- 0
  } else if (distance_method == "cosine") {
    dist_matrix <- as.matrix(qlcMatrix::cosSparse(t(select_for_correlation)))
    if(no_anti_cor)
      dist_matrix[dist_matrix < 0] <- 0
    dist_matrix <- 1 - dist_matrix
  } else if (distance_method == "euclidean") {
    dist_matrix <- as.matrix(dist(select_for_correlation))
  }
  
  print_stat("Distance matrix stats", 
             data = dist_matrix, msg_type = "DEBUG")
  
  # Compute min and max distances.
  # They will be used later to compute weights.
  min_dist <- min(dist_matrix)
  max_dist <- max(dist_matrix)
  
  # Set the rownames / colnames of the distance matrix
  rownames(dist_matrix) <- rownames(select_for_correlation)
  colnames(dist_matrix) <- rownames(select_for_correlation)
  
  # The distance from a gene to itself is 'hidden'
  diag(dist_matrix) <- NA
  
  
  
  # ======================
  #### Compute distance to KNN ####
  # ======================
  
  print_msg(paste0("Computing distances to KNN."), msg_type = "INFO")
  
  # Create a dataframe to store the DKNN values.
  # Gene_id appear both as rownames and column
  # for coding convenience
  df_dknn <- data.frame(
    dknn_values = rep(NA, nrow(dist_matrix)),
    row.names = rownames(dist_matrix),
    gene_id = rownames(dist_matrix)
  )
  
  # This list will be used to store the
  # dknn values
  l_knn <- list()
  
  # Extract the DKNN for each gene
  for (pos in seq_len(nrow(dist_matrix))) {
    gene <- rownames(df_dknn)[pos]
    gene_dist <- dist_matrix[gene, ]
    # Reorder the distance values in increasing order.
    # The distance from a gene to itself (previously set to NA)
    # is placed at the end of the vector.
    gene_dist <- sort(gene_dist, na.last = T)[1:k]
    
    # Add the neigbhors to the list
    l_knn[[gene]] <- gene_dist
    
    # Select the kth distance. 
    # This value corresponds to the DKNN of gene_i
    df_dknn[gene, "dknn_values"] <- gene_dist[k]
  }
  
  print_stat("Observed DKNN stats", 
             data = df_dknn$dknn_values, msg_type = "DEBUG")
  
  
  # ======================
  #### Compute simulated distance to KNN ####
  # ======================
  
  sim_dknn <- vector()
  critical_distance <- vector()
  
  if(!no_dknn_filter){
    
    #################### DKNN simulation
    print_msg(paste0("Computing simulated distances to KNN."), msg_type = "INFO")
    
    nb_selected_genes <- nrow(select_for_correlation)
    
    if (highest == 0) {
      tresh_dknn <- min(df_dknn$dknn_values)
    } else if (highest == 1) {
      tresh_dknn <- max(df_dknn$dknn_values)
    } else {
      tresh_dknn <- stats::quantile(df_dknn$dknn_values, highest)
    }
    
    gene_with_low_dknn <- df_dknn$gene_id[df_dknn$dknn_values > tresh_dknn]
    dist_values_sub <- dist_matrix[gene_with_low_dknn, ]
    dist_values_sub <- dist_values_sub[!is.na(dist_values_sub)]
    
    
    for (sim_nb in 1:nb_selected_genes) {
      # Randomly sample distances for one simulated gene
      dist_sim <- sample(dist_values_sub,
                         size = nb_selected_genes,
                         replace = FALSE)
      
      # Extract the k nearest neighbors of these simulated gene
      dist_sim <- sort(dist_sim)
      sim_dknn[sim_nb] <- dist_sim[k]
    }
    
    print_stat("Simulated DKNN stats", 
               data = sim_dknn, msg_type = "DEBUG")
    
    # The simulated DKNN values follow a normal distribution.
    # Compute the parameters of this distribution
    mean_sim <- mean(sim_dknn)
    sd_sim <- sd(sim_dknn)
    
    
    # ======================
    #### Compute the DKNN threshold (or critical distance) ####
    # ======================
    # Order genes by DKNN values
    print_msg(paste0("Computing threshold of distances to KNN (DKNN threshold)."),
              msg_type = "INFO")
    
    # Order the dknn values from low to high
    df_dknn <- df_dknn[order(df_dknn$dknn_values), ]
    df_dknn[, "FDR"] <- NA
    
    # Compute the FDR
    for (i in 1:nb_selected_genes) {
      gene <- df_dknn$gene_id[i]
      df_dknn[gene, "FDR"] <-
        stats::pnorm(df_dknn[gene, "dknn_values"],
                     mean = mean_sim,
                     sd = sd_sim,
                     lower.tail = T) / (i / nb_selected_genes) * 100
    }
    
    df_dknn$FDR[df_dknn$FDR > 100] <- 100
    

    # ======================
    #### Select genes with a distance value under critical distance ####
    # ======================
    print_msg(paste0("Selecting informative genes."), msg_type = "INFO")
    fdr_tresh_pos <- which(df_dknn$FDR > fdr)[1]
    selected_genes <- df_dknn[1:fdr_tresh_pos,]$gene_id
    critical_distance <-df_dknn$dknn_values[fdr_tresh_pos]
    
  }else{
    selected_genes <- rownames(dist_matrix)
  }
  

  # ======================
    #### Create the ClusterSet object ####
  # ======================
  print_msg(paste0("Create the ClusterSet object."), msg_type = "INFO")
  
  obj <- new("ClusterSet")
  
  if (length(selected_genes) > 0) {
    obj@data <- as.matrix(data[selected_genes, ])
    obj@gene_clusters <- list("1" = rownames(obj@data))
    obj@gene_clusters_metadata <- list("cluster_id" = as.numeric(names(obj@gene_clusters)),
                                       "number" = max(names(obj@gene_clusters)),
                                       "size" = nrow(obj@data))
    
    obs_dknn <- as.vector(df_dknn[, "dknn_values"])
    names(obs_dknn) <- df_dknn[, "gene_id"]
    if(length(selected_genes) > 1){
      center = matrix(apply(obj@data[obj@gene_clusters$`1`, ],
                            2,
                            mean,
                            na.rm = TRUE),
                      nrow = 1)
    }else{
      center = matrix(obj@data[obj@gene_clusters$`1`, ],
                      nrow = 1)
    }
    
    obj@dbf_output <- list("dknn" = obs_dknn,
                           "simulated_dknn" = sim_dknn,
                           "critical_distance" = critical_distance,
                           "fdr" = df_dknn$FDR,
                           "center" = center,
                           "all_gene_expression_matrix" = data,
                           "all_neighbor_distances" = l_knn[selected_genes])
    # obj@normal_model_mean <- mean_sim
    # obj@normal_model_sd <- mean_sd
  }else{
    
    obs_dknn <- as.vector(df_dknn[, "dknn_values"])
    names(obs_dknn) <- df_dknn[, "gene_id"]
    obj@dbf_output <- list("dknn" = obs_dknn,
                           "simulated_dknn" = sim_dknn,
                           "critical_distance" = critical_distance,
                           "fdr" = df_dknn$FDR,
                           "center" = NULL,
                           "all_gene_expression_matrix" = data,
                           "all_neighbor_distances" = NULL)
  }
  
  obj@parameters <- list("distance_method" = distance_method,
                         "k" = k,
                         "highest" = highest,
                         "fdr" = fdr,
                         "row_sum" = row_sum,
                         "no_dknn_filter" = no_dknn_filter,
                         "seed" = seed)
  
  return(obj)
}
