#' Calculate matrix of pair-wise distances between points.
#' Hacked par moi pour réduire les préchargements. Un graph "pré processé" est calculé
#' et utilisé à la place du graph
#' Ca réduit pas mal les coûts fixes.
#' 
#'
#' @param proc_g object representing the network preprocessed par process_graph
#' graph (see Notes)
#' @param from Vector or matrix of points **from** which route distances are to
#' be calculated (see Notes)
#' @param to Vector or matrix of points **to** which route distances are to be
#' calculated (see Notes)
#' @param shortest If `FALSE`, calculate distances along the \emph{fastest}
#' rather than shortest routes (see Notes).
#' @param pairwise If `TRUE`, calculate distances only between the ordered
#' pairs of `from` and `to`.
#' @param heap Type of heap to use in priority queue. Options include
#' Fibonacci Heap (default; `FHeap`), Binary Heap (`BHeap`),
#' `Trinomial Heap (`TriHeap`), Extended Trinomial Heap
#' (`TriHeapExt`, and 2-3 Heap (`Heap23`).
#' @param parallel If `TRUE`, perform routing calculation in parallel (see
#' details)
#' @param quiet If `FALSE`, display progress messages on screen.
#' @param to_from_indices (NULL) si non nul, remplace to/from et prend le racourci
#' @return square matrix of distances between nodes
#'
#'
#' @family distances
#' @export

dodgr_dists_pre <- function (proc_g,
                             from = NULL,
                             to = NULL,
                             shortest = TRUE,
                             pairwise = FALSE,
                             parallel = TRUE,
                             quiet = TRUE,
                             to_from_indices = NULL) {
  # tictoc::tic()
  if(is.null(to_from_indices))
    to_from_indices <- to_from_index_with_tp_pre(proc_g, from, to)
  
  if (proc_g$compound) {
    graph <- proc_g$graph_compound
    vert_map <- proc_g$vert_map_compound
  } else {
    graph <- proc_g$graph
    vert_map <- proc_g$vert_map
  }
  gr_cols <- proc_g$gr_cols
  if (!shortest) {
    if (is.na (gr_cols$time_weighted)) {
      stop (
        "Graph does not contain a weighted time column from ",
        "which to calculate fastest paths."
      )
    }
    graph [[gr_cols$d_weighted]] <- graph [[gr_cols$time_weighted]]
  }
  
  graph <- convert_graph (graph, gr_cols)
  
  if (!quiet) {
    message ("Calculating shortest paths ... ", appendLF = FALSE)
  }
  
  if (parallel && proc_g$heap == "TriHeapExt") {
    if (!quiet) {
      message (
        "Extended TriHeaps can not be calculated in parallel; ",
        "reverting to serial calculation"
      )
    }
    parallel <- FALSE
  }
  # tictoc::toc()
  
  d <- calculate_distmat (
    graph,
    vert_map,
    to_from_indices$from,
    to_from_indices$to,
    proc_g$heap,
    proc_g$is_spatial,
    parallel,
    pairwise
  )
  
  if (!quiet) {
    message ("done.")
  }
  
  return (d)
}

#' Get lists of 'from' and 'to' indices, potentially mapped on to compound
#' junctions for graphs with turn penalties.
#'
#' This function calls the following `fill_to_from_index` to generate the actual
#' values. For graphs with turn penalties, it also returns the modified version
#' of the graph including compound junctions.
#' @noRd
to_from_index_with_tp_pre <- function (fgraph, from, to) {
  vert_map <- fgraph$vert_map
  from_index <-
    fill_to_from_index_pre (fgraph$vertices, vert_map, fgraph$gr_cols, from, from = TRUE)
  to_index <-
    fill_to_from_index_pre (fgraph$vertices, vert_map, fgraph$gr_cols, to, from = FALSE)
  
  if (fgraph$compound) {
    vert_map <- fgraph$vert_map_compound
    from_index <- remap_tf_index_for_tp (from_index, vert_map, from = TRUE)
    to_index <- remap_tf_index_for_tp (to_index, vert_map, from = FALSE)
  }
  
  return (list (from = from_index, to = to_index))
}


#' Get lists of 'from' and 'to' indices, potentially mapped on to compound
#' junctions for graphs with turn penalties.
#'
#' This function calls the following `fill_to_from_index` to generate the actual
#' values. For graphs with turn penalties, it also returns the modified version
#' of the graph including compound junctions.
#' @noRd
to_from_index_with_tp_pre <- function (fgraph, from, to) {
  vert_map <- fgraph$vert_map
  from_index <-
    fill_to_from_index_pre (fgraph$vertices, vert_map, fgraph$gr_cols, from, from = TRUE)
  to_index <-
    fill_to_from_index_pre (fgraph$vertices, vert_map, fgraph$gr_cols, to, from = FALSE)
  
  if (fgraph$compound) {
    vert_map <- fgraph$vert_map_compound
    from_index <- remap_tf_index_for_tp (from_index, vert_map, from = TRUE)
    to_index <- remap_tf_index_for_tp (to_index, vert_map, from = FALSE)
  }
  
  return (list (from = from_index, to = to_index))
}

#' fill_to_from_index
#'
#' @noRd
fill_to_from_index_pre <- function (vertices,
                                    vert_map,
                                    gr_cols,
                                    pts,
                                    from = TRUE) {
  
  id <- NULL
  if (is.null (pts)) {
    index <- seq_len (nrow (vert_map)) - 1L
    if (!is.null (vert_map$vert)) {
      id <- vert_map$vert
    }
  } else {
    index_id <- get_index_id_cols_pre (vertices, gr_cols, vert_map, pts)
    if (any (is.na (index_id$id))) {
      stop ("Unable to match all routing points to graph vertices")
    }
    index <- index_id$index - 1L # 0-based
    id <- index_id$id
  }
  
  list (index = index, id = id)
}

#' Get an index of `pts` matching `vert_map`.
#'
#' Also returns the corresonding names of those `pts`
#'
#' @return list of `index`, which is 0-based for C++, and corresponding
#' `id` values.
#' @noRd
get_index_id_cols_pre <- function (vertices,
                                   gr_cols,
                                   vert_map,
                                   pts) {
  
  index <- -1
  id <- NULL
  if (!missing (pts)) {
    if (is.integer (pts) && is.vector (pts)) {
      index <- pts
    } else if (is.character (pts) ||
               is.numeric (pts) ||
               is.matrix (pts) ||
               is.data.frame (pts)) {
      index <- get_pts_index_pre (vertices, gr_cols, vert_map, pts)
    } else {
      stop (
        "routing points are of unknown form; must be either ",
        "character, matrix, or integer"
      )
    }
    
    if (length (pts == 2) && is.numeric (pts) &&
        ((any (grepl ("x", names (pts), ignore.case = TRUE)) &&
          any (grepl ("y", names (pts), ignore.case = TRUE))) ||
         (any (grepl ("lon", names (pts), ignore.case = TRUE) &&
               (any (grepl ("lat", names (pts), ignore.case = TRUE))))))) {
      names (pts) <- NULL
    }
    
    id <- get_id_cols (pts)
    if (is.null (id)) {
      id <- vert_map$vert [index]
    } # index is 1-based
  }
  list (index = index, id = id)
}

#' Convert `from` or `to` args of \link{dodgr_dists} to indices into
#' \link{dodgr_vertices}.
#'
#' @param graph A dodgr graph
#' @param vert_map Two-column `data.frame` of unique vertices and
#' corresponding IDs, obtained from `make_vert_map`
#' @param gr_cols Returned from `dodgr_graph_cols()`
#' @param pts Either a vector of names, or a matrix or `data.frame` of
#' arbitrary geographical coordinates for which to get index into vertices of
#' graph.
#'
#' @noRd
get_pts_index_pre <- function (vertices,
                               gr_cols,
                               vert_map,
                               pts) {
  
  if (!(is.matrix (pts) || is.data.frame (pts))) {
    if (!is.numeric (pts)) {
      pts <- matrix (pts, ncol = 1)
    } else {
      nms <- names (pts)
      if (is.null (nms)) {
        nms <- c ("x", "y")
      }
      pts <- matrix (pts, nrow = 1) # vector of (x,y) vals
      colnames (pts) <- nms
    }
  }
  
  if (ncol (pts) == 1) {
    
    pts <- get_pts_index_vec (pts, vert_map)
    
  } else {
    
    pts <- get_pts_index_rect_pre (pts, vertices, gr_cols)
  }
  
  pts
}

get_pts_index_rect_pre <- function (pts, vertices, gr_cols) {
  
  nms <- names (pts)
  if (is.null (nms)) {
    nms <- colnames (pts)
  }
  
  ix <- which (grepl ("x", nms, ignore.case = TRUE) |
                 grepl ("lon", nms, ignore.case = TRUE))
  iy <- which (grepl ("y", nms, ignore.case = TRUE) |
                 grepl ("lat", nms, ignore.case = TRUE))
  
  if (length (ix) != 1 || length (iy) != 1) {
    stop (paste0 (
      "Unable to determine geographical ",
      "coordinates in from/to"
    ))
  }
  
  index <- match (c ("xfr", "yfr", "xto", "yto"), names (gr_cols))
  if (any (is.na (gr_cols [index]))) {
    stop (paste0 (
      "Cannot determine geographical coordinates ",
      "against which to match pts"
    ))
  }
  
  if (is.data.frame (pts)) {
    names (pts) [ix] <- "x"
    names (pts) [iy] <- "y"
  } else {
    colnames (pts) [ix] <- "x"
    colnames (pts) [iy] <- "y"
  }
  
  # Result of rcpp_points_index is 0-indexed for C++
  pts <- rcpp_points_index_par (vertices, pts) + 1
  
  return (pts)
}

#' Préproccess le graph pour gagner du temps
#'
#' @param graph le network
#' @param heap le heap
#'
#' @return un graph préprocéssé
#' @export
#'
process_graph <- function(graph, heap = "BHeap") {
  
  graph <- tbl_to_df (graph)
  hps <- get_heap (heap, graph)
  heap <- hps$heap
  graph <- hps$graph
  graph <- preprocess_spatial_cols (graph)
  gr_cols <- dodgr_graph_cols (graph)
  is_spatial <- is_graph_spatial (graph)
  vertices <- dodgr_vertices(graph)
  vert_map <- make_vert_map (graph, gr_cols, is_spatial)
  compound <- (get_turn_penalty (graph) > 0.0)
  vert_map_compound <- graph_compound <- compound_junction_map <- NULL
  
  keep <- c(".vx0", ".vx1", "edge_", 
            ".vx0_x", ".vx0_y",
            ".vx1_x", ".vx1_y", 
            "d", "dz", "time",
            "d_weighted", "time_weighted", "component")
  
  if (compound) {
    if (methods::is (graph, "dodgr_contracted")) {
      warning (
        "graphs with turn penalties should be submitted in full, ",
        "not contracted form;\nsubmitting contracted graphs may ",
        "produce unexpected behaviour."
      )
    }
    res <- create_compound_junctions (graph)
    graph_compound <- res$graph
    compound_junction_map <- res$edge_map
    
    vert_map_compound <- make_vert_map (graph_compound, gr_cols, is_spatial)
    
    graph_compound <- graph_compound[keep]
  } else {
    graph_compound <- graph[keep]
  }
  
  gr_cols <- dodgr_graph_cols (graph_compound)
  
  return(list(
    vertices = vertices, 
    gr_cols = gr_cols,
    heap = heap,
    compound = compound,
    is_spatial = is_spatial,
    graph_compound = graph_compound, 
    d = graph_compound$d,
    compound_junction_map = compound_junction_map,
    vert_map = vert_map,
    vert_map_compound = vert_map_compound))
}
