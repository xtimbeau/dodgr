#include "dgraph.h"
#include "dijkstra.h"
#include "heaps/heap_lib.h"

#include <algorithm> // std::fill, std::reverse

#include <Rcpp.h>

const float INFINITE_FLOAT =  std::numeric_limits<float>::max ();
const float INFINITE_INT =  std::numeric_limits<int>::max ();

template <typename T>
void inst_graph (DGraph *g, unsigned int nedges,
        std::map <std::string, unsigned int> &vert_map,
        std::vector <std::string> &from,
        std::vector <std::string> &to,
        std::vector <T> &dist,
        std::vector <T> &wt)
{
    for (unsigned int i = 0; i < nedges; ++i)
    {
        unsigned int fromi = vert_map [from [i]];
        unsigned int toi = vert_map [to [i]];
        g->addNewEdge (fromi, toi, dist [i], wt [i]);
    }
}
template void inst_graph <float> (DGraph *g, unsigned int nedges,
        std::map <std::string, unsigned int> &vert_map,
        std::vector <std::string> &from,
        std::vector <std::string> &to,
        std::vector <float> &dist,
        std::vector <float> &wt);

Dijkstra * dijkstra_bheap (unsigned int nverts)
{
    HeapD<BHeap> heapD;
    Dijkstra *dijkstra = new Dijkstra (nverts, &heapD);
    return dijkstra;
}

Dijkstra * dijkstra_fheap (unsigned int nverts)
{
    HeapD<FHeap> heapD;
    Dijkstra *dijkstra = new Dijkstra (nverts, &heapD);
    return dijkstra;
}

Dijkstra * dijkstra_heap23 (unsigned int nverts)
{
    HeapD<Heap23> heapD;
    Dijkstra *dijkstra = new Dijkstra (nverts, &heapD);
    return dijkstra;
}

Dijkstra * dijkstra_triheap (unsigned int nverts)
{
    HeapD<TriHeap> heapD;
    Dijkstra *dijkstra = new Dijkstra (nverts, &heapD);
    return dijkstra;
}

Dijkstra * dijkstra_triheapext (unsigned int nverts)
{
    HeapD<TriHeapExt> heapD;
    Dijkstra *dijkstra = new Dijkstra (nverts, &heapD);
    return dijkstra;
}

Dijkstra * dijkstra_radix (unsigned int nverts)
{
    HeapD<RadixHeap> heapD;
    Dijkstra *dijkstra = new Dijkstra (nverts, &heapD);
    return dijkstra;
}

//' rcpp_get_sp_dists
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_get_sp_dists (Rcpp::DataFrame graph,
        Rcpp::DataFrame vert_map_in,
        std::vector <int> fromi,
        std::vector <int> toi,
        std::string heap_type)
{
    if (fromi [0] < 0) // use all vertices
    {
        Rcpp::NumericVector id_vec = vert_map_in ["id"];
        fromi = Rcpp::as <std::vector <int>> (id_vec);
    }
    if (toi [0] < 0) // use all vertices
    {
        Rcpp::NumericVector id_vec = vert_map_in ["id"];
        toi = Rcpp::as <std::vector <int>> (id_vec);
    }
    unsigned int nfrom = fromi.size (), nto = toi.size ();

    std::vector <std::string> from = graph ["from"];
    std::vector <std::string> to = graph ["to"];
    std::vector <float> dist = graph ["d"];
    std::vector <float> wt = graph ["w"];

    unsigned int nedges = graph.nrow ();
    std::map <std::string, unsigned int> vert_map;
    std::vector <std::string> vert_map_id = vert_map_in ["vert"];
    std::vector <unsigned int> vert_map_n = vert_map_in ["id"];
    for (int i = 0; i < vert_map_in.nrow (); ++i)
    {
        vert_map.emplace (vert_map_id [i], vert_map_n [i]);
    }
    unsigned int nverts = vert_map.size ();

    DGraph *g = new DGraph (nverts);
    inst_graph (g, nedges, vert_map, from, to, dist, wt);

    Dijkstra *dijkstra;

    if (heap_type == "FHeap")
        dijkstra = dijkstra_fheap (nverts);
    else if (heap_type == "BHeap")
        dijkstra = dijkstra_bheap (nverts);
    else if (heap_type == "Heap23")
        dijkstra = dijkstra_heap23 (nverts);
    else if (heap_type == "TriHeap")
        dijkstra = dijkstra_triheap (nverts);
    else if (heap_type == "TriHeapExt")
        dijkstra = dijkstra_triheapext (nverts);
    else if (heap_type == "Radix")
        dijkstra = dijkstra_radix (nverts);

    float* w = new float [nverts];
    float* d = new float [nverts];
    int* prev = new int [nverts];

    dijkstra->init (g); // specify the graph

    // initialise dout matrix to NA
    Rcpp::NumericVector na_vec = Rcpp::NumericVector (nfrom * nto,
            Rcpp::NumericVector::get_na ());
    Rcpp::NumericMatrix dout (nfrom, nto, na_vec.begin ());
    for (unsigned int v = 0; v < nfrom; v++)
    {
        std::fill (w, w + nverts, INFINITE_FLOAT);
        std::fill (d, d + nverts, INFINITE_FLOAT);

        dijkstra->run (d, w, prev, fromi [v]);
        for (unsigned int vi = 0; vi < nto; vi++)
            if (w [toi [vi]] < INFINITE_FLOAT)
                dout (v, vi) = d [toi [vi]];
    }

    delete [] d;
    delete [] w;
    delete [] prev;

    delete dijkstra;
    delete g;

    return (dout);
}


//' rcpp_get_paths
//'
//' @param graph The data.frame holding the graph edges
//' @param vert_map_in map from <std::string> vertex ID to (0-indexed) integer
//' index of vertices
//' @param fromi Index into vert_map_in of vertex numbers
//' @param toi Index into vert_map_in of vertex numbers
//'
//' @note The graph is constructed with 0-indexed vertex numbers contained in
//' code{vert_map_in}. Both \code{fromi} and \code{toi} already map directly
//' onto these. The graph has to be constructed by first constructing a
//' \code{std::map} object (\code{vertmap}) for \code{vert_map_in}, then
//' translating all \code{graph["from"/"to"]} values into these indices. This
//' construction is done in \code{inst_graph}.
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::List rcpp_get_paths (Rcpp::DataFrame graph,
        Rcpp::DataFrame vert_map_in,
        std::vector <int> fromi,
        std::vector <int> toi,
        std::string heap_type)
{
    if (fromi [0] < 0) // use all vertices
    {
        Rcpp::NumericVector id_vec = vert_map_in ["id"];
        fromi = Rcpp::as <std::vector <int>> (id_vec);
    }
    if (toi [0] < 0) // use all vertices
    {
        Rcpp::NumericVector id_vec = vert_map_in ["id"];
        toi = Rcpp::as <std::vector <int>> (id_vec);
    }
    unsigned int nfrom = fromi.size (), nto = toi.size ();

    std::vector <std::string> from = graph ["from"];
    std::vector <std::string> to = graph ["to"];
    std::vector <float> dist = graph ["d"];
    std::vector <float> wt = graph ["w"];

    unsigned int nedges = graph.nrow ();
    std::map <std::string, unsigned int> vert_map;
    std::vector <std::string> vert_map_id = vert_map_in ["vert"];
    std::vector <unsigned int> vert_map_n = vert_map_in ["id"];
    for (int i = 0; i < vert_map_in.nrow (); ++i)
    {
        vert_map.emplace (vert_map_id [i], vert_map_n [i]);
    }
    unsigned int nverts = vert_map.size ();

    DGraph *g = new DGraph (nverts);
    inst_graph (g, nedges, vert_map, from, to, dist, wt);

    Dijkstra *dijkstra;

    if (heap_type == "FHeap")
        dijkstra = dijkstra_fheap (nverts);
    else if (heap_type == "BHeap")
        dijkstra = dijkstra_bheap (nverts);
    else if (heap_type == "Heap23")
        dijkstra = dijkstra_heap23 (nverts);
    else if (heap_type == "TriHeap")
        dijkstra = dijkstra_triheap (nverts);
    else if (heap_type == "TriHeapExt")
        dijkstra = dijkstra_triheapext (nverts);
    else if (heap_type == "Radix")
        dijkstra = dijkstra_radix (nverts);

    dijkstra->init (g); // specify the graph

    Rcpp::List res (nfrom);
    float* w = new float [nverts];
    float* d = new float [nverts];
    int* prev = new int [nverts];

    for (unsigned int v = 0; v < nfrom; v++)
    {
        std::fill (w, w + nverts, INFINITE_FLOAT);
        std::fill (d, d + nverts, INFINITE_FLOAT);

        dijkstra->run (d, w, prev, fromi [v]);

        Rcpp::List res1 (nto);
        for (unsigned int vi = 0; vi < nto; vi++)
        {
            std::vector <unsigned int> onePath;
            if (w [toi [vi]] < INFINITE_FLOAT)
            {
                unsigned int target = toi [vi];
                while (target < INFINITE_INT)
                {
                    // Note that targets are all C++ 0-indexed and are converted
                    // directly here to R-style 1-indexes.
                    onePath.push_back (target + 1);
                    target = prev [target];
                }
            }
            std::reverse (onePath.begin (), onePath.end ());
            res1 [vi] = onePath;
        }
        res [v] = res1;
    }


    delete [] d;
    delete [] w;
    delete [] prev;

    delete dijkstra;
    delete g;

    return (res);
}
