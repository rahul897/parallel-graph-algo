# parallel-graph-algo

The parallel algorithm NSD (network similarity decomposition),
an optimized variant of Iso-Rank algorithm, is used to find the similarity matrix on which maximum bipartite matching algorithm
is used to map similar vertices across complex graphs having more than 10k vertices. 

FILES

"serial_both.cpp" is serial implementation of NSD algortihm to find the similarity matrix and serial implementation
of maximum bipartite algorithm, Use "input.txt" as input.

"parallel_nsd_serial_match.cpp" is parallel implementation of NSD algortihm to find the similarity matrix and serial implementation
of maximum bipartite algorithm

"v2_parallel_both.cpp" is parallel implementation of NSD algortihm to find the similarity matrix and parallel implementation
of maximum bipartite algorithm. Use "svdf.txt" as input and "temp_matrix.h" as the header.


REFERENCES

1.Kollias, Georgios; Sathe, Madan; Schenk, Olaf; and Grama, Ananth, "Fast Parallel Algorithms for Graph Similarity and Matching"

2.Madan Sathe: “Parallel Graph Algorithms for Finding Weighted Matchings and Subgraphs in Computational Science”

3.GiorgosKollias, ShahinMohammadi, and AnanthGrama: “Network Similarity Decomposition (NSD): A Fast and Scalable Approach to Network Alignment”
