# kmotif
There are two algorithms implemented in this software to find top-k mmotif
nmotif and kmotif.

To use the software compile it with g++:
g++ Motif.cpp -o motif

Then use it as follow:
./motif file data_size window_length motif_length k naive_or_kmotif query_every_window_of_size_q output_disjoint_motifs
