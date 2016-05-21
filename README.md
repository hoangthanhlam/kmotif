# kmotif
There are two algorithms implemented in this software to find top-k mmotif
nmotif and kmotif.

To use the software compile it with g++:
g++ Motif.cpp -o motif

Then use it as follow:
./motif file data_size window_length motif_length k naive_or_kmotif query_every_window_of_size_q output_disjoint_motifs

where 

file: the input file

data_size: maximum data points 

window_length: the sliding window length

k: the number of expected motifs

naive_or_kmotif: 1 for nmotif and  0  for  kmotif

query_every_window_of_size_q: query frequency

output_disjoint_motifs: both algorithms output non-overlapping  motifs, but  not disjoint motif pairs, so in order to find 
disjoint motif pairs set this value to 1. Only works for  nmotif, not kmotif
