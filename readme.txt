based on seq_min, going to build my tree
first read fasta and convert to sigs by using pbuild


current algorithm:
create ambi node to store NN (neither stay or split), ancestor only union non-ambi nodes.
Ambi nodes will be deleted/ignored for reinsertion, seqs previously in ambi nodes (and other nodes) will find the best suited nodes during reinsertion.

find closest leaves & branch, prioritise leaf, if not proceed with closest branch, prune if distance >= stay threshold for leaf