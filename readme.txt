based on seq_min, going to build my tree
first read fasta and convert to sigs by using pbuild

tree build -f "{folder}/*.fna"
rmb ""! 

tree tree 3 read seqs method, default=0
0: the old ways, read all seqs, then cluster "cap" seqs, only method allow force split root with -f
1: 2D seq vector, only read up to "cap" seq
2: the other methods read seqs into vector first, this method does not store seqs, re-read the file whenever need the seqs again, very slow but doable with large collection