based on seq_min, going to build my tree
first read fasta and convert to sigs by using pbuild

tree build -f "{folder}/*.fna"
rmb ""! 


stats runs:
each runs with "max_seqCount" number of seqs.
Start by randomly sampling  "max_seqCount*runs" number of seq (or until the end)
Using this subset, random sampling again "max_seqCount" from it to run the experiments
* avoid reading the input files multiple times and is faster with large input bin file
* can always read all seqs