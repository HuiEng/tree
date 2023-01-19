based on seq_min, going to build my tree
first read fasta and convert to sigs by using pbuild


current algorithm:
Divide the children of root if the count is create than the limit (5 for now)
use the first child as the first seed, and its furthest sibling as the second seed, group the children to their nearest seed

the seeds are independent subtree, not fair to be treated as normal branch due to its distortion, need to compare them in a different way, especially for the NNs. (haven't fully implement)