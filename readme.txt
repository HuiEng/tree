based on seq_min, going to build my tree
first read fasta and convert to sigs by using pbuild


current algorithm:
Introduce the "ambiguous node" to store elements that are near neighbour to the other siblings.

The ambiguous can be singleton or potentially overlapping with one or more siblings or outliers of the siblings as the node grow. The idea is that the ambiguous elements does not belong/stay to any exising clusters, but not far enough from existing clusters to form independant subtree as well. We want to avoid choosing the ambiguous elements as the centroid of a cluster because the cluster will have large distortion but small cluster size and large overlapping area with some other tight clusters.

Therefore, each node can have a list of child nodes and ONE ambiguous node. The ambiguous node stores any element that "stay" with the node/parent, but near neighbour to some or all of the other siblings. Besides the ambiguous node, the siblings must always "split" among each other to avoid overlaps at the higher level of the tree. 

When an input seq meets a branch/supercluster, it will ignore the ambiguous node for comparison. It will stay in one of the sibling if threshold is met, or create new sibling if mismatches with all siblings or form another supercluster when NN with some of the silbing. The supercluster/branch is the union of all its children. 

Then we need to reinsert the elements in the ambiguous node to see if there is any better candidate. Rotation can happen too, if not maybe the tree will be refined after reinsertion, still thinking about it. The key now is avoiding outlier to be selected as the centroid, the ambiguous node is like a place holder because NNs should be placed under the same subtree anyways, just need to find the best one, we put in the ambiguous node so that they can find a better cluster when there are more candidates later.
