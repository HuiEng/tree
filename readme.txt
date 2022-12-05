based on seq_min, going to build my tree
first read fasta and convert to sigs by using pbuild


current algorithm:
>If distance(seq,centroid)<stay_threshold
- Insert into centroid then move one

>else if distance(seq, all children of a node) > split_threshold
- Spawn new children
* to modify=> if the node is a branch, spawn sibling to node instead

>else if matching to only 1 branch
- traverse again 

>else
- separate matches & mismatches
if node is branch, dont touch matches, promote mismatches to become siblings of node
if node is not a branch, ie is a root, create new branch to merge matches, dont touch mismatches

- if NN only leaves
create new t_parent and merge NN leaves


- find nearest among NN
get the distance between seq and every NN, if NN is a branch, get the distance of its children too.
If NN is a grandchild & grandchild is branch => traverse grandchild
If NN is a grandchild & grandchild is leaf => add level	
if NN is a branch => all children are not closer to seq, create new sibling to NN branch
if NN is a leaf => add level

* Add level means 
turn a leave into a branch, clear matrices & seqIDs
copy the content the leave into a new node, make it children of the leaf
then insert new seq as the second child of the leaf
