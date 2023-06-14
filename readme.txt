based on seq_min, going to build my tree
first read fasta and convert to sigs by using pbuild


########## Read multiple bins and build tree
mkdir data
tree build ../toy/toy.fasta -m -b ./data/ -k 9 -w 100
ls ./data/*.bin >list.txt
tree prim -M -i list.txt -L 0.18 -I 4 >tree.json
###########

tree build -f "{folder}/*.fna"
rmb ""! 

see similarityStatus.xlsx for method, need check
somewhat working now, haven't to check force split
tree properties:
priority: average similarity to the node centroid

Node types: leaf, ambi, branch, super, root(?)

Leaf:
- bottomest level of the tree, priority must be (more or less) greater than the stay threshold

Ambi: 
* can only exist under a branch with at least one leaf
- elements in the ambi node are near neighbour to the centroid of the leaf(s)
- has a chance to be converted into leaf if the priority of the ambi node is greater than the stay threshold
* there can be more than ambi node per branch
- create a new ambi node if the priority of the existing one is less than the split threshold

Branch: 
- Consists of leaf(s) and maybe ambiguous node(s)
- imagine a solid cluster with outliers in outer dotted cluster
*** should it be turned into a super? maybe check the priority?

Super:
- Children can only be leaf, branch and/or super
- priority can be low
- can only have two similarity status; stay or split
- stay if the similarity of the query seq and the node centroid is greater than its priority
* may have overlapping region with sibling node(s)
* priority should consider the priority of its children, not happening at the moment
*** need check. Now dissolve super if the priority is lower than the split threshold or if it is overlapping with more than half of the sibling super.
- dissolve meaning moving the children of the current super to another nearest super (or root)
- if the distance to the nearest super is less than the split threshold, move to grandparent (parent of the ori super) instead 

Root:
*** yet to be check
- designed to handle bounded tree width
- dissimilar cluster forced to placed together due to the children size limit
- distorition can be very big, aka low priority
- should be split/rotate/promote

